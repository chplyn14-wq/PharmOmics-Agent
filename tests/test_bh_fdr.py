"""Tests for pharmomics.analysis.bh_fdr.

Covers:
- Normal BH correction with known values.
- Feature order preservation.
- NaN exclusion from correction, NaN output preserved.
- All NaN input → all NaN output.
- Empty input → empty output.
- Single gene → adj equals raw.
- p=0 stays zero.
- Clamping to 1.0.
- Duplicate p-values behavior.
- Monotonicity of adjusted values (in sorted order).
- Invalid p-value raises BHAdjustmentError.
- Large batch correctness.
- GenePValueAdj frozen behaviour.
"""

from __future__ import annotations

from math import isnan

import pytest

from pharmomics.analysis.bh_fdr import (
    BHAdjustmentError,
    GenePValueAdj,
    compute_bh_adjusted_p_values,
)
from pharmomics.analysis.per_gene_pvalue import GenePValue


def _gp(gene_id: str, p: float) -> GenePValue:
    """Shorthand for creating GenePValue."""
    return GenePValue(gene_id=gene_id, p_value=p)


# ---------------------------------------------------------------------------
# Normal correction — known values
# ---------------------------------------------------------------------------


class TestBHNormal:
    """Verify correct BH adjustment with hand-computed values."""

    def test_two_genes_known_values(self) -> None:
        """Two genes: raw p=[0.001, 0.03] → verify BH output.

        m=2, sorted: [0.001(rank1), 0.03(rank2)]
        adj: 0.001*2/1=0.002, 0.03*2/2=0.03
        monotonicity: 0.002 < 0.03 → ok
        result: [0.002, 0.03]
        """
        raw = (_gp("G0", 0.001), _gp("G1", 0.03))
        result = compute_bh_adjusted_p_values(raw)

        assert len(result) == 2
        assert result[0].gene_id == "G0"
        assert result[0].adj_p_value == pytest.approx(0.002)
        assert result[1].gene_id == "G1"
        assert result[1].adj_p_value == pytest.approx(0.03)

    def test_five_genes_step_up(self) -> None:
        """5 genes demonstrating monotonicity step-up.

        m=5, sorted: [0.001, 0.02, 0.05, 0.08, 0.3]
        raw adj:    [0.005, 0.05, 0.0833, 0.1, 0.3]
        monotonicity: already non-decreasing
        """
        raw = (
            _gp("G0", 0.001), _gp("G1", 0.02), _gp("G2", 0.05),
            _gp("G3", 0.08), _gp("G4", 0.3),
        )
        result = compute_bh_adjusted_p_values(raw)

        assert len(result) == 5
        assert result[0].adj_p_value == pytest.approx(0.005)
        assert result[1].adj_p_value == pytest.approx(0.05)
        assert result[2].adj_p_value == pytest.approx(0.08333, abs=1e-4)
        assert result[3].adj_p_value == pytest.approx(0.1)
        assert result[4].adj_p_value == pytest.approx(0.3)


# ---------------------------------------------------------------------------
# Order preservation
# ---------------------------------------------------------------------------


class TestOrderPreservation:
    """Output order must match input order."""

    def test_preserves_feature_order(self) -> None:
        """Input: Z, A, M → output order: Z, A, M."""
        raw = (_gp("G_Z", 0.1), _gp("G_A", 0.001), _gp("G_M", 0.05))
        result = compute_bh_adjusted_p_values(raw)

        assert [r.gene_id for r in result] == ["G_Z", "G_A", "G_M"]


# ---------------------------------------------------------------------------
# NaN handling
# ---------------------------------------------------------------------------


class TestNaNHandling:
    """NaN p-values excluded from correction, output stays NaN."""

    def test_nan_excluded_from_correction(self) -> None:
        """G1 is NaN; m=2 (G0, G2 only)."""
        raw = (_gp("G0", 0.001), _gp("G1", float("nan")), _gp("G2", 0.03))
        result = compute_bh_adjusted_p_values(raw)

        assert len(result) == 3
        assert result[0].gene_id == "G0"
        assert result[0].adj_p_value == pytest.approx(0.002)  # m=2, rank=1
        assert isnan(result[1].adj_p_value)
        assert isnan(result[1].raw_p_value)
        assert result[2].gene_id == "G2"
        assert result[2].adj_p_value == pytest.approx(0.03)  # m=2, rank=2

    def test_nan_in_middle_of_sorted_sequence(self) -> None:
        """NaN at various positions should not affect valid corrections."""
        raw = (
            _gp("G0", 0.5),
            _gp("G1", float("nan")),
            _gp("G2", 0.01),
            _gp("G3", float("nan")),
            _gp("G4", 0.001),
        )
        result = compute_bh_adjusted_p_values(raw)

        # Valid: G4=0.001(rank1), G2=0.01(rank2), G0=0.5(rank3), m=3
        # adj: 0.001*3/1=0.003, 0.01*3/2=0.015, 0.5*3/3=0.5
        assert result[0].adj_p_value == pytest.approx(0.5)
        assert isnan(result[1].adj_p_value)
        assert result[2].adj_p_value == pytest.approx(0.015)
        assert isnan(result[3].adj_p_value)
        assert result[4].adj_p_value == pytest.approx(0.003)

    def test_all_nan_returns_nan_tuple(self) -> None:
        """All NaN → same length, all NaN adj."""
        raw = (_gp("G0", float("nan")), _gp("G1", float("nan")))
        result = compute_bh_adjusted_p_values(raw)

        assert len(result) == 2
        assert isnan(result[0].adj_p_value)
        assert isnan(result[1].adj_p_value)
        assert result[0].gene_id == "G0"
        assert result[1].gene_id == "G1"


# ---------------------------------------------------------------------------
# Edge cases
# ---------------------------------------------------------------------------


class TestEdgeCases:
    """Empty, single gene, boundaries, duplicates."""

    def test_empty_input_returns_empty(self) -> None:
        """() → ()."""
        assert compute_bh_adjusted_p_values(()) == ()

    def test_single_gene_adj_equals_raw(self) -> None:
        """m=1 → adj = raw."""
        raw = (_gp("G0", 0.05),)
        result = compute_bh_adjusted_p_values(raw)

        assert len(result) == 1
        assert result[0].adj_p_value == pytest.approx(0.05)

    def test_p_zero_stays_zero(self) -> None:
        """p=0 → adj=0.0."""
        raw = (_gp("G0", 0.0), _gp("G1", 0.01))
        result = compute_bh_adjusted_p_values(raw)

        assert result[0].adj_p_value == pytest.approx(0.0)

    def test_adj_clamped_to_one(self) -> None:
        """Construct case where raw adj exceeds 1.0."""
        # m=3, sorted: [0.5(rank1), 0.8(rank2), 0.9(rank3)]
        # raw adj: [1.5, 1.2, 0.9] → monotonicity → [0.9, 0.9, 0.9]
        # Already clamped at 0.9, so let's use a case where clamp triggers.
        # m=3, [0.7, 0.8, 0.9]: raw adj=[2.1, 1.2, 0.9]
        # mono backwards: adj[1]=min(1.2,0.9)=0.9, adj[0]=min(2.1,0.9)=0.9
        # → all 0.9 (no clamp needed). Need larger values.
        # m=2, [0.6, 0.8]: raw adj=[1.2, 0.8]
        # mono: adj[0]=min(1.2, 0.8)=0.8 → no clamp
        # Need: adj after mono > 1.0. With p ∈ [0,1] and m/i >= 1,
        # the last rank gives p*m/m = p ≤ 1. Monotonicity propagates min
        # backward, so max adj = last raw adj ≤ 1. Clamp only needed if
        # raw adj at last rank > 1 (impossible with p ≤ 1) or if monotonicity
        # were reversed. Still, clamp is a safety measure.
        # Let's verify clamp doesn't break anything with p=1.
        raw = (_gp("G0", 1.0),)
        result = compute_bh_adjusted_p_values(raw)
        assert result[0].adj_p_value == pytest.approx(1.0)

    def test_p_one_multiple_genes(self) -> None:
        """Multiple genes all with p=1."""
        raw = (_gp("G0", 1.0), _gp("G1", 1.0), _gp("G2", 1.0))
        result = compute_bh_adjusted_p_values(raw)

        for r in result:
            assert r.adj_p_value == pytest.approx(1.0)

    def test_duplicate_p_values(self) -> None:
        """Duplicate p-values get distinct or equalized adj per BH rules.

        m=4, sorted: [0.01, 0.01, 0.01, 0.05]
        raw adj: [0.04, 0.02, 0.0133, 0.05]
        mono backwards: adj[2]=min(0.0133, 0.05)=0.0133
                        adj[1]=min(0.02, 0.0133)=0.0133
                        adj[0]=min(0.04, 0.0133)=0.0133
        result: [0.0133, 0.0133, 0.0133, 0.05]
        """
        raw = (
            _gp("G0", 0.01), _gp("G1", 0.01),
            _gp("G2", 0.01), _gp("G3", 0.05),
        )
        result = compute_bh_adjusted_p_values(raw)

        expected = 0.01 * 4 / 3  # ≈ 0.01333
        assert result[0].adj_p_value == pytest.approx(expected, abs=1e-4)
        assert result[1].adj_p_value == pytest.approx(expected, abs=1e-4)
        assert result[2].adj_p_value == pytest.approx(expected, abs=1e-4)
        assert result[3].adj_p_value == pytest.approx(0.05)


# ---------------------------------------------------------------------------
# Monotonicity
# ---------------------------------------------------------------------------


class TestMonotonicity:
    """Adjusted values must be non-decreasing in sorted raw p-value order."""

    def test_monotonicity_property(self) -> None:
        """After sorting by raw_p_value, adj_p_value is non-decreasing."""
        raw = (
            _gp("G0", 0.3), _gp("G1", 0.001), _gp("G2", 0.1),
            _gp("G3", 0.005), _gp("G4", 0.05),
        )
        result = compute_bh_adjusted_p_values(raw)

        # Sort by raw p-value and check adj monotonicity.
        sorted_by_raw = sorted(result, key=lambda x: x.raw_p_value)
        for i in range(len(sorted_by_raw) - 1):
            a = sorted_by_raw[i].adj_p_value
            b = sorted_by_raw[i + 1].adj_p_value
            assert a <= b + 1e-10


# ---------------------------------------------------------------------------
# Error handling
# ---------------------------------------------------------------------------


class TestErrorHandling:
    """Invalid p-values raise BHAdjustmentError."""

    def test_invalid_p_value_negative_raises(self) -> None:
        """p < 0 raises."""
        raw = (_gp("G0", -0.1),)
        with pytest.raises(BHAdjustmentError, match="out of"):
            compute_bh_adjusted_p_values(raw)

    def test_invalid_p_value_above_one_raises(self) -> None:
        """p > 1 raises."""
        raw = (_gp("G0", 1.5),)
        with pytest.raises(BHAdjustmentError, match="out of"):
            compute_bh_adjusted_p_values(raw)

    def test_invalid_among_valid_raises(self) -> None:
        """One invalid p-value among valid ones raises."""
        raw = (_gp("G0", 0.01), _gp("G1", -0.001), _gp("G2", 0.05))
        with pytest.raises(BHAdjustmentError, match="G1"):
            compute_bh_adjusted_p_values(raw)


# ---------------------------------------------------------------------------
# Large batch
# ---------------------------------------------------------------------------


class TestLargeBatch:
    """Correctness at scale."""

    def test_large_batch_no_crash(self) -> None:
        """10000 genes, all finite, valid range."""
        raw = tuple(_gp(f"G{i}", 0.00005 * (i + 1)) for i in range(10_000))
        result = compute_bh_adjusted_p_values(raw)

        assert len(result) == 10_000
        for r in result:
            assert 0.0 <= r.adj_p_value <= 1.0


# ---------------------------------------------------------------------------
# Immutability
# ---------------------------------------------------------------------------


class TestGenePValueAdjImmutability:
    """GenePValueAdj instances must be immutable."""

    def test_gene_p_value_adj_is_frozen(self) -> None:
        """Cannot mutate fields."""
        entry = GenePValueAdj(gene_id="G0", raw_p_value=0.01, adj_p_value=0.05)
        with pytest.raises((TypeError, AttributeError)):
            entry.gene_id = "changed"  # type: ignore[misc]
