"""Tests for pharmomics.analysis.log2_fold_change."""

from __future__ import annotations

from math import isnan, log2

import pytest

from pharmomics.analysis.log2_fold_change import (
    GeneLog2FoldChange,
    Log2FoldChangeError,
    compute_log2_fold_change,
)
from pharmomics.analysis.mean_expression import GeneMeanExpression

# ---------------------------------------------------------------------------
# compute_log2_fold_change — normal cases
# ---------------------------------------------------------------------------


class TestComputeLog2FoldChangeNormal:
    """Verify correct log2 fold change computation."""

    def test_simple_upregulation(self) -> None:
        """comparison=200, reference=100 -> log2(2.0) = 1.0."""
        entries = (
            GeneMeanExpression(
                gene_id="G0", comparison_mean=200.0, reference_mean=100.0
            ),
        )
        result = compute_log2_fold_change(entries)
        assert len(result) == 1
        assert result[0].gene_id == "G0"
        assert result[0].log2fc == pytest.approx(1.0)

    def test_simple_downregulation(self) -> None:
        """comparison=50, reference=100 -> log2(0.5) = -1.0."""
        entries = (
            GeneMeanExpression(
                gene_id="G1", comparison_mean=50.0, reference_mean=100.0
            ),
        )
        result = compute_log2_fold_change(entries)
        assert result[0].log2fc == pytest.approx(-1.0)

    def test_no_change(self) -> None:
        """comparison == reference -> log2(1.0) = 0.0."""
        entries = (
            GeneMeanExpression(gene_id="G0", comparison_mean=42.0, reference_mean=42.0),
        )
        result = compute_log2_fold_change(entries)
        assert result[0].log2fc == pytest.approx(0.0)

    def test_multiple_genes(self) -> None:
        """Multiple genes return correct log2FC values."""
        entries = (
            GeneMeanExpression(
                gene_id="EGFR", comparison_mean=210.0, reference_mean=105.0
            ),
            GeneMeanExpression(
                gene_id="TP53", comparison_mean=50.0, reference_mean=55.0
            ),
        )
        result = compute_log2_fold_change(entries)
        assert len(result) == 2
        assert result[0].gene_id == "EGFR"
        assert result[0].log2fc == pytest.approx(log2(210.0 / 105.0))
        assert result[1].gene_id == "TP53"
        assert result[1].log2fc == pytest.approx(log2(50.0 / 55.0))

    def test_preserves_feature_order(self) -> None:
        """Result order must match input order."""
        entries = (
            GeneMeanExpression(gene_id="G_Z", comparison_mean=10.0, reference_mean=5.0),
            GeneMeanExpression(
                gene_id="G_A", comparison_mean=20.0, reference_mean=10.0
            ),
            GeneMeanExpression(
                gene_id="G_M", comparison_mean=30.0, reference_mean=15.0
            ),
        )
        result = compute_log2_fold_change(entries)
        assert [r.gene_id for r in result] == ["G_Z", "G_A", "G_M"]

    def test_returns_tuple(self) -> None:
        """Return type must be tuple."""
        entries = (
            GeneMeanExpression(gene_id="G0", comparison_mean=1.0, reference_mean=2.0),
        )
        result = compute_log2_fold_change(entries)
        assert isinstance(result, tuple)

    def test_single_gene_single_sample_means(self) -> None:
        """Regression-style: (4+6)/2=5 vs (8+12)/2=10 -> log2(0.5)=-1."""
        entries = (
            GeneMeanExpression(gene_id="G0", comparison_mean=5.0, reference_mean=10.0),
        )
        result = compute_log2_fold_change(entries)
        assert result[0].log2fc == pytest.approx(-1.0)


# ---------------------------------------------------------------------------
# compute_log2_fold_change — zero / negative cases
# ---------------------------------------------------------------------------


class TestComputeLog2FoldChangeZeroNegative:
    """Verify NaN is returned for zero or negative means."""

    def test_zero_comparison_mean(self) -> None:
        """comparison=0 -> NaN."""
        entries = (
            GeneMeanExpression(gene_id="G0", comparison_mean=0.0, reference_mean=10.0),
        )
        result = compute_log2_fold_change(entries)
        assert isnan(result[0].log2fc)

    def test_zero_reference_mean(self) -> None:
        """reference=0 -> NaN."""
        entries = (
            GeneMeanExpression(gene_id="G0", comparison_mean=10.0, reference_mean=0.0),
        )
        result = compute_log2_fold_change(entries)
        assert isnan(result[0].log2fc)

    def test_both_zero(self) -> None:
        """Both means zero -> NaN."""
        entries = (
            GeneMeanExpression(gene_id="G0", comparison_mean=0.0, reference_mean=0.0),
        )
        result = compute_log2_fold_change(entries)
        assert isnan(result[0].log2fc)

    def test_negative_comparison_mean(self) -> None:
        """Negative comparison -> NaN."""
        entries = (
            GeneMeanExpression(gene_id="G0", comparison_mean=-5.0, reference_mean=10.0),
        )
        result = compute_log2_fold_change(entries)
        assert isnan(result[0].log2fc)

    def test_negative_reference_mean(self) -> None:
        """Negative reference -> NaN."""
        entries = (
            GeneMeanExpression(gene_id="G0", comparison_mean=10.0, reference_mean=-3.0),
        )
        result = compute_log2_fold_change(entries)
        assert isnan(result[0].log2fc)

    def test_mixed_valid_and_invalid(self) -> None:
        """Valid genes return values; invalid genes return NaN."""
        entries = (
            GeneMeanExpression(
                gene_id="G_VALID", comparison_mean=8.0, reference_mean=4.0
            ),
            GeneMeanExpression(
                gene_id="G_ZERO", comparison_mean=0.0, reference_mean=5.0
            ),
            GeneMeanExpression(
                gene_id="G_NEG", comparison_mean=3.0, reference_mean=-1.0
            ),
        )
        result = compute_log2_fold_change(entries)
        assert len(result) == 3
        assert result[0].log2fc == pytest.approx(1.0)
        assert isnan(result[1].log2fc)
        assert isnan(result[2].log2fc)


# ---------------------------------------------------------------------------
# compute_log2_fold_change — error cases
# ---------------------------------------------------------------------------


class TestComputeLog2FoldChangeErrors:
    """Verify error handling for invalid inputs."""

    def test_raises_empty_input(self) -> None:
        """Empty tuple raises Log2FoldChangeError."""
        with pytest.raises(Log2FoldChangeError, match="No mean expression data"):
            compute_log2_fold_change(())


# ---------------------------------------------------------------------------
# GeneLog2FoldChange model
# ---------------------------------------------------------------------------


class TestGeneLog2FoldChangeModel:
    """Verify the dataclass behavior."""

    def test_is_frozen(self) -> None:
        """GeneLog2FoldChange instances must be immutable."""
        entry = GeneLog2FoldChange(gene_id="G0", log2fc=1.0)
        with pytest.raises((TypeError, AttributeError)):
            entry.gene_id = "changed"  # type: ignore[misc]

    def test_nan_equality_identity(self) -> None:
        """Two entries with NaN log2fc are distinct instances."""
        a = GeneLog2FoldChange(gene_id="G0", log2fc=float("nan"))
        b = GeneLog2FoldChange(gene_id="G0", log2fc=float("nan"))
        # NaN != NaN by IEEE semantics, but frozen dataclass compares fields
        # so a == b will be False for the log2fc field (nan != nan).
        assert a.gene_id == b.gene_id
