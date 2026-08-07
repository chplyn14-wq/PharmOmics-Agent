"""Tests for pharmomics.analysis.assembly.

Covers:
- Single-gene assembly with known values.
- Multi-gene assembly.
- base_mean = (comparison_mean + reference_mean) / 2.
- significant when adj_p_value < fdr_threshold.
- not significant when adj_p_value > fdr_threshold.
- not significant at exact threshold boundary.
- NaN adj_p_value → significant=False.
- NaN p_value and NaN log2fc propagate.
- Custom fdr_threshold from specification.
- Gene ID mismatch raises AssemblyError.
- Different lengths raises AssemblyError.
- Duplicate gene_id raises AssemblyError.
- Empty input returns empty tuple.
- Output order matches mean_expressions order.
"""

from __future__ import annotations

from math import isnan

import pytest

from pharmomics.analysis.assembly import (
    AssemblyError,
    assemble_gene_differentials,
)
from pharmomics.analysis.bh_fdr import GenePValueAdj
from pharmomics.analysis.log2_fold_change import GeneLog2FoldChange
from pharmomics.analysis.mean_expression import GeneMeanExpression
from pharmomics.analysis.per_gene_pvalue import GenePValue
from pharmomics.analysis.schemas import AnalysisSpecification

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _mean(gene_id: str, comp: float, ref: float) -> GeneMeanExpression:
    return GeneMeanExpression(gene_id=gene_id, comparison_mean=comp, reference_mean=ref)


def _log2fc(gene_id: str, lfc: float) -> GeneLog2FoldChange:
    return GeneLog2FoldChange(gene_id=gene_id, log2fc=lfc)


def _pval(gene_id: str, p: float) -> GenePValue:
    return GenePValue(gene_id=gene_id, p_value=p)


def _adj(gene_id: str, raw: float, adj: float) -> GenePValueAdj:
    return GenePValueAdj(gene_id=gene_id, raw_p_value=raw, adj_p_value=adj)


def _spec(threshold: float | None = None) -> AnalysisSpecification:
    params = {}
    if threshold is not None:
        params["fdr_threshold"] = threshold
    return AnalysisSpecification(
        analysis_type="differential_analysis",
        parameters=params,
    )


# ---------------------------------------------------------------------------
# Normal assembly
# ---------------------------------------------------------------------------


class TestAssembleSingleGene:
    """Single gene assembly with known values."""

    def test_single_gene_all_fields(self) -> None:
        """All fields correctly populated."""
        means = (_mean("EGFR", 210.0, 105.0),)
        lfc = (_log2fc("EGFR", 1.0),)
        pval = (_pval("EGFR", 0.001),)
        adj = (_adj("EGFR", 0.001, 0.002),)

        result = assemble_gene_differentials(means, lfc, pval, adj, _spec())

        assert len(result) == 1
        g = result[0]
        assert g.gene_id == "EGFR"
        assert g.log2_fold_change == 1.0
        assert g.p_value == 0.001
        assert g.adj_p_value == 0.002
        assert g.significant is True
        assert g.base_mean == 157.5  # (210 + 105) / 2


class TestAssembleMultipleGenes:
    """Multi-gene assembly."""

    def test_three_genes(self) -> None:
        """Three genes with mixed significance."""
        means = (
            _mean("G1", 100.0, 50.0),
            _mean("G2", 80.0, 75.0),
            _mean("G3", 200.0, 100.0),
        )
        lfc = (
            _log2fc("G1", 1.0),
            _log2fc("G2", 0.09),
            _log2fc("G3", 1.0),
        )
        pval = (
            _pval("G1", 0.001),
            _pval("G2", 0.8),
            _pval("G3", 0.01),
        )
        adj = (
            _adj("G1", 0.001, 0.003),
            _adj("G2", 0.8, 0.8),
            _adj("G3", 0.01, 0.015),
        )

        result = assemble_gene_differentials(means, lfc, pval, adj, _spec())

        assert len(result) == 3
        assert result[0].gene_id == "G1"
        assert result[0].significant is True
        assert result[1].gene_id == "G2"
        assert result[1].significant is False
        assert result[2].gene_id == "G3"
        assert result[2].significant is True

    def test_returns_tuple(self) -> None:
        """Return type must be tuple."""
        means = (_mean("G0", 10.0, 5.0),)
        lfc = (_log2fc("G0", 1.0),)
        pval = (_pval("G0", 0.01),)
        adj = (_adj("G0", 0.01, 0.02),)
        result = assemble_gene_differentials(means, lfc, pval, adj, _spec())
        assert isinstance(result, tuple)


# ---------------------------------------------------------------------------
# base_mean
# ---------------------------------------------------------------------------


class TestBaseMean:
    """base_mean = (comparison_mean + reference_mean) / 2."""

    def test_symmetric_means(self) -> None:
        """comp=100, ref=50 → base_mean=75."""
        means = (_mean("G0", 100.0, 50.0),)
        lfc = (_log2fc("G0", 1.0),)
        pval = (_pval("G0", 0.01),)
        adj = (_adj("G0", 0.01, 0.02),)
        result = assemble_gene_differentials(means, lfc, pval, adj, _spec())
        assert result[0].base_mean == 75.0

    def test_equal_means(self) -> None:
        """comp=42, ref=42 → base_mean=42."""
        means = (_mean("G0", 42.0, 42.0),)
        lfc = (_log2fc("G0", 0.0),)
        pval = (_pval("G0", 0.5),)
        adj = (_adj("G0", 0.5, 0.5),)
        result = assemble_gene_differentials(means, lfc, pval, adj, _spec())
        assert result[0].base_mean == 42.0

    def test_zero_means(self) -> None:
        """comp=0, ref=0 → base_mean=0."""
        means = (_mean("G0", 0.0, 0.0),)
        lfc = (_log2fc("G0", float("nan")),)
        pval = (_pval("G0", float("nan")),)
        adj = (_adj("G0", float("nan"), float("nan")),)
        result = assemble_gene_differentials(means, lfc, pval, adj, _spec())
        assert result[0].base_mean == 0.0


# ---------------------------------------------------------------------------
# significant
# ---------------------------------------------------------------------------


class TestSignificant:
    """significant = adj_p_value < fdr_threshold."""

    def test_below_threshold(self) -> None:
        """adj=0.01, threshold=0.05 → True."""
        means = (_mean("G0", 10.0, 5.0),)
        lfc = (_log2fc("G0", 1.0),)
        pval = (_pval("G0", 0.005),)
        adj = (_adj("G0", 0.005, 0.01),)
        result = assemble_gene_differentials(means, lfc, pval, adj, _spec())
        assert result[0].significant is True

    def test_above_threshold(self) -> None:
        """adj=0.3, threshold=0.05 → False."""
        means = (_mean("G0", 10.0, 5.0),)
        lfc = (_log2fc("G0", 1.0),)
        pval = (_pval("G0", 0.2),)
        adj = (_adj("G0", 0.2, 0.3),)
        result = assemble_gene_differentials(means, lfc, pval, adj, _spec())
        assert result[0].significant is False

    def test_at_threshold_boundary(self) -> None:
        """adj=0.05, threshold=0.05 → False (strict less-than)."""
        means = (_mean("G0", 10.0, 5.0),)
        lfc = (_log2fc("G0", 1.0),)
        pval = (_pval("G0", 0.025),)
        adj = (_adj("G0", 0.025, 0.05),)
        result = assemble_gene_differentials(means, lfc, pval, adj, _spec())
        assert result[0].significant is False

    def test_custom_fdr_threshold(self) -> None:
        """threshold=0.01, adj=0.02 → False (would be True at default 0.05)."""
        means = (_mean("G0", 10.0, 5.0),)
        lfc = (_log2fc("G0", 1.0),)
        pval = (_pval("G0", 0.008),)
        adj = (_adj("G0", 0.008, 0.02),)
        result = assemble_gene_differentials(
            means, lfc, pval, adj, _spec(threshold=0.01)
        )
        assert result[0].significant is False

    def test_custom_stricter_threshold(self) -> None:
        """threshold=0.001, adj=0.0005 → True."""
        means = (_mean("G0", 10.0, 5.0),)
        lfc = (_log2fc("G0", 1.0),)
        pval = (_pval("G0", 0.0002),)
        adj = (_adj("G0", 0.0002, 0.0005),)
        result = assemble_gene_differentials(
            means, lfc, pval, adj, _spec(threshold=0.001)
        )
        assert result[0].significant is True


# ---------------------------------------------------------------------------
# NaN handling
# ---------------------------------------------------------------------------


class TestNaNHandling:
    """NaN values propagate correctly."""

    def test_nan_adj_p_value_not_significant(self) -> None:
        """NaN adj_p_value → significant=False."""
        means = (_mean("G0", 0.0, 10.0),)
        lfc = (_log2fc("G0", float("nan")),)
        pval = (_pval("G0", float("nan")),)
        adj = (_adj("G0", float("nan"), float("nan")),)
        result = assemble_gene_differentials(means, lfc, pval, adj, _spec())
        assert result[0].significant is False
        assert isnan(result[0].log2_fold_change)
        assert isnan(result[0].p_value)
        assert isnan(result[0].adj_p_value)

    def test_nan_p_value_propagated(self) -> None:
        """NaN p_value preserved in output."""
        means = (_mean("G0", 10.0, 5.0),)
        lfc = (_log2fc("G0", 1.0),)
        pval = (_pval("G0", float("nan")),)
        adj = (_adj("G0", float("nan"), float("nan")),)
        result = assemble_gene_differentials(means, lfc, pval, adj, _spec())
        assert isnan(result[0].p_value)
        assert isnan(result[0].adj_p_value)

    def test_nan_log2fc_propagated(self) -> None:
        """NaN log2fc preserved in output."""
        means = (_mean("G0", 0.0, 5.0),)
        lfc = (_log2fc("G0", float("nan")),)
        pval = (_pval("G0", 0.01),)
        adj = (_adj("G0", 0.01, 0.02),)
        result = assemble_gene_differentials(means, lfc, pval, adj, _spec())
        assert isnan(result[0].log2_fold_change)
        # p-value is still valid
        assert result[0].p_value == 0.01
        assert result[0].significant is True

    def test_mixed_valid_and_nan_genes(self) -> None:
        """G0 valid, G1 NaN, G2 valid."""
        means = (
            _mean("G0", 100.0, 50.0),
            _mean("G1", 0.0, 10.0),
            _mean("G2", 200.0, 100.0),
        )
        lfc = (
            _log2fc("G0", 1.0),
            _log2fc("G1", float("nan")),
            _log2fc("G2", 1.0),
        )
        pval = (
            _pval("G0", 0.001),
            _pval("G1", float("nan")),
            _pval("G2", 0.002),
        )
        adj = (
            _adj("G0", 0.001, 0.003),
            _adj("G1", float("nan"), float("nan")),
            _adj("G2", 0.002, 0.004),
        )
        result = assemble_gene_differentials(means, lfc, pval, adj, _spec())
        assert len(result) == 3
        assert result[0].significant is True
        assert result[1].significant is False
        assert isnan(result[1].adj_p_value)
        assert result[2].significant is True


# ---------------------------------------------------------------------------
# Error cases
# ---------------------------------------------------------------------------


class TestGeneIdMismatch:
    """Gene ID sets must match across all inputs."""

    def test_mean_vs_log2fc_mismatch(self) -> None:
        """Different gene in log2fc."""
        means = (_mean("G0", 10.0, 5.0),)
        lfc = (_log2fc("G1", 1.0),)
        pval = (_pval("G0", 0.01),)
        adj = (_adj("G0", 0.01, 0.02),)
        with pytest.raises(AssemblyError, match="mean_expressions vs log2fc_results"):
            assemble_gene_differentials(means, lfc, pval, adj, _spec())

    def test_mean_vs_pval_mismatch(self) -> None:
        """Different gene in pvalue_results."""
        means = (_mean("G0", 10.0, 5.0),)
        lfc = (_log2fc("G0", 1.0),)
        pval = (_pval("G1", 0.01),)
        adj = (_adj("G0", 0.01, 0.02),)
        with pytest.raises(AssemblyError, match="mean_expressions vs pvalue_results"):
            assemble_gene_differentials(means, lfc, pval, adj, _spec())

    def test_mean_vs_adj_mismatch(self) -> None:
        """Different gene in adjusted_results."""
        means = (_mean("G0", 10.0, 5.0),)
        lfc = (_log2fc("G0", 1.0),)
        pval = (_pval("G0", 0.01),)
        adj = (_adj("G1", 0.01, 0.02),)
        with pytest.raises(AssemblyError, match="mean_expressions vs adjusted_results"):
            assemble_gene_differentials(means, lfc, pval, adj, _spec())

    def test_extra_gene_in_one_input(self) -> None:
        """One input has an extra gene."""
        means = (_mean("G0", 10.0, 5.0), _mean("G1", 20.0, 10.0))
        lfc = (_log2fc("G0", 1.0), _log2fc("G1", 1.0))
        pval = (_pval("G0", 0.01), _pval("G1", 0.02))
        adj = (_adj("G0", 0.01, 0.02),)  # Missing G1
        with pytest.raises(AssemblyError, match="missing.*G1"):
            assemble_gene_differentials(means, lfc, pval, adj, _spec())


class TestDuplicateGeneId:
    """Duplicate gene_id within a single input raises."""

    def test_duplicate_in_mean_expressions(self) -> None:
        means = (_mean("G0", 10.0, 5.0), _mean("G0", 20.0, 10.0))
        lfc = (_log2fc("G0", 1.0),)
        pval = (_pval("G0", 0.01),)
        adj = (_adj("G0", 0.01, 0.02),)
        with pytest.raises(AssemblyError, match="Duplicate gene_id.*G0"):
            assemble_gene_differentials(means, lfc, pval, adj, _spec())

    def test_duplicate_in_log2fc(self) -> None:
        means = (_mean("G0", 10.0, 5.0),)
        lfc = (_log2fc("G0", 1.0), _log2fc("G0", 2.0))
        pval = (_pval("G0", 0.01),)
        adj = (_adj("G0", 0.01, 0.02),)
        with pytest.raises(AssemblyError, match="Duplicate gene_id.*G0"):
            assemble_gene_differentials(means, lfc, pval, adj, _spec())


class TestEmptyInput:
    """Empty input returns empty tuple."""

    def test_empty_returns_empty(self) -> None:
        """All empty → empty."""
        result = assemble_gene_differentials((), (), (), (), _spec())
        assert result == ()
        assert isinstance(result, tuple)


class TestOrderPreservation:
    """Output order matches mean_expressions order."""

    def test_output_order_matches_mean_expressions(self) -> None:
        """Input: G_Z, G_A, G_M → output: G_Z, G_A, G_M."""
        means = (
            _mean("G_Z", 30.0, 15.0),
            _mean("G_A", 10.0, 5.0),
            _mean("G_M", 20.0, 10.0),
        )
        lfc = (
            _log2fc("G_Z", 1.0),
            _log2fc("G_A", 1.0),
            _log2fc("G_M", 1.0),
        )
        pval = (
            _pval("G_Z", 0.01),
            _pval("G_A", 0.01),
            _pval("G_M", 0.01),
        )
        adj = (
            _adj("G_Z", 0.01, 0.02),
            _adj("G_A", 0.01, 0.02),
            _adj("G_M", 0.01, 0.02),
        )
        result = assemble_gene_differentials(means, lfc, pval, adj, _spec())
        assert [g.gene_id for g in result] == ["G_Z", "G_A", "G_M"]
