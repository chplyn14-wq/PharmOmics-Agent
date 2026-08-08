"""Tests for pharmomics.analysis.results.

Covers:
- GeneDifferential instantiation and frozen behaviour.
- AnalysisResult minimal and full instances.
- AnalysisResult frozen behaviour.
- n_genes_tested consistency (pass and raise).
- gene_results order preservation.
- dataclasses.asdict() compatibility.
"""

from __future__ import annotations

from dataclasses import asdict

import pytest

from pharmomics.analysis.results import (
    AnalysisConsistencyError,
    AnalysisResult,
    GeneDifferential,
)

# ---------------------------------------------------------------------------
# GeneDifferential
# ---------------------------------------------------------------------------


class TestGeneDifferential:
    """GeneDifferential is a frozen dataclass with the expected fields."""

    def test_minimal_instance(self) -> None:
        g = GeneDifferential(
            gene_id="EGFR",
            log2_fold_change=1.0,
            p_value=0.001,
            adj_p_value=0.01,
            significant=True,
            base_mean=150.0,
        )
        assert g.gene_id == "EGFR"
        assert g.log2_fold_change == 1.0
        assert g.p_value == 0.001
        assert g.adj_p_value == 0.01
        assert g.significant is True
        assert g.base_mean == 150.0

    def test_full_instance(self) -> None:
        """Values at boundaries (zero fold-change, non-significant)."""
        g = GeneDifferential(
            gene_id="TP53",
            log2_fold_change=0.0,
            p_value=0.95,
            adj_p_value=0.95,
            significant=False,
            base_mean=50.0,
        )
        assert g.gene_id == "TP53"
        assert g.log2_fold_change == 0.0
        assert g.p_value == 0.95
        assert g.adj_p_value == 0.95
        assert g.significant is False
        assert g.base_mean == 50.0

    def test_frozen(self) -> None:
        g = GeneDifferential(
            gene_id="EGFR",
            log2_fold_change=1.0,
            p_value=0.001,
            adj_p_value=0.01,
            significant=True,
            base_mean=150.0,
        )
        with pytest.raises((TypeError, AttributeError)):
            g.gene_id = "ERBB2"  # type: ignore[misc]


# ---------------------------------------------------------------------------
# AnalysisResult
# ---------------------------------------------------------------------------


class TestAnalysisResult:
    """AnalysisResult validates internal consistency on construction."""

    _GENES = (
        GeneDifferential(
            gene_id="EGFR",
            log2_fold_change=1.0,
            p_value=0.001,
            adj_p_value=0.01,
            significant=True,
            base_mean=150.0,
        ),
        GeneDifferential(
            gene_id="TP53",
            log2_fold_change=0.1,
            p_value=0.8,
            adj_p_value=0.9,
            significant=False,
            base_mean=50.0,
        ),
    )

    def test_minimal_instance(self) -> None:
        r = AnalysisResult(
            analysis_type="differential_analysis",
            contrast_id="trt_vs_ctrl",
            gene_results=(),
            n_genes_tested=0,
        )
        assert r.analysis_type == "differential_analysis"
        assert r.contrast_id == "trt_vs_ctrl"
        assert r.gene_results == ()
        assert r.n_genes_tested == 0
        assert r.warnings == ()

    def test_full_instance(self) -> None:
        r = AnalysisResult(
            analysis_type="differential_analysis",
            contrast_id="trt_vs_ctrl",
            gene_results=self._GENES,
            n_genes_tested=2,
            warnings=("low replicate count",),
        )
        assert r.analysis_type == "differential_analysis"
        assert r.contrast_id == "trt_vs_ctrl"
        assert len(r.gene_results) == 2
        assert r.gene_results[0].gene_id == "EGFR"
        assert r.gene_results[1].gene_id == "TP53"
        assert r.n_genes_tested == 2
        assert r.warnings == ("low replicate count",)

    def test_frozen(self) -> None:
        r = AnalysisResult(
            analysis_type="differential_analysis",
            contrast_id="trt_vs_ctrl",
            gene_results=(),
            n_genes_tested=0,
        )
        with pytest.raises((TypeError, AttributeError)):
            r.analysis_type = "clustering"  # type: ignore[misc]

    def test_consistent_n_passes(self) -> None:
        """n_genes_tested == len(gene_results) should succeed."""
        r = AnalysisResult(
            analysis_type="differential_analysis",
            contrast_id="trt_vs_ctrl",
            gene_results=self._GENES,
            n_genes_tested=2,
        )
        assert r.n_genes_tested == len(r.gene_results)

    def test_inconsistent_n_raises(self) -> None:
        """n_genes_tested != len(gene_results) should raise."""
        with pytest.raises(AnalysisConsistencyError) as excinfo:
            AnalysisResult(
                analysis_type="differential_analysis",
                contrast_id="trt_vs_ctrl",
                gene_results=self._GENES,
                n_genes_tested=5,
            )
        assert "n_genes_tested" in str(excinfo.value)
        assert "len(gene_results)" in str(excinfo.value)

    def test_gene_results_order(self) -> None:
        """gene_results must preserve input order."""
        genes = (
            GeneDifferential(
                gene_id="C",
                log2_fold_change=0.0,
                p_value=1.0,
                adj_p_value=1.0,
                significant=False,
                base_mean=0.0,
            ),
            GeneDifferential(
                gene_id="A",
                log2_fold_change=0.0,
                p_value=1.0,
                adj_p_value=1.0,
                significant=False,
                base_mean=0.0,
            ),
            GeneDifferential(
                gene_id="B",
                log2_fold_change=0.0,
                p_value=1.0,
                adj_p_value=1.0,
                significant=False,
                base_mean=0.0,
            ),
        )
        r = AnalysisResult(
            analysis_type="differential_analysis",
            contrast_id="x",
            gene_results=genes,
            n_genes_tested=3,
        )
        assert [g.gene_id for g in r.gene_results] == ["C", "A", "B"]

    def test_warnings_preserved(self) -> None:
        """warnings tuple is accessible and immutable in spirit."""
        warnings = ("warning A", "warning B")
        r = AnalysisResult(
            analysis_type="differential_analysis",
            contrast_id="x",
            gene_results=(),
            n_genes_tested=0,
            warnings=warnings,
        )
        assert r.warnings == warnings
        assert isinstance(r.warnings, tuple)

    def test_asdict_works(self) -> None:
        """dataclasses.asdict() returns a JSON-serialisable dict."""
        r = AnalysisResult(
            analysis_type="differential_analysis",
            contrast_id="trt_vs_ctrl",
            gene_results=self._GENES,
            n_genes_tested=2,
            warnings=("low replicate count",),
        )
        d = asdict(r)
        assert d["analysis_type"] == "differential_analysis"
        assert d["contrast_id"] == "trt_vs_ctrl"
        assert d["n_genes_tested"] == 2
        assert isinstance(d["gene_results"], tuple)
        assert len(d["gene_results"]) == 2
        assert d["gene_results"][0]["gene_id"] == "EGFR"
        assert d["gene_results"][1]["gene_id"] == "TP53"
        assert d["warnings"] == ("low replicate count",)
