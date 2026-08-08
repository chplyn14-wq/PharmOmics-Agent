"""Tests for pharmomics.analysis.render TSV export functions.

Covers:
- render_results_tsv() produces correct TSV for all genes.
- render_significant_genes_tsv() produces correct TSV for significant genes.
- Header columns are correct and in fixed order.
- Gene order preserved from AnalysisResult.gene_results.
- Values come from GeneDifferential (no recalculation).
- significant bool exports as True/False string.
- NaN values export as "NaN".
- significant_genes.tsv contains only significant=True genes.
- Relative order of significant genes is preserved.
- Empty significant set produces header-only TSV.
- End-to-end with demo fixture: TSV p-value / adj_p_value are finite or NaN.
"""

from __future__ import annotations

from pharmomics.analysis.render import (
    render_results_tsv,
    render_significant_genes_tsv,
)
from pharmomics.analysis.results import AnalysisResult, GeneDifferential

# ---------------------------------------------------------------------------
# Test fixture helpers
# ---------------------------------------------------------------------------


def _make_result(
    *,
    gene_results: tuple[GeneDifferential, ...] = (),
    n_genes_tested: int = 0,
    warnings: tuple[str, ...] = (),
) -> AnalysisResult:
    """Create a minimal AnalysisResult for testing."""
    return AnalysisResult(
        analysis_type="differential_analysis",
        contrast_id="trt_vs_ctrl",
        gene_results=gene_results,
        n_genes_tested=n_genes_tested,
        warnings=warnings,
    )


def _make_gene(
    *,
    gene_id: str = "EGFR",
    log2_fold_change: float = 1.0,
    p_value: float = 0.001,
    adj_p_value: float = 0.01,
    significant: bool = True,
    base_mean: float = 150.0,
) -> GeneDifferential:
    """Create a single GeneDifferential for testing."""
    return GeneDifferential(
        gene_id=gene_id,
        log2_fold_change=log2_fold_change,
        p_value=p_value,
        adj_p_value=adj_p_value,
        significant=significant,
        base_mean=base_mean,
    )


# ---------------------------------------------------------------------------
# results.tsv header
# ---------------------------------------------------------------------------


class TestResultsTsvHeader:
    """Verify results.tsv header columns."""

    def test_header_columns(self) -> None:
        r = _make_result()
        tsv = render_results_tsv(r)
        header = tsv.split("\n")[0]
        expected = (
            "gene_id\tlog2_fold_change\tp_value\tadj_p_value\tsignificant\tbase_mean"
        )
        assert header == expected

    def test_header_tab_separated(self) -> None:
        r = _make_result()
        tsv = render_results_tsv(r)
        header = tsv.split("\n")[0]
        assert "\t" in header
        assert header.count("\t") == 5

    def test_header_ends_with_newline(self) -> None:
        r = _make_result()
        tsv = render_results_tsv(r)
        assert tsv.endswith("\n")


# ---------------------------------------------------------------------------
# results.tsv content
# ---------------------------------------------------------------------------


class TestResultsTsvContent:
    """Verify results.tsv content correctness."""

    def test_all_genes_present(self) -> None:
        genes = tuple(_make_gene(gene_id=f"G{i}") for i in range(5))
        r = _make_result(gene_results=genes, n_genes_tested=5)
        tsv = render_results_tsv(r)
        lines = [ln for ln in tsv.split("\n") if ln]
        assert len(lines) == 6  # header + 5 genes

    def test_gene_order_preserved(self) -> None:
        genes = (
            _make_gene(gene_id="C", significant=False),
            _make_gene(gene_id="A", significant=False),
            _make_gene(gene_id="B", significant=False),
        )
        r = _make_result(gene_results=genes, n_genes_tested=3)
        tsv = render_results_tsv(r)
        lines = [ln for ln in tsv.split("\n") if ln][1:]  # skip header
        assert lines[0].startswith("C")
        assert lines[1].startswith("A")
        assert lines[2].startswith("B")

    def test_values_from_gene_differential(self) -> None:
        """Values must be the raw GeneDifferential values, not recalculated."""
        g = _make_gene(
            gene_id="TEST1",
            log2_fold_change=-2.5,
            p_value=0.003,
            adj_p_value=0.015,
            base_mean=200.5,
        )
        r = _make_result(gene_results=(g,), n_genes_tested=1)
        tsv = render_results_tsv(r)
        line = tsv.split("\n")[1]
        fields = line.split("\t")
        assert fields[0] == "TEST1"
        assert fields[1] == "-2.5"
        assert fields[2] == "0.003"
        assert fields[3] == "0.015"
        assert fields[5] == "200.5"

    def test_significant_true_exports_true_string(self) -> None:
        g = _make_gene(significant=True)
        r = _make_result(gene_results=(g,), n_genes_tested=1)
        tsv = render_results_tsv(r)
        line = tsv.split("\n")[1]
        fields = line.split("\t")
        assert fields[4] == "True"

    def test_significant_false_exports_false_string(self) -> None:
        g = _make_gene(significant=False)
        r = _make_result(gene_results=(g,), n_genes_tested=1)
        tsv = render_results_tsv(r)
        line = tsv.split("\n")[1]
        fields = line.split("\t")
        assert fields[4] == "False"

    def test_nan_p_value_exports_as_nan_string(self) -> None:
        g = _make_gene(p_value=float("nan"), adj_p_value=float("nan"))
        r = _make_result(gene_results=(g,), n_genes_tested=1)
        tsv = render_results_tsv(r)
        line = tsv.split("\n")[1]
        fields = line.split("\t")
        assert fields[2] == "NaN"
        assert fields[3] == "NaN"

    def test_nan_log2fc_exports_as_nan_string(self) -> None:
        g = _make_gene(log2_fold_change=float("nan"))
        r = _make_result(gene_results=(g,), n_genes_tested=1)
        tsv = render_results_tsv(r)
        line = tsv.split("\n")[1]
        fields = line.split("\t")
        assert fields[1] == "NaN"

    def test_nan_base_mean_exports_as_nan_string(self) -> None:
        g = _make_gene(base_mean=float("nan"))
        r = _make_result(gene_results=(g,), n_genes_tested=1)
        tsv = render_results_tsv(r)
        line = tsv.split("\n")[1]
        fields = line.split("\t")
        assert fields[5] == "NaN"

    def test_deterministic_output(self) -> None:
        genes = tuple(_make_gene(gene_id=f"G{i}", p_value=0.001 * i) for i in range(3))
        r = _make_result(gene_results=genes, n_genes_tested=3)
        assert render_results_tsv(r) == render_results_tsv(r)


# ---------------------------------------------------------------------------
# significant_genes.tsv
# ---------------------------------------------------------------------------


class TestSignificantGenesTsv:
    """Verify significant_genes.tsv correctness."""

    def test_header_columns(self) -> None:
        r = _make_result()
        tsv = render_significant_genes_tsv(r)
        header = tsv.split("\n")[0]
        expected = (
            "gene_id\tlog2_fold_change\tp_value\tadj_p_value\tsignificant\tbase_mean"
        )
        assert header == expected

    def test_only_significant_genes_included(self) -> None:
        genes = (
            _make_gene(gene_id="A", significant=True),
            _make_gene(gene_id="B", significant=False),
            _make_gene(gene_id="C", significant=True),
            _make_gene(gene_id="D", significant=False),
        )
        r = _make_result(gene_results=genes, n_genes_tested=4)
        tsv = render_significant_genes_tsv(r)
        lines = [ln for ln in tsv.split("\n") if ln]
        assert len(lines) == 3  # header + 2 significant
        assert lines[1].startswith("A")
        assert lines[2].startswith("C")

    def test_relative_order_preserved(self) -> None:
        genes = (
            _make_gene(gene_id="C", significant=False),
            _make_gene(gene_id="A", significant=True),
            _make_gene(gene_id="B", significant=True),
            _make_gene(gene_id="D", significant=False),
        )
        r = _make_result(gene_results=genes, n_genes_tested=4)
        tsv = render_significant_genes_tsv(r)
        lines = [ln for ln in tsv.split("\n") if ln][1:]  # skip header
        assert lines[0].startswith("A")
        assert lines[1].startswith("B")

    def test_no_significant_genes_header_only(self) -> None:
        genes = (
            _make_gene(gene_id="A", significant=False),
            _make_gene(gene_id="B", significant=False),
        )
        r = _make_result(gene_results=genes, n_genes_tested=2)
        tsv = render_significant_genes_tsv(r)
        lines = [ln for ln in tsv.split("\n") if ln]
        assert len(lines) == 1  # header only
        # Still a valid TSV with header + trailing newline
        assert tsv.endswith("\n")
        assert "gene_id" in tsv

    def test_empty_gene_results_header_only(self) -> None:
        r = _make_result(gene_results=(), n_genes_tested=0)
        tsv = render_significant_genes_tsv(r)
        lines = [ln for ln in tsv.split("\n") if ln]
        assert len(lines) == 1

    def test_significant_values_not_recomputed(self) -> None:
        """significant column should be the GeneDifferential.significant value."""
        g = _make_gene(gene_id="X", significant=True, p_value=0.001)
        r = _make_result(gene_results=(g,), n_genes_tested=1)
        tsv = render_significant_genes_tsv(r)
        line = tsv.split("\n")[1]
        fields = line.split("\t")
        assert fields[4] == "True"

    def test_deterministic_output(self) -> None:
        genes = (
            _make_gene(gene_id="G1", significant=True),
            _make_gene(gene_id="G2", significant=False),
            _make_gene(gene_id="G3", significant=True),
        )
        r = _make_result(gene_results=genes, n_genes_tested=3)
        assert render_significant_genes_tsv(r) == render_significant_genes_tsv(r)


# ---------------------------------------------------------------------------
# End-to-end with demo data
# ---------------------------------------------------------------------------


class TestEndToEndTsv:
    """Full pipeline: make_demo_inputs → run_analysis → render TSVs."""

    def test_end_to_end_results_tsv(self) -> None:
        from pharmomics.analysis.example_data import make_demo_inputs
        from pharmomics.analysis.run import run_analysis

        omics, design, spec = make_demo_inputs()
        result = run_analysis(spec, design, omics)
        tsv = render_results_tsv(result)

        lines = [ln for ln in tsv.split("\n") if ln]
        assert len(lines) == result.n_genes_tested + 1  # header + genes
        assert "gene_id" in lines[0]

    def test_end_to_end_significant_tsv(self) -> None:
        from pharmomics.analysis.example_data import make_demo_inputs
        from pharmomics.analysis.run import run_analysis

        omics, design, spec = make_demo_inputs()
        result = run_analysis(spec, design, omics)
        tsv = render_significant_genes_tsv(result)

        lines = [ln for ln in tsv.split("\n") if ln]
        sig_count = sum(1 for g in result.gene_results if g.significant)
        assert len(lines) == sig_count + 1  # header + significant

    def test_end_to_end_tsv_finite_values(self) -> None:
        """Demo data TSV should have valid p-value / adj_p_value columns."""
        import math

        from pharmomics.analysis.example_data import make_demo_inputs
        from pharmomics.analysis.run import run_analysis

        omics, design, spec = make_demo_inputs()
        result = run_analysis(spec, design, omics)
        tsv = render_results_tsv(result)

        for line in tsv.split("\n")[1:]:
            if not line:
                continue
            fields = line.split("\t")
            p_val = fields[2]
            adj_p = fields[3]
            # NaN is acceptable; otherwise should be a finite number
            if p_val != "NaN":
                assert math.isfinite(float(p_val)), f"p_value not finite: {p_val}"
            if adj_p != "NaN":
                assert math.isfinite(float(adj_p)), f"adj_p_value not finite: {adj_p}"

    def test_end_to_end_all_genes_in_tsv(self) -> None:
        """All demo genes should appear in results.tsv."""
        from pharmomics.analysis.example_data import make_demo_inputs
        from pharmomics.analysis.run import run_analysis

        omics, design, spec = make_demo_inputs()
        result = run_analysis(spec, design, omics)
        tsv = render_results_tsv(result)

        for gene in ("EGFR", "ERBB2", "TP53", "MYC", "KRAS", "PTEN"):
            assert gene in tsv, f"{gene} not found in results.tsv"
