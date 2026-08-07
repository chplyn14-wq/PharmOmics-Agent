"""Tests for pharmomics.analysis.render.

Covers:
- render_markdown_report produces well-formed Markdown from AnalysisResult.
- Metadata section contains analysis_type, contrast_id, n_genes_tested.
- Gene results table renders all fields with correct formatting.
- NaN / +Inf / -Inf values render as explicit strings.
- Significant bool renders as Yes/No without recalculation.
- Gene results order is preserved (no re-sorting).
- Warnings section shows None when empty, list items when present.
- Empty gene_results produces a table with header but no data rows.
- Numeric formatting is deterministic (.3f for log2FC, .2f for base_mean,
  .4g for p-values).
- End-to-end: make_demo_inputs() → run_analysis() → render_markdown_report().
"""

from __future__ import annotations

from pharmomics.analysis.render import render_markdown_report
from pharmomics.analysis.results import AnalysisResult, GeneDifferential

# ---------------------------------------------------------------------------
# Test fixture helpers
# ---------------------------------------------------------------------------


def _make_result(
    *,
    analysis_type: str = "differential_analysis",
    contrast_id: str = "trt_vs_ctrl",
    gene_results: tuple[GeneDifferential, ...] = (),
    n_genes_tested: int = 0,
    warnings: tuple[str, ...] = (),
) -> AnalysisResult:
    """Create a minimal AnalysisResult for testing."""
    return AnalysisResult(
        analysis_type=analysis_type,
        contrast_id=contrast_id,
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
# Metadata
# ---------------------------------------------------------------------------


class TestMetadata:
    """Verify the summary section renders correct metadata."""

    def test_output_contains_analysis_type(self) -> None:
        r = _make_result(analysis_type="differential_analysis")
        report = render_markdown_report(r)
        assert "| Analysis type | differential_analysis |" in report

    def test_output_contains_contrast_id(self) -> None:
        r = _make_result(contrast_id="treated_vs_control")
        report = render_markdown_report(r)
        assert "| Contrast | treated_vs_control |" in report

    def test_output_contains_n_genes_tested(self) -> None:
        genes = tuple(_make_gene(gene_id=f"G{i}") for i in range(42))
        r = _make_result(n_genes_tested=42, gene_results=genes)
        report = render_markdown_report(r)
        assert "| Genes tested | 42 |" in report

    def test_report_starts_with_title(self) -> None:
        r = _make_result()
        report = render_markdown_report(r)
        assert report.startswith("# Differential Analysis Report")

    def test_result_is_string(self) -> None:
        r = _make_result()
        report = render_markdown_report(r)
        assert isinstance(report, str)


# ---------------------------------------------------------------------------
# Empty gene_results
# ---------------------------------------------------------------------------


class TestEmptyGeneResults:
    """Verify behaviour when no gene results are present."""

    def test_empty_result_produces_valid_report(self) -> None:
        r = _make_result(n_genes_tested=0, gene_results=())
        report = render_markdown_report(r)

        # Summary section present
        assert "| Genes tested | 0 |" in report
        # Table header exists
        assert "| Gene | log2FC |" in report
        assert "|---|---|---|---|---|---|" in report
        # No data rows (the header line is the last gene-related line)
        lines = report.split("\n")
        gene_section_start = next(
            i for i, line in enumerate(lines) if "| Gene | log2FC |" in line
        )
        # Only the header line and separator line follow; no data rows
        table_lines = [
            line
            for line in lines[gene_section_start:]
            if line.startswith("|")
        ]
        assert len(table_lines) == 2  # header + separator


# ---------------------------------------------------------------------------
# Gene table rendering
# ---------------------------------------------------------------------------


class TestGeneTable:
    """Verify gene results table content and formatting."""

    def test_single_gene_report(self) -> None:
        g = _make_gene(gene_id="EGFR", log2_fold_change=1.0, base_mean=150.0)
        r = _make_result(gene_results=(g,), n_genes_tested=1)
        report = render_markdown_report(r)

        assert "EGFR" in report

    def test_gene_id_rendered(self) -> None:
        g = _make_gene(gene_id="TP53")
        r = _make_result(gene_results=(g,), n_genes_tested=1)
        report = render_markdown_report(r)
        assert "| TP53" in report


# ---------------------------------------------------------------------------
# Order preservation
# ---------------------------------------------------------------------------


class TestOrderPreservation:
    """Verify gene_results appear in original input order."""

    def test_multiple_genes_order_preserved(self) -> None:
        genes = (
            _make_gene(
                gene_id="C", log2_fold_change=0.0, base_mean=0.0, significant=False
            ),
            _make_gene(
                gene_id="A", log2_fold_change=0.0, base_mean=0.0, significant=False
            ),
            _make_gene(
                gene_id="B", log2_fold_change=0.0, base_mean=0.0, significant=False
            ),
        )
        r = _make_result(gene_results=genes, n_genes_tested=3)
        report = render_markdown_report(r)

        c_pos = report.index("| C |")
        a_pos = report.index("| A |")
        b_pos = report.index("| B |")
        assert c_pos < a_pos < b_pos


# ---------------------------------------------------------------------------
# Significant rendering
# ---------------------------------------------------------------------------


class TestSignificantRendering:
    """Verify significant bool renders directly without recalculation."""

    def test_significant_true_shows_yes(self) -> None:
        g = _make_gene(significant=True)
        r = _make_result(gene_results=(g,), n_genes_tested=1)
        report = render_markdown_report(r)
        assert "| Yes" in report

    def test_significant_false_shows_no(self) -> None:
        g = _make_gene(significant=False)
        r = _make_result(gene_results=(g,), n_genes_tested=1)
        report = render_markdown_report(r)
        assert "| No" in report


# ---------------------------------------------------------------------------
# Warnings
# ---------------------------------------------------------------------------


class TestWarnings:
    """Verify warnings section rendering."""

    def test_warnings_none_when_empty(self) -> None:
        r = _make_result(warnings=())
        report = render_markdown_report(r)
        assert "## Warnings" in report
        # "None" appears in the warnings section
        warnings_section = report.split("## Warnings")[1].split("## Gene")[0]
        assert "None" in warnings_section
        assert "- " not in warnings_section

    def test_warnings_displayed_when_present(self) -> None:
        r = _make_result(warnings=("Low replicate count", "Gene X has zero expression"))
        report = render_markdown_report(r)

        assert "- Low replicate count" in report
        assert "- Gene X has zero expression" in report


# ---------------------------------------------------------------------------
# NaN / Inf formatting
# ---------------------------------------------------------------------------


class TestNonFiniteFormatting:
    """Verify NaN, +Inf, -Inf render as explicit strings."""

    def test_nan_values_rendered(self) -> None:
        g = _make_gene(
            gene_id="NaN_GENE",
            p_value=float("nan"),
            adj_p_value=float("nan"),
            significant=False,
        )
        r = _make_result(gene_results=(g,), n_genes_tested=1)
        report = render_markdown_report(r)

        # NaN should appear in the gene row
        gene_row = [line for line in report.split("\n") if "NaN_GENE" in line][0]
        assert "NaN" in gene_row

    def test_inf_values_rendered(self) -> None:
        g = _make_gene(
            gene_id="INF_GENE",
            log2_fold_change=float("inf"),
            significant=False,
        )
        r = _make_result(gene_results=(g,), n_genes_tested=1)
        report = render_markdown_report(r)

        gene_row = [line for line in report.split("\n") if "INF_GENE" in line][0]
        assert "+Inf" in gene_row

    def test_negative_inf_rendered(self) -> None:
        g = _make_gene(
            gene_id="NINF_GENE",
            log2_fold_change=float("-inf"),
            significant=False,
        )
        r = _make_result(gene_results=(g,), n_genes_tested=1)
        report = render_markdown_report(r)

        gene_row = [line for line in report.split("\n") if "NINF_GENE" in line][0]
        assert "-Inf" in gene_row


# ---------------------------------------------------------------------------
# Numeric formatting
# ---------------------------------------------------------------------------


class TestNumericFormatting:
    """Verify deterministic formatting for finite values."""

    def test_log2fc_formatting(self) -> None:
        """log2FC uses .3f format."""
        g = _make_gene(gene_id="FMT", log2_fold_change=1.23456, base_mean=100.0)
        r = _make_result(gene_results=(g,), n_genes_tested=1)
        report = render_markdown_report(r)

        gene_row = [line for line in report.split("\n") if "| FMT" in line][0]
        assert "1.235" in gene_row  # .3f rounds

    def test_p_value_formatting_finite(self) -> None:
        """p-value uses .4g format."""
        g = _make_gene(gene_id="PVAL", p_value=0.001234, adj_p_value=0.01234)
        r = _make_result(gene_results=(g,), n_genes_tested=1)
        report = render_markdown_report(r)

        gene_row = [line for line in report.split("\n") if "| PVAL" in line][0]
        # .4g: 0.001234 → "0.001234", 0.01234 → "0.01234"
        assert "0.001234" in gene_row
        assert "0.01234" in gene_row

    def test_base_mean_formatting(self) -> None:
        """base_mean uses .2f format."""
        g = _make_gene(gene_id="BMEAN", base_mean=123.456)
        r = _make_result(gene_results=(g,), n_genes_tested=1)
        report = render_markdown_report(r)

        gene_row = [line for line in report.split("\n") if "| BMEAN" in line][0]
        assert "123.46" in gene_row  # .2f rounds

    def test_small_p_value_scientific_notation(self) -> None:
        """Very small p-values use scientific notation via .4g."""
        g = _make_gene(
            gene_id="TINY",
            p_value=1.234e-10,
            adj_p_value=6.17e-10,
            significant=True,
        )
        r = _make_result(gene_results=(g,), n_genes_tested=1)
        report = render_markdown_report(r)

        gene_row = [line for line in report.split("\n") if "| TINY" in line][0]
        # .4g should produce something like "1.234e-10"
        assert "1.234e-10" in gene_row

    def test_zero_log2fc(self) -> None:
        """Zero log2FC renders as 0.000."""
        g = _make_gene(gene_id="ZERO", log2_fold_change=0.0, base_mean=50.0)
        r = _make_result(gene_results=(g,), n_genes_tested=1)
        report = render_markdown_report(r)

        gene_row = [line for line in report.split("\n") if "| ZERO" in line][0]
        assert "0.000" in gene_row


# ---------------------------------------------------------------------------
# Deterministic output
# ---------------------------------------------------------------------------


class TestDeterministicOutput:
    """Verify the same input always produces the same output."""

    def test_deterministic_output(self) -> None:
        g = _make_gene()
        r = _make_result(gene_results=(g,), n_genes_tested=1)

        report1 = render_markdown_report(r)
        report2 = render_markdown_report(r)

        assert report1 == report2


# ---------------------------------------------------------------------------
# End-to-end with real pipeline
# ---------------------------------------------------------------------------


class TestEndToEndPipeline:
    """Full pipeline: make_demo_inputs → run_analysis → render."""

    def test_end_to_end_render(self) -> None:
        """Render a report from a real analysis run."""
        from pharmomics.analysis.example_data import make_demo_inputs
        from pharmomics.analysis.run import run_analysis

        omics, design, spec = make_demo_inputs()
        result = run_analysis(spec, design, omics)
        report = render_markdown_report(result)

        assert isinstance(report, str)
        assert result.contrast_id in report
        assert result.analysis_type in report
        assert str(result.n_genes_tested) in report

    def test_end_to_end_contains_all_genes(self) -> None:
        """All 6 demo genes appear in the report."""
        from pharmomics.analysis.example_data import make_demo_inputs
        from pharmomics.analysis.run import run_analysis

        omics, design, spec = make_demo_inputs()
        result = run_analysis(spec, design, omics)
        report = render_markdown_report(result)

        for gene_id in ("EGFR", "ERBB2", "TP53", "MYC", "KRAS", "PTEN"):
            assert gene_id in report

    def test_end_to_end_nan_p_values_rendered(self) -> None:
        """Demo data has NaN p-values → rendered as NaN."""
        from pharmomics.analysis.example_data import make_demo_inputs
        from pharmomics.analysis.run import run_analysis

        omics, design, spec = make_demo_inputs()
        result = run_analysis(spec, design, omics)
        report = render_markdown_report(result)

        # All genes have NaN p-values in demo data (zero variance)
        for gene in result.gene_results:
            import math

            if math.isnan(gene.p_value):
                gene_row = [
                    line for line in report.split("\n") if gene.gene_id in line
                ][0]
                assert "NaN" in gene_row

    def test_end_to_end_no_warnings(self) -> None:
        """Successful demo run has no warnings."""
        from pharmomics.analysis.example_data import make_demo_inputs
        from pharmomics.analysis.run import run_analysis

        omics, design, spec = make_demo_inputs()
        result = run_analysis(spec, design, omics)
        report = render_markdown_report(result)

        assert result.warnings == ()
        # Warnings section should show None
        warnings_section = report.split("## Warnings")[1].split("## Gene")[0]
        assert "None" in warnings_section

    def test_end_to_end_gene_order_preserved(self) -> None:
        """Gene order in report matches result.gene_results order."""
        from pharmomics.analysis.example_data import make_demo_inputs
        from pharmomics.analysis.run import run_analysis

        omics, design, spec = make_demo_inputs()
        result = run_analysis(spec, design, omics)
        report = render_markdown_report(result)

        expected_order = ["EGFR", "ERBB2", "TP53", "MYC", "KRAS", "PTEN"]
        positions = [report.index(gid) for gid in expected_order]
        assert positions == sorted(positions)
