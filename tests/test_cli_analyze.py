"""Tests for pharmomics CLI analyze command."""

from __future__ import annotations

from pathlib import Path

from typer.testing import CliRunner

from pharmomics.main import app

runner = CliRunner()

FIXTURES = Path(__file__).parent / "fixtures"

# ---------------------------------------------------------------------------
# Helper: create temporary input files
# ---------------------------------------------------------------------------


def _write_expr_file(tmp_path: Path, name: str = "expression.tsv") -> Path:
    """Write a small expression matrix with non-zero within-group variance."""
    path = tmp_path / name
    path.write_text(
        "gene\tctrl_1\tctrl_2\tctrl_3\ttrt_1\ttrt_2\ttrt_3\n"
        "EGFR\t1000\t1050\t980\t200\t210\t190\n"
        "ERBB2\t500\t520\t480\t800\t790\t810\n"
        "TP53\t300\t310\t290\t305\t295\t300\n"
        "BRCA1\t150\t160\t140\t50\t55\t45\n"
        "MYC\t800\t820\t780\t1200\t1180\t1220\n",
        encoding="utf-8",
    )
    return path


def _write_metadata_json(tmp_path: Path, name: str = "metadata.json") -> Path:
    """Write sample metadata matching the expression file."""
    path = tmp_path / name
    path.write_text(
        '{"samples": {'
        '"ctrl_1": {"cell_line": "PC9", "condition": "DMSO", "replicate": 1},'
        '"ctrl_2": {"cell_line": "PC9", "condition": "DMSO", "replicate": 2},'
        '"ctrl_3": {"cell_line": "PC9", "condition": "DMSO", "replicate": 3},'
        '"trt_1": {"cell_line": "PC9", "condition": "osi_DTP", "replicate": 1},'
        '"trt_2": {"cell_line": "PC9", "condition": "osi_DTP", "replicate": 2},'
        '"trt_3": {"cell_line": "PC9", "condition": "osi_DTP", "replicate": 3}'
        "}}",
        encoding="utf-8",
    )
    return path


# ---------------------------------------------------------------------------
# Success cases
# ---------------------------------------------------------------------------


class TestCLIAnalyzeSuccess:
    """Verify CLI analyze command success cases."""

    def test_end_to_end_with_fixture_files(self, tmp_path: Path) -> None:
        """Real expression + metadata files should produce a report."""
        expr = _write_expr_file(tmp_path)
        meta = _write_metadata_json(tmp_path)
        output = tmp_path / "report.md"
        result = runner.invoke(
            app,
            [
                "analyze",
                "--expression-file",
                str(expr),
                "--metadata-file",
                str(meta),
                "--contrast-control",
                "DMSO",
                "--contrast-treatment",
                "osi_DTP",
                "--output",
                str(output),
            ],
        )
        assert result.exit_code == 0
        assert output.exists()
        content = output.read_text(encoding="utf-8")
        assert "# Differential Analysis Report" in content

    def test_default_output_path(self, tmp_path: Path) -> None:
        """Default output path (report.md) should be created in cwd."""
        expr = _write_expr_file(tmp_path)
        meta = _write_metadata_json(tmp_path)
        output = tmp_path / "report.md"
        result = runner.invoke(
            app,
            [
                "analyze",
                "--expression-file",
                str(expr),
                "--metadata-file",
                str(meta),
                "--contrast-control",
                "DMSO",
                "--contrast-treatment",
                "osi_DTP",
                "--output",
                str(output),
            ],
        )
        assert result.exit_code == 0
        assert output.exists()

    def test_custom_output_path(self, tmp_path: Path) -> None:
        """Custom --output path should be respected."""
        expr = _write_expr_file(tmp_path)
        meta = _write_metadata_json(tmp_path)
        output = tmp_path / "sub" / "my-report.md"
        output.parent.mkdir()
        result = runner.invoke(
            app,
            [
                "analyze",
                "--expression-file",
                str(expr),
                "--metadata-file",
                str(meta),
                "--contrast-control",
                "DMSO",
                "--contrast-treatment",
                "osi_DTP",
                "--output",
                str(output),
            ],
        )
        assert result.exit_code == 0
        assert output.exists()

    def test_output_contains_analysis_type(self, tmp_path: Path) -> None:
        """Report should contain the analysis_type."""
        expr = _write_expr_file(tmp_path)
        meta = _write_metadata_json(tmp_path)
        output = tmp_path / "report.md"
        result = runner.invoke(
            app,
            [
                "analyze",
                "--expression-file",
                str(expr),
                "--metadata-file",
                str(meta),
                "--contrast-control",
                "DMSO",
                "--contrast-treatment",
                "osi_DTP",
                "--output",
                str(output),
            ],
        )
        assert result.exit_code == 0
        content = output.read_text(encoding="utf-8")
        assert "differential_analysis" in content

    def test_output_contains_contrast_id(self, tmp_path: Path) -> None:
        """Report should contain the contrast_id."""
        expr = _write_expr_file(tmp_path)
        meta = _write_metadata_json(tmp_path)
        output = tmp_path / "report.md"
        result = runner.invoke(
            app,
            [
                "analyze",
                "--expression-file",
                str(expr),
                "--metadata-file",
                str(meta),
                "--contrast-control",
                "DMSO",
                "--contrast-treatment",
                "osi_DTP",
                "--output",
                str(output),
            ],
        )
        assert result.exit_code == 0
        content = output.read_text(encoding="utf-8")
        assert "osi_dtp_vs_dmso" in content

    def test_output_contains_gene_results(self, tmp_path: Path) -> None:
        """Report should contain the gene results table."""
        expr = _write_expr_file(tmp_path)
        meta = _write_metadata_json(tmp_path)
        output = tmp_path / "report.md"
        result = runner.invoke(
            app,
            [
                "analyze",
                "--expression-file",
                str(expr),
                "--metadata-file",
                str(meta),
                "--contrast-control",
                "DMSO",
                "--contrast-treatment",
                "osi_DTP",
                "--output",
                str(output),
            ],
        )
        assert result.exit_code == 0
        content = output.read_text(encoding="utf-8")
        assert "| Gene |" in content
        assert "| log2FC |" in content
        assert "EGFR" in content

    def test_report_contains_statistical_values(self, tmp_path: Path) -> None:
        """Report should contain log2FC, p-value, and significance values."""
        expr = _write_expr_file(tmp_path)
        meta = _write_metadata_json(tmp_path)
        output = tmp_path / "report.md"
        result = runner.invoke(
            app,
            [
                "analyze",
                "--expression-file",
                str(expr),
                "--metadata-file",
                str(meta),
                "--contrast-control",
                "DMSO",
                "--contrast-treatment",
                "osi_DTP",
                "--output",
                str(output),
            ],
        )
        assert result.exit_code == 0
        content = output.read_text(encoding="utf-8")
        assert "| p-value |" in content
        assert "| adj p-value |" in content
        assert "| Significant |" in content

    def test_report_written_confirmation(self, tmp_path: Path) -> None:
        """CLI should print confirmation message."""
        expr = _write_expr_file(tmp_path)
        meta = _write_metadata_json(tmp_path)
        output = tmp_path / "report.md"
        result = runner.invoke(
            app,
            [
                "analyze",
                "--expression-file",
                str(expr),
                "--metadata-file",
                str(meta),
                "--contrast-control",
                "DMSO",
                "--contrast-treatment",
                "osi_DTP",
                "--output",
                str(output),
            ],
        )
        assert result.exit_code == 0
        assert "Report written to" in result.output


# ---------------------------------------------------------------------------
# Failure cases
# ---------------------------------------------------------------------------


class TestCLIAnalyzeFailure:
    """Verify CLI analyze command failure cases."""

    def test_missing_expression_file(self, tmp_path: Path) -> None:
        """Should fail when expression file does not exist."""
        meta = _write_metadata_json(tmp_path)
        result = runner.invoke(
            app,
            [
                "analyze",
                "--expression-file",
                str(tmp_path / "nonexistent.tsv"),
                "--metadata-file",
                str(meta),
                "--contrast-control",
                "DMSO",
                "--contrast-treatment",
                "osi_DTP",
            ],
        )
        assert result.exit_code != 0
        assert "not found" in result.output

    def test_missing_metadata_file(self, tmp_path: Path) -> None:
        """Should fail when metadata file does not exist."""
        expr = _write_expr_file(tmp_path)
        result = runner.invoke(
            app,
            [
                "analyze",
                "--expression-file",
                str(expr),
                "--metadata-file",
                str(tmp_path / "nonexistent.json"),
                "--contrast-control",
                "DMSO",
                "--contrast-treatment",
                "osi_DTP",
            ],
        )
        assert result.exit_code != 0
        assert "not found" in result.output

    def test_output_parent_directory_missing(self, tmp_path: Path) -> None:
        """Should fail when output parent directory does not exist."""
        expr = _write_expr_file(tmp_path)
        meta = _write_metadata_json(tmp_path)
        output = tmp_path / "nonexistent" / "dir" / "report.md"
        result = runner.invoke(
            app,
            [
                "analyze",
                "--expression-file",
                str(expr),
                "--metadata-file",
                str(meta),
                "--contrast-control",
                "DMSO",
                "--contrast-treatment",
                "osi_DTP",
                "--output",
                str(output),
            ],
        )
        assert result.exit_code != 0
        assert "does not exist" in result.output

    def test_control_group_not_found(self, tmp_path: Path) -> None:
        """Should fail when control condition has no samples."""
        expr = _write_expr_file(tmp_path)
        # Metadata with no "DMSO" condition
        meta_path = tmp_path / "metadata.json"
        meta_path.write_text(
            '{"samples": {'
            '"ctrl_1": {"condition": "ctrl"},'
            '"ctrl_2": {"condition": "ctrl"},'
            '"ctrl_3": {"condition": "ctrl"},'
            '"trt_1": {"condition": "osi_DTP"},'
            '"trt_2": {"condition": "osi_DTP"},'
            '"trt_3": {"condition": "osi_DTP"}'
            "}}",
            encoding="utf-8",
        )
        result = runner.invoke(
            app,
            [
                "analyze",
                "--expression-file",
                str(expr),
                "--metadata-file",
                str(meta_path),
                "--contrast-control",
                "DMSO",
                "--contrast-treatment",
                "osi_DTP",
            ],
        )
        assert result.exit_code != 0

    def test_treatment_group_not_found(self, tmp_path: Path) -> None:
        """Should fail when treatment condition has no samples."""
        expr = _write_expr_file(tmp_path)
        meta_path = tmp_path / "metadata.json"
        meta_path.write_text(
            '{"samples": {'
            '"ctrl_1": {"condition": "DMSO"},'
            '"ctrl_2": {"condition": "DMSO"},'
            '"ctrl_3": {"condition": "DMSO"},'
            '"trt_1": {"condition": "ctrl"},'
            '"trt_2": {"condition": "ctrl"},'
            '"trt_3": {"condition": "ctrl"}'
            "}}",
            encoding="utf-8",
        )
        result = runner.invoke(
            app,
            [
                "analyze",
                "--expression-file",
                str(expr),
                "--metadata-file",
                str(meta_path),
                "--contrast-control",
                "DMSO",
                "--contrast-treatment",
                "osi_DTP",
            ],
        )
        assert result.exit_code != 0

    def test_sample_id_mismatch(self, tmp_path: Path) -> None:
        """Should fail when metadata sample IDs do not match expression."""
        expr = _write_expr_file(tmp_path)
        meta_path = tmp_path / "metadata.json"
        meta_path.write_text(
            '{"samples": {'
            '"S1": {"condition": "DMSO"},'
            '"S2": {"condition": "DMSO"},'
            '"S3": {"condition": "DMSO"},'
            '"S4": {"condition": "osi_DTP"},'
            '"S5": {"condition": "osi_DTP"},'
            '"S6": {"condition": "osi_DTP"}'
            "}}",
            encoding="utf-8",
        )
        result = runner.invoke(
            app,
            [
                "analyze",
                "--expression-file",
                str(expr),
                "--metadata-file",
                str(meta_path),
                "--contrast-control",
                "DMSO",
                "--contrast-treatment",
                "osi_DTP",
            ],
        )
        assert result.exit_code != 0

    def test_malformed_expression_file(self, tmp_path: Path) -> None:
        """Should fail on malformed expression matrix."""
        expr_path = tmp_path / "bad.tsv"
        expr_path.write_text("just_one_column\n", encoding="utf-8")
        meta = _write_metadata_json(tmp_path)
        result = runner.invoke(
            app,
            [
                "analyze",
                "--expression-file",
                str(expr_path),
                "--metadata-file",
                str(meta),
                "--contrast-control",
                "DMSO",
                "--contrast-treatment",
                "osi_DTP",
            ],
        )
        assert result.exit_code != 0


# ---------------------------------------------------------------------------
# Help
# ---------------------------------------------------------------------------


class TestCLIAnalyzeHelp:
    """Verify CLI analyze --help works."""

    def test_help(self) -> None:
        """--help should show analyze command options."""
        result = runner.invoke(app, ["analyze", "--help"])
        assert result.exit_code == 0
        assert "--expression-file" in result.output
        assert "--metadata-file" in result.output
        assert "--contrast-control" in result.output
        assert "--contrast-treatment" in result.output
        assert "--output" in result.output
