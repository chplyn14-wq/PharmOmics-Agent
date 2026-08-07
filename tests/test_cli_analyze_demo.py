"""Tests for pharmomics CLI analyze-demo command."""

from __future__ import annotations

from pathlib import Path

from typer.testing import CliRunner

from pharmomics.main import app

runner = CliRunner()


class TestCLIAnalyzeDemoSuccess:
    """Verify CLI analyze-demo command success cases."""

    def test_default_output_path(self, tmp_path: Path) -> None:
        """Default output path (report.md) should be created in cwd."""
        cwd = tmp_path
        output = cwd / "report.md"
        result = runner.invoke(
            app,
            ["analyze-demo", "--output", str(output)],
        )
        assert result.exit_code == 0
        assert "Report written to" in result.output
        assert output.exists()
        content = output.read_text(encoding="utf-8")
        assert "# Differential Analysis Report" in content

    def test_custom_output_path(self, tmp_path: Path) -> None:
        """Custom --output path should be respected."""
        output = tmp_path / "sub" / "my-report.md"
        output.parent.mkdir()
        result = runner.invoke(
            app,
            ["analyze-demo", "--output", str(output)],
        )
        assert result.exit_code == 0
        assert output.exists()
        content = output.read_text(encoding="utf-8")
        assert "# Differential Analysis Report" in content

    def test_overwrites_existing_file(self, tmp_path: Path) -> None:
        """Should overwrite an existing report.md without error."""
        output = tmp_path / "report.md"
        output.write_text("old content", encoding="utf-8")
        result = runner.invoke(
            app,
            ["analyze-demo", "--output", str(output)],
        )
        assert result.exit_code == 0
        new_content = output.read_text(encoding="utf-8")
        assert "# Differential Analysis Report" in new_content
        assert "old content" not in new_content

    def test_report_contains_gene_table(self, tmp_path: Path) -> None:
        """Report should contain the gene results table header."""
        output = tmp_path / "report.md"
        result = runner.invoke(
            app,
            ["analyze-demo", "--output", str(output)],
        )
        assert result.exit_code == 0
        content = output.read_text(encoding="utf-8")
        assert "| Gene |" in content
        assert "| log2FC |" in content
        assert "EGFR" in content


class TestCLIAnalyzeDemoFailure:
    """Verify CLI analyze-demo command failure cases."""

    def test_parent_directory_missing(self, tmp_path: Path) -> None:
        """Should fail with non-zero exit when parent dir is missing."""
        output = tmp_path / "nonexistent" / "dir" / "report.md"
        result = runner.invoke(
            app,
            ["analyze-demo", "--output", str(output)],
        )
        assert result.exit_code != 0
        assert "does not exist" in result.output
