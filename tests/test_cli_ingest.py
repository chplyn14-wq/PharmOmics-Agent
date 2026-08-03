"""Tests for pharmomics CLI ingest command."""

from __future__ import annotations

import gzip
import shutil
from pathlib import Path

from typer.testing import CliRunner

from pharmomics.main import app

FIXTURES = Path(__file__).parent / "fixtures"

runner = CliRunner()


def _make_gzip(src: Path, dst: Path) -> Path:
    """Create a gzip-compressed copy of a file."""
    with open(src, "rb") as f_in:
        with gzip.open(dst, "wb") as f_out:
            shutil.copyfileobj(f_in, f_out)
    return dst


class TestCLIIngestSuccess:
    """Verify CLI ingest command success cases."""

    def test_ingest_tsv(self, tmp_path: Path) -> None:
        run_dir = tmp_path / "run"
        result = runner.invoke(
            app,
            [
                "ingest",
                "--expression-file",
                str(FIXTURES / "synthetic_expression.tsv"),
                "--metadata-file",
                str(FIXTURES / "synthetic_metadata.json"),
                "--source-id",
                "GSE_SYNTHETIC",
                "--run-dir",
                str(run_dir),
            ],
        )
        assert result.exit_code == 0
        assert "Ingestion complete" in result.output
        assert "5" in result.output  # n_genes
        assert "6" in result.output  # n_samples

    def test_ingest_gzip(self, tmp_path: Path) -> None:
        gz_path = tmp_path / "expression.tsv.gz"
        _make_gzip(FIXTURES / "synthetic_expression.tsv", gz_path)
        meta_path = tmp_path / "metadata.json"
        meta_path.write_text(
            '{"samples": {'
            '"PC9_DMSO_1": {"condition": "DMSO"},'
            '"PC9_DMSO_2": {"condition": "DMSO"},'
            '"PC9_DMSO_3": {"condition": "DMSO"},'
            '"PC9_osi_DTP_1": {"condition": "osi_DTP"},'
            '"PC9_osi_DTP_2": {"condition": "osi_DTP"},'
            '"PC9_osi_DTP_3": {"condition": "osi_DTP"}}}',
            encoding="utf-8",
        )
        run_dir = tmp_path / "run"
        result = runner.invoke(
            app,
            [
                "ingest",
                "--expression-file",
                str(gz_path),
                "--metadata-file",
                str(meta_path),
                "--source-id",
                "GSE_GZ_TEST",
                "--run-dir",
                str(run_dir),
            ],
        )
        assert result.exit_code == 0
        assert "gzip" in result.output

    def test_ingest_csv(self, tmp_path: Path) -> None:
        run_dir = tmp_path / "run"
        meta_path = tmp_path / "metadata.json"
        meta_path.write_text(
            '{"samples": {'
            '"PC9_DMSO_1": {"condition": "DMSO"},'
            '"PC9_DMSO_2": {"condition": "DMSO"},'
            '"PC9_DMSO_3": {"condition": "DMSO"},'
            '"PC9_osi_DTP_1": {"condition": "osi_DTP"},'
            '"PC9_osi_DTP_2": {"condition": "osi_DTP"},'
            '"PC9_osi_DTP_3": {"condition": "osi_DTP"}}}',
            encoding="utf-8",
        )
        result = runner.invoke(
            app,
            [
                "ingest",
                "--expression-file",
                str(FIXTURES / "synthetic_expression.csv"),
                "--metadata-file",
                str(meta_path),
                "--source-id",
                "GSE_CSV_TEST",
                "--run-dir",
                str(run_dir),
            ],
        )
        assert result.exit_code == 0


class TestCLIIngestFailure:
    """Verify CLI ingest command failure cases."""

    def test_missing_expression_file(self, tmp_path: Path) -> None:
        run_dir = tmp_path / "run"
        result = runner.invoke(
            app,
            [
                "ingest",
                "--expression-file",
                "/nonexistent/file.tsv",
                "--metadata-file",
                str(FIXTURES / "synthetic_metadata.json"),
                "--source-id",
                "GSE_TEST",
                "--run-dir",
                str(run_dir),
            ],
        )
        assert result.exit_code != 0
        assert "not found" in result.output

    def test_missing_metadata_file(self, tmp_path: Path) -> None:
        run_dir = tmp_path / "run"
        result = runner.invoke(
            app,
            [
                "ingest",
                "--expression-file",
                str(FIXTURES / "synthetic_expression.tsv"),
                "--metadata-file",
                "/nonexistent/metadata.json",
                "--source-id",
                "GSE_TEST",
                "--run-dir",
                str(run_dir),
            ],
        )
        assert result.exit_code != 0
        assert "not found" in result.output

    def test_validation_failure(self, tmp_path: Path) -> None:
        """Ingestion should fail with mismatched expression/metadata."""
        run_dir = tmp_path / "run"
        bad_meta = tmp_path / "bad_metadata.json"
        bad_meta.write_text(
            '{"samples": {'
            '"WRONG_ID_1": {"condition": "DMSO"},'
            '"WRONG_ID_2": {"condition": "treated"}}}',
            encoding="utf-8",
        )
        result = runner.invoke(
            app,
            [
                "ingest",
                "--expression-file",
                str(FIXTURES / "synthetic_expression.tsv"),
                "--metadata-file",
                str(bad_meta),
                "--source-id",
                "GSE_TEST",
                "--run-dir",
                str(run_dir),
            ],
        )
        assert result.exit_code != 0
        assert "Ingestion failed" in result.output
