"""Tests for pharmomics.ingestion — expression matrix and metadata loading.

Milestone 1B test coverage includes:
- Valid TSV ingestion
- Valid gzip TSV ingestion
- CSV delimiter detection
- Duplicate sample columns
- Duplicate metadata sample IDs
- Expression/metadata sample mismatch
- Extra metadata samples
- Missing condition
- Empty matrix
- Malformed numeric values
- Integer-count classification
- Non-integer estimated-count classification
- Transformed-value classification
- Unknown classification
- Ensembl version stripping
- HGNC-symbol detection
- Duplicate gene IDs
- Relative persisted paths
- Independent file hashes
- JSON serialization round-trip
- Windows-compatible paths
- CLI success and CLI validation failure
"""

from __future__ import annotations

import gzip
import json
import shutil
from pathlib import Path

import pytest

from pharmomics.ingestion.loader import (
    CompressionType,
    ContrastError,
    ExpressionFileError,
    GeneIdType,
    MetadataFileError,
    ValueType,
    classify_expression_values,
    count_replicates_per_condition,
    ingest,
    inspect_gene_identifiers,
    load_expression_matrix,
    load_sample_metadata,
    validate_contrast,
    write_ingestion_manifest,
)
from pharmomics.run_store import hash_file_sha256

FIXTURES = Path(__file__).parent / "fixtures"

# ---------------------------------------------------------------------------
# Helper: create gzip fixture files
# ---------------------------------------------------------------------------


def _make_gzip(src: Path, dst: Path) -> Path:
    """Create a gzip-compressed copy of a file."""
    with open(src, "rb") as f_in:
        with gzip.open(dst, "wb") as f_out:
            shutil.copyfileobj(f_in, f_out)
    return dst


# ---------------------------------------------------------------------------
# Expression matrix loading
# ---------------------------------------------------------------------------


class TestLoadExpressionMatrixTSV:
    """Verify TSV expression matrix loading."""

    def test_valid_tsv(self, tmp_path: Path) -> None:
        src = FIXTURES / "synthetic_expression.tsv"
        result = load_expression_matrix(src)
        assert result.n_genes == 5
        assert result.n_samples == 6
        assert result.sample_ids == [
            "PC9_DMSO_1",
            "PC9_DMSO_2",
            "PC9_DMSO_3",
            "PC9_osi_DTP_1",
            "PC9_osi_DTP_2",
            "PC9_osi_DTP_3",
        ]
        assert result.gene_ids == ["EGFR", "ERBB2", "TP53", "BRCA1", "MYC"]
        assert result.delimiter == "\t"
        assert result.compression == CompressionType.NONE

    def test_valid_gzip_tsv(self, tmp_path: Path) -> None:
        gz_path = tmp_path / "expression.tsv.gz"
        _make_gzip(FIXTURES / "synthetic_expression.tsv", gz_path)
        result = load_expression_matrix(gz_path)
        assert result.n_genes == 5
        assert result.n_samples == 6
        assert result.compression == CompressionType.GZIP


class TestLoadExpressionMatrixCSV:
    """Verify CSV expression matrix loading."""

    def test_csv_delimiter_detection(self) -> None:
        src = FIXTURES / "synthetic_expression.csv"
        result = load_expression_matrix(src)
        assert result.delimiter == ","
        assert result.n_genes == 5
        assert result.n_samples == 6


class TestLoadExpressionMatrixErrors:
    """Verify expression matrix error handling."""

    def test_file_not_found(self) -> None:
        with pytest.raises(FileNotFoundError):
            load_expression_matrix(Path("/nonexistent/path.tsv"))

    def test_empty_matrix(self, tmp_path: Path) -> None:
        f = tmp_path / "empty.tsv"
        f.write_bytes(b"")
        with pytest.raises(ExpressionFileError, match="empty"):
            load_expression_matrix(f)

    def test_duplicate_sample_columns(self, tmp_path: Path) -> None:
        f = tmp_path / "dup_samples.tsv"
        f.write_text("gene\tS1\tS1\tS2\nG1\t1\t2\t3\n", encoding="utf-8")
        with pytest.raises(ExpressionFileError, match="Duplicate sample"):
            load_expression_matrix(f)

    def test_no_data_rows(self, tmp_path: Path) -> None:
        f = tmp_path / "header_only.tsv"
        f.write_text("gene\tS1\n", encoding="utf-8")
        with pytest.raises(ExpressionFileError, match="no data rows"):
            load_expression_matrix(f)

    def test_malformed_columns(self, tmp_path: Path) -> None:
        f = tmp_path / "malformed.tsv"
        f.write_text("gene\tS1\tS2\nG1\t1\nG2\t2\t3\n", encoding="utf-8")
        with pytest.raises(ExpressionFileError, match="expected 3 columns"):
            load_expression_matrix(f)

    def test_malformed_numeric_values(self, tmp_path: Path) -> None:
        f = tmp_path / "bad_numeric.tsv"
        f.write_text("gene\tS1\tS2\nG1\t1.5\tabc\nG2\t2.0\t3.0\n", encoding="utf-8")
        with pytest.raises(ExpressionFileError):
            load_expression_matrix(f)

    def test_missing_gene_identifiers(self, tmp_path: Path) -> None:
        f = tmp_path / "missing_genes.tsv"
        f.write_text("gene\tS1\tS2\n\t1\t2\nG2\t3\t4\n", encoding="utf-8")
        with pytest.raises(ExpressionFileError, match="Missing gene"):
            load_expression_matrix(f)

    def test_single_column_file(self, tmp_path: Path) -> None:
        f = tmp_path / "single_col.tsv"
        f.write_text("gene\nG1\nG2\n", encoding="utf-8")
        with pytest.raises(ExpressionFileError, match="at least 2 columns"):
            load_expression_matrix(f)


# ---------------------------------------------------------------------------
# Sample metadata loading
# ---------------------------------------------------------------------------


class TestLoadSampleMetadata:
    """Verify sample metadata loading."""

    def test_valid_json(self) -> None:
        src = FIXTURES / "synthetic_metadata.json"
        expr_ids = [
            "PC9_DMSO_1", "PC9_DMSO_2", "PC9_DMSO_3",
            "PC9_osi_DTP_1", "PC9_osi_DTP_2", "PC9_osi_DTP_3",
        ]
        result = load_sample_metadata(src, expression_sample_ids=expr_ids)
        assert result.sample_ids == set(expr_ids)
        assert result.conditions["PC9_DMSO_1"] == "DMSO"
        assert result.conditions["PC9_osi_DTP_1"] == "osi_DTP"

    def test_valid_tsv(self) -> None:
        src = FIXTURES / "synthetic_metadata.tsv"
        expr_ids = ["HCC2935_DMSO_1", "HCC2935_DMSO_2", "HCC2935_DMSO_3"]
        result = load_sample_metadata(src, expression_sample_ids=expr_ids)
        assert len(result.sample_ids) == 3
        assert result.conditions["HCC2935_DMSO_1"] == "DMSO"
        assert result.cell_lines["HCC2935_DMSO_1"] == "HCC2935"
        assert result.replicates["HCC2935_DMSO_1"] == 1
        assert result.batch_values["HCC2935_DMSO_1"] == "batch1"

    def test_file_not_found(self) -> None:
        with pytest.raises(FileNotFoundError):
            load_sample_metadata(Path("/nonexistent/metadata.json"))

    def test_missing_condition(self, tmp_path: Path) -> None:
        f = tmp_path / "no_condition.json"
        f.write_text(
            '{"samples": {"S1": {"cell_line": "PC9"}}}', encoding="utf-8"
        )
        with pytest.raises(MetadataFileError, match="'condition'"):
            load_sample_metadata(f)

    def test_duplicate_metadata_sample_ids(self, tmp_path: Path) -> None:
        f = tmp_path / "dup_meta.json"
        f.write_text(
            '{"samples": {'
            '"S1": {"condition": "DMSO"},'
            '"S1": {"condition": "treated"},'
            '"S2": {"condition": "DMSO"}}}',
            encoding="utf-8",
        )
        # JSON parser will use the last value for key "S1", so no duplicate
        # Instead create TSV with explicit duplicates
        tsv = tmp_path / "dup_meta.tsv"
        tsv.write_text(
            "sample_id\tcondition\nS1\tDMSO\nS1\ttreated\nS2\tDMSO\n",
            encoding="utf-8",
        )
        with pytest.raises(MetadataFileError, match="Duplicate sample"):
            load_sample_metadata(tsv)

    def test_expression_metadata_mismatch(self, tmp_path: Path) -> None:
        f = tmp_path / "mismatch.json"
        f.write_text(
            '{"samples": {'
            '"S1": {"condition": "DMSO"},'
            '"S2": {"condition": "treated"},'
            '"S3": {"condition": "DMSO"}}}',
            encoding="utf-8",
        )
        expr_ids = ["S1", "S2", "S4"]
        with pytest.raises(MetadataFileError, match="missing in metadata"):
            load_sample_metadata(f, expression_sample_ids=expr_ids)

    def test_extra_metadata_samples(self, tmp_path: Path) -> None:
        f = tmp_path / "extra.json"
        f.write_text(
            '{"samples": {'
            '"S1": {"condition": "DMSO"},'
            '"S2": {"condition": "treated"},'
            '"S3": {"condition": "DMSO"}}}',
            encoding="utf-8",
        )
        expr_ids = ["S1", "S2"]
        with pytest.raises(MetadataFileError, match="extra samples not in expression"):
            load_sample_metadata(f, expression_sample_ids=expr_ids)

    def test_missing_batch_recorded_as_unknown(self, tmp_path: Path) -> None:
        f = tmp_path / "no_batch.json"
        f.write_text(
            '{"samples": {'
            '"S1": {"condition": "DMSO"},'
            '"S2": {"condition": "treated"}}}',
            encoding="utf-8",
        )
        result = load_sample_metadata(f, expression_sample_ids=["S1", "S2"])
        assert result.batch_values["S1"] is None
        assert result.batch_values["S2"] is None


# ---------------------------------------------------------------------------
# Expression-value classification
# ---------------------------------------------------------------------------


class TestClassifyExpressionValues:
    """Verify expression-value classification."""

    def _make_df(
        self,
        values: list[list[float]],
        genes: list[str] | None = None,
    ) -> object:
        import pandas as pd
        if genes is None:
            genes = [f"G{i}" for i in range(len(values))]
        return pd.DataFrame(
            [[g] + row for g, row in zip(genes, values)],
            columns=["gene"] + [f"S{i}" for i in range(len(values[0]))],
        )

    def test_integer_count_classification(self) -> None:
        df = self._make_df([
            [100, 200, 300],
            [50, 60, 70],
        ])
        result = classify_expression_values(df)
        assert result.value_type == ValueType.RAW_INTEGER_COUNTS

    def test_non_integer_estimated_classification(self) -> None:
        df = self._make_df([
            [100.5, 200.3, 300.7],
            [50.1, 60.9, 70.2],
        ])
        result = classify_expression_values(df)
        assert result.value_type == ValueType.NON_INTEGER_ESTIMATED_COUNTS

    def test_transformed_value_classification(self) -> None:
        df = self._make_df([
            [-1.5, 2.3, 0.5],
            [1.2, -0.3, 4.0],
        ])
        result = classify_expression_values(df)
        assert result.value_type == ValueType.TRANSFORMED_VALUES

    def test_unknown_classification(self) -> None:
        # Empty dataframe with just gene column
        import pandas as pd
        df = pd.DataFrame({"gene": ["G1"], "S1": [0.0]})
        result = classify_expression_values(df)
        assert result.value_type == ValueType.UNKNOWN

    def test_override_classification(self) -> None:
        df = self._make_df([
            [100, 200, 300],
            [50, 60, 70],
        ])
        result = classify_expression_values(
            df, value_type_override=ValueType.TRANSFORMED_VALUES,
        )
        assert result.value_type == ValueType.TRANSFORMED_VALUES

    def test_mixed_integer_non_integer(self) -> None:
        df = self._make_df([
            [100, 200.5, 300],
            [50, 60, 70],
        ])
        result = classify_expression_values(df)
        assert result.value_type == ValueType.NORMALIZED_NONNEGATIVE_VALUES


# ---------------------------------------------------------------------------
# Gene identifier inspection
# ---------------------------------------------------------------------------


class TestInspectGeneIdentifiers:
    """Verify gene identifier inspection."""

    def test_ensembl_version_stripping(self) -> None:
        ids = ["ENSG00000146648.10", "ENSG00000141510.15", "ENSG00000012048.12"]
        result = inspect_gene_identifiers(ids)
        assert result.id_type == GeneIdType.ENSEMBL
        assert result.normalized_ids == [
            "ENSG00000146648",
            "ENSG00000141510",
            "ENSG00000012048",
        ]
        assert result.ensembl_count == 3

    def test_hgnc_symbol_detection(self) -> None:
        ids = ["EGFR", "TP53", "BRCA1", "MYC", "ZYX"]
        result = inspect_gene_identifiers(ids)
        assert result.id_type == GeneIdType.HGNC_SYMBOLS
        assert result.hgnc_count == 5

    def test_entrez_id_detection(self) -> None:
        ids = ["1956", "7157", "672", "4609"]
        result = inspect_gene_identifiers(ids)
        assert result.id_type == GeneIdType.ENTREZ_IDS
        assert result.entrez_count == 4

    def test_mixed_detection(self) -> None:
        ids = ["ENSG00000146648", "TP53", "1956"]
        result = inspect_gene_identifiers(ids)
        assert result.id_type == GeneIdType.MIXED

    def test_unknown_detection(self) -> None:
        ids = ["!invalid", "@weird", "#bad"]
        result = inspect_gene_identifiers(ids)
        assert result.id_type == GeneIdType.UNKNOWN
        assert result.unknown_count == 3

    def test_duplicate_gene_ids(self) -> None:
        ids = ["EGFR", "TP53", "EGFR", "MYC"]
        result = inspect_gene_identifiers(ids)
        assert "EGFR" in result.duplicate_ids

    def test_override_classification(self) -> None:
        ids = ["EGFR", "TP53"]
        result = inspect_gene_identifiers(
            ids, gene_id_type_override=GeneIdType.ENTREZ_IDS,
        )
        assert result.id_type == GeneIdType.ENTREZ_IDS


# ---------------------------------------------------------------------------
# Replicate counting
# ---------------------------------------------------------------------------


class TestCountReplicatesPerCondition:
    """Verify replicate counting."""

    def test_counts(self) -> None:
        conditions = {
            "S1": "DMSO",
            "S2": "DMSO",
            "S3": "DMSO",
            "S4": "treated",
            "S5": "treated",
        }
        counts = count_replicates_per_condition(conditions)
        assert counts["DMSO"] == 3
        assert counts["treated"] == 2


# ---------------------------------------------------------------------------
# Contrast validation
# ---------------------------------------------------------------------------


class TestValidateContrast:
    """Verify contrast validation."""

    def test_valid_contrast(self) -> None:
        conditions = {
            "S1": "DMSO", "S2": "DMSO", "S3": "DMSO",
            "S4": "treated", "S5": "treated",
        }
        result = validate_contrast(conditions, "DMSO", "treated")
        assert result["valid"] is True
        assert result["control_count"] == 3
        assert result["treatment_count"] == 2

    def test_missing_control(self) -> None:
        conditions = {"S1": "treated", "S2": "treated"}
        with pytest.raises(ContrastError, match="control condition"):
            validate_contrast(conditions, "DMSO", "treated")

    def test_missing_treatment(self) -> None:
        conditions = {"S1": "DMSO", "S2": "DMSO"}
        with pytest.raises(ContrastError, match="treatment condition"):
            validate_contrast(conditions, "DMSO", "treated")


# ---------------------------------------------------------------------------
# Full ingestion pipeline
# ---------------------------------------------------------------------------


class TestFullIngestion:
    """Verify the full ingestion pipeline."""

    def test_valid_tsv_ingestion(self, tmp_path: Path) -> None:
        run_dir = tmp_path / "run"
        result = ingest(
            expression_path=FIXTURES / "synthetic_expression.tsv",
            metadata_path=FIXTURES / "synthetic_metadata.json",
            source_id="GSE_SYNTHETIC",
            run_dir=run_dir,
        )
        assert result.n_genes == 5
        assert result.n_samples == 6
        assert result.source_id == "GSE_SYNTHETIC"

    def test_relative_persisted_paths(self, tmp_path: Path) -> None:
        run_dir = tmp_path / "run"
        result = ingest(
            expression_path=FIXTURES / "synthetic_expression.tsv",
            metadata_path=FIXTURES / "synthetic_metadata.json",
            source_id="GSE_SYNTHETIC",
            run_dir=run_dir,
        )
        # Paths should be relative to run_dir
        assert not Path(result.expression_path).is_absolute()
        assert not Path(result.metadata_path).is_absolute()
        # Stored files should exist
        assert (run_dir / result.expression_path).exists()
        assert (run_dir / result.metadata_path).exists()

    def test_independent_file_hashes(self, tmp_path: Path) -> None:
        run_dir = tmp_path / "run"
        result = ingest(
            expression_path=FIXTURES / "synthetic_expression.tsv",
            metadata_path=FIXTURES / "synthetic_metadata.json",
            source_id="GSE_SYNTHETIC",
            run_dir=run_dir,
        )
        # Hashes should be non-empty and different
        assert len(result.expression_sha256) == 64
        assert len(result.metadata_sha256) == 64
        assert result.expression_sha256 != result.metadata_sha256
        # Verify against direct hashing
        expr_hash = hash_file_sha256(FIXTURES / "synthetic_expression.tsv")
        assert result.expression_sha256 == expr_hash

    def test_json_serialization_roundtrip(self, tmp_path: Path) -> None:
        run_dir = tmp_path / "run"
        result = ingest(
            expression_path=FIXTURES / "synthetic_expression.tsv",
            metadata_path=FIXTURES / "synthetic_metadata.json",
            source_id="GSE_SYNTHETIC",
            run_dir=run_dir,
        )
        d = result.to_dict()
        # Should be valid JSON
        json_str = json.dumps(d)
        parsed = json.loads(json_str)
        assert parsed["source_id"] == "GSE_SYNTHETIC"
        assert parsed["n_genes"] == 5

    def test_windows_compatible_paths(self, tmp_path: Path) -> None:
        """Verify paths work on Windows (forward slashes in JSON)."""
        run_dir = tmp_path / "run"
        result = ingest(
            expression_path=FIXTURES / "synthetic_expression.tsv",
            metadata_path=FIXTURES / "synthetic_metadata.json",
            source_id="GSE_SYNTHETIC",
            run_dir=run_dir,
        )
        # Paths should not contain backslashes in the manifest
        assert "\\" not in result.expression_path
        assert "\\" not in result.metadata_path

    def test_gzip_tsv_ingestion(self, tmp_path: Path) -> None:
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
        result = ingest(
            expression_path=gz_path,
            metadata_path=meta_path,
            source_id="GSE_GZ_TEST",
            run_dir=run_dir,
        )
        assert result.n_genes == 5
        assert result.compression == CompressionType.GZIP

    def test_contrast_validation(self, tmp_path: Path) -> None:
        run_dir = tmp_path / "run"
        result = ingest(
            expression_path=FIXTURES / "synthetic_expression.tsv",
            metadata_path=FIXTURES / "synthetic_metadata.json",
            source_id="GSE_SYNTHETIC",
            run_dir=run_dir,
            contrast_control="DMSO",
            contrast_treatment="osi_DTP",
        )
        assert result.replicate_counts["DMSO"] == 3
        assert result.replicate_counts["osi_DTP"] == 3


# ---------------------------------------------------------------------------
# Manifest writing
# ---------------------------------------------------------------------------


class TestWriteIngestionManifest:
    """Verify manifest writing."""

    def test_manifest_written(self, tmp_path: Path) -> None:
        run_dir = tmp_path / "run"
        result = ingest(
            expression_path=FIXTURES / "synthetic_expression.tsv",
            metadata_path=FIXTURES / "synthetic_metadata.json",
            source_id="GSE_SYNTHETIC",
            run_dir=run_dir,
        )
        manifest_path = write_ingestion_manifest(run_dir, result)
        assert manifest_path.exists()
        assert manifest_path.name == "ingestion_manifest.json"

    def test_manifest_readable(self, tmp_path: Path) -> None:
        run_dir = tmp_path / "run"
        result = ingest(
            expression_path=FIXTURES / "synthetic_expression.tsv",
            metadata_path=FIXTURES / "synthetic_metadata.json",
            source_id="GSE_SYNTHETIC",
            run_dir=run_dir,
        )
        manifest_path = write_ingestion_manifest(run_dir, result)
        with open(manifest_path, encoding="utf-8") as f:
            data = json.load(f)
        assert data["source_id"] == "GSE_SYNTHETIC"
        assert data["n_genes"] == 5
        assert data["n_samples"] == 6
        assert "ingested_at" in data
        assert "software_version" in data
        assert "configuration_hash" in data
