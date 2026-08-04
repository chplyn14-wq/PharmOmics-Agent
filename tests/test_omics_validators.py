"""Tests for pharmomics.omics.validation — OmicsMatrix validators."""

from __future__ import annotations

from datetime import UTC, datetime

import pandas as pd

from pharmomics.omics.schemas import (
    FeatureMetadata,
    OmicsMatrix,
    ProvenanceRecord,
    SampleMetadata,
)
from pharmomics.omics.validation import (
    _check_counts_consistency,
    _check_dataframe_shape,
    _check_feature_metadata_coverage,
    _check_nonempty_ids,
    _check_sample_metadata_coverage,
    _check_schema_version,
    _check_unique_ids,
    _find_duplicates,
    validate,
)

# ---------------------------------------------------------------------------
# Helper: build a valid OmicsMatrix directly (bypass adapter)
# ---------------------------------------------------------------------------


def _make_valid_matrix(
    *,
    n_features: int = 3,
    n_samples: int = 2,
    gene_ids: list[str] | None = None,
    sample_ids: list[str] | None = None,
    schema_version: str = "1.0.0",
    omit_feature_meta: set[str] | None = None,
    omit_sample_meta: set[str] | None = None,
    extra_rows: int = 0,
    extra_cols: int = 0,
    feature_metadata: dict[str, FeatureMetadata] | None = None,
    sample_metadata: dict[str, SampleMetadata] | None = None,
) -> OmicsMatrix:
    """Construct an OmicsMatrix, optionally with deliberate violations."""
    if gene_ids is None:
        gene_ids = [f"G{i}" for i in range(n_features)]
    if sample_ids is None:
        sample_ids = [f"S{i}" for i in range(n_samples)]

    # Build dataframe using actual sample_ids length (not n_samples which
    # may be intentionally wrong for the test).
    actual_rows = n_features + extra_rows
    actual_cols = len(sample_ids) + 1 + extra_cols
    rows = []
    for i in range(actual_rows):
        gid = gene_ids[i] if i < len(gene_ids) else f"G_extra{i}"
        vals = [float(i * actual_cols + j) for j in range(1, actual_cols)]
        rows.append([gid] + vals)
    cols = ["gene"] + sample_ids + [
        f"extra_{j}" for j in range(extra_cols)
    ]
    df = pd.DataFrame(rows, columns=cols)

    # Build feature metadata
    if feature_metadata is not None:
        fm = feature_metadata
    elif omit_feature_meta:
        fm = {}
        for fid in gene_ids:
            if fid not in omit_feature_meta:
                fm[fid] = FeatureMetadata(feature_id=fid)
    else:
        fm = {fid: FeatureMetadata(feature_id=fid) for fid in gene_ids}

    # Build sample metadata
    if sample_metadata is not None:
        sm = sample_metadata
    elif omit_sample_meta:
        sm = {}
        for sid in sample_ids:
            if sid not in omit_sample_meta:
                sm[sid] = SampleMetadata(sample_id=sid)
    else:
        sm = {sid: SampleMetadata(sample_id=sid) for sid in sample_ids}

    now = datetime.now(UTC).strftime("%Y-%m-%dT%H:%M:%SZ")

    return OmicsMatrix(
        matrix_id="test-mx-001",
        schema_version=schema_version,
        modality="transcriptomics",
        feature_type="gene",
        measurement_type="raw_counts",
        normalization_status="raw",
        n_features=n_features,
        n_samples=n_samples,
        feature_ids=gene_ids,
        sample_ids=sample_ids,
        dataframe=df,
        feature_metadata=fm,
        sample_metadata=sm,
        provenance=[
            ProvenanceRecord(
                source_id="GSE_TEST",
                source_file="expr.tsv",
                sha256="abc123",
                ingested_at=now,
                software_version="0.1.0",
            )
        ],
        created_at=now,
    )


# ---------------------------------------------------------------------------
# _find_duplicates utility
# ---------------------------------------------------------------------------


class TestFindDuplicates:
    """Verify the _find_duplicates helper."""

    def test_no_duplicates(self) -> None:
        assert _find_duplicates(["A", "B", "C"]) == set()

    def test_one_duplicate(self) -> None:
        assert _find_duplicates(["A", "B", "A"]) == {"A"}

    def test_multiple_duplicates(self) -> None:
        assert _find_duplicates(["A", "B", "A", "B", "C"]) == {"A", "B"}

    def test_empty_list(self) -> None:
        assert _find_duplicates([]) == set()

    def test_all_same(self) -> None:
        assert _find_duplicates(["X", "X", "X"]) == {"X"}


# ---------------------------------------------------------------------------
# _check_schema_version
# ---------------------------------------------------------------------------


class TestCheckSchemaVersion:
    """Verify schema version validation."""

    def test_valid_version(self) -> None:
        matrix = _make_valid_matrix()
        assert _check_schema_version(matrix) == []

    def test_unknown_version(self) -> None:
        matrix = _make_valid_matrix(schema_version="9.9.9")
        violations = _check_schema_version(matrix)
        assert len(violations) == 1
        assert "9.9.9" in violations[0]

    def test_empty_version(self) -> None:
        matrix = _make_valid_matrix(schema_version="")
        violations = _check_schema_version(matrix)
        assert len(violations) == 1


# ---------------------------------------------------------------------------
# _check_counts_consistency
# ---------------------------------------------------------------------------


class TestCheckCountsConsistency:
    """Verify n_features/n_samples vs id list length."""

    def test_consistent(self) -> None:
        matrix = _make_valid_matrix()
        assert _check_counts_consistency(matrix) == []

    def test_n_features_mismatch(self) -> None:
        matrix = _make_valid_matrix(n_features=5, gene_ids=["G0", "G1"])
        violations = _check_counts_consistency(matrix)
        assert len(violations) == 1
        assert "n_features" in violations[0]

    def test_n_samples_mismatch(self) -> None:
        matrix = _make_valid_matrix(n_samples=4, sample_ids=["S0", "S1"])
        violations = _check_counts_consistency(matrix)
        assert len(violations) == 1
        assert "n_samples" in violations[0]


# ---------------------------------------------------------------------------
# _check_feature_metadata_coverage
# ---------------------------------------------------------------------------


class TestCheckFeatureMetadataCoverage:
    """Verify every feature_id has metadata."""

    def test_full_coverage(self) -> None:
        matrix = _make_valid_matrix()
        assert _check_feature_metadata_coverage(matrix) == []

    def test_missing_one(self) -> None:
        matrix = _make_valid_matrix(
            gene_ids=["EGFR", "TP53", "MYC"],
            omit_feature_meta={"TP53"},
        )
        violations = _check_feature_metadata_coverage(matrix)
        assert len(violations) == 1
        assert "TP53" in violations[0]

    def test_missing_multiple(self) -> None:
        matrix = _make_valid_matrix(
            gene_ids=["G0", "G1", "G2", "G3", "G4"],
            omit_feature_meta={"G1", "G3"},
        )
        violations = _check_feature_metadata_coverage(matrix)
        assert len(violations) == 1
        assert "2" in violations[0]


# ---------------------------------------------------------------------------
# _check_sample_metadata_coverage
# ---------------------------------------------------------------------------


class TestCheckSampleMetadataCoverage:
    """Verify every sample_id has metadata."""

    def test_full_coverage(self) -> None:
        matrix = _make_valid_matrix()
        assert _check_sample_metadata_coverage(matrix) == []

    def test_missing_one(self) -> None:
        matrix = _make_valid_matrix(
            sample_ids=["S0", "S1", "S2"],
            omit_sample_meta={"S1"},
        )
        violations = _check_sample_metadata_coverage(matrix)
        assert len(violations) == 1
        assert "S1" in violations[0]


# ---------------------------------------------------------------------------
# _check_dataframe_shape
# ---------------------------------------------------------------------------


class TestCheckDataFrameShape:
    """Verify DataFrame dimensions match n_features/n_samples."""

    def test_correct_shape(self) -> None:
        matrix = _make_valid_matrix()
        assert _check_dataframe_shape(matrix) == []

    def test_extra_rows(self) -> None:
        matrix = _make_valid_matrix(extra_rows=2)
        violations = _check_dataframe_shape(matrix)
        assert len(violations) == 1
        assert "rows" in violations[0]

    def test_extra_cols(self) -> None:
        matrix = _make_valid_matrix(extra_cols=3)
        violations = _check_dataframe_shape(matrix)
        assert len(violations) == 1
        assert "columns" in violations[0]


# ---------------------------------------------------------------------------
# _check_unique_ids
# ---------------------------------------------------------------------------


class TestCheckUniqueIds:
    """Verify no duplicate IDs."""

    def test_all_unique(self) -> None:
        matrix = _make_valid_matrix()
        assert _check_unique_ids(matrix) == []

    def test_duplicate_feature_ids(self) -> None:
        matrix = _make_valid_matrix(
            n_features=3,
            gene_ids=["EGFR", "EGFR", "TP53"],
        )
        violations = _check_unique_ids(matrix)
        assert len(violations) == 1
        assert "EGFR" in violations[0]

    def test_duplicate_sample_ids(self) -> None:
        matrix = _make_valid_matrix(
            n_samples=3,
            sample_ids=["S0", "S0", "S1"],
        )
        violations = _check_unique_ids(matrix)
        assert len(violations) == 1
        assert "S0" in violations[0]


# ---------------------------------------------------------------------------
# _check_nonempty_ids
# ---------------------------------------------------------------------------


class TestCheckNonEmptyIds:
    """Verify no empty strings in ID lists."""

    def test_no_empty(self) -> None:
        matrix = _make_valid_matrix()
        assert _check_nonempty_ids(matrix) == []

    def test_empty_feature_id(self) -> None:
        matrix = _make_valid_matrix(
            gene_ids=["EGFR", "", "TP53"],
        )
        violations = _check_nonempty_ids(matrix)
        assert len(violations) == 1
        assert "feature_id" in violations[0]

    def test_empty_sample_id(self) -> None:
        matrix = _make_valid_matrix(
            sample_ids=["S0", " ", "S2"],
        )
        violations = _check_nonempty_ids(matrix)
        assert len(violations) == 1
        assert "sample_id" in violations[0]


# ---------------------------------------------------------------------------
# Full validate() integration
# ---------------------------------------------------------------------------


class TestValidateFull:
    """Verify the combined validate() function."""

    def test_valid_matrix_no_violations(self) -> None:
        matrix = _make_valid_matrix()
        assert validate(matrix) == []

    def test_fixture_matrix_valid(self) -> None:
        from pathlib import Path

        from pharmomics.ingestion.loader import (
            load_expression_matrix,
            load_sample_metadata,
        )
        from pharmomics.omics.adapter import from_load_results

        fixtures = Path(__file__).parent / "fixtures"
        expr = load_expression_matrix(
            fixtures / "synthetic_expression.tsv"
        )
        meta = load_sample_metadata(
            fixtures / "synthetic_metadata.json",
            expression_sample_ids=expr.sample_ids,
        )
        matrix = from_load_results(
            expr, meta,
            source_id="GSE_SYNTHETIC",
            sha256="test-hash",
        )
        assert validate(matrix) == []
