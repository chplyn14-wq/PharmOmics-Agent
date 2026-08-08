"""Tests for pharmomics.omics.adapter — ingestion-to-omics conversion."""

from __future__ import annotations

import json
import re
from pathlib import Path

import pandas as pd

from pharmomics.ingestion.loader import (
    CompressionType,
    ExpressionLoadResult,
    GeneIdInspectionResult,
    MetadataLoadResult,
    ValueClassificationResult,
    ValueType,
)
from pharmomics.omics.adapter import from_load_results
from pharmomics.omics.enums import MeasurementType, NormalizationStatus
from pharmomics.omics.schemas import OmicsMatrix

FIXTURES = Path(__file__).parent / "fixtures"


# ---------------------------------------------------------------------------
# Test fixture helpers
# ---------------------------------------------------------------------------


def _make_expr_result(
    dataframe: pd.DataFrame | None = None,
    *,
    n_genes: int = 3,
    n_samples: int = 2,
    gene_ids: list[str] | None = None,
    sample_ids: list[str] | None = None,
    original_filename: str = "expr.tsv",
) -> ExpressionLoadResult:
    """Create a minimal ExpressionLoadResult for testing."""
    if dataframe is None:
        if gene_ids is None:
            gene_ids = [f"G{i}" for i in range(n_genes)]
        if sample_ids is None:
            sample_ids = [f"S{i}" for i in range(n_samples)]
        # Build rows: each gene is a row with numeric sample values
        rows = []
        for i, g in enumerate(gene_ids):
            row = [g] + [float(i * len(sample_ids) + j) for j in range(len(sample_ids))]
            rows.append(row)
        columns = ["gene"] + sample_ids
        dataframe = pd.DataFrame(rows, columns=columns)
    elif gene_ids is None:
        gene_ids = dataframe.iloc[:, 0].tolist()
    if sample_ids is None:
        sample_ids = [c for c in dataframe.columns if c != dataframe.columns[0]]
    return ExpressionLoadResult(
        dataframe=dataframe,
        n_genes=len(gene_ids),
        n_samples=len(sample_ids),
        sample_ids=sample_ids,
        gene_ids=gene_ids,
        delimiter="\t",
        compression=CompressionType.NONE,
        original_filename=original_filename,
    )


def _make_meta_result(
    *,
    sample_ids: list[str] | None = None,
    conditions: dict[str, str] | None = None,
    original_filename: str = "meta.json",
) -> MetadataLoadResult:
    """Create a minimal MetadataLoadResult for testing."""
    if sample_ids is None:
        sample_ids = ["S0", "S1"]
    if conditions is None:
        conditions = {sid: "control" for sid in sample_ids}
    conds_list = [conditions.get(sid, "control") for sid in sample_ids]
    return MetadataLoadResult(
        dataframe=pd.DataFrame(
            {
                "sample_id": sample_ids,
                "condition": conds_list,
            }
        ),
        sample_ids=set(sample_ids),
        conditions=conditions,
        cell_lines={sid: None for sid in sample_ids},
        replicates={sid: None for sid in sample_ids},
        batch_values={sid: None for sid in sample_ids},
        original_filename=original_filename,
    )


def _make_value_class(
    value_type: ValueType = ValueType.RAW_INTEGER_COUNTS,
    total_values: int = 6,
    integer_count: int = 6,
) -> ValueClassificationResult:
    """Create a minimal ValueClassificationResult."""
    return ValueClassificationResult(
        value_type=value_type,
        total_values=total_values,
        integer_count=integer_count,
        non_integer_count=total_values - integer_count,
        zero_count=0,
        negative_count=0,
        min_value=0.0,
        max_value=float(total_values),
    )


def _make_gene_inspection(
    gene_ids: list[str] | None = None,
) -> GeneIdInspectionResult:
    """Create a minimal GeneIdInspectionResult."""
    if gene_ids is None:
        gene_ids = ["EGFR", "TP53", "MYC"]
    return GeneIdInspectionResult(
        id_type=gene_ids[0].startswith("ENS") and "ensembl_ids" or "hgnc_symbols",
        original_ids=gene_ids,
        normalized_ids=[re.sub(r"\.[0-9]+$", "", g) for g in gene_ids],
        duplicate_ids=[],
        missing_ids=[],
        ensembl_count=sum(1 for g in gene_ids if g.startswith("ENS")),
        hgnc_count=sum(1 for g in gene_ids if not g.startswith("ENS")),
        entrez_count=0,
        unknown_count=0,
    )


# ---------------------------------------------------------------------------
# Basic construction
# ---------------------------------------------------------------------------


class TestFromLoadResultsBasic:
    """Verify core OmicsMatrix construction from load results."""

    def test_returns_omics_matrix(self) -> None:
        expr = _make_expr_result()
        meta = _make_meta_result()
        result = from_load_results(expr, meta)
        assert isinstance(result, OmicsMatrix)

    def test_n_features_matches(self) -> None:
        expr = _make_expr_result(n_genes=5, n_samples=2)
        meta = _make_meta_result(sample_ids=["S0", "S1"])
        result = from_load_results(expr, meta)
        assert result.n_features == 5

    def test_n_samples_matches(self) -> None:
        expr = _make_expr_result(n_genes=3, n_samples=4)
        meta = _make_meta_result(sample_ids=["S0", "S1", "S2", "S3"])
        result = from_load_results(expr, meta)
        assert result.n_samples == 4

    def test_feature_ids_order(self) -> None:
        expr = _make_expr_result(gene_ids=["EGFR", "TP53", "BRCA1"])
        meta = _make_meta_result(sample_ids=["S0", "S1"])
        result = from_load_results(expr, meta)
        assert result.feature_ids == ["EGFR", "TP53", "BRCA1"]

    def test_sample_ids_order(self) -> None:
        expr = _make_expr_result(
            gene_ids=["G0"],
            sample_ids=["PC9_DMSO_1", "PC9_osi_DTP_1"],
        )
        meta = _make_meta_result(
            sample_ids=["PC9_DMSO_1", "PC9_osi_DTP_1"],
            conditions={"PC9_DMSO_1": "DMSO", "PC9_osi_DTP_1": "osi_DTP"},
        )
        result = from_load_results(expr, meta)
        assert result.sample_ids == ["PC9_DMSO_1", "PC9_osi_DTP_1"]

    def test_dataframe_carried_through(self) -> None:
        expr = _make_expr_result()
        meta = _make_meta_result()
        result = from_load_results(expr, meta)
        assert isinstance(result.dataframe, pd.DataFrame)
        assert result.dataframe.shape == (3, 3)  # gene col + 2 samples

    def test_modality_is_transcriptomics(self) -> None:
        expr = _make_expr_result()
        meta = _make_meta_result()
        result = from_load_results(expr, meta)
        assert result.modality == "transcriptomics"

    def test_feature_type_is_gene(self) -> None:
        expr = _make_expr_result()
        meta = _make_meta_result()
        result = from_load_results(expr, meta)
        assert result.feature_type == "gene"

    def test_schema_version(self) -> None:
        expr = _make_expr_result()
        meta = _make_meta_result()
        result = from_load_results(expr, meta)
        assert result.schema_version == "1.0.0"


# ---------------------------------------------------------------------------
# Auto-generated IDs and timestamps
# ---------------------------------------------------------------------------


class TestAutoGeneratedFields:
    """Verify auto-generated matrix_id and created_at."""

    def test_auto_matrix_id_format(self) -> None:
        expr = _make_expr_result()
        meta = _make_meta_result()
        result = from_load_results(expr, meta, source_id="GSE_TEST")
        assert result.matrix_id.startswith("mx-GSE_TEST-")
        # Should have 8 hex chars after the dash
        parts = result.matrix_id.split("-")
        assert len(parts[-1]) == 8

    def test_matrix_ids_are_unique(self) -> None:
        expr = _make_expr_result()
        meta = _make_meta_result()
        r1 = from_load_results(expr, meta, source_id="GSE_A")
        r2 = from_load_results(expr, meta, source_id="GSE_A")
        assert r1.matrix_id != r2.matrix_id

    def test_custom_matrix_id(self) -> None:
        expr = _make_expr_result()
        meta = _make_meta_result()
        result = from_load_results(expr, meta, matrix_id="custom-mx-001")
        assert result.matrix_id == "custom-mx-001"

    def test_auto_created_at_iso_format(self) -> None:
        expr = _make_expr_result()
        meta = _make_meta_result()
        result = from_load_results(expr, meta)
        # Should match ISO-8601 pattern
        assert re.match(r"\d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2}Z", result.created_at)

    def test_custom_created_at(self) -> None:
        expr = _make_expr_result()
        meta = _make_meta_result()
        result = from_load_results(expr, meta, created_at="2025-06-01T00:00:00Z")
        assert result.created_at == "2025-06-01T00:00:00Z"


# ---------------------------------------------------------------------------
# Value classification mapping
# ---------------------------------------------------------------------------


class TestValueClassificationMapping:
    """Verify MeasurementType and NormalizationStatus inference."""

    def test_raw_integer_counts(self) -> None:
        expr = _make_expr_result()
        meta = _make_meta_result()
        vc = _make_value_class(ValueType.RAW_INTEGER_COUNTS)
        result = from_load_results(expr, meta, value_class=vc)
        assert result.measurement_type == MeasurementType.RAW_COUNTS.value
        assert result.normalization_status == NormalizationStatus.RAW.value

    def test_non_integer_estimated_counts(self) -> None:
        expr = _make_expr_result()
        meta = _make_meta_result()
        vc = _make_value_class(ValueType.NON_INTEGER_ESTIMATED_COUNTS)
        result = from_load_results(expr, meta, value_class=vc)
        assert result.measurement_type == MeasurementType.ESTIMATED_COUNTS.value
        assert result.normalization_status == NormalizationStatus.UNKNOWN.value

    def test_transformed_values(self) -> None:
        expr = _make_expr_result()
        meta = _make_meta_result()
        vc = _make_value_class(ValueType.TRANSFORMED_VALUES)
        result = from_load_results(expr, meta, value_class=vc)
        assert result.measurement_type == MeasurementType.UNKNOWN.value
        assert result.normalization_status == NormalizationStatus.TRANSFORMED.value

    def test_normalized_nonnegative_values(self) -> None:
        expr = _make_expr_result()
        meta = _make_meta_result()
        vc = _make_value_class(ValueType.NORMALIZED_NONNEGATIVE_VALUES)
        result = from_load_results(expr, meta, value_class=vc)
        assert result.measurement_type == MeasurementType.UNKNOWN.value
        assert result.normalization_status == NormalizationStatus.UNKNOWN.value

    def test_unknown_value_type(self) -> None:
        expr = _make_expr_result()
        meta = _make_meta_result()
        vc = _make_value_class(ValueType.UNKNOWN)
        result = from_load_results(expr, meta, value_class=vc)
        assert result.measurement_type == MeasurementType.UNKNOWN.value
        assert result.normalization_status == NormalizationStatus.UNKNOWN.value

    def test_no_value_class_defaults_to_unknown(self) -> None:
        expr = _make_expr_result()
        meta = _make_meta_result()
        result = from_load_results(expr, meta)
        assert result.measurement_type == MeasurementType.UNKNOWN.value
        assert result.normalization_status == NormalizationStatus.UNKNOWN.value


# ---------------------------------------------------------------------------
# Sample metadata propagation
# ---------------------------------------------------------------------------


class TestSampleMetadataPropagation:
    """Verify sample metadata flows into OmicsMatrix."""

    def test_sample_count(self) -> None:
        sids = ["S0", "S1", "S2"]
        expr = _make_expr_result(sample_ids=sids)
        meta = _make_meta_result(sample_ids=sids)
        result = from_load_results(expr, meta)
        assert len(result.sample_metadata) == 3

    def test_condition_propagates(self) -> None:
        sids = ["S0", "S1"]
        conds = {"S0": "DMSO", "S1": "treated"}
        expr = _make_expr_result(sample_ids=sids)
        meta = _make_meta_result(sample_ids=sids, conditions=conds)
        result = from_load_results(expr, meta)
        assert result.sample_metadata["S0"].condition == "DMSO"
        assert result.sample_metadata["S1"].condition == "treated"

    def test_cell_line_in_annotations(self) -> None:
        sids = ["PC9_1", "PC9_2"]
        expr = _make_expr_result(sample_ids=sids)
        meta = MetadataLoadResult(
            dataframe=pd.DataFrame({"sample_id": sids, "condition": ["DMSO"] * 2}),
            sample_ids=set(sids),
            conditions={"PC9_1": "DMSO", "PC9_2": "DMSO"},
            cell_lines={"PC9_1": "PC9", "PC9_2": "PC9"},
            replicates={"PC9_1": 1, "PC9_2": 2},
            batch_values={"PC9_1": "batch1", "PC9_2": "batch1"},
            original_filename="meta.json",
        )
        result = from_load_results(expr, meta)
        assert result.sample_metadata["PC9_1"].annotations["cell_line"] == "PC9"

    def test_replicate_in_annotations(self) -> None:
        sids = ["S0", "S1"]
        expr = _make_expr_result(sample_ids=sids)
        meta = MetadataLoadResult(
            dataframe=pd.DataFrame({"sample_id": sids, "condition": ["ctrl"] * 2}),
            sample_ids=set(sids),
            conditions={"S0": "ctrl", "S1": "ctrl"},
            cell_lines={"S0": None, "S1": None},
            replicates={"S0": 1, "S1": 2},
            batch_values={"S0": None, "S1": None},
            original_filename="meta.json",
        )
        result = from_load_results(expr, meta)
        # At least one sample has replicate in annotations
        has_replicate = any(
            "replicate" in sm.annotations for sm in result.sample_metadata.values()
        )
        assert has_replicate

    def test_batch_in_annotations(self) -> None:
        sids = ["S0", "S1"]
        expr = _make_expr_result(sample_ids=sids)
        meta = MetadataLoadResult(
            dataframe=pd.DataFrame({"sample_id": sids, "condition": ["ctrl"] * 2}),
            sample_ids=set(sids),
            conditions={"S0": "ctrl", "S1": "ctrl"},
            cell_lines={"S0": None, "S1": None},
            replicates={"S0": None, "S1": None},
            batch_values={"S0": "B1", "S1": "B2"},
            original_filename="meta.json",
        )
        result = from_load_results(expr, meta)
        assert result.sample_metadata["S0"].annotations.get("batch") == "B1"
        assert result.sample_metadata["S1"].annotations.get("batch") == "B2"

    def test_null_fields_excluded_from_annotations(self) -> None:
        """None values should not appear in annotations."""
        sids = ["S0"]
        expr = _make_expr_result(sample_ids=sids)
        meta = _make_meta_result(sample_ids=sids)
        result = from_load_results(expr, meta)
        ann = result.sample_metadata["S0"].annotations
        assert "cell_line" not in ann
        assert "replicate" not in ann
        assert "batch" not in ann


# ---------------------------------------------------------------------------
# Feature metadata propagation
# ---------------------------------------------------------------------------


class TestFeatureMetadataPropagation:
    """Verify feature metadata flows into OmicsMatrix."""

    def test_feature_count(self) -> None:
        expr = _make_expr_result(gene_ids=["G0", "G1", "G2"])
        meta = _make_meta_result()
        result = from_load_results(expr, meta)
        assert len(result.feature_metadata) == 3

    def test_feature_ids_match_keys(self) -> None:
        expr = _make_expr_result(gene_ids=["EGFR", "TP53"])
        meta = _make_meta_result()
        result = from_load_results(expr, meta)
        assert set(result.feature_metadata.keys()) == {"EGFR", "TP53"}

    def test_without_gene_inspection(self) -> None:
        expr = _make_expr_result(gene_ids=["G0", "G1"])
        meta = _make_meta_result()
        result = from_load_results(expr, meta)
        # Metadata exists but minimal
        assert result.feature_metadata["G0"].feature_id == "G0"
        assert result.feature_metadata["G0"].normalized_id is None

    def test_with_gene_inspection_ensembl_stripping(self) -> None:
        gene_ids = ["ENSG00000146648.10", "ENSG00000141510.15"]
        expr = _make_expr_result(gene_ids=gene_ids)
        meta = _make_meta_result()
        gi = _make_gene_inspection(gene_ids)
        result = from_load_results(expr, meta, gene_inspection=gi)
        fm = result.feature_metadata["ENSG00000146648.10"]
        assert fm.normalized_id == "ENSG00000146648"

    def test_with_gene_inspection_hgnc_no_stripping(self) -> None:
        gene_ids = ["EGFR", "TP53"]
        expr = _make_expr_result(gene_ids=gene_ids)
        meta = _make_meta_result()
        gi = _make_gene_inspection(gene_ids)
        result = from_load_results(expr, meta, gene_inspection=gi)
        # HGNC symbols don't get normalized
        assert result.feature_metadata["EGFR"].normalized_id is None


# ---------------------------------------------------------------------------
# Provenance
# ---------------------------------------------------------------------------


class TestProvenance:
    """Verify provenance record construction."""

    def test_provenance_created_with_source_id(self) -> None:
        expr = _make_expr_result()
        meta = _make_meta_result()
        result = from_load_results(expr, meta, source_id="GSE193258", sha256="abc123")
        assert len(result.provenance) == 1
        assert result.provenance[0].source_id == "GSE193258"
        assert result.provenance[0].sha256 == "abc123"
        assert result.provenance[0].source_file == "expr.tsv"

    def test_no_provenance_without_source_id(self) -> None:
        expr = _make_expr_result()
        meta = _make_meta_result()
        result = from_load_results(expr, meta)
        assert result.provenance == []

    def test_software_version_in_provenance(self) -> None:
        expr = _make_expr_result()
        meta = _make_meta_result()
        result = from_load_results(
            expr,
            meta,
            source_id="GSE_X",
            software_version="0.2.0",
        )
        assert result.provenance[0].software_version == "0.2.0"


# ---------------------------------------------------------------------------
# model_dump and descriptor
# ---------------------------------------------------------------------------


class TestSerialization:
    """Verify model_dump excludes dataframe and descriptor property works."""

    def test_model_dump_excludes_dataframe(self) -> None:
        expr = _make_expr_result()
        meta = _make_meta_result()
        result = from_load_results(expr, meta)
        d = result.model_dump()
        assert "dataframe" not in d

    def test_model_dump_is_json_serializable(self) -> None:
        expr = _make_expr_result()
        meta = _make_meta_result()
        result = from_load_results(expr, meta)
        d = result.model_dump()
        # Should not raise
        json.dumps(d)

    def test_model_dump_preserves_provenance(self) -> None:
        expr = _make_expr_result()
        meta = _make_meta_result()
        result = from_load_results(expr, meta, source_id="GSE_TEST")
        d = result.model_dump()
        assert d["provenance"][0]["source_id"] == "GSE_TEST"

    def test_descriptor_property(self) -> None:
        expr = _make_expr_result()
        meta = _make_meta_result()
        result = from_load_results(expr, meta)
        desc = result.descriptor
        assert desc.matrix_id == result.matrix_id
        assert desc.n_features == result.n_features
        assert desc.n_samples == result.n_samples

    def test_descriptor_is_json_serializable(self) -> None:
        expr = _make_expr_result()
        meta = _make_meta_result()
        result = from_load_results(expr, meta)
        desc = result.descriptor
        # Should not raise
        json.dumps(desc.model_dump())


# ---------------------------------------------------------------------------
# Integration: real fixture data
# ---------------------------------------------------------------------------


class TestIntegrationWithFixtures:
    """Verify adapter works with real ingestion fixtures."""

    def test_from_fixture_expression_and_metadata(self) -> None:
        """Load real fixtures and construct an OmicsMatrix."""
        from pharmomics.ingestion.loader import (
            load_expression_matrix,
            load_sample_metadata,
        )

        expr = load_expression_matrix(FIXTURES / "synthetic_expression.tsv")
        meta = load_sample_metadata(
            FIXTURES / "synthetic_metadata.json",
            expression_sample_ids=expr.sample_ids,
        )
        result = from_load_results(
            expr, meta, source_id="GSE_SYNTHETIC", sha256="test-hash"
        )
        assert result.n_features == 5
        assert result.n_samples == 6
        assert result.feature_ids == ["EGFR", "ERBB2", "TP53", "BRCA1", "MYC"]
        assert "PC9_DMSO_1" in result.sample_metadata
        assert result.sample_metadata["PC9_DMSO_1"].condition == "DMSO"
        assert len(result.provenance) == 1
