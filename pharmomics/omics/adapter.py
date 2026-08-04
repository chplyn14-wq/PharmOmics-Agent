"""Adapter: construct an ``OmicsMatrix`` from ingestion-layer objects.

This module bridges the ingestion I/O layer (``ExpressionLoadResult``,
``MetadataLoadResult``) and the domain layer (``OmicsMatrix``).  It
performs pure in-memory conversion — no file I/O occurs here.
"""

from __future__ import annotations

import secrets
from datetime import UTC, datetime
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from pharmomics.ingestion.loader import (
        ExpressionLoadResult,
        GeneIdInspectionResult,
        MetadataLoadResult,
        ValueClassificationResult,
    )
    from pharmomics.omics.enums import MeasurementType, NormalizationStatus

from pharmomics.omics.schemas import (
    FeatureMetadata,
    OmicsMatrix,
    ProvenanceRecord,
    SampleMetadata,
)

# ---------------------------------------------------------------------------
# Private helpers
# ---------------------------------------------------------------------------


def _infer_measurement_type(
    value_class: ValueClassificationResult,
) -> MeasurementType:
    """Derive a MeasurementType from the ingestion value classification."""
    from pharmomics.omics.enums import MeasurementType

    vt = value_class.value_type.value
    if vt == "raw_integer_counts":
        return MeasurementType.RAW_COUNTS
    if vt == "non_integer_estimated_counts":
        return MeasurementType.ESTIMATED_COUNTS
    if vt == "transformed_values":
        return MeasurementType.UNKNOWN
    if vt == "normalized_nonnegative_values":
        return MeasurementType.UNKNOWN
    return MeasurementType.UNKNOWN


def _infer_normalization_status(
    value_class: ValueClassificationResult,
) -> NormalizationStatus:
    """Derive a NormalizationStatus from the ingestion value classification."""
    from pharmomics.omics.enums import NormalizationStatus

    vt = value_class.value_type.value
    if vt == "raw_integer_counts":
        return NormalizationStatus.RAW
    if vt == "non_integer_estimated_counts":
        return NormalizationStatus.UNKNOWN
    if vt == "transformed_values":
        return NormalizationStatus.TRANSFORMED
    if vt == "normalized_nonnegative_values":
        return NormalizationStatus.UNKNOWN
    return NormalizationStatus.UNKNOWN


def _build_feature_metadata(
    gene_ids: list[str],
    gene_inspection: GeneIdInspectionResult | None,
) -> dict[str, FeatureMetadata]:
    """Create FeatureMetadata entries for every gene/feature."""
    meta: dict[str, FeatureMetadata] = {}
    for gid in gene_ids:
        gid = str(gid).strip()
        norm_id: str | None = None
        annotations: dict[str, str] = {}

        if gene_inspection is not None:
            orig = gene_inspection.original_ids
            idx = orig.index(gid) if gid in orig else None
            if idx is not None and idx < len(gene_inspection.normalized_ids):
                norm_val = gene_inspection.normalized_ids[idx]
                if norm_val != gid:
                    norm_id = norm_val

            if gid in gene_inspection.original_ids:
                g_idx = gene_inspection.original_ids.index(gid)
                # Determine annotation key based on what type this gene is
                if g_idx < len(gene_inspection.original_ids):
                    annotations["raw_id"] = gid

        meta[gid] = FeatureMetadata(
            feature_id=gid,
            annotations=annotations,
            normalized_id=norm_id,
        )
    return meta


def _build_sample_metadata(
    meta_result: MetadataLoadResult,
) -> dict[str, SampleMetadata]:
    """Create SampleMetadata entries for every sample."""
    meta: dict[str, SampleMetadata] = {}
    for sid in meta_result.sample_ids:
        annotations: dict[str, object] = {}
        cl = meta_result.cell_lines.get(sid)
        if cl is not None:
            annotations["cell_line"] = cl
        rep = meta_result.replicates.get(sid)
        if rep is not None:
            annotations["replicate"] = rep
        batch = meta_result.batch_values.get(sid)
        if batch is not None:
            annotations["batch"] = batch

        meta[sid] = SampleMetadata(
            sample_id=sid,
            condition=meta_result.conditions.get(sid),
            annotations=annotations,
        )
    return meta


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------


def from_load_results(
    expr_result: ExpressionLoadResult,
    meta_result: MetadataLoadResult,
    *,
    value_class: ValueClassificationResult | None = None,
    gene_inspection: GeneIdInspectionResult | None = None,
    source_id: str = "",
    sha256: str = "",
    software_version: str = "",
    matrix_id: str | None = None,
    created_at: str | None = None,
) -> OmicsMatrix:
    """Construct an ``OmicsMatrix`` from ingestion-layer load results.

    This is a **pure in-memory conversion** — it reads no files and
    performs no I/O.  All data comes from objects already resident in
    memory after a successful ingestion run.

    Parameters
    ----------
    expr_result : ExpressionLoadResult
        The expression matrix load result (contains the dataframe).
    meta_result : MetadataLoadResult
        The sample metadata load result.
    value_class : ValueClassificationResult, optional
        Expression value classification.  If omitted, a conservative
        unknown mapping is used.
    gene_inspection : GeneIdInspectionResult, optional
        Gene identifier inspection result.  If omitted, feature metadata
        contains only the raw feature_id.
    source_id : str, optional
        Source identifier (e.g. GEO accession).  Used for provenance.
    sha256 : str, optional
        SHA-256 hash of the source file.  Used for provenance.
    software_version : str, optional
        Software version string.  Used for provenance.
    matrix_id : str, optional
        Unique matrix identifier.  Auto-generated if omitted.
    created_at : str, optional
        ISO-8601 creation timestamp.  Auto-generated if omitted.

    Returns
    -------
    OmicsMatrix
        A fully populated domain object ready for downstream analysis.
    """
    if matrix_id is None:
        matrix_id = f"mx-{source_id}-{secrets.token_hex(4)}"

    if created_at is None:
        created_at = datetime.now(UTC).strftime("%Y-%m-%dT%H:%M:%SZ")

    # Infer measurement type and normalization status
    if value_class is not None:
        measurement_type = _infer_measurement_type(value_class)
        normalization_status = _infer_normalization_status(value_class)
    else:
        from pharmomics.omics.enums import MeasurementType, NormalizationStatus

        measurement_type = MeasurementType.UNKNOWN
        normalization_status = NormalizationStatus.UNKNOWN

    # Build metadata dicts
    feature_meta = _build_feature_metadata(
        expr_result.gene_ids, gene_inspection,
    )
    sample_meta = _build_sample_metadata(meta_result)

    # Build provenance
    provenance: list[ProvenanceRecord] = []
    if source_id:
        provenance.append(
            ProvenanceRecord(
                source_id=source_id,
                source_file=expr_result.original_filename,
                sha256=sha256,
                ingested_at=created_at,
                software_version=software_version,
            )
        )

    return OmicsMatrix(
        matrix_id=matrix_id,
        schema_version="1.0.0",
        modality="transcriptomics",
        feature_type="gene",
        measurement_type=measurement_type.value,
        normalization_status=normalization_status.value,
        n_features=expr_result.n_genes,
        n_samples=expr_result.n_samples,
        feature_ids=expr_result.gene_ids,
        sample_ids=expr_result.sample_ids,
        dataframe=expr_result.dataframe,
        feature_metadata=feature_meta,
        sample_metadata=sample_meta,
        provenance=provenance,
        created_at=created_at,
    )
