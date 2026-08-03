"""PharmOmics ingestion — file loading, validation, and provenance tracking.

Schemas for gene mapping and ingested dataset metadata are in
``pharmomics.ingestion.schemas``.
"""

from pharmomics.ingestion.loader import (
    CompressionType,
    GeneIdInspectionResult,
    GeneIdType,
    IngestionError,
    IngestionResult,
    MetadataLoadResult,
    ValueClassificationResult,
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

__all__ = [
    "CompressionType",
    "GeneIdInspectionResult",
    "GeneIdType",
    "IngestionError",
    "IngestionResult",
    "MetadataLoadResult",
    "ValueClassificationResult",
    "ValueType",
    "classify_expression_values",
    "count_replicates_per_condition",
    "ingest",
    "inspect_gene_identifiers",
    "load_expression_matrix",
    "load_sample_metadata",
    "validate_contrast",
    "write_ingestion_manifest",
]
