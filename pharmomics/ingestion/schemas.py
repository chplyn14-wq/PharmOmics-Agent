"""Ingestion-domain Pydantic schemas for PharmOmics.

Provides schemas for gene identifier mapping and ingested dataset metadata.
These are separate from run provenance (run_store) and agent/LLM logging
(agents.schemas).
"""

from __future__ import annotations

from pydantic import BaseModel, ConfigDict, Field


class GeneMappingRecord(BaseModel):
    """Mapping entry for a single gene identifier."""

    model_config = ConfigDict(frozen=True)

    gene_symbol: str
    ensembl_id: str | None = None
    entrez_id: str | None = None
    organism: str = "Homo sapiens"
    source: str = "HGNC"


class GeneMappingSet(BaseModel):
    """Collection of gene identifier mappings."""

    model_config = ConfigDict(frozen=True)

    mappings: list[GeneMappingRecord] = Field(default_factory=list)


class IngestedDataset(BaseModel):
    """Metadata for a successfully ingested expression dataset."""

    model_config = ConfigDict(frozen=True)

    gse_accession: str
    expression_file: str
    n_genes: int
    n_samples: int
    gene_id_type: str  # e.g. "HGNC gene symbols", "Ensembl IDs"
    value_type: str  # e.g. "estimated counts", "TPM", "raw counts"
    ingested_at: str  # ISO-8601 datetime string
    sha256: str
