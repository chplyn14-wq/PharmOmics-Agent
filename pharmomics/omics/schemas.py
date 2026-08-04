"""Pydantic schemas for the PharmOmics omics core.

Defines the domain-layer data objects used throughout Phase 2:

- ``OmicsDescriptor`` — lightweight metadata-only summary of a matrix
  (serialisable to JSON without carrying numeric data).
- ``OmicsMatrix`` — the full domain object including the numeric
  dataframe, feature/sample metadata, and provenance.

Both models are frozen (``frozen=True``) so that once constructed their
contents cannot be mutated.  ``OmicsMatrix`` allows arbitrary types
so it can hold a ``pandas.DataFrame``.
"""

from __future__ import annotations

from typing import Any

import pandas as pd
from pydantic import BaseModel, ConfigDict, Field

# ---------------------------------------------------------------------------
# Sub-models
# ---------------------------------------------------------------------------


class FeatureMetadata(BaseModel):
    """Metadata for a single feature (gene, protein, metabolite, …)."""

    model_config = ConfigDict(frozen=True)

    feature_id: str
    annotations: dict[str, Any] = Field(default_factory=dict)
    """e.g. ``{"hgnc_symbol": "EGFR", "chromosome": "7"}``"""

    normalized_id: str | None = None
    """Stripped Ensembl version or mapped synonym."""


class SampleMetadata(BaseModel):
    """Metadata for a single sample."""

    model_config = ConfigDict(frozen=True)

    sample_id: str
    condition: str | None = None
    annotations: dict[str, Any] = Field(default_factory=dict)
    """e.g. ``{"cell_line": "PC9", "replicate": 1, "batch": "B1"}``"""


class Transformation(BaseModel):
    """Applied transformation on the matrix values."""

    model_config = ConfigDict(frozen=True)

    method: str
    """e.g. ``"log2"``, ``"log1p"``, ``"zscore"``."""

    parameters: dict[str, Any] = Field(default_factory=dict)
    """e.g. ``{"base": 2, "pseudocount": 1}``."""

    applied_by: str = ""
    """Who applied it: ``"user"`` | ``"pipeline"`` | ``""``."""


class ProvenanceRecord(BaseModel):
    """Source provenance for one input file that contributed to the matrix."""

    model_config = ConfigDict(frozen=True)

    source_id: str
    """e.g. GEO accession or internal identifier."""

    source_file: str
    """Original filename."""

    sha256: str
    """SHA-256 hex digest of the source file."""

    ingested_at: str
    """ISO-8601 timestamp."""

    software_version: str = ""


# ---------------------------------------------------------------------------
# OmicsDescriptor — metadata-only, always JSON-serialisable
# ---------------------------------------------------------------------------


class OmicsDescriptor(BaseModel):
    """Lightweight, JSON-serialisable summary of an omics matrix.

    Contains every field of ``OmicsMatrix`` *except* the numeric dataframe.
    Useful for cataloguing, discovery, and manifest entries where carrying
    the full data in memory is unnecessary.
    """

    model_config = ConfigDict(frozen=True)

    matrix_id: str
    schema_version: str = "1.0.0"

    # Import deferred to avoid circular import at module level.
    # These fields use the enums defined in the sibling ``enums`` module.
    modality: str
    feature_type: str
    measurement_type: str
    normalization_status: str

    n_features: int
    n_samples: int
    feature_ids: list[str]
    sample_ids: list[str]

    transformation: Transformation | None = None

    feature_metadata: dict[str, FeatureMetadata] = Field(default_factory=dict)
    sample_metadata: dict[str, SampleMetadata] = Field(default_factory=dict)
    provenance: list[ProvenanceRecord] = Field(default_factory=list)
    warnings: list[str] = Field(default_factory=list)

    created_at: str


# ---------------------------------------------------------------------------
# OmicsMatrix — full domain object with dataframe
# ---------------------------------------------------------------------------


class OmicsMatrix(BaseModel):
    """Core domain object: a multi-omics compatible data matrix.

    Holds the numeric dataframe (with features as rows, samples as columns),
    all associated metadata, and full provenance.  Constructed from
    ingestion results via the adapter (added separately).

    The dataframe is excluded from ``model_dump()`` so that the descriptor
    portion always serialises to JSON without carrying bulk numeric data.
    """

    model_config = ConfigDict(frozen=True, arbitrary_types_allowed=True)

    matrix_id: str
    schema_version: str = "1.0.0"

    modality: str
    feature_type: str
    measurement_type: str
    normalization_status: str

    n_features: int
    n_samples: int
    feature_ids: list[str]
    sample_ids: list[str]

    dataframe: pd.DataFrame
    """Feature-by-sample numeric matrix.
    Rows correspond to ``feature_ids`` (in order);
    columns correspond to ``sample_ids`` (in order).
    """

    transformation: Transformation | None = None

    feature_metadata: dict[str, FeatureMetadata] = Field(default_factory=dict)
    sample_metadata: dict[str, SampleMetadata] = Field(default_factory=dict)
    provenance: list[ProvenanceRecord] = Field(default_factory=list)
    warnings: list[str] = Field(default_factory=list)

    created_at: str

    def model_dump(self, **kwargs: Any) -> dict[str, Any]:
        """Return a JSON-serialisable dict, excluding the dataframe."""
        d = super().model_dump(exclude={"dataframe"}, **kwargs)
        return d

    @property
    def descriptor(self) -> OmicsDescriptor:
        """Return the metadata-only ``OmicsDescriptor`` for this matrix."""
        return OmicsDescriptor(
            matrix_id=self.matrix_id,
            schema_version=self.schema_version,
            modality=self.modality,
            feature_type=self.feature_type,
            measurement_type=self.measurement_type,
            normalization_status=self.normalization_status,
            n_features=self.n_features,
            n_samples=self.n_samples,
            feature_ids=self.feature_ids,
            sample_ids=self.sample_ids,
            transformation=self.transformation,
            feature_metadata=self.feature_metadata,
            sample_metadata=self.sample_metadata,
            provenance=self.provenance,
            warnings=self.warnings,
            created_at=self.created_at,
        )
