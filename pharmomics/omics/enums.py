"""Controlled vocabularies for the PharmOmics omics core.

These enums describe *what* an omics matrix contains (modality, feature
type, measurement type, normalisation state) in a way that is independent
of any particular file format or quantification pipeline.
"""

from __future__ import annotations

from enum import StrEnum

# ---------------------------------------------------------------------------
# Modality — which omics technology produced the data?
# ---------------------------------------------------------------------------


class Modality(StrEnum):
    """Omics data modality."""

    TRANSCRIPTOMICS = "transcriptomics"
    PROTEOMICS = "proteomics"
    METABOLOMICS = "metabolomics"
    EPIGENOMICS = "epigenomics"


# ---------------------------------------------------------------------------
# FeatureType — what does each feature row represent?
# ---------------------------------------------------------------------------


class FeatureType(StrEnum):
    """What a feature row represents."""

    GENE = "gene"
    TRANSCRIPT = "transcript"
    PROTEIN = "protein"
    METABOLITE = "metabolite"
    PEAK = "peak"
    UNKNOWN = "unknown"


# ---------------------------------------------------------------------------
# MeasurementType — what do the numeric values measure?
# ---------------------------------------------------------------------------


class MeasurementType(StrEnum):
    """What the numeric values in the matrix measure."""

    RAW_COUNTS = "raw_counts"
    ESTIMATED_COUNTS = "estimated_counts"
    TPM = "tpm"
    FPKM = "fpkm"
    CPM = "cpm"
    INTENSITY = "intensity"
    ABUNDANCE = "abundance"
    LOG2FC = "log2fc"
    ZSCORE = "zscore"
    BETA = "beta"
    UNKNOWN = "unknown"


# ---------------------------------------------------------------------------
# NormalizationStatus — how far has the data been normalised?
# ---------------------------------------------------------------------------


class NormalizationStatus(StrEnum):
    """Normalisation state of the matrix values."""

    RAW = "raw"
    WITHIN_SAMPLE = "within_sample"
    BETWEEN_SAMPLE = "between_sample"
    TRANSFORMED = "transformed"
    UNKNOWN = "unknown"
