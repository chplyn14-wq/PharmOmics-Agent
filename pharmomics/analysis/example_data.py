"""Deterministic example data for PharmOmics development and testing.

Provides ``make_demo_inputs()``, a single convenience function that returns
a small, self-contained ``(OmicsMatrix, ExperimentDesign, AnalysisSpecification)``
triple.  All values are hard-coded — no randomness, no file I/O, no network
access.

The example models a simple two-condition experiment:

- 6 genes × 6 samples (3 control + 3 treatment)
- Treatment values are set to approximately 2× control (log₂FC ≈ 1)
- One contrast: ``treated_vs_control``
- One analysis specification: ``differential_analysis``
"""

from __future__ import annotations

import pandas as pd

from pharmomics.analysis.schemas import AnalysisSpecification
from pharmomics.experiment.enums import (
    FactorType,
    GroupRole,
)
from pharmomics.experiment.schemas import (
    Contrast,
    DesignSample,
    ExperimentalFactor,
    ExperimentalGroup,
    ExperimentDesign,
)
from pharmomics.omics.enums import (
    MeasurementType,
    Modality,
    NormalizationStatus,
)
from pharmomics.omics.schemas import OmicsMatrix

# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------


def make_demo_inputs() -> tuple[OmicsMatrix, ExperimentDesign, AnalysisSpecification]:
    """Return a deterministic, small-scale example triple.

    Returns
    -------
    OmicsMatrix
        6 genes × 6 samples expression matrix (3 control + 3 treatment).
        Treatment values are ~2× control values (log₂FC ≈ 1).
    ExperimentDesign
        Two groups, one factor (``condition``), one contrast
        (``treated_vs_control``).
    AnalysisSpecification
        ``analysis_type="differential_analysis"`` with
        ``contrast_references=["treated_vs_control"]``.

    Notes
    -----
    Pure function — no I/O, no randomness, no external dependencies.
    """
    return (
        _build_demo_omics(),
        _build_demo_design(),
        _build_demo_spec(),
    )


# ---------------------------------------------------------------------------
# Private helpers
# ---------------------------------------------------------------------------

# Control baseline expression values for 6 genes (arbitrary counts)
_CONTROL_VALUES = [100.0, 200.0, 50.0, 500.0, 300.0, 150.0]
# Treatment values ≈ 2× control (log₂FC ≈ 1 for all genes)
_TREATMENT_VALUES = [200.0, 400.0, 100.0, 1000.0, 600.0, 300.0]

_DEMO_GENES = [
    "EGFR", "ERBB2", "TP53", "MYC", "KRAS", "PTEN",
]
_DEMO_SAMPLES = [
    "ctrl_1", "ctrl_2", "ctrl_3",
    "trt_1", "trt_2", "trt_3",
]


def _build_demo_omics() -> OmicsMatrix:
    """Construct a 6-gene × 6-sample OmicsMatrix with deterministic values."""
    df = pd.DataFrame(
        {"feature_id": _DEMO_GENES}
        | {
            "ctrl_1": [100.0, 200.0, 50.0, 500.0, 300.0, 150.0],
            "ctrl_2": [100.0, 200.0, 50.0, 500.0, 300.0, 150.0],
            "ctrl_3": [100.0, 200.0, 50.0, 500.0, 300.0, 150.0],
            "trt_1": [200.0, 400.0, 100.0, 1000.0, 600.0, 300.0],
            "trt_2": [200.0, 400.0, 100.0, 1000.0, 600.0, 300.0],
            "trt_3": [200.0, 400.0, 100.0, 1000.0, 600.0, 300.0],
        },
    )
    return OmicsMatrix(
        matrix_id="mx-demo",
        schema_version="1.0.0",
        modality=Modality.TRANSCRIPTOMICS,
        feature_type="gene",
        measurement_type=MeasurementType.ESTIMATED_COUNTS,
        normalization_status=NormalizationStatus.RAW,
        n_features=len(_DEMO_GENES),
        n_samples=len(_DEMO_SAMPLES),
        feature_ids=list(_DEMO_GENES),
        sample_ids=list(_DEMO_SAMPLES),
        dataframe=df,
        created_at="2026-01-01T00:00:00Z",
    )


def _build_demo_design() -> ExperimentDesign:
    """Construct an ExperimentDesign with two groups and one contrast."""
    return ExperimentDesign(
        experiment_id="exp-demo",
        description="Demo: control vs treatment differential expression",
        samples=[
            DesignSample(
                sample_id="ctrl_1", group_id="ctrl",
                factor_values={"condition": "control"},
            ),
            DesignSample(
                sample_id="ctrl_2", group_id="ctrl",
                factor_values={"condition": "control"},
            ),
            DesignSample(
                sample_id="ctrl_3", group_id="ctrl",
                factor_values={"condition": "control"},
            ),
            DesignSample(
                sample_id="trt_1", group_id="trt",
                factor_values={"condition": "treatment"},
            ),
            DesignSample(
                sample_id="trt_2", group_id="trt",
                factor_values={"condition": "treatment"},
            ),
            DesignSample(
                sample_id="trt_3", group_id="trt",
                factor_values={"condition": "treatment"},
            ),
        ],
        groups=[
            ExperimentalGroup(
                group_id="ctrl", label="Control",
                role=GroupRole.CONTROL,
            ),
            ExperimentalGroup(
                group_id="trt", label="Treatment",
                role=GroupRole.TREATMENT,
            ),
        ],
        factors=[
            ExperimentalFactor(
                factor_id="condition",
                factor_type=FactorType.CATEGORICAL,
                description="Treatment condition",
                levels=["control", "treatment"],
            ),
        ],
        contrasts=[
            Contrast(
                contrast_id="treated_vs_control",
                comparison_group_id="trt",
                reference_group_id="ctrl",
                description="Treatment vs control differential expression",
            ),
        ],
    )


def _build_demo_spec() -> AnalysisSpecification:
    """Construct an AnalysisSpecification for differential analysis."""
    return AnalysisSpecification(
        analysis_type="differential_analysis",
        factor_references=["condition"],
        contrast_references=["treated_vs_control"],
    )
