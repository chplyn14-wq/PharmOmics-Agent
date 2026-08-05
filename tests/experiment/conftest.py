"""Shared fixtures for experiment domain tests."""

from __future__ import annotations

import pytest

from pharmomics.experiment.enums import (
    CovariateRole,
    CovariateValueType,
    FactorType,
    GroupRole,
)
from pharmomics.experiment.schemas import (
    Contrast,
    CovariateDefinition,
    DesignSample,
    ExperimentalFactor,
    ExperimentalGroup,
    ExperimentDesign,
    PairingDefinition,
    Quantity,
    Treatment,
)

# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture
def minimal_valid_design() -> ExperimentDesign:
    """Smallest valid ExperimentDesign: 1 sample, 1 group, 1 factor."""
    return ExperimentDesign(
        schema_version="1.0.0",
        experiment_id="exp-minimal",
        samples=[
            DesignSample(
                sample_id="s1",
                group_id="grp1",
                factor_values={"drug": "A"},
            ),
        ],
        groups=[
            ExperimentalGroup(
                group_id="grp1",
                label="Group 1",
                role=GroupRole.TREATMENT,
            ),
        ],
        factors=[
            ExperimentalFactor(
                factor_id="drug",
                factor_type=FactorType.CATEGORICAL,
                levels=["A", "B"],
            ),
        ],
    )


@pytest.fixture
def simple_drug_control_design() -> ExperimentDesign:
    """2 groups (treatment/control), 3 samples each, 1 contrast, 1 covariate."""
    return ExperimentDesign(
        schema_version="1.0.0",
        experiment_id="exp-drug",
        samples=[
            DesignSample(
                sample_id="ctrl_1",
                group_id="control",
                factor_values={"drug": "vehicle"},
            ),
            DesignSample(
                sample_id="ctrl_2",
                group_id="control",
                factor_values={"drug": "vehicle"},
            ),
            DesignSample(
                sample_id="ctrl_3",
                group_id="control",
                factor_values={"drug": "vehicle"},
            ),
            DesignSample(
                sample_id="trt_1",
                group_id="treatment",
                factor_values={"drug": "compound_x"},
                treatment=Treatment(
                    compound="compound_x",
                    dose=Quantity(value=500, unit="nM"),
                    duration="24h",
                ),
            ),
            DesignSample(
                sample_id="trt_2",
                group_id="treatment",
                factor_values={"drug": "compound_x"},
                treatment=Treatment(
                    compound="compound_x",
                    dose=Quantity(value=500, unit="nM"),
                    duration="24h",
                ),
            ),
            DesignSample(
                sample_id="trt_3",
                group_id="treatment",
                factor_values={"drug": "compound_x"},
                treatment=Treatment(
                    compound="compound_x",
                    dose=Quantity(value=500, unit="nM"),
                    duration="24h",
                ),
            ),
        ],
        groups=[
            ExperimentalGroup(
                group_id="control",
                label="Vehicle control",
                role=GroupRole.CONTROL,
            ),
            ExperimentalGroup(
                group_id="treatment",
                label="Compound X",
                role=GroupRole.TREATMENT,
            ),
        ],
        factors=[
            ExperimentalFactor(
                factor_id="drug",
                factor_type=FactorType.CATEGORICAL,
                levels=["vehicle", "compound_x"],
            ),
        ],
        contrasts=[
            Contrast(
                contrast_id="trt_vs_ctrl",
                comparison_group_id="treatment",
                reference_group_id="control",
                description="Compound X vs vehicle",
            ),
        ],
        covariates=[
            CovariateDefinition(
                covariate_id="batch",
                role=CovariateRole.BATCH,
                value_type=CovariateValueType.CATEGORICAL,
            ),
        ],
    )


@pytest.fixture
def paired_before_after_design() -> ExperimentDesign:
    """3 pairs, each with before/after samples."""
    return ExperimentDesign(
        schema_version="1.0.0",
        experiment_id="exp-paired",
        samples=[
            DesignSample(
                sample_id="subj1_before",
                group_id="baseline",
                factor_values={"time": "before"},
                pair_id="pair1",
            ),
            DesignSample(
                sample_id="subj1_after",
                group_id="post",
                factor_values={"time": "after"},
                pair_id="pair1",
            ),
            DesignSample(
                sample_id="subj2_before",
                group_id="baseline",
                factor_values={"time": "before"},
                pair_id="pair2",
            ),
            DesignSample(
                sample_id="subj2_after",
                group_id="post",
                factor_values={"time": "after"},
                pair_id="pair2",
            ),
            DesignSample(
                sample_id="subj3_before",
                group_id="baseline",
                factor_values={"time": "before"},
                pair_id="pair3",
            ),
            DesignSample(
                sample_id="subj3_after",
                group_id="post",
                factor_values={"time": "after"},
                pair_id="pair3",
            ),
        ],
        groups=[
            ExperimentalGroup(
                group_id="baseline",
                label="Baseline",
                role=GroupRole.OBSERVATIONAL,
            ),
            ExperimentalGroup(
                group_id="post",
                label="Post treatment",
                role=GroupRole.TREATMENT,
            ),
        ],
        factors=[
            ExperimentalFactor(
                factor_id="time",
                factor_type=FactorType.CATEGORICAL,
                levels=["before", "after"],
            ),
        ],
        pairing=PairingDefinition(
            pairing_type="before_after",
            description="Before/after on same subject",
        ),
    )


@pytest.fixture
def multi_factor_design() -> ExperimentDesign:
    """2 factors (drug categorical, time continuous), 4 groups."""
    return ExperimentDesign(
        schema_version="1.0.0",
        experiment_id="exp-multifact",
        samples=[
            DesignSample(
                sample_id="s1",
                group_id="g1",
                factor_values={"drug": "A", "time": 0},
            ),
            DesignSample(
                sample_id="s2",
                group_id="g2",
                factor_values={"drug": "A", "time": 24},
            ),
            DesignSample(
                sample_id="s3",
                group_id="g3",
                factor_values={"drug": "B", "time": 0},
            ),
            DesignSample(
                sample_id="s4",
                group_id="g4",
                factor_values={"drug": "B", "time": 24},
            ),
        ],
        groups=[
            ExperimentalGroup(
                group_id="g1",
                label="Drug A, T=0",
                role=GroupRole.CONTROL,
            ),
            ExperimentalGroup(
                group_id="g2",
                label="Drug A, T=24",
                role=GroupRole.TREATMENT,
            ),
            ExperimentalGroup(
                group_id="g3",
                label="Drug B, T=0",
                role=GroupRole.CONTROL,
            ),
            ExperimentalGroup(
                group_id="g4",
                label="Drug B, T=24",
                role=GroupRole.TREATMENT,
            ),
        ],
        factors=[
            ExperimentalFactor(
                factor_id="drug",
                factor_type=FactorType.CATEGORICAL,
                levels=["A", "B"],
            ),
            ExperimentalFactor(
                factor_id="time",
                factor_type=FactorType.CONTINUOUS,
                description="Time in hours",
            ),
        ],
        covariates=[
            CovariateDefinition(
                covariate_id="age",
                role=CovariateRole.CLINICAL,
                value_type=CovariateValueType.CONTINUOUS,
            ),
        ],
    )


@pytest.fixture
def batch_adjusted_design() -> ExperimentDesign:
    """2 batches, samples distributed across batches."""
    return ExperimentDesign(
        schema_version="1.0.0",
        experiment_id="exp-batch",
        samples=[
            DesignSample(
                sample_id="s1",
                group_id="grp1",
                covariate_values={"batch": "B1"},
            ),
            DesignSample(
                sample_id="s2",
                group_id="grp1",
                covariate_values={"batch": "B1"},
            ),
            DesignSample(
                sample_id="s3",
                group_id="grp1",
                covariate_values={"batch": "B2"},
            ),
            DesignSample(
                sample_id="s4",
                group_id="grp1",
                covariate_values={"batch": "B2"},
            ),
        ],
        groups=[
            ExperimentalGroup(
                group_id="grp1",
                label="All samples",
                role=GroupRole.OBSERVATIONAL,
            ),
        ],
        covariates=[
            CovariateDefinition(
                covariate_id="batch",
                role=CovariateRole.BATCH,
                value_type=CovariateValueType.CATEGORICAL,
            ),
        ],
    )


@pytest.fixture
def invalid_reference_design() -> ExperimentDesign:
    """Design with broken references (unknown group_id, factor, covariate)."""
    return ExperimentDesign(
        schema_version="1.0.0",
        experiment_id="exp-broken",
        samples=[
            DesignSample(
                sample_id="s1",
                group_id="nonexistent_group",
                factor_values={"unknown_factor": "x"},
                covariate_values={"unknown_cov": "y"},
            ),
        ],
        groups=[
            ExperimentalGroup(
                group_id="grp1",
                label="Real group",
                role=GroupRole.CONTROL,
            ),
        ],
        factors=[
            ExperimentalFactor(
                factor_id="known_factor",
                factor_type=FactorType.CATEGORICAL,
                levels=["x"],
            ),
        ],
        covariates=[
            CovariateDefinition(
                covariate_id="known_cov",
                role=CovariateRole.CLINICAL,
                value_type=CovariateValueType.CATEGORICAL,
            ),
        ],
        contrasts=[
            Contrast(
                contrast_id="bad_contrast",
                comparison_group_id="missing_a",
                reference_group_id="missing_b",
            ),
        ],
    )
