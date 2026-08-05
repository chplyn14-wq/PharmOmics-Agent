"""Tests for covariate special value validation (CV-01 … CV-02)."""

from __future__ import annotations

from pharmomics.experiment.enums import CovariateRole, CovariateValueType, GroupRole
from pharmomics.experiment.schemas import (
    CovariateDefinition,
    DesignSample,
    ExperimentalGroup,
    ExperimentDesign,
)
from pharmomics.experiment.validation import validate


class TestContinuousCovariateNaN:
    """CV-01: NaN rejected for continuous covariate."""

    def test_nan_rejected(self) -> None:
        d = ExperimentDesign(
            experiment_id="exp",
            groups=[
                ExperimentalGroup(group_id="g1", label="G", role=GroupRole.CONTROL)
            ],
            samples=[
                DesignSample(
                    sample_id="s1",
                    group_id="g1",
                    covariate_values={"age": float("nan")},
                )
            ],
            covariates=[
                CovariateDefinition(
                    covariate_id="age",
                    role=CovariateRole.CLINICAL,
                    value_type=CovariateValueType.CONTINUOUS,
                )
            ],
        )
        errs = validate(d)
        assert any("NaN" in e for e in errs)


class TestContinuousCovariateInfinity:
    """CV-02: infinity rejected for continuous covariate."""

    def test_inf_rejected(self) -> None:
        d = ExperimentDesign(
            experiment_id="exp",
            groups=[
                ExperimentalGroup(group_id="g1", label="G", role=GroupRole.CONTROL)
            ],
            samples=[
                DesignSample(
                    sample_id="s1",
                    group_id="g1",
                    covariate_values={"age": float("inf")},
                )
            ],
            covariates=[
                CovariateDefinition(
                    covariate_id="age",
                    role=CovariateRole.CLINICAL,
                    value_type=CovariateValueType.CONTINUOUS,
                )
            ],
        )
        errs = validate(d)
        assert any("infinity" in e for e in errs)


class TestBatchAsCovariate:
    """Batch is represented as CovariateDefinition(role="batch")."""

    def test_batch_as_covariate(self) -> None:
        d = ExperimentDesign(
            experiment_id="exp",
            groups=[
                ExperimentalGroup(group_id="g1", label="G", role=GroupRole.CONTROL)
            ],
            samples=[
                DesignSample(
                    sample_id="s1",
                    group_id="g1",
                    covariate_values={"batch": "B1"},
                )
            ],
            covariates=[
                CovariateDefinition(
                    covariate_id="batch",
                    role=CovariateRole.BATCH,
                    value_type=CovariateValueType.CATEGORICAL,
                )
            ],
        )
        errs = validate(d)
        assert not any("batch" in e.lower() for e in errs)

    def test_batch_categorical_accepts_string(self) -> None:
        d = ExperimentDesign(
            experiment_id="exp",
            groups=[
                ExperimentalGroup(group_id="g1", label="G", role=GroupRole.CONTROL)
            ],
            samples=[
                DesignSample(
                    sample_id="s1", group_id="g1", covariate_values={"batch": "B1"}
                ),
                DesignSample(
                    sample_id="s2", group_id="g1", covariate_values={"batch": "B2"}
                ),
            ],
            covariates=[
                CovariateDefinition(
                    covariate_id="batch",
                    role=CovariateRole.BATCH,
                    value_type=CovariateValueType.CATEGORICAL,
                )
            ],
        )
        errs = validate(d)
        assert not any("covariate" in e.lower() for e in errs)
