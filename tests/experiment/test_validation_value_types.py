"""Tests for value type consistency validation (V-01 … V-05)."""

from __future__ import annotations

from pharmomics.experiment.enums import (
    CovariateRole,
    CovariateValueType,
    FactorType,
    GroupRole,
)
from pharmomics.experiment.schemas import (
    CovariateDefinition,
    DesignSample,
    ExperimentalFactor,
    ExperimentalGroup,
    ExperimentDesign,
)
from pharmomics.experiment.validation import validate


class TestCategoricalFactorLevels:
    """V-01: categorical factor value must be in levels."""

    def test_valid_level(self) -> None:
        d = ExperimentDesign(
            experiment_id="exp",
            samples=[
                DesignSample(
                    sample_id="s1",
                    group_id="g1",
                    factor_values={"drug": "A"},
                )
            ],
            groups=[
                ExperimentalGroup(group_id="g1", label="G", role=GroupRole.CONTROL)
            ],
            factors=[
                ExperimentalFactor(
                    factor_id="drug",
                    factor_type=FactorType.CATEGORICAL,
                    levels=["A", "B"],
                )
            ],
        )
        errs = validate(d)
        assert not any("not in levels" in e for e in errs)

    def test_invalid_level(self) -> None:
        d = ExperimentDesign(
            experiment_id="exp",
            samples=[
                DesignSample(
                    sample_id="s1",
                    group_id="g1",
                    factor_values={"drug": "X"},
                )
            ],
            groups=[
                ExperimentalGroup(group_id="g1", label="G", role=GroupRole.CONTROL)
            ],
            factors=[
                ExperimentalFactor(
                    factor_id="drug",
                    factor_type=FactorType.CATEGORICAL,
                    levels=["A", "B"],
                )
            ],
        )
        errs = validate(d)
        assert any("not in levels" in e for e in errs)

    def test_non_string_for_categorical(self) -> None:
        d = ExperimentDesign(
            experiment_id="exp",
            samples=[
                DesignSample(
                    sample_id="s1",
                    group_id="g1",
                    factor_values={"drug": 1},
                )
            ],
            groups=[
                ExperimentalGroup(group_id="g1", label="G", role=GroupRole.CONTROL)
            ],
            factors=[
                ExperimentalFactor(
                    factor_id="drug",
                    factor_type=FactorType.CATEGORICAL,
                    levels=["A", "B"],
                )
            ],
        )
        errs = validate(d)
        assert any("is not a string for type categorical" in e for e in errs)


class TestContinuousFactor:
    """V-02: continuous factor value must be numeric."""

    def test_valid_numeric(self) -> None:
        d = ExperimentDesign(
            experiment_id="exp",
            samples=[
                DesignSample(
                    sample_id="s1",
                    group_id="g1",
                    factor_values={"time": 24.0},
                )
            ],
            groups=[
                ExperimentalGroup(group_id="g1", label="G", role=GroupRole.CONTROL)
            ],
            factors=[
                ExperimentalFactor(
                    factor_id="time",
                    factor_type=FactorType.CONTINUOUS,
                )
            ],
        )
        errs = validate(d)
        assert not any("is not numeric" in e for e in errs)

    def test_invalid_string(self) -> None:
        d = ExperimentDesign(
            experiment_id="exp",
            samples=[
                DesignSample(
                    sample_id="s1",
                    group_id="g1",
                    factor_values={"time": "fast"},
                )
            ],
            groups=[
                ExperimentalGroup(group_id="g1", label="G", role=GroupRole.CONTROL)
            ],
            factors=[
                ExperimentalFactor(
                    factor_id="time",
                    factor_type=FactorType.CONTINUOUS,
                )
            ],
        )
        errs = validate(d)
        assert any("is not numeric for type continuous" in e for e in errs)

    def test_bool_rejected(self) -> None:
        d = ExperimentDesign(
            experiment_id="exp",
            samples=[
                DesignSample(
                    sample_id="s1",
                    group_id="g1",
                    factor_values={"time": True},
                )
            ],
            groups=[
                ExperimentalGroup(group_id="g1", label="G", role=GroupRole.CONTROL)
            ],
            factors=[
                ExperimentalFactor(
                    factor_id="time",
                    factor_type=FactorType.CONTINUOUS,
                )
            ],
        )
        errs = validate(d)
        assert any("is not numeric for type continuous" in e for e in errs)


class TestCovariateTypes:
    """V-03, V-04: covariate value type checks."""

    def test_categorical_covariate_string(self) -> None:
        d = ExperimentDesign(
            experiment_id="exp",
            samples=[
                DesignSample(
                    sample_id="s1",
                    group_id="g1",
                    covariate_values={"batch": "B1"},
                )
            ],
            groups=[
                ExperimentalGroup(group_id="g1", label="G", role=GroupRole.CONTROL)
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
        assert not any("is not categorical" in e for e in errs)

    def test_continuous_covariate_numeric(self) -> None:
        d = ExperimentDesign(
            experiment_id="exp",
            samples=[
                DesignSample(
                    sample_id="s1",
                    group_id="g1",
                    covariate_values={"age": 45},
                )
            ],
            groups=[
                ExperimentalGroup(group_id="g1", label="G", role=GroupRole.CONTROL)
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
        assert not any("is not numeric" in e for e in errs)

    def test_continuous_covariate_string_rejected(self) -> None:
        d = ExperimentDesign(
            experiment_id="exp",
            samples=[
                DesignSample(
                    sample_id="s1",
                    group_id="g1",
                    covariate_values={"age": "old"},
                )
            ],
            groups=[
                ExperimentalGroup(group_id="g1", label="G", role=GroupRole.CONTROL)
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
        assert any("is not numeric for type continuous" in e for e in errs)

    def test_integer_accepted_for_continuous(self) -> None:
        d = ExperimentDesign(
            experiment_id="exp",
            samples=[
                DesignSample(
                    sample_id="s1",
                    group_id="g1",
                    covariate_values={"age": 45},
                )
            ],
            groups=[
                ExperimentalGroup(group_id="g1", label="G", role=GroupRole.CONTROL)
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
        assert not any("is not numeric" in e for e in errs)


class TestCategoricalFactorMissingLevels:
    """V-05: categorical factor without levels is a violation."""

    def test_categorical_no_levels(self) -> None:
        d = ExperimentDesign(
            experiment_id="exp",
            factors=[
                ExperimentalFactor(
                    factor_id="drug",
                    factor_type=FactorType.CATEGORICAL,
                )
            ],
        )
        errs = validate(d)
        assert any("categorical but has no levels" in e for e in errs)
