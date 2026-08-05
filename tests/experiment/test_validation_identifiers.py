"""Tests for identifier uniqueness validation (I-01 … I-12)."""

from __future__ import annotations

from pharmomics.experiment.enums import FactorType, GroupRole
from pharmomics.experiment.schemas import (
    Contrast,
    CovariateDefinition,
    DesignSample,
    ExperimentalFactor,
    ExperimentalGroup,
    ExperimentDesign,
)
from pharmomics.experiment.validation import validate


class TestIdentifierEmpty:
    """I-01, I-02, I-08, I-09, I-10, I-11, I-12: empty/whitespace IDs."""

    def test_empty_sample_id(self) -> None:
        d = ExperimentDesign(
            experiment_id="exp",
            samples=[DesignSample(sample_id="", group_id="g1")],
            groups=[
                ExperimentalGroup(group_id="g1", label="G", role=GroupRole.CONTROL)
            ],
        )
        errs = validate(d)
        assert any("Empty sample_id" in e for e in errs)

    def test_whitespace_only_sample_id(self) -> None:
        d = ExperimentDesign(
            experiment_id="exp",
            samples=[DesignSample(sample_id="   ", group_id="g1")],
            groups=[
                ExperimentalGroup(group_id="g1", label="G", role=GroupRole.CONTROL)
            ],
        )
        errs = validate(d)
        assert any("Whitespace-only sample_id" in e for e in errs)

    def test_empty_group_id(self) -> None:
        d = ExperimentDesign(
            experiment_id="exp",
            groups=[ExperimentalGroup(group_id="", label="G", role=GroupRole.CONTROL)],
        )
        errs = validate(d)
        assert any("Empty group_id" in e for e in errs)

    def test_empty_factor_id(self) -> None:
        d = ExperimentDesign(
            experiment_id="exp",
            factors=[
                ExperimentalFactor(
                    factor_id="",
                    factor_type=FactorType.CATEGORICAL,
                    levels=["A"],
                )
            ],
        )
        errs = validate(d)
        assert any("Empty factor_id" in e for e in errs)

    def test_empty_covariate_id(self) -> None:
        from pharmomics.experiment.enums import CovariateRole, CovariateValueType

        d = ExperimentDesign(
            experiment_id="exp",
            covariates=[
                CovariateDefinition(
                    covariate_id="  ",
                    role=CovariateRole.BATCH,
                    value_type=CovariateValueType.CATEGORICAL,
                )
            ],
        )
        errs = validate(d)
        assert any("Empty covariate_id" in e for e in errs)

    def test_empty_contrast_id(self) -> None:
        d = ExperimentDesign(
            experiment_id="exp",
            contrasts=[
                Contrast(
                    contrast_id="",
                    comparison_group_id="A",
                    reference_group_id="B",
                )
            ],
            groups=[
                ExperimentalGroup(group_id="A", label="A", role=GroupRole.TREATMENT),
                ExperimentalGroup(group_id="B", label="B", role=GroupRole.CONTROL),
            ],
        )
        errs = validate(d)
        assert any("Empty contrast_id" in e for e in errs)

    def test_empty_experiment_id(self) -> None:
        d = ExperimentDesign(experiment_id="  ")
        errs = validate(d)
        assert any("Empty experiment_id" in e for e in errs)


class TestIdentifierDuplicates:
    """I-03, I-04, I-05, I-06, I-07: duplicate IDs."""

    def test_duplicate_sample_id(self) -> None:
        d = ExperimentDesign(
            experiment_id="exp",
            samples=[
                DesignSample(sample_id="s1", group_id="g1"),
                DesignSample(sample_id="s1", group_id="g1"),
            ],
            groups=[
                ExperimentalGroup(group_id="g1", label="G", role=GroupRole.CONTROL)
            ],
        )
        errs = validate(d)
        assert any("Duplicate sample_id" in e for e in errs)

    def test_duplicate_group_id(self) -> None:
        d = ExperimentDesign(
            experiment_id="exp",
            groups=[
                ExperimentalGroup(group_id="g1", label="A", role=GroupRole.CONTROL),
                ExperimentalGroup(group_id="g1", label="B", role=GroupRole.TREATMENT),
            ],
        )
        errs = validate(d)
        assert any("Duplicate group_id" in e for e in errs)

    def test_duplicate_factor_id(self) -> None:
        d = ExperimentDesign(
            experiment_id="exp",
            factors=[
                ExperimentalFactor(
                    factor_id="f", factor_type=FactorType.CATEGORICAL, levels=["A"]
                ),
                ExperimentalFactor(factor_id="f", factor_type=FactorType.CONTINUOUS),
            ],
        )
        errs = validate(d)
        assert any("Duplicate factor_id" in e for e in errs)

    def test_duplicate_covariate_id(self) -> None:
        from pharmomics.experiment.enums import CovariateRole, CovariateValueType

        d = ExperimentDesign(
            experiment_id="exp",
            covariates=[
                CovariateDefinition(
                    covariate_id="c",
                    role=CovariateRole.BATCH,
                    value_type=CovariateValueType.CATEGORICAL,
                ),
                CovariateDefinition(
                    covariate_id="c",
                    role=CovariateRole.CLINICAL,
                    value_type=CovariateValueType.CONTINUOUS,
                ),
            ],
        )
        errs = validate(d)
        assert any("Duplicate covariate_id" in e for e in errs)

    def test_duplicate_contrast_id(self) -> None:
        d = ExperimentDesign(
            experiment_id="exp",
            contrasts=[
                Contrast(
                    contrast_id="c", comparison_group_id="A", reference_group_id="B"
                ),
                Contrast(
                    contrast_id="c", comparison_group_id="C", reference_group_id="D"
                ),
            ],
            groups=[
                ExperimentalGroup(group_id="A", label="A", role=GroupRole.TREATMENT),
                ExperimentalGroup(group_id="B", label="B", role=GroupRole.CONTROL),
                ExperimentalGroup(group_id="C", label="C", role=GroupRole.TREATMENT),
                ExperimentalGroup(group_id="D", label="D", role=GroupRole.CONTROL),
            ],
        )
        errs = validate(d)
        assert any("Duplicate contrast_id" in e for e in errs)


class TestIdentifierCaseSensitive:
    """IDs are case-sensitive: 'A' ≠ 'a'."""

    def test_case_sensitive_sample_ids(self) -> None:
        d = ExperimentDesign(
            experiment_id="exp",
            samples=[
                DesignSample(sample_id="S1", group_id="g1"),
                DesignSample(sample_id="s1", group_id="g1"),
            ],
            groups=[
                ExperimentalGroup(group_id="g1", label="G", role=GroupRole.CONTROL)
            ],
        )
        errs = validate(d)
        assert not any("Duplicate sample_id" in e for e in errs)
