"""Tests for reference integrity validation (R-01 … R-05)."""

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


class TestUnknownGroupReference:
    """R-01: sample → group."""

    def test_unknown_group_id(self) -> None:
        d = ExperimentDesign(
            experiment_id="exp",
            samples=[DesignSample(sample_id="s1", group_id="ghost")],
            groups=[
                ExperimentalGroup(group_id="g1", label="G", role=GroupRole.CONTROL)
            ],
        )
        errs = validate(d)
        assert any("unknown group_id" in e for e in errs)


class TestUnknownContrastGroups:
    """R-02, R-03: contrast → groups."""

    def test_unknown_comparison_group(self) -> None:
        d = ExperimentDesign(
            experiment_id="exp",
            contrasts=[
                Contrast(
                    contrast_id="c1",
                    comparison_group_id="missing",
                    reference_group_id="g1",
                )
            ],
            groups=[
                ExperimentalGroup(group_id="g1", label="G", role=GroupRole.CONTROL)
            ],
        )
        errs = validate(d)
        assert any("unknown comparison_group_id" in e for e in errs)

    def test_unknown_reference_group(self) -> None:
        d = ExperimentDesign(
            experiment_id="exp",
            contrasts=[
                Contrast(
                    contrast_id="c1",
                    comparison_group_id="g1",
                    reference_group_id="missing",
                )
            ],
            groups=[
                ExperimentalGroup(group_id="g1", label="G", role=GroupRole.CONTROL)
            ],
        )
        errs = validate(d)
        assert any("unknown reference_group_id" in e for e in errs)


class TestUnknownFactorKey:
    """R-04: sample.factor_values keys → factor_ids."""

    def test_unknown_factor_key(self) -> None:
        d = ExperimentDesign(
            experiment_id="exp",
            samples=[
                DesignSample(
                    sample_id="s1",
                    group_id="g1",
                    factor_values={"phantom": "x"},
                )
            ],
            groups=[
                ExperimentalGroup(group_id="g1", label="G", role=GroupRole.CONTROL)
            ],
            factors=[
                ExperimentalFactor(
                    factor_id="real",
                    factor_type=FactorType.CATEGORICAL,
                    levels=["x"],
                )
            ],
        )
        errs = validate(d)
        assert any(
            "factor_values key" in e and "not in defined factors" in e for e in errs
        )


class TestUnknownCovariateKey:
    """R-05: sample.covariate_values keys → covariate_ids."""

    def test_unknown_covariate_key(self) -> None:
        from pharmomics.experiment.enums import CovariateRole, CovariateValueType

        d = ExperimentDesign(
            experiment_id="exp",
            samples=[
                DesignSample(
                    sample_id="s1",
                    group_id="g1",
                    covariate_values={"phantom": "x"},
                )
            ],
            groups=[
                ExperimentalGroup(group_id="g1", label="G", role=GroupRole.CONTROL)
            ],
            covariates=[
                CovariateDefinition(
                    covariate_id="real",
                    role=CovariateRole.CLINICAL,
                    value_type=CovariateValueType.CATEGORICAL,
                )
            ],
        )
        errs = validate(d)
        assert any(
            "covariate_values key" in e and "not in defined covariates" in e
            for e in errs
        )


class TestValidReferences:
    """All references resolve — no violations."""

    def test_all_refs_valid(self, simple_drug_control_design: ExperimentDesign) -> None:
        errs = validate(simple_drug_control_design)
        assert not any("unknown" in e for e in errs)
        assert not any("not in defined" in e for e in errs)
