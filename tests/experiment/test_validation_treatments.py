"""Tests for treatment validation (T-01 … T-07)."""

from __future__ import annotations

from pharmomics.experiment.enums import GroupRole
from pharmomics.experiment.schemas import (
    DesignSample,
    ExperimentalGroup,
    ExperimentDesign,
    Quantity,
    Treatment,
)
from pharmomics.experiment.validation import validate


class TestEmptyCompound:
    """T-01: empty compound string."""

    def test_empty_compound_rejected(self) -> None:
        d = ExperimentDesign(
            experiment_id="exp",
            groups=[
                ExperimentalGroup(group_id="g1", label="G", role=GroupRole.CONTROL)
            ],
            samples=[
                DesignSample(
                    sample_id="s1",
                    group_id="g1",
                    treatment=Treatment(compound=""),
                )
            ],
        )
        errs = validate(d)
        assert any("empty compound" in e for e in errs)

    def test_whitespace_compound_rejected(self) -> None:
        d = ExperimentDesign(
            experiment_id="exp",
            groups=[
                ExperimentalGroup(group_id="g1", label="G", role=GroupRole.CONTROL)
            ],
            samples=[
                DesignSample(
                    sample_id="s1",
                    group_id="g1",
                    treatment=Treatment(compound="   "),
                )
            ],
        )
        errs = validate(d)
        assert any("empty compound" in e for e in errs)


class TestDoseValue:
    """T-02, T-03, T-04, T-05: dose value checks."""

    def test_negative_dose_rejected(self) -> None:
        d = ExperimentDesign(
            experiment_id="exp",
            groups=[
                ExperimentalGroup(group_id="g1", label="G", role=GroupRole.CONTROL)
            ],
            samples=[
                DesignSample(
                    sample_id="s1",
                    group_id="g1",
                    treatment=Treatment(
                        compound="X", dose=Quantity(value=-5.0, unit="nM")
                    ),
                )
            ],
        )
        errs = validate(d)
        assert any("is negative" in e for e in errs)

    def test_zero_dose_rejected(self) -> None:
        d = ExperimentDesign(
            experiment_id="exp",
            groups=[
                ExperimentalGroup(group_id="g1", label="G", role=GroupRole.CONTROL)
            ],
            samples=[
                DesignSample(
                    sample_id="s1",
                    group_id="g1",
                    treatment=Treatment(
                        compound="X", dose=Quantity(value=0, unit="nM")
                    ),
                )
            ],
        )
        errs = validate(d)
        assert any("dose value is zero" in e for e in errs)

    def test_nan_dose_rejected(self) -> None:
        d = ExperimentDesign(
            experiment_id="exp",
            groups=[
                ExperimentalGroup(group_id="g1", label="G", role=GroupRole.CONTROL)
            ],
            samples=[
                DesignSample(
                    sample_id="s1",
                    group_id="g1",
                    treatment=Treatment(
                        compound="X", dose=Quantity(value=float("nan"), unit="nM")
                    ),
                )
            ],
        )
        errs = validate(d)
        assert any("dose is NaN" in e for e in errs)

    def test_inf_dose_rejected(self) -> None:
        d = ExperimentDesign(
            experiment_id="exp",
            groups=[
                ExperimentalGroup(group_id="g1", label="G", role=GroupRole.CONTROL)
            ],
            samples=[
                DesignSample(
                    sample_id="s1",
                    group_id="g1",
                    treatment=Treatment(
                        compound="X", dose=Quantity(value=float("inf"), unit="nM")
                    ),
                )
            ],
        )
        errs = validate(d)
        assert any("dose is infinity" in e for e in errs)


class TestDoseUnit:
    """T-06, T-07: dose unit checks."""

    def test_empty_unit_rejected(self) -> None:
        d = ExperimentDesign(
            experiment_id="exp",
            groups=[
                ExperimentalGroup(group_id="g1", label="G", role=GroupRole.CONTROL)
            ],
            samples=[
                DesignSample(
                    sample_id="s1",
                    group_id="g1",
                    treatment=Treatment(
                        compound="X", dose=Quantity(value=100, unit="")
                    ),
                )
            ],
        )
        errs = validate(d)
        assert any("empty unit" in e for e in errs)

    def test_whitespace_unit_rejected(self) -> None:
        d = ExperimentDesign(
            experiment_id="exp",
            groups=[
                ExperimentalGroup(group_id="g1", label="G", role=GroupRole.CONTROL)
            ],
            samples=[
                DesignSample(
                    sample_id="s1",
                    group_id="g1",
                    treatment=Treatment(
                        compound="X", dose=Quantity(value=100, unit="  ")
                    ),
                )
            ],
        )
        errs = validate(d)
        assert any("whitespace-only" in e for e in errs)


class TestValidTreatment:
    """Valid treatment produces no violations."""

    def test_valid_treatment(self) -> None:
        d = ExperimentDesign(
            experiment_id="exp",
            groups=[
                ExperimentalGroup(group_id="g1", label="G", role=GroupRole.CONTROL)
            ],
            samples=[
                DesignSample(
                    sample_id="s1",
                    group_id="g1",
                    treatment=Treatment(
                        compound="erlotinib",
                        dose=Quantity(value=500, unit="nM"),
                        duration="24h",
                    ),
                )
            ],
        )
        errs = validate(d)
        assert not any("treatment" in e.lower() for e in errs)

    def test_compound_only(self) -> None:
        """Treatment with compound only (no dose/duration) is valid."""
        d = ExperimentDesign(
            experiment_id="exp",
            groups=[
                ExperimentalGroup(group_id="g1", label="G", role=GroupRole.CONTROL)
            ],
            samples=[
                DesignSample(
                    sample_id="s1",
                    group_id="g1",
                    treatment=Treatment(compound="DMSO"),
                )
            ],
        )
        errs = validate(d)
        assert not any("treatment" in e.lower() for e in errs)


class TestNoUnitConversion:
    """No unit conversion or equivalence is performed."""

    def test_arbitrary_unit_accepted(self) -> None:
        """Any unit string is accepted as-is."""
        d = ExperimentDesign(
            experiment_id="exp",
            groups=[
                ExperimentalGroup(group_id="g1", label="G", role=GroupRole.CONTROL)
            ],
            samples=[
                DesignSample(
                    sample_id="s1",
                    group_id="g1",
                    treatment=Treatment(
                        compound="X", dose=Quantity(value=1, unit="banana_units")
                    ),
                )
            ],
        )
        errs = validate(d)
        assert not any("unit" in e.lower() for e in errs)
