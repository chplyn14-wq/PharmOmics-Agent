"""Tests for experiment domain schema construction and serialization."""

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
# Quantity
# ---------------------------------------------------------------------------


class TestQuantity:
    def test_valid(self) -> None:
        q = Quantity(value=500.0, unit="nM")
        assert q.value == 500.0
        assert q.unit == "nM"

    def test_zero_value(self) -> None:
        q = Quantity(value=0, unit="mg")
        assert q.value == 0

    def test_frozen(self) -> None:
        q = Quantity(value=1, unit="nM")
        with pytest.raises(Exception):
            q.value = 2


# ---------------------------------------------------------------------------
# Treatment
# ---------------------------------------------------------------------------


class TestTreatment:
    def test_valid_full(self) -> None:
        t = Treatment(
            compound="erlotinib",
            dose=Quantity(value=500, unit="nM"),
            duration="24h",
            description="EGFR inhibitor",
        )
        assert t.compound == "erlotinib"
        assert t.dose is not None
        assert t.dose.value == 500
        assert t.duration == "24h"

    def test_valid_minimal(self) -> None:
        t = Treatment(compound="DMSO")
        assert t.compound == "DMSO"
        assert t.dose is None
        assert t.duration is None

    def test_frozen(self) -> None:
        t = Treatment(compound="X")
        with pytest.raises(Exception):
            t.compound = "Y"


# ---------------------------------------------------------------------------
# ExperimentalGroup
# ---------------------------------------------------------------------------


class TestExperimentalGroup:
    def test_valid(self) -> None:
        g = ExperimentalGroup(
            group_id="ctrl",
            label="Control",
            description="Vehicle control group",
            role=GroupRole.CONTROL,
        )
        assert g.group_id == "ctrl"
        assert g.role == GroupRole.CONTROL

    def test_all_roles(self) -> None:
        for role in GroupRole:
            g = ExperimentalGroup(
                group_id="g",
                label="L",
                role=role,
            )
            assert g.role == role

    def test_frozen(self) -> None:
        g = ExperimentalGroup(group_id="g", label="L", role=GroupRole.OTHER)
        with pytest.raises(Exception):
            g.label = "X"


# ---------------------------------------------------------------------------
# ExperimentalFactor
# ---------------------------------------------------------------------------


class TestExperimentalFactor:
    def test_valid_categorical(self) -> None:
        f = ExperimentalFactor(
            factor_id="drug",
            factor_type=FactorType.CATEGORICAL,
            levels=["A", "B", "C"],
        )
        assert f.levels == ["A", "B", "C"]

    def test_valid_continuous(self) -> None:
        f = ExperimentalFactor(
            factor_id="time",
            factor_type=FactorType.CONTINUOUS,
            description="Time in hours",
        )
        assert f.levels is None

    def test_frozen(self) -> None:
        f = ExperimentalFactor(
            factor_id="f",
            factor_type=FactorType.CATEGORICAL,
        )
        with pytest.raises(Exception):
            f.factor_id = "x"


# ---------------------------------------------------------------------------
# CovariateDefinition
# ---------------------------------------------------------------------------


class TestCovariateDefinition:
    def test_valid_batch(self) -> None:
        c = CovariateDefinition(
            covariate_id="batch",
            role=CovariateRole.BATCH,
            value_type=CovariateValueType.CATEGORICAL,
        )
        assert c.role == CovariateRole.BATCH

    def test_valid_clinical(self) -> None:
        c = CovariateDefinition(
            covariate_id="age",
            role=CovariateRole.CLINICAL,
            value_type=CovariateValueType.CONTINUOUS,
            description="Subject age",
        )
        assert c.value_type == CovariateValueType.CONTINUOUS

    def test_frozen(self) -> None:
        c = CovariateDefinition(
            covariate_id="c",
            role=CovariateRole.OTHER,
            value_type=CovariateValueType.CATEGORICAL,
        )
        with pytest.raises(Exception):
            c.covariate_id = "x"


# ---------------------------------------------------------------------------
# DesignSample
# ---------------------------------------------------------------------------


class TestDesignSample:
    def test_valid_minimal(self) -> None:
        s = DesignSample(sample_id="s1", group_id="g1")
        assert s.factor_values == {}
        assert s.covariate_values == {}
        assert s.treatment is None

    def test_valid_full(self) -> None:
        s = DesignSample(
            sample_id="s1",
            group_id="g1",
            factor_values={"drug": "A"},
            treatment=Treatment(compound="X", dose=Quantity(value=100, unit="nM")),
            biological_replicate="R1",
            covariate_values={"batch": "B1"},
            pair_id="P1",
        )
        assert s.factor_values["drug"] == "A"
        assert s.treatment is not None
        assert s.pair_id == "P1"

    def test_mutable_default_isolation(self) -> None:
        s1 = DesignSample(sample_id="s1", group_id="g1")
        s2 = DesignSample(sample_id="s2", group_id="g1")
        assert s1.factor_values is not s2.factor_values
        assert s1.covariate_values is not s2.covariate_values

    def test_frozen(self) -> None:
        s = DesignSample(sample_id="s1", group_id="g1")
        with pytest.raises(Exception):
            s.sample_id = "s2"


# ---------------------------------------------------------------------------
# Contrast
# ---------------------------------------------------------------------------


class TestContrast:
    def test_valid(self) -> None:
        c = Contrast(
            contrast_id="c1",
            comparison_group_id="treatment",
            reference_group_id="control",
            description="Drug vs vehicle",
        )
        assert c.comparison_group_id == "treatment"
        assert c.reference_group_id == "control"

    def test_minimal(self) -> None:
        c = Contrast(
            contrast_id="c1",
            comparison_group_id="A",
            reference_group_id="B",
        )
        assert c.description is None

    def test_frozen(self) -> None:
        c = Contrast(
            contrast_id="c",
            comparison_group_id="A",
            reference_group_id="B",
        )
        with pytest.raises(Exception):
            c.comparison_group_id = "C"


# ---------------------------------------------------------------------------
# PairingDefinition
# ---------------------------------------------------------------------------


class TestPairingDefinition:
    def test_valid(self) -> None:
        p = PairingDefinition(
            pairing_type="before_after",
            description="Before and after treatment",
        )
        assert p.pairing_type == "before_after"

    def test_frozen(self) -> None:
        p = PairingDefinition(pairing_type="matched")
        with pytest.raises(Exception):
            p.pairing_type = "other"


# ---------------------------------------------------------------------------
# ExperimentDesign
# ---------------------------------------------------------------------------


class TestExperimentDesign:
    def test_valid_minimal(self) -> None:
        d = ExperimentDesign(
            schema_version="1.0.0",
            experiment_id="exp-1",
        )
        assert d.samples == []
        assert d.contrasts == []
        assert d.pairing is None

    def test_valid_full(
        self,
        simple_drug_control_design: ExperimentDesign,
    ) -> None:
        d = simple_drug_control_design
        assert len(d.samples) == 6
        assert len(d.groups) == 2
        assert len(d.factors) == 1
        assert len(d.contrasts) == 1
        assert d.schema_version == "1.0.0"

    def test_mutable_default_isolation(self) -> None:
        d1 = ExperimentDesign(
            schema_version="1.0.0",
            experiment_id="exp-1",
        )
        d2 = ExperimentDesign(
            schema_version="1.0.0",
            experiment_id="exp-2",
        )
        assert d1.samples is not d2.samples
        assert d1.groups is not d2.groups
        assert d1.contrasts is not d2.contrasts

    def test_frozen(self) -> None:
        d = ExperimentDesign(
            schema_version="1.0.0",
            experiment_id="exp-1",
        )
        with pytest.raises(Exception):
            d.experiment_id = "exp-2"

    def test_json_roundtrip(
        self,
        simple_drug_control_design: ExperimentDesign,
    ) -> None:
        d = simple_drug_control_design
        data = d.model_dump()
        roundtrip = ExperimentDesign(**data)
        assert roundtrip == d

    def test_schema_version_default(self) -> None:
        d = ExperimentDesign(
            experiment_id="exp",
        )
        assert d.schema_version == "1.0.0"
