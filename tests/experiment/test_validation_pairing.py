"""Tests for pairing validation (P-01 … P-03)."""

from __future__ import annotations

from pharmomics.experiment.enums import GroupRole
from pharmomics.experiment.schemas import (
    DesignSample,
    ExperimentalGroup,
    ExperimentDesign,
    PairingDefinition,
)
from pharmomics.experiment.validation import validate


class TestPairingDefinedNoSamples:
    """P-01: pairing defined but no samples have pair_id."""

    def test_pairing_without_samples_flagged(self) -> None:
        d = ExperimentDesign(
            experiment_id="exp",
            groups=[
                ExperimentalGroup(group_id="g1", label="G", role=GroupRole.CONTROL)
            ],
            samples=[DesignSample(sample_id="s1", group_id="g1")],
            pairing=PairingDefinition(pairing_type="before_after"),
        )
        errs = validate(d)
        assert any("no samples have a pair_id" in e for e in errs)


class TestSingleSamplePair:
    """P-02: pair with only one sample is flagged."""

    def test_one_sample_pair_flagged(self) -> None:
        d = ExperimentDesign(
            experiment_id="exp",
            groups=[
                ExperimentalGroup(group_id="g1", label="G", role=GroupRole.CONTROL)
            ],
            samples=[
                DesignSample(sample_id="s1", group_id="g1", pair_id="lonely"),
            ],
            pairing=PairingDefinition(pairing_type="before_after"),
        )
        errs = validate(d)
        assert any("has only 1 sample" in e for e in errs)


class TestValidPairing:
    """Valid pairing produces no violations."""

    def test_valid_pairs(self, paired_before_after_design: ExperimentDesign) -> None:
        errs = validate(paired_before_after_design)
        assert not any("Pair" in e for e in errs)

    def test_multiple_pairs(self) -> None:
        d = ExperimentDesign(
            experiment_id="exp",
            groups=[
                ExperimentalGroup(group_id="g1", label="G", role=GroupRole.CONTROL)
            ],
            samples=[
                DesignSample(sample_id="s1", group_id="g1", pair_id="p1"),
                DesignSample(sample_id="s2", group_id="g1", pair_id="p1"),
                DesignSample(sample_id="s3", group_id="g1", pair_id="p2"),
                DesignSample(sample_id="s4", group_id="g1", pair_id="p2"),
            ],
            pairing=PairingDefinition(pairing_type="matched"),
        )
        errs = validate(d)
        assert not any("Pair" in e for e in errs)


class TestUnpairedSamplesAllowed:
    """P-03: missing pair_id is valid."""

    def test_mixed_paired_unpaired(self) -> None:
        d = ExperimentDesign(
            experiment_id="exp",
            groups=[
                ExperimentalGroup(group_id="g1", label="G", role=GroupRole.CONTROL)
            ],
            samples=[
                DesignSample(sample_id="s1", group_id="g1", pair_id="p1"),
                DesignSample(sample_id="s2", group_id="g1", pair_id="p1"),
                DesignSample(sample_id="s3", group_id="g1"),  # unpaired
            ],
            pairing=PairingDefinition(pairing_type="before_after"),
        )
        errs = validate(d)
        assert not any("unpaired" in e.lower() for e in errs)

    def test_no_pairing_no_pair_ids_valid(self) -> None:
        """No pairing defined and no pair_ids is valid."""
        d = ExperimentDesign(
            experiment_id="exp",
            groups=[
                ExperimentalGroup(group_id="g1", label="G", role=GroupRole.CONTROL)
            ],
            samples=[
                DesignSample(sample_id="s1", group_id="g1"),
                DesignSample(sample_id="s2", group_id="g1"),
            ],
        )
        errs = validate(d)
        assert not any("pair" in e.lower() for e in errs)
