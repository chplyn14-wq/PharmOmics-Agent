"""Tests for validation purity and determinism."""

from __future__ import annotations

import copy

from pharmomics.experiment.enums import GroupRole
from pharmomics.experiment.schemas import (
    DesignSample,
    ExperimentalGroup,
    ExperimentDesign,
)
from pharmomics.experiment.validation import validate


class TestNoMutation:
    """validate() must not mutate the input ExperimentDesign."""

    def test_validate_does_not_mutate_design(
        self,
        simple_drug_control_design: ExperimentDesign,
    ) -> None:
        original = copy.deepcopy(simple_drug_control_design)
        validate(simple_drug_control_design)
        assert simple_drug_control_design == original


class TestDeterminism:
    """Same input always produces the same ordered error list."""

    def test_repeated_calls_identical_result(
        self,
        invalid_reference_design: ExperimentDesign,
    ) -> None:
        r1 = validate(invalid_reference_design)
        r2 = validate(invalid_reference_design)
        assert r1 == r2

    def test_error_order_stable(
        self,
        invalid_reference_design: ExperimentDesign,
    ) -> None:
        r1 = validate(invalid_reference_design)
        r2 = validate(invalid_reference_design)
        assert len(r1) == len(r2)
        for a, b in zip(r1, r2):
            assert a == b


class TestValidDesignEmptyErrors:
    """Valid designs produce an empty error list."""

    def test_minimal_valid_no_errors(
        self,
        minimal_valid_design: ExperimentDesign,
    ) -> None:
        assert validate(minimal_valid_design) == []

    def test_simple_drug_control_no_errors(
        self,
        simple_drug_control_design: ExperimentDesign,
    ) -> None:
        assert validate(simple_drug_control_design) == []

    def test_paired_no_errors(
        self,
        paired_before_after_design: ExperimentDesign,
    ) -> None:
        assert validate(paired_before_after_design) == []

    def test_multi_factor_no_errors(
        self,
        multi_factor_design: ExperimentDesign,
    ) -> None:
        assert validate(multi_factor_design) == []

    def test_batch_adjusted_no_errors(
        self,
        batch_adjusted_design: ExperimentDesign,
    ) -> None:
        assert validate(batch_adjusted_design) == []


class TestInvalidDesignHasErrors:
    """Invalid designs produce at least one error."""

    def test_invalid_reference_has_errors(
        self,
        invalid_reference_design: ExperimentDesign,
    ) -> None:
        assert len(validate(invalid_reference_design)) > 0


class TestReplicateNoEnforcement:
    """RP-01…RP-03: replicate checks do not produce violations."""

    def test_empty_replicate_id_allowed(
        self,
        minimal_valid_design: ExperimentDesign,
    ) -> None:
        # Rebuild with empty replicate — model is frozen so we construct fresh
        d = ExperimentDesign(
            experiment_id="exp-rep",
            samples=[
                DesignSample(
                    sample_id="s1",
                    group_id="g1",
                    biological_replicate="",
                ),
            ],
            groups=[
                ExperimentalGroup(
                    group_id="g1",
                    label="G",
                    role=GroupRole.CONTROL,
                ),
            ],
        )
        # Just verify no validation error from empty replicate
        errs = validate(d)
        assert not any("replicate" in e.lower() for e in errs)

    def test_no_replicate_id_allowed(
        self,
        minimal_valid_design: ExperimentDesign,
    ) -> None:
        # replicate is already None by default — should not error
        assert validate(minimal_valid_design) == []
