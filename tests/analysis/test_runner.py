"""Tests for pharmomics.analysis.runner.

Covers:
- validate_analysis_inputs passes on valid demo inputs.
- design validation failure raises AnalysisValidationError.
- design/omics compatibility failure raises AnalysisValidationError.
- spec/design compatibility failure raises AnalysisValidationError.
- Multiple violations are preserved in the exception message.
"""

from __future__ import annotations

import pytest

from pharmomics.analysis.example_data import make_demo_inputs
from pharmomics.analysis.runner import (
    AnalysisValidationError,
    validate_analysis_inputs,
)


def _make_invalid_design(**overrides: object):
    """Return a broken ExperimentDesign for testing."""
    from pharmomics.experiment.enums import GroupRole
    from pharmomics.experiment.schemas import (
        Contrast,
        DesignSample,
        ExperimentalFactor,
        ExperimentalGroup,
        ExperimentDesign,
    )

    return ExperimentDesign(
        experiment_id=overrides.get("experiment_id", ""),
        description=None,
        samples=overrides.get(
            "samples",
            [DesignSample(sample_id="s1", group_id="g1")],
        ),
        groups=overrides.get(
            "groups",
            [
                ExperimentalGroup(
                    group_id="g1", label="Group 1", role=GroupRole.CONTROL
                ),
            ],
        ),
        factors=overrides.get(
            "factors",
            [
                ExperimentalFactor(
                    factor_id="condition",
                    factor_type="categorical",
                    levels=["control", "treatment"],
                ),
            ],
        ),
        contrasts=overrides.get(
            "contrasts",
            [
                Contrast(
                    contrast_id="treated_vs_control",
                    comparison_group_id="g1",
                    reference_group_id="g1",
                ),
            ],
        ),
    )


def _make_omics_with_samples(sample_ids: list[str]):
    """Return an OmicsMatrix with the given sample_ids."""
    from datetime import UTC, datetime

    import pandas as pd

    from pharmomics.omics.enums import MeasurementType, NormalizationStatus
    from pharmomics.omics.schemas import OmicsMatrix

    df = pd.DataFrame(
        {"feature_id": ["G1", "G2"]}
        | {sid: [1.0, 2.0] for sid in sample_ids},
    )
    return OmicsMatrix(
        matrix_id="mx-test",
        modality="transcriptomics",
        feature_type="gene",
        measurement_type=MeasurementType.ESTIMATED_COUNTS,
        normalization_status=NormalizationStatus.RAW,
        n_features=2,
        n_samples=len(sample_ids),
        feature_ids=["G1", "G2"],
        sample_ids=list(sample_ids),
        dataframe=df,
        created_at=datetime.now(UTC).strftime("%Y-%m-%dT%H:%M:%SZ"),
    )


# ---------------------------------------------------------------------------
# Valid inputs
# ---------------------------------------------------------------------------


class TestValidInputs:
    """validate_analysis_inputs succeeds on a well-formed triple."""

    def test_valid_demo_inputs_pass(self) -> None:
        omics, design, spec = make_demo_inputs()
        # Should not raise
        validate_analysis_inputs(spec, design, omics)


# ---------------------------------------------------------------------------
# Design validation failure
# ---------------------------------------------------------------------------


class TestDesignValidationFailure:
    """Empty experiment_id violates validate(design)."""

    def test_invalid_design_raises(self) -> None:
        omics, _, spec = make_demo_inputs()
        bad_design = _make_invalid_design(experiment_id="")

        with pytest.raises(AnalysisValidationError) as excinfo:
            validate_analysis_inputs(spec, bad_design, omics)

        assert "Empty experiment_id" in str(excinfo.value)


# ---------------------------------------------------------------------------
# Design / omics compatibility failure
# ---------------------------------------------------------------------------


class TestOmicsCompatibilityFailure:
    """OmicsMatrix missing samples referenced by the design."""

    def test_missing_samples_raises(self) -> None:
        omics, design, spec = make_demo_inputs()
        # Replace omics with one that has none of the design's samples
        bad_omics = _make_omics_with_samples(["unknown_1", "unknown_2"])

        with pytest.raises(AnalysisValidationError) as excinfo:
            validate_analysis_inputs(spec, design, bad_omics)

        assert "No overlap" in str(excinfo.value)


# ---------------------------------------------------------------------------
# Spec / design compatibility failure
# ---------------------------------------------------------------------------


class TestSpecDesignCompatibilityFailure:
    """AnalysisSpecification references a contrast not in the design."""

    def test_unknown_contrast_raises(self) -> None:
        from pharmomics.analysis.schemas import AnalysisSpecification

        omics, design, _ = make_demo_inputs()
        bad_spec = AnalysisSpecification(
            analysis_type="differential_analysis",
            contrast_references=["nonexistent_contrast"],
        )

        with pytest.raises(AnalysisValidationError) as excinfo:
            validate_analysis_inputs(bad_spec, design, omics)

        assert "nonexistent_contrast" in str(excinfo.value)
        assert "contrast_reference" in str(excinfo.value)


# ---------------------------------------------------------------------------
# Multiple violations preserved
# ---------------------------------------------------------------------------


class TestMultipleViolationsPreserved:
    """All violations from all checks appear in the exception message."""

    def test_all_violations_in_message(self) -> None:
        """Design has empty experiment_id AND contrast references itself
        (same group on both sides) AND spec references unknown factor."""
        from pharmomics.analysis.schemas import AnalysisSpecification
        from pharmomics.experiment.enums import GroupRole
        from pharmomics.experiment.schemas import (
            Contrast,
            DesignSample,
            ExperimentalFactor,
            ExperimentalGroup,
            ExperimentDesign,
        )

        bad_design = ExperimentDesign(
            experiment_id="",  # violation I-12
            samples=[
                DesignSample(sample_id="s1", group_id="g1"),
            ],
            groups=[
                ExperimentalGroup(
                    group_id="g1", label="G1", role=GroupRole.CONTROL,
                ),
            ],
            factors=[
                ExperimentalFactor(
                    factor_id="drug",
                    factor_type="categorical",
                    levels=["a", "b"],
                ),
            ],
            contrasts=[
                Contrast(
                    contrast_id="c1",
                    comparison_group_id="g1",
                    reference_group_id="g1",  # violation C-01
                ),
            ],
        )
        bad_spec = AnalysisSpecification(
            analysis_type="differential_analysis",
            factor_references=["nonexistent_factor"],  # violation
        )
        bad_omics = _make_omics_with_samples(["s1"])

        with pytest.raises(AnalysisValidationError) as excinfo:
            validate_analysis_inputs(bad_spec, bad_design, bad_omics)

        msg = str(excinfo.value)
        # All three violation sources should appear
        assert "Empty experiment_id" in msg
        assert "nonexistent_factor" in msg
        # C-01: same group on both sides
        assert "same group on both sides" in msg
