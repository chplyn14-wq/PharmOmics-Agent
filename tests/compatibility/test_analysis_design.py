"""Tests for AnalysisSpecification ↔ ExperimentDesign compatibility checks."""

from __future__ import annotations

import copy

from pharmomics.analysis.schemas import AnalysisSpecification
from pharmomics.compatibility.analysis_design import (
    check_analysis_design_compatibility,
)
from pharmomics.experiment.enums import GroupRole
from pharmomics.experiment.schemas import (
    Contrast,
    DesignSample,
    ExperimentalFactor,
    ExperimentalGroup,
    ExperimentDesign,
)

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

_FACTOR_TYPE_CATEGORICAL = "categorical"


def _make_design(
    *,
    factor_ids: list[str] | None = None,
    contrast_ids: list[str] | None = None,
) -> ExperimentDesign:
    """Construct a minimal ExperimentDesign with the given IDs."""
    factor_ids = factor_ids or ["drug"]
    contrast_ids = contrast_ids or ["treated_vs_vehicle"]

    return ExperimentDesign(
        experiment_id="exp-test",
        samples=[
            DesignSample(sample_id="s1", group_id="g1"),
        ],
        groups=[
            ExperimentalGroup(
                group_id="g1",
                label="Group 1",
                role=GroupRole.CONTROL,
            ),
            ExperimentalGroup(
                group_id="g2",
                label="Group 2",
                role=GroupRole.TREATMENT,
            ),
        ],
        factors=[
            ExperimentalFactor(
                factor_id=fid,
                factor_type=_FACTOR_TYPE_CATEGORICAL,
            )
            for fid in factor_ids
        ],
        contrasts=[
            Contrast(
                contrast_id=cid,
                comparison_group_id="g2",
                reference_group_id="g1",
            )
            for cid in contrast_ids
        ],
    )


# ---------------------------------------------------------------------------
# All references present
# ---------------------------------------------------------------------------


class TestAllReferencesPresent:
    """All factor_references and contrast_references exist in the design."""

    def test_all_references_match(self) -> None:
        design = _make_design(
            factor_ids=["drug", "time"],
            contrast_ids=["treated_vs_vehicle", "high_vs_low"],
        )
        spec = AnalysisSpecification(
            analysis_type="differential_analysis",
            factor_references=["drug", "time"],
            contrast_references=["treated_vs_vehicle", "high_vs_low"],
        )
        assert check_analysis_design_compatibility(spec, design) == []

    def test_subset_of_references_is_valid(self) -> None:
        design = _make_design(
            factor_ids=["drug", "time", "dose"],
            contrast_ids=["a_vs_b", "c_vs_d"],
        )
        spec = AnalysisSpecification(
            analysis_type="differential_analysis",
            factor_references=["drug"],
            contrast_references=["a_vs_b"],
        )
        assert check_analysis_design_compatibility(spec, design) == []


# ---------------------------------------------------------------------------
# Missing factor references
# ---------------------------------------------------------------------------


class TestMissingFactorReferences:
    """factor_reference(s) not found in ExperimentDesign."""

    def test_single_missing_factor(self) -> None:
        design = _make_design(factor_ids=["drug"])
        spec = AnalysisSpecification(
            analysis_type="differential_analysis",
            factor_references=["drug", "nonexistent"],
        )
        errs = check_analysis_design_compatibility(spec, design)
        assert any("factor_reference" in e for e in errs)
        assert any("nonexistent" in e for e in errs)

    def test_all_missing_factors(self) -> None:
        design = _make_design(factor_ids=["drug"])
        spec = AnalysisSpecification(
            analysis_type="differential_analysis",
            factor_references=["foo", "bar"],
        )
        errs = check_analysis_design_compatibility(spec, design)
        assert any("factor_reference" in e for e in errs)
        assert any("foo" in e for e in errs)
        assert any("bar" in e for e in errs)


# ---------------------------------------------------------------------------
# Missing contrast references
# ---------------------------------------------------------------------------


class TestMissingContrastReferences:
    """contrast_reference(s) not found in ExperimentDesign."""

    def test_single_missing_contrast(self) -> None:
        design = _make_design(contrast_ids=["treated_vs_vehicle"])
        spec = AnalysisSpecification(
            analysis_type="differential_analysis",
            contrast_references=["treated_vs_vehicle", "ghost_contrast"],
        )
        errs = check_analysis_design_compatibility(spec, design)
        assert any("contrast_reference" in e for e in errs)
        assert any("ghost_contrast" in e for e in errs)

    def test_all_missing_contrasts(self) -> None:
        design = _make_design(contrast_ids=["treated_vs_vehicle"])
        spec = AnalysisSpecification(
            analysis_type="differential_analysis",
            contrast_references=["x", "y"],
        )
        errs = check_analysis_design_compatibility(spec, design)
        assert any("contrast_reference" in e for e in errs)
        assert any("x" in e for e in errs)
        assert any("y" in e for e in errs)


# ---------------------------------------------------------------------------
# Both missing
# ---------------------------------------------------------------------------


class TestBothMissing:
    """Both factor and contrast references have missing entries."""

    def test_both_types_missing_reported(self) -> None:
        design = _make_design(factor_ids=["drug"], contrast_ids=["a_vs_b"])
        spec = AnalysisSpecification(
            analysis_type="differential_analysis",
            factor_references=["no_factor"],
            contrast_references=["no_contrast"],
        )
        errs = check_analysis_design_compatibility(spec, design)
        assert any("factor_reference" in e for e in errs)
        assert any("contrast_reference" in e for e in errs)
        assert len(errs) == 2


# ---------------------------------------------------------------------------
# Empty reference lists (valid)
# ---------------------------------------------------------------------------


class TestEmptyReferences:
    """Empty factor_references or contrast_references are valid."""

    def test_empty_factor_references(self) -> None:
        design = _make_design()
        spec = AnalysisSpecification(
            analysis_type="clustering",
            factor_references=[],
            contrast_references=["treated_vs_vehicle"],
        )
        assert check_analysis_design_compatibility(spec, design) == []

    def test_empty_contrast_references(self) -> None:
        design = _make_design()
        spec = AnalysisSpecification(
            analysis_type="pca",
            factor_references=["drug"],
            contrast_references=[],
        )
        assert check_analysis_design_compatibility(spec, design) == []

    def test_both_empty_references(self) -> None:
        design = _make_design()
        spec = AnalysisSpecification(
            analysis_type="custom",
            factor_references=[],
            contrast_references=[],
        )
        assert check_analysis_design_compatibility(spec, design) == []


# ---------------------------------------------------------------------------
# Purity
# ---------------------------------------------------------------------------


class TestPurity:
    """Compatibility check does not mutate inputs."""

    def test_does_not_mutate_specification(self) -> None:
        design = _make_design()
        spec = AnalysisSpecification(
            analysis_type="differential_analysis",
            factor_references=["drug"],
            contrast_references=["treated_vs_vehicle"],
            parameters={"alpha": 0.05},
        )
        original = copy.deepcopy(spec)
        check_analysis_design_compatibility(spec, design)
        assert spec == original

    def test_does_not_mutate_design(self) -> None:
        design = _make_design()
        spec = AnalysisSpecification(
            analysis_type="differential_analysis",
            factor_references=["drug"],
        )
        original = copy.deepcopy(design)
        check_analysis_design_compatibility(spec, design)
        assert design == original
