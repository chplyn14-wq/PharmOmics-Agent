"""Tests for pharmomics.analysis.example_data.

Covers:
- make_demo_inputs() returns three objects of the correct types.
- OmicsMatrix shape and determinism.
- ExperimentDesign has two groups and one contrast.
- AnalysisSpecification references match design entities.
- Cross-validation: validate(design), check_compatibility,
  check_analysis_design_compatibility.
- Determinism: consecutive calls produce equal data.
"""

from __future__ import annotations

import pandas as pd

from pharmomics.analysis.example_data import make_demo_inputs
from pharmomics.analysis.schemas import AnalysisSpecification
from pharmomics.compatibility.analysis_design import (
    check_analysis_design_compatibility,
)
from pharmomics.compatibility.omics_design import check_compatibility
from pharmomics.experiment.schemas import ExperimentDesign
from pharmomics.experiment.validation import validate
from pharmomics.omics.schemas import OmicsMatrix


class TestMakeDemoInputs:
    """make_demo_inputs returns three correctly typed objects."""

    def test_returns_three_objects(self) -> None:
        omics, design, spec = make_demo_inputs()
        assert isinstance(omics, OmicsMatrix)
        assert isinstance(design, ExperimentDesign)
        assert isinstance(spec, AnalysisSpecification)

    def test_omics_has_expected_shape(self) -> None:
        omics, _, _ = make_demo_inputs()
        assert omics.n_features == 6
        assert omics.n_samples == 6
        assert omics.dataframe.shape == (6, 7)  # feature_id + 6 sample cols

    def test_omics_feature_ids_are_correct(self) -> None:
        omics, _, _ = make_demo_inputs()
        assert omics.feature_ids == [
            "EGFR",
            "ERBB2",
            "TP53",
            "MYC",
            "KRAS",
            "PTEN",
        ]

    def test_omics_sample_ids_are_correct(self) -> None:
        omics, _, _ = make_demo_inputs()
        assert omics.sample_ids == [
            "ctrl_1",
            "ctrl_2",
            "ctrl_3",
            "trt_1",
            "trt_2",
            "trt_3",
        ]

    def test_design_has_two_groups(self) -> None:
        _, design, _ = make_demo_inputs()
        assert len(design.groups) == 2
        group_ids = {g.group_id for g in design.groups}
        assert group_ids == {"ctrl", "trt"}

    def test_design_has_six_samples(self) -> None:
        _, design, _ = make_demo_inputs()
        assert len(design.samples) == 6

    def test_design_has_one_contrast(self) -> None:
        _, design, _ = make_demo_inputs()
        assert len(design.contrasts) == 1
        assert design.contrasts[0].contrast_id == "treated_vs_control"

    def test_design_has_one_factor(self) -> None:
        _, design, _ = make_demo_inputs()
        assert len(design.factors) == 1
        assert design.factors[0].factor_id == "condition"

    def test_spec_is_differential_analysis(self) -> None:
        _, _, spec = make_demo_inputs()
        assert spec.analysis_type == "differential_analysis"

    def test_spec_references_match(self) -> None:
        _, _, spec = make_demo_inputs()
        assert spec.contrast_references == ["treated_vs_control"]
        assert spec.factor_references == ["condition"]


class TestCrossValidation:
    """The triple passes all existing cross-domain validators."""

    def test_design_passes_validate(self) -> None:
        _, design, _ = make_demo_inputs()
        assert validate(design) == []

    def test_omics_and_design_are_compatible(self) -> None:
        omics, design, _ = make_demo_inputs()
        assert check_compatibility(design, omics) == []

    def test_spec_and_design_are_compatible(self) -> None:
        _, design, spec = make_demo_inputs()
        assert check_analysis_design_compatibility(spec, design) == []


class TestDeterminism:
    """Consecutive calls produce structurally equal results."""

    def test_two_calls_produce_equal_omics(self) -> None:
        omics1, _, _ = make_demo_inputs()
        omics2, _, _ = make_demo_inputs()
        assert omics1.feature_ids == omics2.feature_ids
        assert omics1.sample_ids == omics2.sample_ids
        assert omics1.dataframe.equals(omics2.dataframe)
        assert omics1.n_features == omics2.n_features
        assert omics1.n_samples == omics2.n_samples

    def test_two_calls_produce_equal_design(self) -> None:
        _, design1, _ = make_demo_inputs()
        _, design2, _ = make_demo_inputs()
        assert design1 == design2

    def test_two_calls_produce_equal_spec(self) -> None:
        _, _, spec1 = make_demo_inputs()
        _, _, spec2 = make_demo_inputs()
        assert spec1.analysis_type == spec2.analysis_type
        assert spec1.factor_references == spec2.factor_references
        assert spec1.contrast_references == spec2.contrast_references
        assert spec1.parameters == spec2.parameters

    def test_dataframe_values_are_not_random(self) -> None:
        """Treatment values should be approximately 2× control values."""
        omics, _, _ = make_demo_inputs()
        df = omics.dataframe
        ctrl_cols = ["ctrl_1", "ctrl_2", "ctrl_3"]
        trt_cols = ["trt_1", "trt_2", "trt_3"]
        ctrl_mean = df[ctrl_cols].mean(axis=1)
        trt_mean = df[trt_cols].mean(axis=1)
        # log2FC ≈ 1 means trt/ctrl ≈ 2
        ratios = trt_mean / ctrl_mean
        pd.testing.assert_series_equal(
            ratios,
            pd.Series([2.0] * 6, index=range(6)),
            check_names=False,
        )
