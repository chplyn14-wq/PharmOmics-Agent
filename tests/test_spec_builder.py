"""Tests for pharmomics.cli.analyze — build_analysis_spec."""

from __future__ import annotations

from pharmomics.analysis.schemas import AnalysisSpecification
from pharmomics.cli.analyze import build_analysis_spec


class TestBuildAnalysisSpec:
    """Verify AnalysisSpecification construction."""

    def test_differential_analysis_type(self) -> None:
        spec = build_analysis_spec(contrast_id="trt_vs_ctrl")
        assert spec.analysis_type == "differential_analysis"

    def test_contrast_reference_set(self) -> None:
        spec = build_analysis_spec(contrast_id="osi_dtp_vs_dmso")
        assert spec.contrast_references == ["osi_dtp_vs_dmso"]

    def test_factor_reference_set(self) -> None:
        spec = build_analysis_spec(contrast_id="trt_vs_ctrl")
        assert spec.factor_references == ["condition"]

    def test_default_fdr_threshold(self) -> None:
        spec = build_analysis_spec(contrast_id="trt_vs_ctrl")
        assert spec.parameters["fdr_threshold"] == 0.05

    def test_custom_fdr_threshold(self) -> None:
        spec = build_analysis_spec(
            contrast_id="trt_vs_ctrl",
            fdr_threshold=0.01,
        )
        assert spec.parameters["fdr_threshold"] == 0.01

    def test_custom_analysis_type(self) -> None:
        spec = build_analysis_spec(
            analysis_type="some_other_type",
            contrast_id="trt_vs_ctrl",
        )
        assert spec.analysis_type == "some_other_type"

    def test_returns_analysis_specification(self) -> None:
        spec = build_analysis_spec(contrast_id="trt_vs_ctrl")
        assert isinstance(spec, AnalysisSpecification)

    def test_empty_contrast_id(self) -> None:
        spec = build_analysis_spec(contrast_id="")
        assert spec.contrast_references == [""]
