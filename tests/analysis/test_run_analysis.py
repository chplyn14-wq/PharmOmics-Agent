"""Tests for pharmomics.analysis.run.

Covers:
- run_analysis returns a valid AnalysisResult on valid inputs.
- Result has the expected warning placeholder.
- Result structure is internally consistent.
- Unsupported analysis_type raises AnalysisValidationError.
- Validation failures propagate correctly.
"""

from __future__ import annotations

import pytest

from pharmomics.analysis.example_data import make_demo_inputs
from pharmomics.analysis.results import AnalysisResult
from pharmomics.analysis.run import run_analysis
from pharmomics.analysis.runner import AnalysisValidationError


class TestRunAnalysisSuccess:
    """run_analysis succeeds and returns a well-formed AnalysisResult."""

    def test_returns_result_on_valid_inputs(self) -> None:
        omics, design, spec = make_demo_inputs()
        result = run_analysis(spec, design, omics)
        assert isinstance(result, AnalysisResult)
        assert result.analysis_type == "differential_analysis"
        assert result.contrast_id == "treated_vs_control"

    def test_result_has_warning(self) -> None:
        omics, design, spec = make_demo_inputs()
        result = run_analysis(spec, design, omics)
        assert "not implemented" in result.warnings[0].lower()

    def test_result_structure_valid(self) -> None:
        omics, design, spec = make_demo_inputs()
        result = run_analysis(spec, design, omics)
        # n_genes_tested must equal len(gene_results) — enforced by __post_init__
        assert result.n_genes_tested == len(result.gene_results)


class TestUnsupportedAnalysisType:
    """run_analysis rejects analysis types it does not support."""

    def test_unsupported_type_raises(self) -> None:
        from pharmomics.analysis.schemas import AnalysisSpecification

        omics, design, _ = make_demo_inputs()
        unsupported_spec = AnalysisSpecification(
            analysis_type="clustering",
        )

        with pytest.raises(AnalysisValidationError) as excinfo:
            run_analysis(unsupported_spec, design, omics)

        assert "Unsupported analysis_type" in str(excinfo.value)
        assert "clustering" in str(excinfo.value)


class TestRunAnalysisValidationGate:
    """Validation failures propagate before dispatch."""

    def test_validation_failure_propagates(self) -> None:
        omics, _, spec = make_demo_inputs()
        # Create a design with an empty experiment_id
        from pharmomics.experiment.enums import GroupRole
        from pharmomics.experiment.schemas import (
            DesignSample,
            ExperimentalGroup,
            ExperimentDesign,
        )

        bad_design = ExperimentDesign(
            experiment_id="",
            samples=[DesignSample(sample_id="s1", group_id="g1")],
            groups=[
                ExperimentalGroup(
                    group_id="g1", label="G1", role=GroupRole.CONTROL,
                ),
            ],
        )

        with pytest.raises(AnalysisValidationError) as excinfo:
            run_analysis(spec, bad_design, omics)

        # The error should come from validate(design), not dispatch
        assert "Empty experiment_id" in str(excinfo.value)
