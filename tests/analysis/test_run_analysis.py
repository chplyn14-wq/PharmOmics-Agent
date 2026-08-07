"""Tests for pharmomics.analysis.run.

Covers:
- run_analysis returns a valid AnalysisResult on valid inputs.
- End-to-end differential analysis produces real statistical results.
- Result structure is internally consistent.
- Unsupported analysis_type raises AnalysisValidationError.
- Validation failures propagate correctly.
"""

from __future__ import annotations

import math

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
        assert result.n_genes_tested == 6
        assert len(result.gene_results) == 6

    def test_no_placeholder_warning(self) -> None:
        omics, design, spec = make_demo_inputs()
        result = run_analysis(spec, design, omics)
        assert result.warnings == ()

    def test_result_structure_valid(self) -> None:
        omics, design, spec = make_demo_inputs()
        result = run_analysis(spec, design, omics)
        # n_genes_tested must equal len(gene_results) — enforced by __post_init__
        assert result.n_genes_tested == len(result.gene_results)


class TestEndToEndDifferentialAnalysis:
    """End-to-end deterministic differential analysis verification."""

    def test_full_pipeline_results(self) -> None:
        """make_demo_inputs() → run_analysis() → real statistical results."""
        omics, design, spec = make_demo_inputs()
        result = run_analysis(spec, design, omics)

        # Non-empty results
        assert len(result.gene_results) == 6
        assert result.n_genes_tested == 6
        assert result.contrast_id == "treated_vs_control"
        assert result.analysis_type == "differential_analysis"

        # All 6 demo genes present
        gene_ids = {g.gene_id for g in result.gene_results}
        assert gene_ids == {"EGFR", "ERBB2", "TP53", "MYC", "KRAS", "PTEN"}

        # Demo data: trt ≈ 2× ctrl → log2FC ≈ 1.0 for all genes.
        # Ctrl values are repeated (zero variance within group).
        # With zero variance, Welch's t-test returns NaN p-values.
        for g in result.gene_results:
            assert g.log2_fold_change == pytest.approx(1.0, abs=1e-9)
            # Zero variance within groups → NaN p-value, NaN adj_p_value,
            # significant=False
            assert math.isnan(g.p_value)
            assert math.isnan(g.adj_p_value)
            assert g.significant is False

        # Base mean: (ctrl_mean + trt_mean) / 2
        # EGFR: ctrl=100, trt=200 → base_mean=150
        egfr = next(g for g in result.gene_results if g.gene_id == "EGFR")
        assert egfr.base_mean == pytest.approx(150.0, abs=1e-9)

        # ERBB2: ctrl=200, trt=400 → base_mean=300
        erbb2 = next(g for g in result.gene_results if g.gene_id == "ERBB2")
        assert erbb2.base_mean == pytest.approx(300.0, abs=1e-9)

        # MYC: ctrl=500, trt=1000 → base_mean=750
        myc = next(g for g in result.gene_results if g.gene_id == "MYC")
        assert myc.base_mean == pytest.approx(750.0, abs=1e-9)

        # No warnings in successful run
        assert result.warnings == ()

    def test_gene_results_order_preserved(self) -> None:
        """Gene results appear in original feature-id order."""
        omics, design, spec = make_demo_inputs()
        result = run_analysis(spec, design, omics)
        expected_order = ["EGFR", "ERBB2", "TP53", "MYC", "KRAS", "PTEN"]
        assert [g.gene_id for g in result.gene_results] == expected_order


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
                    group_id="g1",
                    label="G1",
                    role=GroupRole.CONTROL,
                ),
            ],
        )

        with pytest.raises(AnalysisValidationError) as excinfo:
            run_analysis(spec, bad_design, omics)

        # The error should come from validate(design), not dispatch
        assert "Empty experiment_id" in str(excinfo.value)
