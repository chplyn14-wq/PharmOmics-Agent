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
        """Demo data has zero within-group variance → NaN p-values → warning."""
        omics, design, spec = make_demo_inputs()
        result = run_analysis(spec, design, omics)
        # Zero variance within groups produces NaN p-values; the warning
        # must be surfaced so the user is not silently misled.
        assert len(result.warnings) == 1
        assert "undefined p-values" in result.warnings[0]
        assert "zero within-group variance" in result.warnings[0]
        # All 6 demo genes are affected
        assert "6 gene(s)" in result.warnings[0]

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

        # Demo data has zero within-group variance → NaN p-values → warning
        assert len(result.warnings) == 1
        assert "6 gene(s)" in result.warnings[0]

    def test_gene_results_order_preserved(self) -> None:
        """Gene results appear in original feature-id order."""
        omics, design, spec = make_demo_inputs()
        result = run_analysis(spec, design, omics)
        expected_order = ["EGFR", "ERBB2", "TP53", "MYC", "KRAS", "PTEN"]
        assert [g.gene_id for g in result.gene_results] == expected_order


class TestInsufficientReplicates:
    """run_analysis rejects contrasts with < 2 samples per group."""

    def test_single_sample_per_group_raises(self) -> None:
        """One sample in each group → AnalysisValidationError before inference."""
        import pandas as pd

        from pharmomics.analysis.schemas import AnalysisSpecification
        from pharmomics.experiment.enums import FactorType, GroupRole
        from pharmomics.experiment.schemas import (
            Contrast,
            DesignSample,
            ExperimentalFactor,
            ExperimentalGroup,
            ExperimentDesign,
        )
        from pharmomics.omics.enums import (
            MeasurementType,
            Modality,
            NormalizationStatus,
        )
        from pharmomics.omics.schemas import OmicsMatrix

        df = pd.DataFrame(
            {"feature_id": ["GeneA"]}
            | {
                "trt_1": [100.0],
                "ref_1": [200.0],
            },
        )

        omics = OmicsMatrix(
            matrix_id="mx-single",
            schema_version="1.0.0",
            modality=Modality.TRANSCRIPTOMICS,
            feature_type="gene",
            measurement_type=MeasurementType.ESTIMATED_COUNTS,
            normalization_status=NormalizationStatus.RAW,
            n_features=1,
            n_samples=2,
            feature_ids=["GeneA"],
            sample_ids=["trt_1", "ref_1"],
            dataframe=df,
            created_at="2026-01-01T00:00:00Z",
        )

        design = ExperimentDesign(
            experiment_id="exp-single",
            description="Single replicate per group",
            samples=[
                DesignSample(
                    sample_id="trt_1",
                    group_id="trt",
                    factor_values={"condition": "treatment"},
                ),
                DesignSample(
                    sample_id="ref_1",
                    group_id="ref",
                    factor_values={"condition": "reference"},
                ),
            ],
            groups=[
                ExperimentalGroup(
                    group_id="trt", label="Treatment", role=GroupRole.TREATMENT
                ),
                ExperimentalGroup(
                    group_id="ref", label="Reference", role=GroupRole.CONTROL
                ),
            ],
            factors=[
                ExperimentalFactor(
                    factor_id="condition",
                    factor_type=FactorType.CATEGORICAL,
                    levels=["treatment", "reference"],
                ),
            ],
            contrasts=[
                Contrast(
                    contrast_id="trt_vs_ref",
                    comparison_group_id="trt",
                    reference_group_id="ref",
                    description="Treatment vs Reference",
                ),
            ],
        )

        spec = AnalysisSpecification(
            analysis_type="differential_analysis",
            factor_references=["condition"],
            contrast_references=["trt_vs_ref"],
        )

        with pytest.raises(AnalysisValidationError) as excinfo:
            run_analysis(spec, design, omics)

        assert "Insufficient replicates" in str(excinfo.value)
        assert "1 sample(s)" in str(excinfo.value)

    def test_two_vs_one_raises(self) -> None:
        """Two samples in comparison but only one in reference → rejected."""
        import pandas as pd

        from pharmomics.analysis.schemas import AnalysisSpecification
        from pharmomics.experiment.enums import FactorType, GroupRole
        from pharmomics.experiment.schemas import (
            Contrast,
            DesignSample,
            ExperimentalFactor,
            ExperimentalGroup,
            ExperimentDesign,
        )
        from pharmomics.omics.enums import (
            MeasurementType,
            Modality,
            NormalizationStatus,
        )
        from pharmomics.omics.schemas import OmicsMatrix

        df = pd.DataFrame(
            {"feature_id": ["GeneA"]}
            | {
                "trt_1": [100.0],
                "trt_2": [110.0],
                "ref_1": [200.0],
            },
        )

        omics = OmicsMatrix(
            matrix_id="mx-asymmetric",
            schema_version="1.0.0",
            modality=Modality.TRANSCRIPTOMICS,
            feature_type="gene",
            measurement_type=MeasurementType.ESTIMATED_COUNTS,
            normalization_status=NormalizationStatus.RAW,
            n_features=1,
            n_samples=3,
            feature_ids=["GeneA"],
            sample_ids=["trt_1", "trt_2", "ref_1"],
            dataframe=df,
            created_at="2026-01-01T00:00:00Z",
        )

        design = ExperimentDesign(
            experiment_id="exp-asymmetric",
            description="2 vs 1 replicates",
            samples=[
                DesignSample(
                    sample_id="trt_1",
                    group_id="trt",
                    factor_values={"condition": "treatment"},
                ),
                DesignSample(
                    sample_id="trt_2",
                    group_id="trt",
                    factor_values={"condition": "treatment"},
                ),
                DesignSample(
                    sample_id="ref_1",
                    group_id="ref",
                    factor_values={"condition": "reference"},
                ),
            ],
            groups=[
                ExperimentalGroup(
                    group_id="trt", label="Treatment", role=GroupRole.TREATMENT
                ),
                ExperimentalGroup(
                    group_id="ref", label="Reference", role=GroupRole.CONTROL
                ),
            ],
            factors=[
                ExperimentalFactor(
                    factor_id="condition",
                    factor_type=FactorType.CATEGORICAL,
                    levels=["treatment", "reference"],
                ),
            ],
            contrasts=[
                Contrast(
                    contrast_id="trt_vs_ref",
                    comparison_group_id="trt",
                    reference_group_id="ref",
                    description="Treatment vs Reference",
                ),
            ],
        )

        spec = AnalysisSpecification(
            analysis_type="differential_analysis",
            factor_references=["condition"],
            contrast_references=["trt_vs_ref"],
        )

        with pytest.raises(AnalysisValidationError) as excinfo:
            run_analysis(spec, design, omics)

        assert "Insufficient replicates" in str(excinfo.value)
        assert "reference group has 1 sample(s)" in str(excinfo.value)


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


class TestEndToEndFiniteStatistics:
    """End-to-end deterministic test exercising finite p-value → BH-FDR →
    significant/non-significant classification through run_analysis().

    Uses a dedicated fixture with 5 genes × 3 replicates per group,
    constructed so that group-internal variance is non-zero and group
    means differ by design.  This forces the full statistical path to
    produce finite (non-NaN) results.
    """

    @staticmethod
    def _make_finite_fixture():
        """Build (OmicsMatrix, ExperimentDesign, AnalysisSpecification) with
        non-zero within-group variance for all genes."""
        import pandas as pd

        from pharmomics.analysis.schemas import AnalysisSpecification
        from pharmomics.experiment.enums import FactorType, GroupRole
        from pharmomics.experiment.schemas import (
            Contrast,
            DesignSample,
            ExperimentalFactor,
            ExperimentalGroup,
            ExperimentDesign,
        )
        from pharmomics.omics.enums import (
            MeasurementType,
            Modality,
            NormalizationStatus,
        )
        from pharmomics.omics.schemas import OmicsMatrix

        genes = ["GeneA", "GeneB", "GeneC", "GeneD", "GeneE"]
        comp_samples = ["trt_a1", "trt_a2", "trt_a3"]
        ref_samples = ["ref_b1", "ref_b2", "ref_b3"]
        all_samples = comp_samples + ref_samples

        # Deterministic expression matrix — every group has non-zero variance.
        # Rows (by gene index): GeneA, GeneB, GeneC, GeneD, GeneE
        #
        # Design rationale:
        #   GeneA: large comp vs ref difference → very small p, significant
        #   GeneB: identical comp/ref values → equal means, equal var → p=1.0
        #   GeneC: different per-sample values but equal group means (12 each)
        #     → t=0 → p=1.0
        #   GeneD: small comp/ref difference → moderate p, non-significant
        #   GeneE: 2× fold-change (220 vs 110) → log2fc=1.0, borderline significant
        df = pd.DataFrame(
            {"feature_id": genes}
            | {
                # Comparison group (trt_a) samples
                "trt_a1": [10.0, 10.0, 9.0, 4.0, 200.0],
                "trt_a2": [12.0, 12.0, 12.0, 5.0, 220.0],
                "trt_a3": [14.0, 14.0, 15.0, 6.0, 240.0],
                # Reference group (ref_b) samples
                "ref_b1": [100.0, 10.0, 10.0, 6.0, 100.0],
                "ref_b2": [110.0, 12.0, 11.0, 7.0, 110.0],
                "ref_b3": [120.0, 14.0, 15.0, 8.0, 120.0],
            },
        )

        omics = OmicsMatrix(
            matrix_id="mx-finite",
            schema_version="1.0.0",
            modality=Modality.TRANSCRIPTOMICS,
            feature_type="gene",
            measurement_type=MeasurementType.ESTIMATED_COUNTS,
            normalization_status=NormalizationStatus.RAW,
            n_features=len(genes),
            n_samples=len(all_samples),
            feature_ids=list(genes),
            sample_ids=list(all_samples),
            dataframe=df,
            created_at="2026-01-01T00:00:00Z",
        )

        design = ExperimentDesign(
            experiment_id="exp-finite-stats",
            description="Deterministic finite-statistics verification",
            samples=[
                DesignSample(sample_id="trt_a1", group_id="trt_a"),
                DesignSample(sample_id="trt_a2", group_id="trt_a"),
                DesignSample(sample_id="trt_a3", group_id="trt_a"),
                DesignSample(sample_id="ref_b1", group_id="ref_b"),
                DesignSample(sample_id="ref_b2", group_id="ref_b"),
                DesignSample(sample_id="ref_b3", group_id="ref_b"),
            ],
            groups=[
                ExperimentalGroup(
                    group_id="trt_a", label="Treatment A", role=GroupRole.TREATMENT
                ),
                ExperimentalGroup(
                    group_id="ref_b", label="Reference B", role=GroupRole.CONTROL
                ),
            ],
            factors=[
                ExperimentalFactor(
                    factor_id="condition",
                    factor_type=FactorType.CATEGORICAL,
                    levels=["treatment_a", "reference_b"],
                ),
            ],
            contrasts=[
                Contrast(
                    contrast_id="a_vs_b",
                    comparison_group_id="trt_a",
                    reference_group_id="ref_b",
                    description="Treatment A vs Reference B",
                ),
            ],
        )

        spec = AnalysisSpecification(
            analysis_type="differential_analysis",
            factor_references=["condition"],
            contrast_references=["a_vs_b"],
        )

        return omics, design, spec

    def test_all_p_values_finite(self) -> None:
        """Every gene receives a finite (non-NaN) raw p-value."""
        omics, design, spec = self._make_finite_fixture()
        result = run_analysis(spec, design, omics)
        for g in result.gene_results:
            assert math.isfinite(g.p_value), f"{g.gene_id} p_value is not finite"

    def test_all_adj_p_values_finite(self) -> None:
        """Every gene receives a finite (non-NaN) BH-adjusted p-value."""
        omics, design, spec = self._make_finite_fixture()
        result = run_analysis(spec, design, omics)
        for g in result.gene_results:
            assert math.isfinite(g.adj_p_value), (
                f"{g.gene_id} adj_p_value is not finite"
            )

    def test_significant_genes(self) -> None:
        """GeneA and GeneE are classified as significant (adj_p < 0.05)."""
        omics, design, spec = self._make_finite_fixture()
        result = run_analysis(spec, design, omics)

        gene_a = next(g for g in result.gene_results if g.gene_id == "GeneA")
        assert gene_a.significant is True
        assert gene_a.adj_p_value < 0.05

        gene_e = next(g for g in result.gene_results if g.gene_id == "GeneE")
        assert gene_e.significant is True
        assert gene_e.adj_p_value < 0.05

    def test_non_significant_genes(self) -> None:
        """GeneB, GeneC, GeneD are classified as non-significant."""
        omics, design, spec = self._make_finite_fixture()
        result = run_analysis(spec, design, omics)

        for gid in ("GeneB", "GeneC", "GeneD"):
            gene = next(g for g in result.gene_results if g.gene_id == gid)
            assert gene.significant is False, f"{gid} should not be significant"

    def test_n_genes_tested(self) -> None:
        """Result reports exactly 5 genes tested."""
        omics, design, spec = self._make_finite_fixture()
        result = run_analysis(spec, design, omics)
        assert result.n_genes_tested == 5
        assert len(result.gene_results) == 5

    def test_gene_order_preserved(self) -> None:
        """Gene results appear in original feature_id order."""
        omics, design, spec = self._make_finite_fixture()
        result = run_analysis(spec, design, omics)
        expected = ["GeneA", "GeneB", "GeneC", "GeneD", "GeneE"]
        assert [g.gene_id for g in result.gene_results] == expected

    def test_log2_fold_change_values(self) -> None:
        """Log2FC values match deterministic arithmetic."""
        omics, design, spec = self._make_finite_fixture()
        result = run_analysis(spec, design, omics)

        # GeneE: comp_mean=220, ref_mean=110 → log2(2)=1.0 (exact)
        gene_e = next(g for g in result.gene_results if g.gene_id == "GeneE")
        assert gene_e.log2_fold_change == pytest.approx(1.0, abs=1e-9)

        # GeneC: comp_mean=12, ref_mean=12 → log2(1)=0.0 (equal means from
        #   [9,12,15] vs [10,11,15], different per-sample values but same mean)
        gene_c = next(g for g in result.gene_results if g.gene_id == "GeneC")
        assert gene_c.log2_fold_change == pytest.approx(0.0, abs=1e-9)

        # GeneB: comp=[10,12,14], ref=[10,12,14] → identical → mean=12, base_mean=12
        gene_b = next(g for g in result.gene_results if g.gene_id == "GeneB")
        assert gene_b.log2_fold_change == pytest.approx(0.0, abs=1e-9)
        assert gene_b.base_mean == pytest.approx(12.0, abs=1e-9)

        # GeneA: comp_mean=12, ref_mean=110 → log2(12/110) ≈ -3.1964
        gene_a = next(g for g in result.gene_results if g.gene_id == "GeneA")
        assert gene_a.log2_fold_change == pytest.approx(-3.196397, abs=1e-4)
        assert gene_a.log2_fold_change < 0  # down-regulated (comp < ref)

        # GeneD: comp_mean=5, ref_mean=7 → log2(5/7) ≈ -0.485
        gene_d = next(g for g in result.gene_results if g.gene_id == "GeneD")
        assert gene_d.log2_fold_change == pytest.approx(-0.485427, abs=1e-4)

    def test_base_mean_values(self) -> None:
        """Base mean = (comp_mean + ref_mean) / 2 matches arithmetic."""
        omics, design, spec = self._make_finite_fixture()
        result = run_analysis(spec, design, omics)

        # GeneA: (12 + 110) / 2 = 61
        gene_a = next(g for g in result.gene_results if g.gene_id == "GeneA")
        assert gene_a.base_mean == pytest.approx(61.0, abs=1e-9)

        # GeneE: (220 + 110) / 2 = 165
        gene_e = next(g for g in result.gene_results if g.gene_id == "GeneE")
        assert gene_e.base_mean == pytest.approx(165.0, abs=1e-9)

        # GeneC: (12 + 12) / 2 = 12 (equal means from [9,12,15] vs [10,11,15])
        gene_c = next(g for g in result.gene_results if g.gene_id == "GeneC")
        assert gene_c.base_mean == pytest.approx(12.0, abs=1e-9)

        # GeneB: (12 + 12) / 2 = 12
        gene_b = next(g for g in result.gene_results if g.gene_id == "GeneB")
        assert gene_b.base_mean == pytest.approx(12.0, abs=1e-9)

        # GeneD: (5 + 7) / 2 = 6
        gene_d = next(g for g in result.gene_results if g.gene_id == "GeneD")
        assert gene_d.base_mean == pytest.approx(6.0, abs=1e-9)

    def test_p_value_magnitude_gene_a(self) -> None:
        """GeneA has a very small p-value due to huge fold-change."""
        omics, design, spec = self._make_finite_fixture()
        result = run_analysis(spec, design, omics)
        gene_a = next(g for g in result.gene_results if g.gene_id == "GeneA")
        # t ≈ 10–15 with low df → p in the 1e-3 range
        assert gene_a.p_value < 0.01
        assert math.isfinite(gene_a.p_value)

    def test_gene_c_p_value_is_one(self) -> None:
        """GeneC has identical group means → t=0 → p=1.0 (exact)."""
        omics, design, spec = self._make_finite_fixture()
        result = run_analysis(spec, design, omics)
        gene_c = next(g for g in result.gene_results if g.gene_id == "GeneC")
        assert gene_c.p_value == pytest.approx(1.0, abs=1e-9)
        assert gene_c.adj_p_value == pytest.approx(1.0, abs=1e-9)

    def test_gene_b_p_value_is_one(self) -> None:
        """GeneB has identical group means (same values, different order)
        → t=0 → p=1.0 (exact)."""
        omics, design, spec = self._make_finite_fixture()
        result = run_analysis(spec, design, omics)
        gene_b = next(g for g in result.gene_results if g.gene_id == "GeneB")
        assert gene_b.p_value == pytest.approx(1.0, abs=1e-9)
        assert gene_b.adj_p_value == pytest.approx(1.0, abs=1e-9)

    def test_no_warnings(self) -> None:
        """Successful finite-statistics run produces no warnings."""
        omics, design, spec = self._make_finite_fixture()
        result = run_analysis(spec, design, omics)
        assert result.warnings == ()

    def test_full_results_summary(self) -> None:
        """Print all gene results for manual inspection of the full pipeline."""
        omics, design, spec = self._make_finite_fixture()
        result = run_analysis(spec, design, omics)

        for g in result.gene_results:
            assert math.isfinite(g.log2_fold_change)
            assert math.isfinite(g.base_mean)
            assert math.isfinite(g.p_value)
            assert math.isfinite(g.adj_p_value)


class TestPartialNaNWarning:
    """Regression: when one gene has zero variance among otherwise finite
    genes, the analysis succeeds but a warning is emitted listing the
    affected gene(s)."""

    @staticmethod
    def _make_partial_nan_fixture():
        """Build (OmicsMatrix, ExperimentDesign, AnalysisSpecification) where
        GeneB has zero within-group variance (identical replicates) but
        all other genes have non-zero variance."""
        import pandas as pd

        from pharmomics.analysis.schemas import AnalysisSpecification
        from pharmomics.experiment.enums import FactorType, GroupRole
        from pharmomics.experiment.schemas import (
            Contrast,
            DesignSample,
            ExperimentalFactor,
            ExperimentalGroup,
            ExperimentDesign,
        )
        from pharmomics.omics.enums import (
            MeasurementType,
            Modality,
            NormalizationStatus,
        )
        from pharmomics.omics.schemas import OmicsMatrix

        genes = ["GeneA", "GeneB", "GeneC"]
        comp_samples = ["trt_1", "trt_2", "trt_3"]
        ref_samples = ["ref_1", "ref_2", "ref_3"]
        all_samples = comp_samples + ref_samples

        # GeneA: different values per replicate, large diff → finite p
        # GeneB: IDENTICAL replicates in both groups → zero variance → NaN
        # GeneC: different values per replicate, moderate diff → finite p
        df = pd.DataFrame(
            {"feature_id": genes}
            | {
                "trt_1": [10.0, 100.0, 50.0],
                "trt_2": [12.0, 100.0, 55.0],
                "trt_3": [14.0, 100.0, 45.0],
                "ref_1": [100.0, 200.0, 20.0],
                "ref_2": [110.0, 200.0, 25.0],
                "ref_3": [120.0, 200.0, 15.0],
            },
        )

        omics = OmicsMatrix(
            matrix_id="mx-partial-nan",
            schema_version="1.0.0",
            modality=Modality.TRANSCRIPTOMICS,
            feature_type="gene",
            measurement_type=MeasurementType.ESTIMATED_COUNTS,
            normalization_status=NormalizationStatus.RAW,
            n_features=len(genes),
            n_samples=len(all_samples),
            feature_ids=list(genes),
            sample_ids=list(all_samples),
            dataframe=df,
            created_at="2026-01-01T00:00:00Z",
        )

        design = ExperimentDesign(
            experiment_id="exp-partial-nan",
            description="Partial NaN verification: one gene with zero variance",
            samples=[
                DesignSample(sample_id="trt_1", group_id="trt"),
                DesignSample(sample_id="trt_2", group_id="trt"),
                DesignSample(sample_id="trt_3", group_id="trt"),
                DesignSample(sample_id="ref_1", group_id="ref"),
                DesignSample(sample_id="ref_2", group_id="ref"),
                DesignSample(sample_id="ref_3", group_id="ref"),
            ],
            groups=[
                ExperimentalGroup(
                    group_id="trt", label="Treatment", role=GroupRole.TREATMENT
                ),
                ExperimentalGroup(
                    group_id="ref", label="Reference", role=GroupRole.CONTROL
                ),
            ],
            factors=[
                ExperimentalFactor(
                    factor_id="condition",
                    factor_type=FactorType.CATEGORICAL,
                    levels=["treatment", "reference"],
                ),
            ],
            contrasts=[
                Contrast(
                    contrast_id="trt_vs_ref",
                    comparison_group_id="trt",
                    reference_group_id="ref",
                    description="Treatment vs Reference",
                ),
            ],
        )

        spec = AnalysisSpecification(
            analysis_type="differential_analysis",
            factor_references=["condition"],
            contrast_references=["trt_vs_ref"],
        )

        return omics, design, spec

    def test_analysis_succeeds_with_warning(self) -> None:
        """Analysis completes; a warning lists the NaN gene."""
        omics, design, spec = self._make_partial_nan_fixture()
        result = run_analysis(spec, design, omics)

        assert result.n_genes_tested == 3
        assert len(result.warnings) == 1
        assert "1 gene(s)" in result.warnings[0]
        assert "GeneB" in result.warnings[0]

    def test_nan_gene_is_not_significant(self) -> None:
        """The zero-variance gene is marked non-significant."""
        omics, design, spec = self._make_partial_nan_fixture()
        result = run_analysis(spec, design, omics)

        gene_b = next(g for g in result.gene_results if g.gene_id == "GeneB")
        assert math.isnan(gene_b.p_value)
        assert math.isnan(gene_b.adj_p_value)
        assert gene_b.significant is False

    def test_finite_genes_unchanged(self) -> None:
        """Non-zero-variance genes receive finite p-values."""
        omics, design, spec = self._make_partial_nan_fixture()
        result = run_analysis(spec, design, omics)

        for gid in ("GeneA", "GeneC"):
            gene = next(g for g in result.gene_results if g.gene_id == gid)
            assert math.isfinite(gene.p_value), f"{gid} p_value is not finite"
            assert math.isfinite(gene.adj_p_value), f"{gid} adj_p_value is not finite"
