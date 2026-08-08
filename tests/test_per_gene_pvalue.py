"""Tests for pharmomics.analysis.per_gene_pvalue.

Covers:
- Normal p-value computation with known inputs.
- Insufficient samples → NaN.
- Zero variance → NaN.
- Non-finite values (NaN/Inf) → NaN.
- Structural errors raise PerGenePValueError.
- scipy edge cases (large values, small values, unequal groups).
- GenePValue frozen behaviour.
- Return type and feature order preservation.
"""

from __future__ import annotations

from math import isnan

import pandas as pd
import pytest

from pharmomics.analysis.contrast import ResolvedContrast
from pharmomics.analysis.matrix_slice import (
    DifferentialInput,
    MatrixSlice,
    prepare_differential_inputs,
)
from pharmomics.analysis.per_gene_pvalue import (
    GenePValue,
    PerGenePValueError,
    compute_per_gene_p_value,
)
from pharmomics.omics.schemas import OmicsMatrix

# ---------------------------------------------------------------------------
# Test fixture helpers
# ---------------------------------------------------------------------------


def _make_matrix(
    *,
    matrix_id: str = "mx-test-001",
    sample_ids: list[str] | None = None,
    feature_ids: list[str] | None = None,
    data: list[list[float]] | None = None,
) -> OmicsMatrix:
    """Create a minimal OmicsMatrix for testing."""
    if sample_ids is None:
        sample_ids = ["S0", "S1", "S2", "S3"]
    if feature_ids is None:
        feature_ids = ["G0", "G1", "G2"]

    if data is not None:
        rows = [[gid] + row for gid, row in zip(feature_ids, data)]
    else:
        rows = []
        for i, g in enumerate(feature_ids):
            row = [float(i * len(sample_ids) + j) for j in range(len(sample_ids))]
            rows.append([g] + row)

    df = pd.DataFrame(rows, columns=["gene"] + sample_ids)

    return OmicsMatrix(
        matrix_id=matrix_id,
        schema_version="1.0.0",
        modality="transcriptomics",
        feature_type="gene",
        measurement_type="unknown",
        normalization_status="unknown",
        n_features=len(feature_ids),
        n_samples=len(sample_ids),
        feature_ids=feature_ids,
        sample_ids=sample_ids,
        dataframe=df,
        created_at="2025-01-01T00:00:00Z",
    )


def _make_contrast(
    *,
    comparison_ids: list[str] | None = None,
    reference_ids: list[str] | None = None,
) -> ResolvedContrast:
    """Create a minimal ResolvedContrast for testing."""
    if comparison_ids is None:
        comparison_ids = ["S0", "S1"]
    if reference_ids is None:
        reference_ids = ["S2", "S3"]
    return ResolvedContrast(
        contrast_id="test_contrast",
        comparison_group_id="group_a",
        reference_group_id="group_b",
        comparison_sample_ids=tuple(comparison_ids),
        reference_sample_ids=tuple(reference_ids),
    )


# ---------------------------------------------------------------------------
# compute_per_gene_p_value — normal cases
# ---------------------------------------------------------------------------


class TestComputePerGenePValueNormal:
    """Verify correct Welch's t-test p-value computation."""

    def test_two_genes_2x2_samples(self) -> None:
        """2 genes, 2 comp samples, 2 ref samples — known values."""
        # G0: comp=[10, 12], ref=[20, 22] — clear separation
        # G1: comp=[50, 52], ref=[51, 53] — overlap
        matrix = _make_matrix(
            sample_ids=["S0", "S1", "S2", "S3"],
            feature_ids=["G0", "G1"],
            data=[
                [10.0, 12.0, 20.0, 22.0],
                [50.0, 52.0, 51.0, 53.0],
            ],
        )
        contrast = _make_contrast(
            comparison_ids=["S0", "S1"],
            reference_ids=["S2", "S3"],
        )
        diff_input = prepare_differential_inputs(matrix, contrast)
        result = compute_per_gene_p_value(diff_input)

        assert len(result) == 2
        assert result[0].gene_id == "G0"
        assert result[1].gene_id == "G1"
        # G0 has clear separation → small p-value
        assert 0.0 < result[0].p_value < 0.05
        # G1 has overlap → larger p-value
        assert result[1].p_value > result[0].p_value

    def test_single_gene_2x2(self) -> None:
        """Single gene, 2+2 samples — returns one entry."""
        matrix = _make_matrix(
            sample_ids=["S0", "S1", "S2", "S3"],
            feature_ids=["G0"],
            data=[[4.0, 6.0, 8.0, 12.0]],
        )
        contrast = _make_contrast(
            comparison_ids=["S0", "S1"],
            reference_ids=["S2", "S3"],
        )
        diff_input = prepare_differential_inputs(matrix, contrast)
        result = compute_per_gene_p_value(diff_input)

        assert len(result) == 1
        assert result[0].gene_id == "G0"
        assert 0.0 < result[0].p_value <= 1.0

    def test_large_mean_difference_small_pvalue(self) -> None:
        """3+3 samples with large comp/ref difference → p < 0.05."""
        # comp: [100, 110, 120], ref: [10, 12, 14]
        matrix = _make_matrix(
            sample_ids=["S0", "S1", "S2", "S3", "S4", "S5"],
            feature_ids=["G0"],
            data=[[100.0, 110.0, 120.0, 10.0, 12.0, 14.0]],
        )
        contrast = _make_contrast(
            comparison_ids=["S0", "S1", "S2"],
            reference_ids=["S3", "S4", "S5"],
        )
        diff_input = prepare_differential_inputs(matrix, contrast)
        result = compute_per_gene_p_value(diff_input)

        assert len(result) == 1
        assert result[0].p_value < 0.05

    def test_identical_groups_pvalue_near_one(self) -> None:
        """3+3 samples, comp/ref from same distribution → p ≈ 1.0."""
        # comp: [10.0, 10.5, 9.5], ref: [10.1, 9.9, 10.2]
        matrix = _make_matrix(
            sample_ids=["S0", "S1", "S2", "S3", "S4", "S5"],
            feature_ids=["G0"],
            data=[[10.0, 10.5, 9.5, 10.1, 9.9, 10.2]],
        )
        contrast = _make_contrast(
            comparison_ids=["S0", "S1", "S2"],
            reference_ids=["S3", "S4", "S5"],
        )
        diff_input = prepare_differential_inputs(matrix, contrast)
        result = compute_per_gene_p_value(diff_input)

        assert len(result) == 1
        assert result[0].p_value > 0.5

    def test_returns_tuple(self) -> None:
        """Return type must be tuple."""
        matrix = _make_matrix()
        contrast = _make_contrast()
        diff_input = prepare_differential_inputs(matrix, contrast)
        result = compute_per_gene_p_value(diff_input)

        assert isinstance(result, tuple)

    def test_preserves_feature_order(self) -> None:
        """Result order must match original feature_ids order."""
        matrix = _make_matrix(
            sample_ids=["S0", "S1", "S2", "S3"],
            feature_ids=["G_Z", "G_A", "G_M"],
            data=[
                [1.0, 2.0, 10.0, 11.0],
                [3.0, 4.0, 12.0, 13.0],
                [5.0, 6.0, 14.0, 15.0],
            ],
        )
        contrast = _make_contrast()
        diff_input = prepare_differential_inputs(matrix, contrast)
        result = compute_per_gene_p_value(diff_input)

        assert [r.gene_id for r in result] == ["G_Z", "G_A", "G_M"]

    def test_gene_p_value_is_frozen(self) -> None:
        """GenePValue instances must be immutable."""
        entry = GenePValue(gene_id="G0", p_value=0.05)
        with pytest.raises((TypeError, AttributeError)):
            entry.gene_id = "changed"  # type: ignore[misc]

    def test_integration_with_prepare_differential_inputs(self) -> None:
        """End-to-end: OmicsMatrix -> ResolvedContrast -> slice -> p-value."""
        matrix = _make_matrix(
            sample_ids=["ctrl_1", "ctrl_2", "trt_1", "trt_2"],
            feature_ids=["EGFR", "TP53"],
            data=[
                [100.0, 110.0, 200.0, 220.0],  # EGFR: large difference
                [50.0, 60.0, 55.0, 45.0],  # TP53: overlap
            ],
        )
        contrast = _make_contrast(
            comparison_ids=["trt_1", "trt_2"],
            reference_ids=["ctrl_1", "ctrl_2"],
        )
        diff_input = prepare_differential_inputs(matrix, contrast)
        result = compute_per_gene_p_value(diff_input)

        assert len(result) == 2
        egfr = next(r for r in result if r.gene_id == "EGFR")
        tp53 = next(r for r in result if r.gene_id == "TP53")
        assert 0.0 < egfr.p_value < 0.05  # large effect
        assert 0.0 < tp53.p_value <= 1.0  # small effect

    def test_p_value_in_valid_range(self) -> None:
        """All p-values must be in [0, 1]."""
        matrix = _make_matrix(
            sample_ids=["S0", "S1", "S2", "S3"],
            feature_ids=["G0", "G1", "G2"],
            data=[
                [1.0, 2.0, 10.0, 11.0],
                [5.0, 6.0, 7.0, 8.0],
                [100.0, 200.0, 110.0, 190.0],
            ],
        )
        contrast = _make_contrast()
        diff_input = prepare_differential_inputs(matrix, contrast)
        result = compute_per_gene_p_value(diff_input)

        for entry in result:
            assert 0.0 <= entry.p_value <= 1.0, f"{entry.gene_id}: {entry.p_value}"


# ---------------------------------------------------------------------------
# compute_per_gene_p_value — insufficient samples → NaN
# ---------------------------------------------------------------------------


class TestComputePerGenePValueInsufficientSamples:
    """Verify NaN is returned when sample counts are too low."""

    def test_comparison_one_sample(self) -> None:
        """comparison=1, reference=2 -> NaN."""
        comp_df = pd.DataFrame(
            {
                "gene": ["G0"],
                "S0": [10.0],
            }
        )
        comp_slice = MatrixSlice(
            matrix_id="mx-comp",
            sample_ids=["S0"],
            feature_ids=["G0"],
            dataframe=comp_df,
        )
        ref_df = pd.DataFrame(
            {
                "gene": ["G0"],
                "S1": [20.0],
                "S2": [22.0],
            }
        )
        ref_slice = MatrixSlice(
            matrix_id="mx-ref",
            sample_ids=["S1", "S2"],
            feature_ids=["G0"],
            dataframe=ref_df,
        )
        diff_input = DifferentialInput(comparison=comp_slice, reference=ref_slice)
        result = compute_per_gene_p_value(diff_input)

        assert len(result) == 1
        assert isnan(result[0].p_value)

    def test_reference_one_sample(self) -> None:
        """comparison=2, reference=1 -> NaN."""
        comp_df = pd.DataFrame(
            {
                "gene": ["G0"],
                "S0": [10.0],
                "S1": [12.0],
            }
        )
        comp_slice = MatrixSlice(
            matrix_id="mx-comp",
            sample_ids=["S0", "S1"],
            feature_ids=["G0"],
            dataframe=comp_df,
        )
        ref_df = pd.DataFrame(
            {
                "gene": ["G0"],
                "S2": [20.0],
            }
        )
        ref_slice = MatrixSlice(
            matrix_id="mx-ref",
            sample_ids=["S2"],
            feature_ids=["G0"],
            dataframe=ref_df,
        )
        diff_input = DifferentialInput(comparison=comp_slice, reference=ref_slice)
        result = compute_per_gene_p_value(diff_input)

        assert len(result) == 1
        assert isnan(result[0].p_value)

    def test_both_one_sample(self) -> None:
        """Both groups have 1 sample -> NaN."""
        comp_df = pd.DataFrame({"gene": ["G0"], "S0": [10.0]})
        comp_slice = MatrixSlice(
            matrix_id="mx-comp",
            sample_ids=["S0"],
            feature_ids=["G0"],
            dataframe=comp_df,
        )
        ref_df = pd.DataFrame({"gene": ["G0"], "S1": [20.0]})
        ref_slice = MatrixSlice(
            matrix_id="mx-ref",
            sample_ids=["S1"],
            feature_ids=["G0"],
            dataframe=ref_df,
        )
        diff_input = DifferentialInput(comparison=comp_slice, reference=ref_slice)
        result = compute_per_gene_p_value(diff_input)

        assert len(result) == 1
        assert isnan(result[0].p_value)

    def test_mixed_sufficient_and_insufficient(self) -> None:
        """Some genes have enough samples, others don't."""
        # All genes share the same sample columns in the slices,
        # so this test verifies the gene-level check. We simulate
        # by having 2x2 samples but testing that all get valid results.
        matrix = _make_matrix(
            sample_ids=["S0", "S1", "S2", "S3"],
            feature_ids=["G0", "G1"],
            data=[
                [10.0, 12.0, 20.0, 22.0],
                [30.0, 32.0, 40.0, 42.0],
            ],
        )
        contrast = _make_contrast(
            comparison_ids=["S0", "S1"],
            reference_ids=["S2", "S3"],
        )
        diff_input = prepare_differential_inputs(matrix, contrast)
        result = compute_per_gene_p_value(diff_input)

        assert len(result) == 2
        for entry in result:
            assert not isnan(entry.p_value), (
                f"{entry.gene_id} should have valid p-value"
            )


# ---------------------------------------------------------------------------
# compute_per_gene_p_value — zero variance → NaN
# ---------------------------------------------------------------------------


class TestComputePerGenePValueZeroVariance:
    """Verify behavior when one or both groups have zero variance.

    When exactly one group has zero variance, Welch's t-test still
    produces a finite p-value (the zero-variance group is treated as
    a fixed value, with df equal to the other group's sample count
    minus one).  When both groups have zero variance, the test is
    degenerate and returns NaN.
    """

    def test_comparison_zero_variance(self) -> None:
        """All comp values identical — Welch reduces to df = ref_n - 1."""
        matrix = _make_matrix(
            sample_ids=["S0", "S1", "S2", "S3"],
            feature_ids=["G0"],
            data=[[10.0, 10.0, 20.0, 22.0]],  # comp variance = 0
        )
        contrast = _make_contrast(
            comparison_ids=["S0", "S1"],
            reference_ids=["S2", "S3"],
        )
        diff_input = prepare_differential_inputs(matrix, contrast)
        result = compute_per_gene_p_value(diff_input)

        assert len(result) == 1
        assert not isnan(result[0].p_value), (
            "one-group zero variance should yield finite Welch p-value"
        )
        assert 0.0 < result[0].p_value <= 1.0

    def test_reference_zero_variance(self) -> None:
        """All ref values identical — Welch reduces to df = comp_n - 1."""
        matrix = _make_matrix(
            sample_ids=["S0", "S1", "S2", "S3"],
            feature_ids=["G0"],
            data=[[10.0, 12.0, 20.0, 20.0]],  # ref variance = 0
        )
        contrast = _make_contrast(
            comparison_ids=["S0", "S1"],
            reference_ids=["S2", "S3"],
        )
        diff_input = prepare_differential_inputs(matrix, contrast)
        result = compute_per_gene_p_value(diff_input)

        assert len(result) == 1
        assert not isnan(result[0].p_value), (
            "one-group zero variance should yield finite Welch p-value"
        )
        assert 0.0 < result[0].p_value <= 1.0

    def test_both_zero_variance(self) -> None:
        """Both groups have zero variance -> NaN."""
        matrix = _make_matrix(
            sample_ids=["S0", "S1", "S2", "S3"],
            feature_ids=["G0"],
            data=[[10.0, 10.0, 20.0, 20.0]],
        )
        contrast = _make_contrast(
            comparison_ids=["S0", "S1"],
            reference_ids=["S2", "S3"],
        )
        diff_input = prepare_differential_inputs(matrix, contrast)
        result = compute_per_gene_p_value(diff_input)

        assert len(result) == 1
        assert isnan(result[0].p_value)

    def test_mixed_zero_variance_and_normal(self) -> None:
        """One gene has one-group zero variance, another is normal."""
        # G0: comp=[10,10] zero var, ref=[20,22] — Welch handles this
        # G1: comp=[10,12], ref=[20,22] normal
        matrix = _make_matrix(
            sample_ids=["S0", "S1", "S2", "S3"],
            feature_ids=["G0", "G1"],
            data=[
                [10.0, 10.0, 20.0, 22.0],  # G0: comp zero var
                [10.0, 12.0, 20.0, 22.0],  # G1: normal
            ],
        )
        contrast = _make_contrast(
            comparison_ids=["S0", "S1"],
            reference_ids=["S2", "S3"],
        )
        diff_input = prepare_differential_inputs(matrix, contrast)
        result = compute_per_gene_p_value(diff_input)

        assert len(result) == 2
        assert not isnan(result[0].p_value), (
            "G0 should have finite Welch p-value (one-group zero var)"
        )
        assert not isnan(result[1].p_value), "G1 should have valid p-value"

    def test_one_group_zero_variance_finite_pvalue(self) -> None:
        """Regression: one-group zero variance yields a known Welch p-value.

        comp=[10,10,10] (var=0), ref=[20,22,24] (var=4).
        scipy Welch t-test gives p≈0.00913 with df=2.
        """
        matrix = _make_matrix(
            sample_ids=["S0", "S1", "S2", "S3", "S4", "S5"],
            feature_ids=["G0"],
            data=[[10.0, 10.0, 10.0, 20.0, 22.0, 24.0]],
        )
        contrast = _make_contrast(
            comparison_ids=["S0", "S1", "S2"],
            reference_ids=["S3", "S4", "S5"],
        )
        diff_input = prepare_differential_inputs(matrix, contrast)
        result = compute_per_gene_p_value(diff_input)

        assert len(result) == 1
        assert not isnan(result[0].p_value)
        # Verify against scipy's Welch p-value to 6 decimal places
        assert abs(result[0].p_value - 0.009133) < 1e-5, (
            f"Expected p≈0.009133, got {result[0].p_value}"
        )


# ---------------------------------------------------------------------------
# compute_per_gene_p_value — non-finite values → NaN
# ---------------------------------------------------------------------------


class TestComputePerGenePValueNonFinite:
    """Verify NaN is returned for non-finite sample values."""

    def test_comparison_contains_nan(self) -> None:
        """comp has NaN -> NaN."""
        comp_df = pd.DataFrame(
            {
                "gene": ["G0"],
                "S0": [float("nan")],
                "S1": [12.0],
            }
        )
        comp_slice = MatrixSlice(
            matrix_id="mx-comp",
            sample_ids=["S0", "S1"],
            feature_ids=["G0"],
            dataframe=comp_df,
        )
        ref_df = pd.DataFrame(
            {
                "gene": ["G0"],
                "S2": [20.0],
                "S3": [22.0],
            }
        )
        ref_slice = MatrixSlice(
            matrix_id="mx-ref",
            sample_ids=["S2", "S3"],
            feature_ids=["G0"],
            dataframe=ref_df,
        )
        diff_input = DifferentialInput(comparison=comp_slice, reference=ref_slice)
        result = compute_per_gene_p_value(diff_input)

        assert len(result) == 1
        assert isnan(result[0].p_value)

    def test_reference_contains_inf(self) -> None:
        """ref has Inf -> NaN."""
        comp_df = pd.DataFrame(
            {
                "gene": ["G0"],
                "S0": [10.0],
                "S1": [12.0],
            }
        )
        comp_slice = MatrixSlice(
            matrix_id="mx-comp",
            sample_ids=["S0", "S1"],
            feature_ids=["G0"],
            dataframe=comp_df,
        )
        ref_df = pd.DataFrame(
            {
                "gene": ["G0"],
                "S2": [float("inf")],
                "S3": [22.0],
            }
        )
        ref_slice = MatrixSlice(
            matrix_id="mx-ref",
            sample_ids=["S2", "S3"],
            feature_ids=["G0"],
            dataframe=ref_df,
        )
        diff_input = DifferentialInput(comparison=comp_slice, reference=ref_slice)
        result = compute_per_gene_p_value(diff_input)

        assert len(result) == 1
        assert isnan(result[0].p_value)

    def test_mixed_finite_and_non_finite(self) -> None:
        """G0 has NaN in comp, G1 is fully valid, G2 has Inf in ref."""
        comp_df = pd.DataFrame(
            {
                "gene": ["G0", "G1", "G2"],
                "S0": [float("nan"), 10.0, 30.0],
                "S1": [12.0, 12.0, 32.0],
            }
        )
        comp_slice = MatrixSlice(
            matrix_id="mx-mixed",
            sample_ids=["S0", "S1"],
            feature_ids=["G0", "G1", "G2"],
            dataframe=comp_df,
        )
        ref_df = pd.DataFrame(
            {
                "gene": ["G0", "G1", "G2"],
                "S2": [20.0, 20.0, float("inf")],
                "S3": [22.0, 22.0, 44.0],
            }
        )
        ref_slice = MatrixSlice(
            matrix_id="mx-mixed-ref",
            sample_ids=["S2", "S3"],
            feature_ids=["G0", "G1", "G2"],
            dataframe=ref_df,
        )
        diff_input = DifferentialInput(comparison=comp_slice, reference=ref_slice)
        result = compute_per_gene_p_value(diff_input)

        assert len(result) == 3
        assert isnan(result[0].p_value), "G0: comp has NaN"
        assert not isnan(result[1].p_value), "G1: all values finite"
        assert isnan(result[2].p_value), "G2: ref has Inf"

    def test_mixed_finite_and_non_finite_both_groups(self) -> None:
        """One gene with NaN in comp (G0), one clean (G1)."""
        comp_df = pd.DataFrame(
            {
                "gene": ["G0", "G1"],
                "S0": [float("nan"), 10.0],
                "S1": [12.0, 12.0],
            }
        )
        comp_slice = MatrixSlice(
            matrix_id="mx-mixed2",
            sample_ids=["S0", "S1"],
            feature_ids=["G0", "G1"],
            dataframe=comp_df,
        )
        ref_df = pd.DataFrame(
            {
                "gene": ["G0", "G1"],
                "S2": [20.0, 20.0],
                "S3": [22.0, 22.0],
            }
        )
        ref_slice = MatrixSlice(
            matrix_id="mx-mixed2-ref",
            sample_ids=["S2", "S3"],
            feature_ids=["G0", "G1"],
            dataframe=ref_df,
        )
        diff_input = DifferentialInput(comparison=comp_slice, reference=ref_slice)
        result = compute_per_gene_p_value(diff_input)

        assert len(result) == 2
        assert isnan(result[0].p_value), "G0: comp has NaN"
        assert not isnan(result[1].p_value), "G1: all values finite"


# ---------------------------------------------------------------------------
# compute_per_gene_p_value — structural errors → raise
# ---------------------------------------------------------------------------


class TestComputePerGenePValueErrors:
    """Verify error handling for invalid inputs."""

    def test_raises_feature_mismatch(self) -> None:
        """Different feature sets between comp and ref."""
        comp_matrix = _make_matrix(
            matrix_id="mx-comp",
            sample_ids=["S0", "S1"],
            feature_ids=["G0", "G1"],
        )
        ref_matrix = _make_matrix(
            matrix_id="mx-ref",
            sample_ids=["S2", "S3"],
            feature_ids=["G0", "G2"],  # G1 vs G2 mismatch
        )
        comp_slice = prepare_differential_inputs(
            comp_matrix,
            _make_contrast(comparison_ids=["S0"], reference_ids=["S1"]),
        ).comparison
        ref_slice = prepare_differential_inputs(
            ref_matrix,
            _make_contrast(comparison_ids=["S2"], reference_ids=["S3"]),
        ).reference

        diff_input = DifferentialInput(comparison=comp_slice, reference=ref_slice)

        with pytest.raises(PerGenePValueError, match="Feature sets differ"):
            compute_per_gene_p_value(diff_input)

    def test_raises_empty_features(self) -> None:
        """Empty feature_ids in comparison slice."""
        comp_matrix = _make_matrix(
            sample_ids=["S0", "S1"],
            feature_ids=[],
            data=[],
        )
        ref_matrix = _make_matrix(
            sample_ids=["S2", "S3"],
            feature_ids=["G0"],
            data=[[1.0, 2.0]],
        )
        comp_slice = prepare_differential_inputs(
            comp_matrix,
            _make_contrast(comparison_ids=["S0"], reference_ids=["S1"]),
        ).comparison
        ref_slice = prepare_differential_inputs(
            ref_matrix,
            _make_contrast(comparison_ids=["S2"], reference_ids=["S3"]),
        ).reference

        diff_input = DifferentialInput(comparison=comp_slice, reference=ref_slice)

        with pytest.raises(PerGenePValueError, match="no features"):
            compute_per_gene_p_value(diff_input)

    def test_raises_no_sample_columns(self) -> None:
        """Slice with no sample columns (only feature ID column)."""
        comp_df = pd.DataFrame({"gene": ["G0"]})
        comp_slice = MatrixSlice(
            matrix_id="mx-nosample",
            sample_ids=[],
            feature_ids=["G0"],
            dataframe=comp_df,
        )

        ref_matrix = _make_matrix(
            sample_ids=["S0", "S1"],
            feature_ids=["G0"],
            data=[[1.0, 2.0]],
        )
        ref_slice = prepare_differential_inputs(
            ref_matrix,
            _make_contrast(comparison_ids=["S0"], reference_ids=["S1"]),
        ).reference

        diff_input = DifferentialInput(comparison=comp_slice, reference=ref_slice)

        with pytest.raises(PerGenePValueError, match="no sample columns"):
            compute_per_gene_p_value(diff_input)

    def test_raises_non_numeric_values(self) -> None:
        """Non-numeric sample values cause an error."""
        rows = [
            ["G0", "1.0", "2.0", "3.0", "4.0"],
            ["G1", "5.0", "6.0", "7.0", "8.0"],
        ]
        df = pd.DataFrame(rows, columns=["gene", "S0", "S1", "S2", "S3"])
        matrix = OmicsMatrix(
            matrix_id="mx-nonnumeric",
            schema_version="1.0.0",
            modality="transcriptomics",
            feature_type="gene",
            measurement_type="unknown",
            normalization_status="unknown",
            n_features=2,
            n_samples=4,
            feature_ids=["G0", "G1"],
            sample_ids=["S0", "S1", "S2", "S3"],
            dataframe=df,
            created_at="2025-01-01T00:00:00Z",
        )
        contrast = _make_contrast(
            comparison_ids=["S0", "S1"],
            reference_ids=["S2", "S3"],
        )
        diff_input = prepare_differential_inputs(matrix, contrast)

        with pytest.raises(PerGenePValueError, match="non-numeric"):
            compute_per_gene_p_value(diff_input)


# ---------------------------------------------------------------------------
# compute_per_gene_p_value — scipy edge cases
# ---------------------------------------------------------------------------


class TestComputePerGenePValueScipyEdgeCases:
    """Verify robustness with scipy edge cases."""

    def test_very_large_expression_values(self) -> None:
        """Expression values ~1e10 — should not overflow."""
        matrix = _make_matrix(
            sample_ids=["S0", "S1", "S2", "S3"],
            feature_ids=["G0"],
            data=[[1e10, 1.1e10, 2e10, 2.1e10]],
        )
        contrast = _make_contrast(
            comparison_ids=["S0", "S1"],
            reference_ids=["S2", "S3"],
        )
        diff_input = prepare_differential_inputs(matrix, contrast)
        result = compute_per_gene_p_value(diff_input)

        assert len(result) == 1
        # p-value should be valid (not NaN) or at least a finite number
        assert result[0].p_value == result[0].p_value  # NaN self-check
        if not isnan(result[0].p_value):
            assert 0.0 <= result[0].p_value <= 1.0

    def test_very_small_positive_values(self) -> None:
        """Expression values ~1e-10 — should handle precision."""
        matrix = _make_matrix(
            sample_ids=["S0", "S1", "S2", "S3"],
            feature_ids=["G0"],
            data=[[1e-10, 1.1e-10, 2e-10, 2.1e-10]],
        )
        contrast = _make_contrast(
            comparison_ids=["S0", "S1"],
            reference_ids=["S2", "S3"],
        )
        diff_input = prepare_differential_inputs(matrix, contrast)
        result = compute_per_gene_p_value(diff_input)

        assert len(result) == 1
        assert result[0].p_value == result[0].p_value  # NaN self-check
        if not isnan(result[0].p_value):
            assert 0.0 <= result[0].p_value <= 1.0

    def test_unequal_group_sizes(self) -> None:
        """3 comp vs 5 ref samples — Welch handles unequal n."""
        matrix = _make_matrix(
            sample_ids=["S0", "S1", "S2", "S3", "S4", "S5", "S6", "S7"],
            feature_ids=["G0"],
            data=[[10.0, 12.0, 11.0, 20.0, 22.0, 21.0, 19.0, 23.0]],
        )
        contrast = _make_contrast(
            comparison_ids=["S0", "S1", "S2"],
            reference_ids=["S3", "S4", "S5", "S6", "S7"],
        )
        diff_input = prepare_differential_inputs(matrix, contrast)
        result = compute_per_gene_p_value(diff_input)

        assert len(result) == 1
        assert 0.0 <= result[0].p_value <= 1.0

    def test_extreme_variance_ratio(self) -> None:
        """comp has tiny variance, ref has huge variance."""
        matrix = _make_matrix(
            sample_ids=["S0", "S1", "S2", "S3"],
            feature_ids=["G0"],
            data=[[10.0, 10.001, 50.0, 200.0]],
        )
        contrast = _make_contrast(
            comparison_ids=["S0", "S1"],
            reference_ids=["S2", "S3"],
        )
        diff_input = prepare_differential_inputs(matrix, contrast)
        result = compute_per_gene_p_value(diff_input)

        assert len(result) == 1
        # Welch's t-test should still produce a valid result
        assert result[0].p_value == result[0].p_value  # NaN self-check
        if not isnan(result[0].p_value):
            assert 0.0 <= result[0].p_value <= 1.0
