"""Tests for pharmomics.analysis.mean_expression."""

from __future__ import annotations

import pandas as pd
import pytest

from pharmomics.analysis.contrast import ResolvedContrast
from pharmomics.analysis.matrix_slice import (
    DifferentialInput,
    prepare_differential_inputs,
)
from pharmomics.analysis.mean_expression import (
    GeneMeanExpression,
    MeanExpressionError,
    compute_mean_expression,
)
from pharmomics.omics.schemas import OmicsMatrix

# ---------------------------------------------------------------------------
# Test fixture helpers (reuse pattern from test_matrix_slice)
# ---------------------------------------------------------------------------


def _make_matrix(
    *,
    matrix_id: str = "mx-test-001",
    sample_ids: list[str] | None = None,
    feature_ids: list[str] | None = None,
    data: list[list[float]] | None = None,
) -> OmicsMatrix:
    """Create a minimal OmicsMatrix for testing.

    If ``data`` is provided, it must be a list of rows (one per feature),
    each row having exactly ``len(sample_ids)`` numeric values.
    Otherwise, deterministic defaults are generated.
    """
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
# compute_mean_expression — normal cases
# ---------------------------------------------------------------------------


class TestComputeMeanExpressionNormal:
    """Verify correct arithmetic mean computation."""

    def test_mean_expression_2x2(self) -> None:
        """2 genes, 2 comp samples, 2 ref samples — known values."""
        # G0: comp=[10.0, 20.0] -> mean=15.0; ref=[30.0, 40.0] -> mean=35.0
        # G1: comp=[50.0, 60.0] -> mean=55.0; ref=[70.0, 80.0] -> mean=75.0
        matrix = _make_matrix(
            sample_ids=["S0", "S1", "S2", "S3"],
            feature_ids=["G0", "G1"],
            data=[
                [10.0, 20.0, 30.0, 40.0],
                [50.0, 60.0, 70.0, 80.0],
            ],
        )
        contrast = _make_contrast(
            comparison_ids=["S0", "S1"],
            reference_ids=["S2", "S3"],
        )
        diff_input = prepare_differential_inputs(matrix, contrast)
        result = compute_mean_expression(diff_input)

        assert len(result) == 2
        assert result[0] == GeneMeanExpression(
            gene_id="G0", comparison_mean=15.0, reference_mean=35.0
        )
        assert result[1] == GeneMeanExpression(
            gene_id="G1", comparison_mean=55.0, reference_mean=75.0
        )

    def test_mean_expression_single_gene(self) -> None:
        """Single feature — returns one entry."""
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
        result = compute_mean_expression(diff_input)

        assert len(result) == 1
        assert result[0].gene_id == "G0"
        assert result[0].comparison_mean == 5.0  # (4+6)/2
        assert result[0].reference_mean == 10.0  # (8+12)/2

    def test_mean_expression_single_sample_per_group(self) -> None:
        """1 comp sample, 1 ref sample — mean equals the single value."""
        matrix = _make_matrix(
            sample_ids=["S0", "S1"],
            feature_ids=["G0", "G1"],
            data=[
                [3.5, 7.5],
                [1.0, 2.0],
            ],
        )
        contrast = _make_contrast(
            comparison_ids=["S0"],
            reference_ids=["S1"],
        )
        diff_input = prepare_differential_inputs(matrix, contrast)
        result = compute_mean_expression(diff_input)

        assert result[0].comparison_mean == 3.5
        assert result[0].reference_mean == 7.5
        assert result[1].comparison_mean == 1.0
        assert result[1].reference_mean == 2.0

    def test_mean_expression_unequal_group_sizes(self) -> None:
        """comp=2 samples, ref=3 samples — each group's mean is independent."""
        # G0: comp=[2.0, 4.0] -> 3.0; ref=[6.0, 8.0, 10.0] -> 8.0
        matrix = _make_matrix(
            sample_ids=["S0", "S1", "S2", "S3", "S4"],
            feature_ids=["G0"],
            data=[[2.0, 4.0, 6.0, 8.0, 10.0]],
        )
        contrast = _make_contrast(
            comparison_ids=["S0", "S1"],
            reference_ids=["S2", "S3", "S4"],
        )
        diff_input = prepare_differential_inputs(matrix, contrast)
        result = compute_mean_expression(diff_input)

        assert result[0].comparison_mean == 3.0
        assert result[0].reference_mean == 8.0

    def test_preserves_feature_order(self) -> None:
        """Result order must match original feature_ids order."""
        matrix = _make_matrix(
            sample_ids=["S0", "S1", "S2", "S3"],
            feature_ids=["G_Z", "G_A", "G_M"],
            data=[
                [1.0, 1.0, 1.0, 1.0],
                [2.0, 2.0, 2.0, 2.0],
                [3.0, 3.0, 3.0, 3.0],
            ],
        )
        contrast = _make_contrast()
        diff_input = prepare_differential_inputs(matrix, contrast)
        result = compute_mean_expression(diff_input)

        assert [r.gene_id for r in result] == ["G_Z", "G_A", "G_M"]

    def test_returns_tuple(self) -> None:
        """Return type must be tuple."""
        matrix = _make_matrix()
        contrast = _make_contrast()
        diff_input = prepare_differential_inputs(matrix, contrast)
        result = compute_mean_expression(diff_input)

        assert isinstance(result, tuple)

    def test_gene_mean_expression_is_frozen(self) -> None:
        """GeneMeanExpression instances must be immutable."""
        entry = GeneMeanExpression(
            gene_id="G0", comparison_mean=1.0, reference_mean=2.0
        )
        with pytest.raises((TypeError, AttributeError)):
            entry.gene_id = "changed"  # type: ignore[misc]

    def test_integration_with_prepare_differential_inputs(self) -> None:
        """End-to-end: OmicsMatrix -> ResolvedContrast -> slice -> mean."""
        matrix = _make_matrix(
            sample_ids=["ctrl_1", "ctrl_2", "trt_1", "trt_2"],
            feature_ids=["EGFR", "TP53"],
            data=[
                [100.0, 110.0, 200.0, 220.0],  # EGFR: ctrl=105, trt=210
                [50.0, 60.0, 45.0, 55.0],  # TP53: ctrl=55, trt=50
            ],
        )
        contrast = _make_contrast(
            comparison_ids=["trt_1", "trt_2"],
            reference_ids=["ctrl_1", "ctrl_2"],
        )
        diff_input = prepare_differential_inputs(matrix, contrast)
        result = compute_mean_expression(diff_input)

        assert len(result) == 2
        egfr = next(r for r in result if r.gene_id == "EGFR")
        tp53 = next(r for r in result if r.gene_id == "TP53")
        assert egfr.comparison_mean == 210.0
        assert egfr.reference_mean == 105.0
        assert tp53.comparison_mean == 50.0
        assert tp53.reference_mean == 55.0


# ---------------------------------------------------------------------------
# compute_mean_expression — error cases
# ---------------------------------------------------------------------------


class TestComputeMeanExpressionErrors:
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

        with pytest.raises(MeanExpressionError, match="Feature sets differ"):
            compute_mean_expression(diff_input)

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

        with pytest.raises(MeanExpressionError, match="no features"):
            compute_mean_expression(diff_input)

    def test_raises_no_sample_columns(self) -> None:
        """Slice with no sample columns (only feature ID column)."""
        from pharmomics.analysis.matrix_slice import MatrixSlice

        # Manually construct a slice with no sample columns
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

        with pytest.raises(MeanExpressionError, match="no sample columns"):
            compute_mean_expression(diff_input)

    def test_raises_non_numeric_values(self) -> None:
        """Non-numeric sample values cause an error."""
        # Build a matrix with string values in sample columns
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

        with pytest.raises(MeanExpressionError, match="non-numeric"):
            compute_mean_expression(diff_input)
