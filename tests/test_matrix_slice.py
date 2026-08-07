"""Tests for pharmomics.analysis.matrix_slice."""

from __future__ import annotations

import pandas as pd
import pytest

from pharmomics.analysis.contrast import ResolvedContrast
from pharmomics.analysis.matrix_slice import (
    DifferentialInput,
    MatrixSlice,
    MatrixSliceError,
    _extract_slice,
    prepare_differential_inputs,
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
) -> OmicsMatrix:
    """Create a minimal OmicsMatrix for testing."""
    if sample_ids is None:
        sample_ids = ["S0", "S1", "S2", "S3"]
    if feature_ids is None:
        feature_ids = ["G0", "G1", "G2"]

    rows = []
    for i, g in enumerate(feature_ids):
        row = [g] + [float(i * len(sample_ids) + j) for j in range(len(sample_ids))]
        rows.append(row)

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
# _extract_slice — normal cases
# ---------------------------------------------------------------------------


class TestExtractSliceNormal:
    """Verify correct extraction of sample sub-matrices."""

    def test_extract_two_samples(self) -> None:
        matrix = _make_matrix(sample_ids=["S0", "S1", "S2", "S3"])
        result = _extract_slice(matrix, ["S0", "S2"])

        assert result.sample_ids == ["S0", "S2"]
        assert result.dataframe.shape == (3, 3)  # 3 features, gene + 2 samples
        assert list(result.dataframe.columns) == ["gene", "S0", "S2"]

    def test_extract_single_sample(self) -> None:
        matrix = _make_matrix(sample_ids=["S0", "S1", "S2"])
        result = _extract_slice(matrix, ["S1"])

        assert result.sample_ids == ["S1"]
        assert result.dataframe.shape == (3, 2)  # 3 features, gene + 1 sample
        assert list(result.dataframe.columns) == ["gene", "S1"]

    def test_extract_all_samples(self) -> None:
        matrix = _make_matrix(sample_ids=["A", "B", "C"])
        result = _extract_slice(matrix, ["A", "B", "C"])

        assert result.sample_ids == ["A", "B", "C"]
        assert result.dataframe.shape == (3, 4)  # 3 features, gene + 3 samples
        assert list(result.dataframe.columns) == ["gene", "A", "B", "C"]

    def test_preserves_column_order(self) -> None:
        """Requested order ['S2', 'S0'] should return matrix order ['S0', 'S2']."""
        matrix = _make_matrix(sample_ids=["S0", "S1", "S2", "S3"])
        result = _extract_slice(matrix, ["S2", "S0"])

        assert result.sample_ids == ["S0", "S2"]
        assert list(result.dataframe.columns) == ["gene", "S0", "S2"]

    def test_preserves_all_features(self) -> None:
        matrix = _make_matrix(feature_ids=["G0", "G1", "G2"])
        result = _extract_slice(matrix, ["S0"])

        assert result.feature_ids == ["G0", "G1", "G2"]
        assert len(result.dataframe) == 3
        assert result.dataframe.iloc[:, 0].tolist() == ["G0", "G1", "G2"]

    def test_preserves_feature_id_column(self) -> None:
        """Column 0 must be the feature ID column."""
        matrix = _make_matrix()
        result = _extract_slice(matrix, ["S0"])

        assert result.dataframe.columns[0] == "gene"

    def test_numeric_values_preserved(self) -> None:
        """Extracted values must match the original matrix."""
        matrix = _make_matrix(sample_ids=["S0", "S1"], feature_ids=["G0", "G1"])
        # G0: S0=0.0, S1=1.0
        # G1: S0=2.0, S1=3.0
        result = _extract_slice(matrix, ["S1"])

        assert result.dataframe.iloc[0, 1] == 1.0
        assert result.dataframe.iloc[1, 1] == 3.0

    def test_matrix_id_propagates(self) -> None:
        matrix = _make_matrix(matrix_id="mx-abc-000")
        result = _extract_slice(matrix, ["S0"])

        assert result.matrix_id == "mx-abc-000"


# ---------------------------------------------------------------------------
# _extract_slice — error cases
# ---------------------------------------------------------------------------


class TestExtractSliceErrors:
    """Verify error handling for invalid slice requests."""

    def test_raises_unknown_sample_id(self) -> None:
        matrix = _make_matrix(sample_ids=["S0", "S1"])
        with pytest.raises(MatrixSliceError, match="not found"):
            _extract_slice(matrix, ["S0", "S99"])

    def test_raises_empty_sample_ids(self) -> None:
        matrix = _make_matrix(sample_ids=["S0", "S1"])
        with pytest.raises(MatrixSliceError, match="No sample IDs"):
            _extract_slice(matrix, [])

    def test_raises_all_missing(self) -> None:
        matrix = _make_matrix(sample_ids=["S0", "S1"])
        with pytest.raises(MatrixSliceError, match="not found"):
            _extract_slice(matrix, ["X1", "X2"])

    def test_error_message_contains_missing_ids(self) -> None:
        matrix = _make_matrix(sample_ids=["S0", "S1"])
        with pytest.raises(MatrixSliceError) as exc_info:
            _extract_slice(matrix, ["S0", "BAD1", "BAD2"])

        assert "BAD1" in str(exc_info.value)
        assert "BAD2" in str(exc_info.value)


# ---------------------------------------------------------------------------
# prepare_differential_inputs — integration
# ---------------------------------------------------------------------------


class TestPrepareDifferentialInputs:
    """Verify end-to-end slicing with ResolvedContrast."""

    def test_returns_differential_input(self) -> None:
        matrix = _make_matrix(sample_ids=["S0", "S1", "S2", "S3"])
        contrast = _make_contrast(
            comparison_ids=["S0", "S1"],
            reference_ids=["S2", "S3"],
        )

        result = prepare_differential_inputs(matrix, contrast)

        assert isinstance(result, DifferentialInput)
        assert isinstance(result.comparison, MatrixSlice)
        assert isinstance(result.reference, MatrixSlice)

    def test_comparison_has_correct_samples(self) -> None:
        matrix = _make_matrix(sample_ids=["S0", "S1", "S2", "S3"])
        contrast = _make_contrast(
            comparison_ids=["S0", "S2"],
            reference_ids=["S1", "S3"],
        )

        result = prepare_differential_inputs(matrix, contrast)

        assert result.comparison.sample_ids == ["S0", "S2"]
        assert result.comparison.dataframe.shape[1] == 3  # gene + 2 samples

    def test_reference_has_correct_samples(self) -> None:
        matrix = _make_matrix(sample_ids=["S0", "S1", "S2", "S3"])
        contrast = _make_contrast(
            comparison_ids=["S0"],
            reference_ids=["S1", "S2", "S3"],
        )

        result = prepare_differential_inputs(matrix, contrast)

        assert result.reference.sample_ids == ["S1", "S2", "S3"]
        assert result.reference.dataframe.shape[1] == 4  # gene + 3 samples

    def test_order_preserved_from_matrix(self) -> None:
        """Slices must use matrix column order, not contrast tuple order."""
        matrix = _make_matrix(sample_ids=["S0", "S1", "S2", "S3"])
        # Contrast requests in reverse order
        contrast = _make_contrast(
            comparison_ids=["S3", "S1"],
            reference_ids=["S2", "S0"],
        )

        result = prepare_differential_inputs(matrix, contrast)

        assert result.comparison.sample_ids == ["S1", "S3"]
        assert result.reference.sample_ids == ["S0", "S2"]

    def test_disjoint_groups(self) -> None:
        """Comparison and reference slices must be independent."""
        matrix = _make_matrix(sample_ids=["S0", "S1", "S2"])
        contrast = _make_contrast(
            comparison_ids=["S0"],
            reference_ids=["S1"],
        )

        result = prepare_differential_inputs(matrix, contrast)

        comp_vals = set(result.comparison.dataframe.iloc[:, 1])
        ref_vals = set(result.reference.dataframe.iloc[:, 1])
        assert comp_vals != ref_vals

    def test_raises_if_comparison_sample_missing(self) -> None:
        matrix = _make_matrix(sample_ids=["S0", "S1"])
        contrast = _make_contrast(
            comparison_ids=["S0", "MISSING"],
            reference_ids=["S1"],
        )

        with pytest.raises(MatrixSliceError):
            prepare_differential_inputs(matrix, contrast)

    def test_raises_if_reference_sample_missing(self) -> None:
        matrix = _make_matrix(sample_ids=["S0", "S1"])
        contrast = _make_contrast(
            comparison_ids=["S0"],
            reference_ids=["S1", "MISSING"],
        )

        with pytest.raises(MatrixSliceError):
            prepare_differential_inputs(matrix, contrast)

    def test_shared_matrix_id(self) -> None:
        matrix = _make_matrix()
        contrast = _make_contrast()

        result = prepare_differential_inputs(matrix, contrast)

        assert result.comparison.matrix_id == matrix.matrix_id
        assert result.reference.matrix_id == matrix.matrix_id


# ---------------------------------------------------------------------------
# MatrixSlice immutability
# ---------------------------------------------------------------------------


class TestMatrixSliceImmutability:
    """Verify MatrixSlice is frozen."""

    def test_slice_is_frozen(self) -> None:
        matrix = _make_matrix()
        result = _extract_slice(matrix, ["S0"])

        with pytest.raises((TypeError, AttributeError)):
            result.matrix_id = "changed"

    def test_differential_input_is_frozen(self) -> None:
        matrix = _make_matrix()
        contrast = _make_contrast()
        result = prepare_differential_inputs(matrix, contrast)

        with pytest.raises((TypeError, AttributeError)):
            result.comparison = None
