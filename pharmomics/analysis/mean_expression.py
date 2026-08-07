"""Compute per-feature arithmetic mean expression from a DifferentialInput.

Takes a ``DifferentialInput`` (comparison and reference ``MatrixSlice``
objects) and returns the arithmetic mean of each feature across the
sample columns for both groups.

No statistical tests, no log2FC, no p-values.  Pure descriptive
summaries of the numeric values in each slice.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from pharmomics.analysis.matrix_slice import DifferentialInput


class MeanExpressionError(ValueError):
    """Raised when mean expression cannot be computed."""


@dataclass(frozen=True)
class GeneMeanExpression:
    """Arithmetic mean expression for a single feature."""

    gene_id: str
    comparison_mean: float
    reference_mean: float


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------


def compute_mean_expression(
    differential_input: DifferentialInput,
) -> tuple[GeneMeanExpression, ...]:
    """Compute per-feature arithmetic mean for comparison and reference groups.

    Parameters
    ----------
    differential_input : DifferentialInput
        Comparison and reference MatrixSlices from ``prepare_differential_inputs()``.

    Returns
    -------
    tuple[GeneMeanExpression, ...]
        One entry per feature, in the original feature-id order.

    Raises
    ------
    MeanExpressionError
        If feature sets are inconsistent, a slice has no features,
        a slice has no sample columns, or non-numeric values are found.
    """
    comp = differential_input.comparison
    ref = differential_input.reference

    _validate_slices(comp, ref)

    # Feature ID column is column 0; sample columns are the rest.
    comp_sample_cols = comp.dataframe.columns[1:]
    ref_sample_cols = ref.dataframe.columns[1:]

    # Compute arithmetic mean across sample columns per feature row.
    comp_means = comp.dataframe[comp_sample_cols].mean(axis=1)
    ref_means = ref.dataframe[ref_sample_cols].mean(axis=1)

    return tuple(
        GeneMeanExpression(
            gene_id=str(comp.dataframe.iloc[i, 0]),
            comparison_mean=float(comp_means.iloc[i]),
            reference_mean=float(ref_means.iloc[i]),
        )
        for i in range(len(comp.feature_ids))
    )


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------


def _validate_slices(comp: object, ref: object) -> None:
    """Validate that comparison and reference slices are compatible.

    Parameters
    ----------
    comp : MatrixSlice
        The comparison slice.
    ref : MatrixSlice
        The reference slice.

    Raises
    ------
    MeanExpressionError
        If feature sets are inconsistent or slices are empty.
    """
    comp_features = tuple(comp.feature_ids)  # type: ignore[union-attr]
    ref_features = tuple(ref.feature_ids)  # type: ignore[union-attr]

    if not comp_features:
        raise MeanExpressionError("Comparison slice has no features")
    if not ref_features:
        raise MeanExpressionError("Reference slice has no features")
    if comp_features != ref_features:
        raise MeanExpressionError(
            "Feature sets differ between comparison and reference slices"
        )

    # Check for sample columns (dataframe must have more than just the feature column).
    comp_n_cols = comp.dataframe.shape[1]  # type: ignore[union-attr]
    ref_n_cols = ref.dataframe.shape[1]  # type: ignore[union-attr]

    if comp_n_cols < 2:
        raise MeanExpressionError("Comparison slice has no sample columns")
    if ref_n_cols < 2:
        raise MeanExpressionError("Reference slice has no sample columns")

    # Verify numeric data in sample columns.
    comp_sample_cols = comp.dataframe.columns[1:]  # type: ignore[union-attr]
    ref_sample_cols = ref.dataframe.columns[1:]  # type: ignore[union-attr]

    comp_values = comp.dataframe[comp_sample_cols]  # type: ignore[union-attr]
    ref_values = ref.dataframe[ref_sample_cols]  # type: ignore[union-attr]

    if not _is_numeric(comp_values):
        raise MeanExpressionError("Comparison slice contains non-numeric sample values")
    if not _is_numeric(ref_values):
        raise MeanExpressionError("Reference slice contains non-numeric sample values")


def _is_numeric(df: object) -> bool:
    """Check if all columns in a DataFrame are numeric.

    Parameters
    ----------
    df : pd.DataFrame
        The DataFrame to check.

    Returns
    -------
    bool
        True if all columns are numeric dtype.
    """
    import pandas as pd

    if not isinstance(df, pd.DataFrame):
        return False
    non_numeric = df.select_dtypes(exclude="number")
    return non_numeric.empty
