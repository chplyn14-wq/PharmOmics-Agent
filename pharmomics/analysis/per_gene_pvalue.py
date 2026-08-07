"""Compute per-gene Welch's two-sample t-test p-values from a DifferentialInput.

Takes a ``DifferentialInput`` (comparison and reference ``MatrixSlice``
objects) and returns the Welch's t-test p-value for each gene, computed
from raw sample-level expression values.

No FDR correction, no fold-change computation.  Pure per-gene p-values
derived from deterministic arithmetic on the underlying sample data.
"""

from __future__ import annotations

from dataclasses import dataclass
from math import isinf, isnan
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from pharmomics.analysis.matrix_slice import DifferentialInput


class PerGenePValueError(ValueError):
    """Raised when per-gene p-value cannot be computed due to structural issues."""


@dataclass(frozen=True)
class GenePValue:
    """Welch's t-test p-value for a single gene."""

    gene_id: str
    p_value: float


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------


def compute_per_gene_p_value(
    differential_input: DifferentialInput,
) -> tuple[GenePValue, ...]:
    """Compute per-gene Welch's two-sample t-test p-values.

    For each gene, performs an independent Welch's t-test comparing
    the comparison-group sample expression values against the
    reference-group values.

    Parameters
    ----------
    differential_input : DifferentialInput
        Comparison and reference MatrixSlices from
        ``prepare_differential_inputs()``.  Sample columns must be
        numeric.

    Returns
    -------
    tuple[GenePValue, ...]
        One entry per gene, in original feature-id order.  Genes that
        cannot be tested (insufficient samples, zero variance, or
        non-finite values) receive ``p_value = NaN``.

    Raises
    ------
    PerGenePValueError
        If feature sets differ between groups, a slice has no genes,
        a group has no sample columns, or sample data is non-numeric.
    """
    comp = differential_input.comparison
    ref = differential_input.reference

    _validate_slices(comp, ref)

    comp_sample_cols = comp.dataframe.columns[1:]
    ref_sample_cols = ref.dataframe.columns[1:]

    comp_values = comp.dataframe[comp_sample_cols]
    ref_values = ref.dataframe[ref_sample_cols]

    return tuple(
        _gene_p_value(
            gene_id=str(comp.dataframe.iloc[i, 0]),
            comp_row=comp_values.iloc[i],
            ref_row=ref_values.iloc[i],
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
    PerGenePValueError
        If feature sets are inconsistent or slices are empty.
    """
    comp_features = tuple(comp.feature_ids)  # type: ignore[union-attr]
    ref_features = tuple(ref.feature_ids)  # type: ignore[union-attr]

    if not comp_features:
        raise PerGenePValueError("Comparison slice has no features")
    if not ref_features:
        raise PerGenePValueError("Reference slice has no features")
    if comp_features != ref_features:
        raise PerGenePValueError(
            "Feature sets differ between comparison and reference slices"
        )

    # Check for sample columns (dataframe must have more than just the feature column).
    comp_n_cols = comp.dataframe.shape[1]  # type: ignore[union-attr]
    ref_n_cols = ref.dataframe.shape[1]  # type: ignore[union-attr]

    if comp_n_cols < 2:
        raise PerGenePValueError("Comparison slice has no sample columns")
    if ref_n_cols < 2:
        raise PerGenePValueError("Reference slice has no sample columns")

    # Verify numeric data in sample columns.
    comp_sample_cols = comp.dataframe.columns[1:]  # type: ignore[union-attr]
    ref_sample_cols = ref.dataframe.columns[1:]  # type: ignore[union-attr]

    comp_values = comp.dataframe[comp_sample_cols]  # type: ignore[union-attr]
    ref_values = ref.dataframe[ref_sample_cols]  # type: ignore[union-attr]

    if not _is_numeric(comp_values):
        raise PerGenePValueError(
            "Comparison slice contains non-numeric sample values"
        )
    if not _is_numeric(ref_values):
        raise PerGenePValueError("Reference slice contains non-numeric sample values")


def _is_numeric(df: object) -> bool:
    """Check if all columns in a DataFrame are numeric dtype.

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


def _gene_p_value(
    gene_id: str,
    comp_row: object,
    ref_row: object,
) -> GenePValue:
    """Compute Welch's t-test p-value for a single gene.

    Parameters
    ----------
    gene_id : str
        The gene identifier.
    comp_row : pd.Series
        Sample expression values for the comparison group.
    ref_row : pd.Series
        Sample expression values for the reference group.

    Returns
    -------
    GenePValue
        The gene's p-value, or ``NaN`` when the test cannot be
        computed (insufficient samples, zero variance, non-finite
        values, or scipy returns NaN).
    """
    from scipy.stats import ttest_ind

    comp_vals = list(comp_row)  # type: ignore[union-attr]
    ref_vals = list(ref_row)  # type: ignore[union-attr]

    # Insufficient samples: need >= 2 per group to estimate variance.
    if len(comp_vals) < 2 or len(ref_vals) < 2:
        return GenePValue(gene_id=gene_id, p_value=float("nan"))

    # Non-finite values in either group.
    for v in comp_vals + ref_vals:
        if not isinstance(v, (int, float)) or isnan(v) or isinf(v):
            return GenePValue(gene_id=gene_id, p_value=float("nan"))

    # Zero variance in either group.
    comp_var = _variance(comp_vals)
    ref_var = _variance(ref_vals)
    if comp_var == 0.0 or ref_var == 0.0:
        return GenePValue(gene_id=gene_id, p_value=float("nan"))

    # Run Welch's t-test.
    try:
        result = ttest_ind(comp_vals, ref_vals, equal_var=False)
        p = float(result.pvalue)
        if isnan(p):
            return GenePValue(gene_id=gene_id, p_value=float("nan"))
        return GenePValue(gene_id=gene_id, p_value=p)
    except Exception:
        return GenePValue(gene_id=gene_id, p_value=float("nan"))


def _variance(values: list[float]) -> float:
    """Compute sample variance (Bessel-corrected, ddof=1).

    Parameters
    ----------
    values : list[float]
        Numeric values.

    Returns
    -------
    float
        Sample variance.
    """
    n = len(values)
    mean = sum(values) / n
    return sum((x - mean) ** 2 for x in values) / (n - 1)
