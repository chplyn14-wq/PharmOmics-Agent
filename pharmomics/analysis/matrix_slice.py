"""Extract comparison/reference sub-matrices from an OmicsMatrix for a contrast.

Provides ``prepare_differential_inputs()``, which takes a resolved contrast
and an omics matrix, then returns two ``MatrixSlice`` objects — one for the
comparison group and one for the reference group.

No statistical computation is performed here.  This layer only isolates
the numeric data needed by downstream differential analysis.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    import pandas as pd

    from pharmomics.analysis.contrast import ResolvedContrast
    from pharmomics.omics.schemas import OmicsMatrix


class MatrixSliceError(ValueError):
    """Raised when a matrix slice cannot be extracted."""


@dataclass(frozen=True)
class MatrixSlice:
    """A feature-by-sample sub-matrix extracted from an OmicsMatrix.

    Contains only the columns for the requested sample IDs, plus
    the feature ID column.  All feature rows are preserved.
    """

    matrix_id: str
    """Source OmicsMatrix matrix_id."""

    sample_ids: list[str] = field(default_factory=list)
    """The sample IDs in this slice, in original matrix order."""

    feature_ids: list[str] = field(default_factory=list)
    """All feature IDs (full row set, unchanged)."""

    dataframe: pd.DataFrame = field(default=None, repr=False)
    """Feature-by-sample sub-matrix with feature ID column.

    Same layout as ``OmicsMatrix.dataframe``: column 0 = feature IDs,
    remaining columns = samples in ``sample_ids`` order.
    """


@dataclass(frozen=True)
class DifferentialInput:
    """Comparison and reference sub-matrices for differential analysis."""

    comparison: MatrixSlice
    """Sub-matrix for the comparison (numerator) group."""

    reference: MatrixSlice
    """Sub-matrix for the reference (denominator) group."""


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------


def prepare_differential_inputs(
    matrix: OmicsMatrix,
    contrast: ResolvedContrast,
) -> DifferentialInput:
    """Extract comparison and reference sub-matrices for a resolved contrast.

    Parameters
    ----------
    matrix : OmicsMatrix
        The source matrix to slice.
    contrast : ResolvedContrast
        A resolved contrast with populated sample ID tuples.

    Returns
    -------
    DifferentialInput
        A frozen pair of ``MatrixSlice`` objects, one per group.

    Raises
    ------
    MatrixSliceError
        If any sample ID from the contrast is not found in
        ``matrix.sample_ids``, or if a group has no samples.
    """
    comparison_slice = _extract_slice(matrix, list(contrast.comparison_sample_ids))
    reference_slice = _extract_slice(matrix, list(contrast.reference_sample_ids))

    return DifferentialInput(
        comparison=comparison_slice,
        reference=reference_slice,
    )


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------


def _extract_slice(
    matrix: OmicsMatrix,
    sample_ids: list[str],
) -> MatrixSlice:
    """Extract a sub-matrix for the given sample IDs.

    Parameters
    ----------
    matrix : OmicsMatrix
        The source matrix.
    sample_ids : list[str]
        Sample IDs to include.  Must all exist in ``matrix.sample_ids``.

    Returns
    -------
    MatrixSlice
        A frozen sub-matrix containing all features and only the
        requested sample columns.

    Raises
    ------
    MatrixSliceError
        If ``sample_ids`` is empty, or if any ID is not found in
        ``matrix.sample_ids``.
    """
    if not sample_ids:
        raise MatrixSliceError(
            f"No sample IDs provided for slice of matrix '{matrix.matrix_id}'"
        )

    matrix_id_set = set(matrix.sample_ids)
    missing = [sid for sid in sample_ids if sid not in matrix_id_set]
    if missing:
        raise MatrixSliceError(
            f"Sample(s) {sorted(missing)} not found in OmicsMatrix '{matrix.matrix_id}'"
        )

    # Preserve original matrix column order
    ordered = [sid for sid in matrix.sample_ids if sid in set(sample_ids)]

    # Build sub-dataframe: feature ID column (index 0) + requested sample columns
    feature_col = matrix.dataframe.columns[0]
    sub_df = matrix.dataframe[[feature_col] + ordered].copy()

    return MatrixSlice(
        matrix_id=matrix.matrix_id,
        sample_ids=ordered,
        feature_ids=list(matrix.feature_ids),
        dataframe=sub_df,
    )
