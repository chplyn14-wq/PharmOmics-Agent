"""Compute per-feature log2 fold change from mean expression results.

Takes a sequence of ``GeneMeanExpression`` objects (produced by
``compute_mean_expression()``) and returns the log2 fold change for
each feature as ``comparison_mean / reference_mean``.

No statistical tests, no p-values, no FDR.  Pure arithmetic derived
from deterministic upstream means.
"""

from __future__ import annotations

from dataclasses import dataclass
from math import log2
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from pharmomics.analysis.mean_expression import GeneMeanExpression


class Log2FoldChangeError(ValueError):
    """Raised when log2 fold change cannot be computed."""


@dataclass(frozen=True)
class GeneLog2FoldChange:
    """Log2 fold change for a single feature.

    ``log2fc`` is ``NaN`` (``float("nan")``) when either
    ``comparison_mean`` or ``reference_mean`` is zero or negative,
    because a log-ratio is undefined in that case.
    """

    gene_id: str
    log2fc: float


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------


def compute_log2_fold_change(
    mean_expressions: tuple[GeneMeanExpression, ...],
) -> tuple[GeneLog2FoldChange, ...]:
    """Compute log2 fold change for each feature from pre-computed means.

    For each feature, log2FC = log2(comparison_mean / reference_mean).
    When ``comparison_mean <= 0`` or ``reference_mean <= 0``, the
    result for that feature is ``NaN``.

    Parameters
    ----------
    mean_expressions : tuple[GeneMeanExpression, ...]
        Output of ``compute_mean_expression()``.

    Returns
    -------
    tuple[GeneLog2FoldChange, ...]
        One entry per feature, in the original feature-id order.

    Raises
    ------
    Log2FoldChangeError
        If the input tuple is empty.
    """
    if not mean_expressions:
        raise Log2FoldChangeError("No mean expression data provided")

    return tuple(_gene_log2fc(entry) for entry in mean_expressions)


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------


def _gene_log2fc(entry: GeneMeanExpression) -> GeneLog2FoldChange:
    """Compute log2FC for a single gene, returning NaN for invalid ratios.

    Parameters
    ----------
    entry : GeneMeanExpression
        A single feature's comparison and reference means.

    Returns
    -------
    GeneLog2FoldChange
        The feature's log2 fold change, or ``NaN`` when the ratio
        is undefined (zero or negative means).
    """
    if entry.comparison_mean <= 0 or entry.reference_mean <= 0:
        return GeneLog2FoldChange(
            gene_id=entry.gene_id,
            log2fc=float("nan"),
        )

    return GeneLog2FoldChange(
        gene_id=entry.gene_id,
        log2fc=log2(entry.comparison_mean / entry.reference_mean),
    )
