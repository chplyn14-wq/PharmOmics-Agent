"""Benjamini-Hochberg FDR correction for per-gene raw p-values.

Takes the output of ``compute_per_gene_p_value()`` and produces
BH-adjusted p-values.  NaN raw p-values are excluded from the
correction procedure and receive NaN as their adjusted value.

The multiplicity parameter ``m`` equals the number of valid (finite,
non-NaN) p-values, not the total feature count.
"""

from __future__ import annotations

from dataclasses import dataclass
from math import isfinite, isnan
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from pharmomics.analysis.per_gene_pvalue import GenePValue


# ---------------------------------------------------------------------------
# Result model
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class GenePValueAdj:
    """Single gene with its raw and BH-adjusted p-value."""

    gene_id: str
    raw_p_value: float
    adj_p_value: float


# ---------------------------------------------------------------------------
# Error
# ---------------------------------------------------------------------------


class BHAdjustmentError(ValueError):
    """Raised when BH adjustment cannot proceed due to invalid input."""


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------


def compute_bh_adjusted_p_values(
    raw_p_values: tuple[GenePValue, ...],
) -> tuple[GenePValueAdj, ...]:
    """Apply Benjamini-Hochberg FDR correction to per-gene raw p-values.

    Parameters
    ----------
    raw_p_values : tuple[GenePValue, ...]
        Output from ``compute_per_gene_p_value()``.  Order is preserved.

    Returns
    -------
    tuple[GenePValueAdj, ...]
        One entry per gene, in the same order as input.  Genes with
        ``NaN`` raw p-value receive ``NaN`` for ``adj_p_value``.
        All ``adj_p_value`` entries are clamped to [0, 1].

    Raises
    ------
    BHAdjustmentError
        If any finite raw p-value falls outside [0, 1].
    """
    if not raw_p_values:
        return ()

    # Validate all finite p-values are in [0, 1].
    for entry in raw_p_values:
        p = entry.p_value
        if isfinite(p) and not (0.0 <= p <= 1.0):
            raise BHAdjustmentError(
                f"p-value for gene {entry.gene_id!r} is out of [0, 1]: {p}"
            )

    # Identify valid (non-NaN) entries with their original indices.
    valid_indices = [
        i for i, e in enumerate(raw_p_values) if not isnan(e.p_value)
    ]

    # All NaN — return all-NaN result.
    if not valid_indices:
        return tuple(
            GenePValueAdj(
                gene_id=e.gene_id,
                raw_p_value=float("nan"),
                adj_p_value=float("nan"),
            )
            for e in raw_p_values
        )

    m = len(valid_indices)

    # Build (original_index, raw_p) pairs and sort by raw p-value.
    sorted_entries = sorted(
        ((i, raw_p_values[i].p_value) for i in valid_indices),
        key=lambda x: x[1],
    )

    # Step-up BH: adj[i] = p[i] * m / rank (rank is 1-based).
    raw_adjusted = [p * m / rank for rank, (_, p) in enumerate(sorted_entries, 1)]

    # Monotonicity: walk backwards, ensuring non-decreasing adjusted values.
    for i in range(m - 2, -1, -1):
        if raw_adjusted[i] > raw_adjusted[i + 1]:
            raw_adjusted[i] = raw_adjusted[i + 1]

    # Clamp to [0, 1].
    adj_values = [max(0.0, min(1.0, v)) for v in raw_adjusted]

    # Map back to original indices.
    adj_map = {sorted_entries[i][0]: adj_values[i] for i in range(m)}

    return tuple(
        GenePValueAdj(
            gene_id=e.gene_id,
            raw_p_value=e.p_value if i in adj_map else float("nan"),
            adj_p_value=adj_map[i] if i in adj_map else float("nan"),
        )
        for i, e in enumerate(raw_p_values)
    )
