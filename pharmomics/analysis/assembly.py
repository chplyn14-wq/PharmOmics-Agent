"""Assemble GeneDifferential results from intermediate analysis outputs.

Takes the outputs of the four analysis stages (mean expression, log2 fold
change, raw p-values, BH-adjusted p-values) and assembles them into
``GeneDifferential`` tuples, aligned by ``gene_id``.

No statistical computation is performed here.  This is purely a structural
assembly step that consumes already-computed intermediate results.
"""

from __future__ import annotations

import math
from typing import TYPE_CHECKING

from pharmomics.analysis.results import GeneDifferential

if TYPE_CHECKING:
    from pharmomics.analysis.bh_fdr import GenePValueAdj
    from pharmomics.analysis.log2_fold_change import GeneLog2FoldChange
    from pharmomics.analysis.mean_expression import GeneMeanExpression
    from pharmomics.analysis.per_gene_pvalue import GenePValue
    from pharmomics.analysis.schemas import AnalysisSpecification


# ---------------------------------------------------------------------------
# Error
# ---------------------------------------------------------------------------


class AssemblyError(ValueError):
    """Raised when intermediate results cannot be assembled."""


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

_DEFAULT_FDR_THRESHOLD = 0.05


def assemble_gene_differentials(
    mean_expressions: tuple[GeneMeanExpression, ...],
    log2fc_results: tuple[GeneLog2FoldChange, ...],
    pvalue_results: tuple[GenePValue, ...],
    adjusted_results: tuple[GenePValueAdj, ...],
    specification: AnalysisSpecification,
) -> tuple[GeneDifferential, ...]:
    """Assemble per-gene differential analysis results from intermediate outputs.

    Parameters
    ----------
    mean_expressions : tuple[GeneMeanExpression, ...]
        Output of ``compute_mean_expression()``.
    log2fc_results : tuple[GeneLog2FoldChange, ...]
        Output of ``compute_log2_fold_change()``.
    pvalue_results : tuple[GenePValue, ...]
        Output of ``compute_per_gene_p_value()``.
    adjusted_results : tuple[GenePValueAdj, ...]
        Output of ``compute_bh_adjusted_p_values()``.
    specification : AnalysisSpecification
        Analysis specification carrying the ``fdr_threshold`` parameter
        (default 0.05 if not specified).

    Returns
    -------
    tuple[GeneDifferential, ...]
        One entry per gene, in the same order as ``mean_expressions``.

    Raises
    ------
    AssemblyError
        If any input is empty, if gene_id sets are inconsistent across
        the four intermediate result sets, or if duplicate gene_ids
        are found within any single input.
    """
    if not mean_expressions:
        return ()

    # Validate no duplicates within each input.
    _check_no_duplicates(mean_expressions)
    _check_no_duplicates(log2fc_results)
    _check_no_duplicates(pvalue_results)
    _check_no_duplicates(adjusted_results)

    # Resolve gene_id sets.
    base_ids = _gene_ids(mean_expressions)
    log2fc_ids = _gene_ids(log2fc_results)
    pval_ids = _gene_ids(pvalue_results)
    adj_ids = _gene_ids(adjusted_results)

    if base_ids != log2fc_ids:
        raise AssemblyError(
            f"gene_id mismatch: mean_expressions vs log2fc_results: "
            f"missing={sorted(base_ids - log2fc_ids)}, "
            f"extra={sorted(log2fc_ids - base_ids)}"
        )
    if base_ids != pval_ids:
        raise AssemblyError(
            f"gene_id mismatch: mean_expressions vs pvalue_results: "
            f"missing={sorted(base_ids - pval_ids)}, "
            f"extra={sorted(pval_ids - base_ids)}"
        )
    if base_ids != adj_ids:
        raise AssemblyError(
            f"gene_id mismatch: mean_expressions vs adjusted_results: "
            f"missing={sorted(base_ids - adj_ids)}, "
            f"extra={sorted(adj_ids - base_ids)}"
        )

    # Build lookup maps keyed by gene_id.
    log2fc_map = {f.gene_id: f.log2fc for f in log2fc_results}
    pval_map = {p.gene_id: p.p_value for p in pvalue_results}
    adj_map = {a.gene_id: a.adj_p_value for a in adjusted_results}

    # Read threshold from specification parameters.
    fdr_threshold = float(
        specification.parameters.get("fdr_threshold", _DEFAULT_FDR_THRESHOLD)
    )

    return tuple(
        _build_differential(m, log2fc_map, pval_map, adj_map, fdr_threshold)
        for m in mean_expressions
    )


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------


def _gene_ids(seq: tuple[object, ...]) -> set[str]:
    """Extract the set of gene_ids from a sequence of named-tuple-like objects."""
    return {getattr(item, "gene_id") for item in seq}


def _check_no_duplicates(seq: tuple[object, ...]) -> None:
    """Raise AssemblyError if any gene_id appears more than once."""
    seen: set[str] = set()
    for item in seq:
        gid = getattr(item, "gene_id")
        if gid in seen:
            raise AssemblyError(f"Duplicate gene_id: {gid!r}")
        seen.add(gid)


def _build_differential(
    mean: object,
    log2fc_map: dict[str, float],
    pval_map: dict[str, float],
    adj_map: dict[str, float],
    fdr_threshold: float,
) -> GeneDifferential:
    """Build a single GeneDifferential from lookup maps."""
    gid = getattr(mean, "gene_id")
    comp_mean = getattr(mean, "comparison_mean")
    ref_mean = getattr(mean, "reference_mean")

    base_mean = (comp_mean + ref_mean) / 2.0
    log2fc = log2fc_map[gid]
    p_value = pval_map[gid]
    adj_p_value = adj_map[gid]

    if math.isnan(adj_p_value):
        significant = False
    else:
        significant = adj_p_value < fdr_threshold

    return GeneDifferential(
        gene_id=gid,
        log2_fold_change=log2fc,
        p_value=p_value,
        adj_p_value=adj_p_value,
        significant=significant,
        base_mean=base_mean,
    )
