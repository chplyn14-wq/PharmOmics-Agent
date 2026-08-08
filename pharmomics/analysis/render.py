"""Render an ``AnalysisResult`` as a deterministic Markdown report and TSV exports.

Provides ``render_markdown_report()``, a pure function that consumes an
``AnalysisResult`` and returns a Markdown string suitable for display,
export, or attachment.  Also provides ``render_results_tsv()`` and
``render_significant_genes_tsv()`` for tab-separated export of all
genes and significant-only genes, respectively.  No statistical values
are computed or modified.
"""

from __future__ import annotations

import math
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from pharmomics.analysis.results import AnalysisResult

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

_LOG2FC_FMT = ".3f"
_BASE_MEAN_FMT = ".2f"
_PVALUE_FMT = ".4g"

# TSV column order for both results files
_TSV_COLUMNS = (
    "gene_id",
    "log2_fold_change",
    "p_value",
    "adj_p_value",
    "significant",
    "base_mean",
)

# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------


def render_markdown_report(result: AnalysisResult) -> str:
    """Render an ``AnalysisResult`` as a deterministic Markdown report.

    Parameters
    ----------
    result : AnalysisResult
        A completed analysis result.  All values are consumed as-is;
        no statistical computation is performed.

    Returns
    -------
    str
        A Markdown document containing: a metadata header, a warnings
        section (or ``None``), and a gene results table.

    Notes
    -----
    Output is fully deterministic for a given ``AnalysisResult`` input.
    Gene results appear in the original ``result.gene_results`` order.
    """
    sections: list[str] = []

    # Title
    sections.append("# Differential Analysis Report\n")

    # Summary
    sections.append("## Summary\n")
    sections.append("| Field | Value |")
    sections.append("|---|---|")
    sections.append(f"| Analysis type | {result.analysis_type} |")
    sections.append(f"| Contrast | {result.contrast_id} |")
    sections.append(f"| Genes tested | {result.n_genes_tested} |")
    sections.append("")

    # Warnings
    sections.append("## Warnings\n")
    if result.warnings:
        for w in result.warnings:
            sections.append(f"- {w}")
    else:
        sections.append("None")
    sections.append("")

    # Gene results table
    sections.append("## Gene Results\n")
    sections.append(
        "| Gene | log2FC | p-value | adj p-value | Significant | Base mean |"
    )
    sections.append("|---|---|---|---|---|---|")

    for g in result.gene_results:
        sig = "Yes" if g.significant else "No"
        sections.append(
            f"| {g.gene_id}"
            f" | {_format_float(g.log2_fold_change, _LOG2FC_FMT)}"
            f" | {_format_float(g.p_value, _PVALUE_FMT)}"
            f" | {_format_float(g.adj_p_value, _PVALUE_FMT)}"
            f" | {sig}"
            f" | {_format_float(g.base_mean, _BASE_MEAN_FMT)} |"
        )

    return "\n".join(sections)


def render_results_tsv(result: AnalysisResult) -> str:
    """Render all gene results as a tab-separated string.

    Parameters
    ----------
    result : AnalysisResult
        A completed analysis result.  All values are consumed as-is.

    Returns
    -------
    str
        A TSV document with a header row and one row per gene in
        ``result.gene_results`` order.  Non-finite floats render as
        ``NaN``, ``+Inf`` or ``-Inf``.

    Notes
    -----
    Output is fully deterministic for a given ``AnalysisResult`` input.
    Gene order matches ``result.gene_results`` exactly; no sorting is
    applied.
    """
    lines: list[str] = []
    lines.append("\t".join(_TSV_COLUMNS))
    for g in result.gene_results:
        lines.append(
            "\t".join(
                [
                    g.gene_id,
                    _tsv_float(g.log2_fold_change),
                    _tsv_float(g.p_value),
                    _tsv_float(g.adj_p_value),
                    str(g.significant),
                    _tsv_float(g.base_mean),
                ]
            )
        )
    return "\n".join(lines) + "\n"


def render_significant_genes_tsv(result: AnalysisResult) -> str:
    """Render only significant genes as a tab-separated string.

    Parameters
    ----------
    result : AnalysisResult
        A completed analysis result.

    Returns
    -------
    str
        A TSV document with a header row followed by rows for genes
        where ``gene.significant is True``.  Relative order matches
        their position in ``result.gene_results``.  If no genes are
        significant, only the header row is returned.
    """
    lines: list[str] = []
    lines.append("\t".join(_TSV_COLUMNS))
    for g in result.gene_results:
        if g.significant:
            lines.append(
                "\t".join(
                    [
                        g.gene_id,
                        _tsv_float(g.log2_fold_change),
                        _tsv_float(g.p_value),
                        _tsv_float(g.adj_p_value),
                        str(g.significant),
                        _tsv_float(g.base_mean),
                    ]
                )
            )
    return "\n".join(lines) + "\n"


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------


def _format_float(value: float, fmt: str) -> str:
    """Format a float value, handling NaN and Inf deterministically.

    Parameters
    ----------
    value : float
        The numeric value to format.
    fmt : str
        A Python format specifier (e.g. ``".3f"``).

    Returns
    -------
    str
        Formatted string, or ``NaN``, ``+Inf``, ``-Inf`` for non-finite
        values.
    """
    if math.isnan(value):
        return "NaN"
    if math.isinf(value):
        return "+Inf" if value > 0 else "-Inf"
    return f"{value:{fmt}}"


def _tsv_float(value: float) -> str:
    """Format a float for TSV output, preserving NaN as explicit string.

    Unlike ``_format_float``, this uses the raw ``repr`` for finite
    values so the TSV carries the full-precision number without
    rounding.  Non-finite values render as ``NaN`` / ``+Inf`` / ``-Inf``.
    """
    if math.isnan(value):
        return "NaN"
    if math.isinf(value):
        return "+Inf" if value > 0 else "-Inf"
    return repr(value)
