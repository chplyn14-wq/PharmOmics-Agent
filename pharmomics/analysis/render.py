"""Render an AnalysisResult as a deterministic Markdown report.

Provides ``render_markdown_report()``, a pure function that consumes an
``AnalysisResult`` and returns a Markdown string suitable for display,
export, or attachment.  No statistical values are computed or modified.
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
