"""Structured result models for PharmOmics analysis execution.

Provides ``GeneDifferential`` and ``AnalysisResult`` — frozen dataclasses
that capture the output of a differential analysis run without carrying
any statistical logic or interpretation.
"""

from __future__ import annotations

from dataclasses import dataclass


class AnalysisConsistencyError(ValueError):
    """Raised when AnalysisResult fields are internally inconsistent."""


@dataclass(frozen=True)
class GeneDifferential:
    """Statistical result for a single gene.

    All values are computed by the analysis backend; this class is a
    pure data container with no methods or validation logic.
    """

    gene_id: str
    log2_fold_change: float
    p_value: float
    adj_p_value: float
    significant: bool
    base_mean: float


@dataclass(frozen=True)
class AnalysisResult:
    """Complete result of a single analysis execution.

    Parameters
    ----------
    analysis_type : str
        The type of analysis performed (e.g. ``"differential_analysis"``).
    contrast_id : str
        Identifier of the contrast this result corresponds to.
    gene_results : tuple[GeneDifferential, ...]
        Per-gene statistical results, in feature-id order.
    n_genes_tested : int
        Number of genes included in the test.  Must equal
        ``len(gene_results)``.
    warnings : tuple[str, ...]
        Non-fatal warnings from the analysis pipeline.

    Raises
    ------
    AnalysisConsistencyError
        If ``n_genes_tested`` does not equal ``len(gene_results)``.
    """

    analysis_type: str
    contrast_id: str
    gene_results: tuple[GeneDifferential, ...]
    n_genes_tested: int
    warnings: tuple[str, ...] = ()

    def __post_init__(self) -> None:
        if self.n_genes_tested != len(self.gene_results):
            raise AnalysisConsistencyError(
                f"n_genes_tested ({self.n_genes_tested}) != "
                f"len(gene_results) ({len(self.gene_results)})"
            )
