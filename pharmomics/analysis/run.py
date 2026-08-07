"""Analysis execution entry point for PharmOmics.

Provides ``run_analysis()``, the top-level orchestration function that
validates inputs, dispatches by analysis type, and returns a structured
``AnalysisResult``.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

from pharmomics.analysis.assembly import assemble_gene_differentials
from pharmomics.analysis.bh_fdr import compute_bh_adjusted_p_values
from pharmomics.analysis.contrast import resolve_contrast
from pharmomics.analysis.log2_fold_change import compute_log2_fold_change
from pharmomics.analysis.matrix_slice import prepare_differential_inputs
from pharmomics.analysis.mean_expression import compute_mean_expression
from pharmomics.analysis.per_gene_pvalue import compute_per_gene_p_value
from pharmomics.analysis.results import AnalysisResult
from pharmomics.analysis.runner import (
    AnalysisValidationError,
    validate_analysis_inputs,
)

if TYPE_CHECKING:
    from pharmomics.analysis.schemas import AnalysisSpecification
    from pharmomics.experiment.schemas import ExperimentDesign
    from pharmomics.omics.schemas import OmicsMatrix


def run_analysis(
    specification: AnalysisSpecification,
    design: ExperimentDesign,
    omics: OmicsMatrix,
) -> AnalysisResult:
    """Run an analysis pipeline end-to-end.

    Parameters
    ----------
    specification : AnalysisSpecification
        The analysis intent (type, factor/contrast references, parameters).
    design : ExperimentDesign
        The experimental design to run against.
    omics : OmicsMatrix
        The omics data matrix to analyse.

    Returns
    -------
    AnalysisResult
        Structured result from the analysis.

    Raises
    ------
    AnalysisValidationError
        If input validation fails or the analysis type is unsupported.
    """
    validate_analysis_inputs(specification, design, omics)

    if specification.analysis_type == "differential_analysis":
        return _execute_differential_analysis(specification, design, omics)

    raise AnalysisValidationError(
        f"Unsupported analysis_type: {specification.analysis_type}"
    )


def _execute_differential_analysis(
    specification: AnalysisSpecification,
    design: ExperimentDesign,
    omics: OmicsMatrix,
) -> AnalysisResult:
    """Execute a differential analysis pipeline end-to-end.

    Resolves the contrast, extracts matrix slices, computes per-gene
    statistics (mean expression, log2FC, Welch's p-value, BH-FDR),
    assembles ``GeneDifferential`` results, and returns an
    ``AnalysisResult``.

    Parameters
    ----------
    specification : AnalysisSpecification
        The analysis intent carrying the contrast reference and
        optional ``fdr_threshold`` parameter.
    design : ExperimentDesign
        The experimental design for contrast resolution.
    omics : OmicsMatrix
        The omics data matrix to analyse.

    Returns
    -------
    AnalysisResult
        Structured result with per-gene ``GeneDifferential`` tuples.

    Raises
    ------
    ContrastResolutionError
        If the contrast cannot be resolved in the design.
    MatrixSliceError
        If the matrix cannot be sliced for the resolved contrast.
    MeanExpressionError
        If mean expression cannot be computed from the slices.
    Log2FoldChangeError
        If log2 fold change cannot be computed.
    PerGenePValueError
        If per-gene p-values cannot be computed.
    BHAdjustmentError
        If BH-FDR adjustment fails.
    AssemblyError
        If intermediate results are inconsistent.
    """
    contrast_id = specification.contrast_references[0]

    resolved = resolve_contrast(design, contrast_id)
    diff_input = prepare_differential_inputs(omics, resolved)
    mean_expr = compute_mean_expression(diff_input)
    log2fc = compute_log2_fold_change(mean_expr)
    pvals = compute_per_gene_p_value(diff_input)
    adj = compute_bh_adjusted_p_values(pvals)
    gene_results = assemble_gene_differentials(
        mean_expr, log2fc, pvals, adj, specification
    )

    return AnalysisResult(
        analysis_type="differential_analysis",
        contrast_id=resolved.contrast_id,
        gene_results=gene_results,
        n_genes_tested=len(gene_results),
        warnings=(),
    )
