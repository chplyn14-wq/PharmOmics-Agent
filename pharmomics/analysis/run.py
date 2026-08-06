"""Analysis execution entry point for PharmOmics.

Provides ``run_analysis()``, the top-level orchestration function that
validates inputs, dispatches by analysis type, and returns a structured
``AnalysisResult``.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

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
    """MVP placeholder for differential analysis execution.

    Returns a structurally valid ``AnalysisResult`` with no gene-level
    data.  All numerical fields are deferred to a future statistical
    backend.

    Parameters are intentionally unused — this function is a stub.
    """
    _ = design, omics  # unused in MVP

    return AnalysisResult(
        analysis_type="differential_analysis",
        contrast_id=specification.contrast_references[0],
        gene_results=(),
        n_genes_tested=0,
        warnings=("Differential analysis not implemented yet.",),
    )
