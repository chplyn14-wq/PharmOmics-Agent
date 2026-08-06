"""Orchestration runner for PharmOmics analysis pipelines.

Provides ``validate_analysis_inputs()``, which runs the full cross-domain
validation chain before analysis execution.  All violations are collected
and reported together — no early return on the first failure.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from pharmomics.analysis.schemas import AnalysisSpecification
    from pharmomics.experiment.schemas import ExperimentDesign
    from pharmomics.omics.schemas import OmicsMatrix


class AnalysisValidationError(ValueError):
    """Raised when analysis input validation fails."""


def validate_analysis_inputs(
    specification: AnalysisSpecification,
    design: ExperimentDesign,
    omics: OmicsMatrix,
) -> None:
    """Run the full validation chain for an analysis run.

    Checks performed (in order):

    1. ``validate(design)`` — internal consistency of the experiment design.
    2. ``check_compatibility(design, omics)`` — design samples align with
       omics sample IDs.
    3. ``check_analysis_design_compatibility(specification, design)`` —
       analysis specification references resolve to design entities.

    All violations are collected before raising.  If any exist, they are
    joined into a single ``AnalysisValidationError``.

    Parameters
    ----------
    specification : AnalysisSpecification
        The analysis intent to validate.
    design : ExperimentDesign
        The experimental design to validate.
    omics : OmicsMatrix
        The omics matrix to validate against.

    Raises
    ------
    AnalysisValidationError
        If any validation check returns violations.  The exception message
        contains all violations, one per line.

    Returns
    -------
    None
    """
    from pharmomics.compatibility.analysis_design import (
        check_analysis_design_compatibility,
    )
    from pharmomics.compatibility.omics_design import check_compatibility
    from pharmomics.experiment.validation import validate

    violations: list[str] = []
    violations.extend(validate(design))
    violations.extend(check_compatibility(design, omics))
    violations.extend(check_analysis_design_compatibility(specification, design))

    if violations:
        message = "Validation failed:\n" + "\n".join(
            f"- {v}" for v in violations
        )
        raise AnalysisValidationError(message)
