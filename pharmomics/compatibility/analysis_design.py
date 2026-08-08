"""Compatibility checks between AnalysisSpecification and ExperimentDesign.

Ensures that every factor_reference and contrast_reference declared in an
AnalysisSpecification corresponds to an actual entity in the ExperimentDesign.
All checks are read-only — neither the specification nor the design is
mutated.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from pharmomics.analysis.schemas import AnalysisSpecification
    from pharmomics.experiment.schemas import ExperimentDesign


def check_analysis_design_compatibility(
    specification: AnalysisSpecification,
    design: ExperimentDesign,
) -> list[str]:
    """Check whether *specification* references exist in *design*.

    Parameters
    ----------
    specification : AnalysisSpecification
        The analysis specification whose references to validate.
    design : ExperimentDesign
        The experimental design to resolve references against.

    Returns
    -------
    list[str]
        Empty list if all references resolve;
        otherwise a list of human-readable violation descriptions.

    Notes
    -----
    Pure function — does not mutate either input.
    Does not perform statistical validation or design matrix generation.
    """
    violations: list[str] = []
    violations.extend(
        _check_factor_references(specification, design),
    )
    violations.extend(
        _check_contrast_references(specification, design),
    )
    return violations


# ---------------------------------------------------------------------------
# Individual checks
# ---------------------------------------------------------------------------


def _check_factor_references(
    specification: AnalysisSpecification,
    design: ExperimentDesign,
) -> list[str]:
    """Ensure every factor_reference matches a declared factor_id."""
    valid_ids = {f.factor_id for f in design.factors}
    missing = [ref for ref in specification.factor_references if ref not in valid_ids]

    if not missing:
        return []

    examples = missing[:5]
    return [
        f"AnalysisSpecification factor_reference(s) not found in ExperimentDesign: "
        f"{examples}",
    ]


def _check_contrast_references(
    specification: AnalysisSpecification,
    design: ExperimentDesign,
) -> list[str]:
    """Ensure every contrast_reference matches a declared contrast_id."""
    valid_ids = {c.contrast_id for c in design.contrasts}
    missing = [ref for ref in specification.contrast_references if ref not in valid_ids]

    if not missing:
        return []

    examples = missing[:5]
    return [
        f"AnalysisSpecification contrast_reference(s) not found in ExperimentDesign: "
        f"{examples}",
    ]
