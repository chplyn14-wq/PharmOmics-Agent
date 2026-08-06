"""Validation for AnalysisSpecification.

All checks are read-only — the specification is never mutated (schema
contract: "Validation MUST NOT modify AnalysisSpecification").
"""

from .schemas import AnalysisSpecification


def validate_analysis_specification(
    specification: AnalysisSpecification,
) -> list[str]:
    """Return a list of validation error strings (empty if valid)."""
    errors: list[str] = []

    if not specification.analysis_type or not specification.analysis_type.strip():
        errors.append("analysis_type must be provided")

    for ref in specification.factor_references:
        if not ref or not ref.strip():
            errors.append("factor_references must not contain empty strings")
            break

    for ref in specification.contrast_references:
        if not ref or not ref.strip():
            errors.append("contrast_references must not contain empty strings")
            break

    return errors
