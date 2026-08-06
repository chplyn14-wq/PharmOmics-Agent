from .schemas import AnalysisSpecification


def validate_analysis_specification(
    specification: AnalysisSpecification,
) -> list[str]:
    errors: list[str] = []

    if not specification.analysis_type:
        errors.append("analysis_type must be provided")

    return errors
