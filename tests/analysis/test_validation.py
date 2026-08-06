from pharmomics.analysis.schemas import AnalysisSpecification
from pharmomics.analysis.validation import validate_analysis_specification


def test_valid_analysis_specification():
    specification = AnalysisSpecification(
        analysis_type="differential_analysis"
    )

    assert validate_analysis_specification(specification) == []


def test_invalid_analysis_specification():
    specification = AnalysisSpecification(
        analysis_type=""
    )

    assert validate_analysis_specification(specification) == [
        "analysis_type must be provided"
    ]
