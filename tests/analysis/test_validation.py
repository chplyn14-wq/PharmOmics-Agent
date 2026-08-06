"""Tests for AnalysisSpecification validation.

Covers:
- Valid specs with minimal and fully-populated fields.
- Empty / whitespace-only analysis_type.
- Empty reference lists (valid).
- Empty strings inside reference lists (invalid).
- Arbitrary parameters dict.
- Immutability (frozen dataclass).
- Validation purity (no mutation).
"""

from __future__ import annotations

import copy

from pharmomics.analysis.schemas import AnalysisSpecification
from pharmomics.analysis.validation import validate_analysis_specification

# ---------------------------------------------------------------------------
# Valid specifications
# ---------------------------------------------------------------------------


class TestValidSpecification:
    """Minimal and fully-populated valid AnalysisSpecification."""

    def test_minimal_valid(self) -> None:
        spec = AnalysisSpecification(analysis_type="differential_analysis")
        assert validate_analysis_specification(spec) == []

    def test_with_all_fields(self) -> None:
        spec = AnalysisSpecification(
            analysis_type="differential_analysis",
            factor_references=["drug", "time"],
            contrast_references=["treated_vs_vehicle", "high_vs_low"],
            parameters={"alpha": 0.01, "method": "limma"},
        )
        assert validate_analysis_specification(spec) == []

    def test_empty_factor_references_is_valid(self) -> None:
        spec = AnalysisSpecification(
            analysis_type="pathway_enrichment",
            factor_references=[],
        )
        assert validate_analysis_specification(spec) == []

    def test_empty_contrast_references_is_valid(self) -> None:
        spec = AnalysisSpecification(
            analysis_type="clustering",
            contrast_references=[],
        )
        assert validate_analysis_specification(spec) == []

    def test_arbitrary_parameters_accepted(self) -> None:
        spec = AnalysisSpecification(
            analysis_type="custom",
            parameters={
                "nested": {"key": "value"},
                "list": [1, 2, 3],
                "flag": True,
            },
        )
        assert validate_analysis_specification(spec) == []


# ---------------------------------------------------------------------------
# Invalid analysis_type
# ---------------------------------------------------------------------------


class TestInvalidAnalysisType:
    """Empty or whitespace-only analysis_type must be rejected."""

    def test_empty_analysis_type(self) -> None:
        spec = AnalysisSpecification(analysis_type="")
        errors = validate_analysis_specification(spec)
        assert any("analysis_type must be provided" in e for e in errors)

    def test_whitespace_only_analysis_type(self) -> None:
        spec = AnalysisSpecification(analysis_type="   ")
        errors = validate_analysis_specification(spec)
        assert any("analysis_type must be provided" in e for e in errors)


# ---------------------------------------------------------------------------
# Invalid reference lists
# ---------------------------------------------------------------------------


class TestInvalidReferences:
    """Empty strings inside reference lists must be rejected."""

    def test_empty_string_in_factor_references(self) -> None:
        spec = AnalysisSpecification(
            analysis_type="differential_analysis",
            factor_references=["drug", ""],
        )
        errors = validate_analysis_specification(spec)
        assert any(
            "factor_references must not contain empty strings" in e for e in errors
        )

    def test_whitespace_in_factor_references(self) -> None:
        spec = AnalysisSpecification(
            analysis_type="differential_analysis",
            factor_references=["drug", "  "],
        )
        errors = validate_analysis_specification(spec)
        assert any(
            "factor_references must not contain empty strings" in e for e in errors
        )

    def test_empty_string_in_contrast_references(self) -> None:
        spec = AnalysisSpecification(
            analysis_type="differential_analysis",
            contrast_references=["treated_vs_vehicle", ""],
        )
        errors = validate_analysis_specification(spec)
        assert any(
            "contrast_references must not contain empty strings" in e for e in errors
        )

    def test_whitespace_in_contrast_references(self) -> None:
        spec = AnalysisSpecification(
            analysis_type="differential_analysis",
            contrast_references=["treated_vs_vehicle", "  "],
        )
        errors = validate_analysis_specification(spec)
        assert any(
            "contrast_references must not contain empty strings" in e for e in errors
        )


# ---------------------------------------------------------------------------
# Immutability
# ---------------------------------------------------------------------------


class TestImmutability:
    """AnalysisSpecification is a frozen dataclass."""

    def test_frozen_dataclass_cannot_mutate(self) -> None:
        spec = AnalysisSpecification(analysis_type="differential_analysis")
        try:
            spec.analysis_type = "something_else"  # type: ignore[misc]
            assert False, "Expected frozen dataclass to raise"
        except (AttributeError, TypeError):
            pass  # Expected — the dataclass is frozen


# ---------------------------------------------------------------------------
# Validation purity
# ---------------------------------------------------------------------------


class TestValidationPurity:
    """Validation MUST NOT modify the specification (schema contract)."""

    def test_does_not_mutate_specification(self) -> None:
        spec = AnalysisSpecification(
            analysis_type="differential_analysis",
            factor_references=["drug"],
            contrast_references=["treated_vs_vehicle"],
            parameters={"alpha": 0.05},
        )
        original = copy.deepcopy(spec)
        validate_analysis_specification(spec)
        assert spec == original

    def test_does_not_mutate_on_invalid(self) -> None:
        spec = AnalysisSpecification(
            analysis_type="",
            factor_references=["", "drug"],
            contrast_references=["  "],
        )
        original = copy.deepcopy(spec)
        validate_analysis_specification(spec)
        assert spec == original
