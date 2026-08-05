"""Tests for ExperimentDesign ↔ OmicsMatrix compatibility checks."""

from __future__ import annotations

from datetime import UTC, datetime

import pandas as pd

from pharmomics.compatibility.omics_design import check_compatibility
from pharmomics.experiment.enums import GroupRole
from pharmomics.experiment.schemas import (
    DesignSample,
    ExperimentalGroup,
    ExperimentDesign,
)
from pharmomics.omics.schemas import OmicsMatrix

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _make_omics(sample_ids: list[str]) -> OmicsMatrix:
    """Construct a minimal OmicsMatrix with the given sample_ids."""
    n_features = 3
    n_samples = len(sample_ids)
    df = pd.DataFrame(
        {"feature_id": ["G1", "G2", "G3"]}
        | {sid: [1.0, 2.0, 3.0] for sid in sample_ids},
    )
    return OmicsMatrix(
        matrix_id="mx-test",
        modality="transcriptomics",
        feature_type="gene",
        measurement_type="estimated_counts",
        normalization_status="raw",
        n_features=n_features,
        n_samples=n_samples,
        feature_ids=["G1", "G2", "G3"],
        sample_ids=sample_ids,
        dataframe=df,
        created_at=datetime.now(UTC).strftime("%Y-%m-%dT%H:%M:%SZ"),
    )


def _make_design(sample_ids: list[str]) -> ExperimentDesign:
    """Construct a minimal ExperimentDesign with the given sample_ids."""
    return ExperimentDesign(
        experiment_id="exp-test",
        samples=[DesignSample(sample_id=sid, group_id="g1") for sid in sample_ids],
        groups=[
            ExperimentalGroup(
                group_id="g1",
                label="Group 1",
                role=GroupRole.CONTROL,
            ),
        ],
    )


# ---------------------------------------------------------------------------
# Exact match
# ---------------------------------------------------------------------------


class TestExactMatch:
    """Design samples exactly match OmicsMatrix samples."""

    def test_exact_match_is_compatible(self) -> None:
        design = _make_design(["s1", "s2", "s3"])
        omics = _make_omics(["s1", "s2", "s3"])
        assert check_compatibility(design, omics) == []


# ---------------------------------------------------------------------------
# Design subset of OmicsMatrix
# ---------------------------------------------------------------------------


class TestDesignSubset:
    """Design samples are a subset of OmicsMatrix samples (valid)."""

    def test_subset_is_compatible(self) -> None:
        design = _make_design(["s1", "s2"])
        omics = _make_omics(["s1", "s2", "s3", "s4"])
        assert check_compatibility(design, omics) == []


# ---------------------------------------------------------------------------
# OmicsMatrix subset of design
# ---------------------------------------------------------------------------


class TestOmicsSubset:
    """OmicsMatrix samples are a subset of design samples (violation)."""

    def test_missing_design_samples(self) -> None:
        design = _make_design(["s1", "s2", "s3"])
        omics = _make_omics(["s1", "s2"])
        errs = check_compatibility(design, omics)
        assert any("not found in OmicsMatrix" in e for e in errs)


# ---------------------------------------------------------------------------
# Partial overlap
# ---------------------------------------------------------------------------


class TestPartialOverlap:
    """Some samples overlap, some do not."""

    def test_partial_overlap_reports_missing(self) -> None:
        design = _make_design(["s1", "s2", "s3"])
        omics = _make_omics(["s2", "s3", "s4"])
        errs = check_compatibility(design, omics)
        assert any("not found in OmicsMatrix" in e for e in errs)
        # s1 is missing, s2/s3 overlap
        assert any("s1" in e for e in errs)


# ---------------------------------------------------------------------------
# No overlap
# ---------------------------------------------------------------------------


class TestNoOverlap:
    """No samples in common."""

    def test_no_overlap_violation(self) -> None:
        design = _make_design(["a1", "a2"])
        omics = _make_omics(["b1", "b2"])
        errs = check_compatibility(design, omics)
        assert any("No overlap" in e for e in errs)


# ---------------------------------------------------------------------------
# Duplicate sample IDs in OmicsMatrix
# ---------------------------------------------------------------------------


class TestDuplicateOmicsSamples:
    """OmicsMatrix contains duplicate sample_ids."""

    def test_duplicate_omics_detected(self) -> None:
        design = _make_design(["s1"])
        omics = _make_omics(["s1", "s1", "s2"])
        errs = check_compatibility(design, omics)
        assert any("duplicate sample_id" in e.lower() for e in errs)


# ---------------------------------------------------------------------------
# Multi-modality
# ---------------------------------------------------------------------------


class TestMultiModality:
    """One design applied to multiple modalities."""

    def test_design_compatible_with_transcriptomics(self) -> None:
        design = _make_design(["s1", "s2"])
        tx = _make_omics(["s1", "s2", "s3"])
        assert check_compatibility(design, tx) == []

    def test_design_compatible_with_proteomics(self) -> None:
        design = _make_design(["s1", "s2"])
        px = _make_omics(["s1", "s2", "s4"])
        assert check_compatibility(design, px) == []


# ---------------------------------------------------------------------------
# Purity
# ---------------------------------------------------------------------------


class TestPurity:
    """Compatibility check does not mutate inputs."""

    def test_does_not_mutate_design(
        self,
    ) -> None:
        import copy

        design = _make_design(["s1", "s2"])
        omics = _make_omics(["s1", "s2"])
        original = copy.deepcopy(design)
        check_compatibility(design, omics)
        assert design == original

    def test_does_not_mutate_omics(
        self,
    ) -> None:

        design = _make_design(["s1"])
        omics = _make_omics(["s1", "s2"])
        original_ids = list(omics.sample_ids)
        check_compatibility(design, omics)
        assert omics.sample_ids == original_ids


# ---------------------------------------------------------------------------
# Empty cases
# ---------------------------------------------------------------------------


class TestEmptyCases:
    """Edge cases with empty designs or matrices."""

    def test_empty_design_compatible(self) -> None:
        design = _make_design([])
        omics = _make_omics(["s1", "s2"])
        assert check_compatibility(design, omics) == []

    def test_empty_omics_violation(self) -> None:
        design = _make_design(["s1"])
        omics = _make_omics([])
        errs = check_compatibility(design, omics)
        assert len(errs) > 0
