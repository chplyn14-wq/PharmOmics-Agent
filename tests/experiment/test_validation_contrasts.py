"""Tests for contrast validation (C-01 … C-03)."""

from __future__ import annotations

from pharmomics.experiment.enums import GroupRole
from pharmomics.experiment.schemas import (
    Contrast,
    DesignSample,
    ExperimentalGroup,
    ExperimentDesign,
)
from pharmomics.experiment.validation import validate


class TestSameGroupBothSides:
    """C-01: comparison_group == reference_group is invalid."""

    def test_same_group_rejected(self) -> None:
        d = ExperimentDesign(
            experiment_id="exp",
            groups=[
                ExperimentalGroup(group_id="g1", label="G", role=GroupRole.CONTROL)
            ],
            contrasts=[
                Contrast(
                    contrast_id="c1",
                    comparison_group_id="g1",
                    reference_group_id="g1",
                )
            ],
        )
        errs = validate(d)
        assert any("same group on both sides" in e for e in errs)


class TestEmptyReferencedGroup:
    """C-02: referenced group has no samples."""

    def test_empty_comparison_group(self) -> None:
        d = ExperimentDesign(
            experiment_id="exp",
            groups=[
                ExperimentalGroup(
                    group_id="empty", label="E", role=GroupRole.TREATMENT
                ),
                ExperimentalGroup(group_id="g1", label="G", role=GroupRole.CONTROL),
            ],
            samples=[DesignSample(sample_id="s1", group_id="g1")],
            contrasts=[
                Contrast(
                    contrast_id="c1",
                    comparison_group_id="empty",
                    reference_group_id="g1",
                )
            ],
        )
        errs = validate(d)
        assert any("with no samples" in e for e in errs)

    def test_empty_reference_group(self) -> None:
        d = ExperimentDesign(
            experiment_id="exp",
            groups=[
                ExperimentalGroup(group_id="g1", label="G", role=GroupRole.TREATMENT),
                ExperimentalGroup(group_id="empty", label="E", role=GroupRole.CONTROL),
            ],
            samples=[DesignSample(sample_id="s1", group_id="g1")],
            contrasts=[
                Contrast(
                    contrast_id="c1",
                    comparison_group_id="g1",
                    reference_group_id="empty",
                )
            ],
        )
        errs = validate(d)
        assert any("with no samples" in e for e in errs)


class TestDuplicateContrast:
    """C-03: duplicate same-direction contrast."""

    def test_duplicate_same_direction(self) -> None:
        d = ExperimentDesign(
            experiment_id="exp",
            groups=[
                ExperimentalGroup(group_id="A", label="A", role=GroupRole.TREATMENT),
                ExperimentalGroup(group_id="B", label="B", role=GroupRole.CONTROL),
            ],
            contrasts=[
                Contrast(
                    contrast_id="c1", comparison_group_id="A", reference_group_id="B"
                ),
                Contrast(
                    contrast_id="c2", comparison_group_id="A", reference_group_id="B"
                ),
            ],
        )
        errs = validate(d)
        assert any("Duplicate contrast" in e for e in errs)

    def test_opposite_direction_allowed(self) -> None:
        """A vs B and B vs A are both allowed."""
        d = ExperimentDesign(
            experiment_id="exp",
            groups=[
                ExperimentalGroup(group_id="A", label="A", role=GroupRole.TREATMENT),
                ExperimentalGroup(group_id="B", label="B", role=GroupRole.CONTROL),
            ],
            contrasts=[
                Contrast(
                    contrast_id="c1", comparison_group_id="A", reference_group_id="B"
                ),
                Contrast(
                    contrast_id="c2", comparison_group_id="B", reference_group_id="A"
                ),
            ],
        )
        errs = validate(d)
        assert not any("Duplicate contrast" in e for e in errs)
