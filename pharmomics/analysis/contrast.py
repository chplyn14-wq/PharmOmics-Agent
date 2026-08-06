"""Resolve a contrast to its comparison and reference sample IDs.

Provides a minimal resolver that takes an ``ExperimentDesign`` and a
``contrast_id``, then returns the corresponding sample IDs for the
comparison and reference groups.

No OmicsMatrix is read.  No statistics are computed.  The
``ExperimentDesign`` is not modified.
"""

from __future__ import annotations

from dataclasses import dataclass, field

from pharmomics.experiment.schemas import ExperimentDesign


class ContrastResolutionError(ValueError):
    """Raised when a contrast_id cannot be resolved to sample groups."""


@dataclass(frozen=True)
class ResolvedContrast:
    """The resolved sample sets for one contrast."""

    contrast_id: str
    """The contrast that was resolved."""

    comparison_group_id: str
    """The group ID used as the comparison (numerator) group."""

    reference_group_id: str
    """The group ID used as the reference (denominator) group."""

    comparison_sample_ids: tuple[str, ...] = field(default_factory=tuple)
    """Sample IDs whose ``group_id`` matches ``comparison_group_id``."""

    reference_sample_ids: tuple[str, ...] = field(default_factory=tuple)
    """Sample IDs whose ``group_id`` matches ``reference_group_id``."""


def resolve_contrast(
    design: ExperimentDesign,
    contrast_id: str,
) -> ResolvedContrast:
    """Resolve a single contrast_id to comparison and reference sample IDs.

    Parameters
    ----------
    design : ExperimentDesign
        The experiment design containing samples, groups, and contrasts.
    contrast_id : str
        The contrast to resolve.  Must match exactly one
        ``Contrast.contrast_id`` in ``design.contrasts``.

    Returns
    -------
    ResolvedContrast
        Frozen result with group IDs and sample ID tuples.

    Raises
    ------
    ContrastResolutionError
        When the contrast, group, or samples cannot be resolved, or
        when comparison and reference sample sets overlap.
    """
    exp_id = design.experiment_id

    # --- Find the contrast ------------------------------------------------
    contrast = None
    for c in design.contrasts:
        if c.contrast_id == contrast_id:
            contrast = c
            break
    if contrast is None:
        raise ContrastResolutionError(
            f"Contrast '{contrast_id}' not found in experiment '{exp_id}'"
        )

    # --- Find groups ------------------------------------------------------
    group_ids = {g.group_id for g in design.groups}

    if contrast.comparison_group_id not in group_ids:
        raise ContrastResolutionError(
            f"Comparison group '{contrast.comparison_group_id}' not found "
            f"in experiment '{exp_id}'"
        )

    if contrast.reference_group_id not in group_ids:
        raise ContrastResolutionError(
            f"Reference group '{contrast.reference_group_id}' not found "
            f"in experiment '{exp_id}'"
        )

    # --- Collect samples --------------------------------------------------
    comparison_ids: list[str] = []
    reference_ids: list[str] = []

    for sample in design.samples:
        if sample.group_id == contrast.comparison_group_id:
            comparison_ids.append(sample.sample_id)
        elif sample.group_id == contrast.reference_group_id:
            reference_ids.append(sample.sample_id)

    if not comparison_ids:
        raise ContrastResolutionError(
            f"No samples found for comparison group "
            f"'{contrast.comparison_group_id}' in experiment '{exp_id}'"
        )

    if not reference_ids:
        raise ContrastResolutionError(
            f"No samples found for reference group "
            f"'{contrast.reference_group_id}' in experiment '{exp_id}'"
        )

    # --- Check overlap (defensive) ----------------------------------------
    comparison_set = set(comparison_ids)
    reference_set = set(reference_ids)
    overlap = comparison_set & reference_set
    if overlap:
        raise ContrastResolutionError(
            f"Sample(s) {sorted(overlap)} appear in both comparison and "
            f"reference groups for contrast '{contrast_id}' in experiment "
            f"'{exp_id}'"
        )

    return ResolvedContrast(
        contrast_id=contrast.contrast_id,
        comparison_group_id=contrast.comparison_group_id,
        reference_group_id=contrast.reference_group_id,
        comparison_sample_ids=tuple(comparison_ids),
        reference_sample_ids=tuple(reference_ids),
    )
