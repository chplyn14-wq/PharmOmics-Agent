"""Compatibility checks between ExperimentDesign and OmicsMatrix.

Ensures that an experiment design's declared samples are compatible
with the samples present in an omics data matrix.  All checks use
``OmicsMatrix.sample_ids`` as the source of truth — they do not
inspect the underlying dataframe columns.

See ADR 0003 for the rationale: design and matrix are independent
domain objects aligned only through stable sample_id strings.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from pharmomics.experiment.schemas import ExperimentDesign
    from pharmomics.omics.schemas import OmicsMatrix


def check_compatibility(
    design: ExperimentDesign,
    omics: OmicsMatrix,
) -> list[str]:
    """Check whether *design* samples are compatible with *omics* samples.

    Parameters
    ----------
    design : ExperimentDesign
        The experimental design to check.
    omics : OmicsMatrix
        The omics matrix to check against.  ``omics.sample_ids`` is
        the source of truth for sample identity.

    Returns
    -------
    list[str]
        Empty list if the design is compatible with the matrix;
        otherwise a list of human-readable violation descriptions.

    Notes
    -----
    Pure function — does not mutate either input.
    Does not import or depend on ``pharmomics.experiment.validation``.
    """
    violations: list[str] = []
    violations.extend(_check_omics_duplicate_samples(omics))
    violations.extend(_check_design_samples_in_omics(design, omics))
    return violations


# ---------------------------------------------------------------------------
# Individual checks
# ---------------------------------------------------------------------------


def _check_omics_duplicate_samples(omics: OmicsMatrix) -> list[str]:
    """Detect duplicate sample_ids within the OmicsMatrix."""
    seen: set[str] = set()
    dups: set[str] = set()
    for sid in omics.sample_ids:
        if sid in seen:
            dups.add(sid)
        seen.add(sid)

    if not dups:
        return []
    examples = sorted(dups)[:5]
    return [f"OmicsMatrix contains duplicate sample_id: {examples}"]


def _check_design_samples_in_omics(
    design: ExperimentDesign,
    omics: OmicsMatrix,
) -> list[str]:
    """Ensure every design sample_id exists in omics.sample_ids."""
    omics_ids = set(omics.sample_ids)
    design_ids = {s.sample_id for s in design.samples}

    missing = design_ids - omics_ids
    if not missing:
        return []

    # Special case: completely no overlap
    overlap = design_ids & omics_ids
    if not overlap and omics_ids:
        return [
            "No overlap between ExperimentDesign and OmicsMatrix samples"
        ]

    examples = sorted(missing)[:5]
    return [
        f"Design sample(s) not found in OmicsMatrix: {examples}"
    ]
