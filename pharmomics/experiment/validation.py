"""Validation rules for ``ExperimentDesign`` domain objects.

All validators operate purely on in-memory ``ExperimentDesign`` objects —
no file I/O, no network calls.  Each function returns a list of
human-readable violation strings (empty list = valid).

See the Phase 2B.0 schema contract and ADR 0006 for the purity contract:
``validate()`` does not mutate its input and performs no side effects.
"""

from __future__ import annotations

import math
from typing import TYPE_CHECKING

from pharmomics.experiment.enums import CovariateValueType, FactorType

if TYPE_CHECKING:
    from pharmomics.experiment.schemas import (
        CovariateDefinition,
        DesignSample,
        ExperimentalFactor,
        ExperimentDesign,
    )


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------


def validate(design: ExperimentDesign) -> list[str]:
    """Validate internal consistency of an *ExperimentDesign*.

    Parameters
    ----------
    design : ExperimentDesign
        The design to validate.

    Returns
    -------
    list[str]
        Empty list if the design is internally consistent; otherwise a
        list of human-readable violation descriptions.

    Notes
    -----
    Pure in-memory validation.  Does not mutate *design*.
    Performs no file I/O, no network calls.
    Does not check compatibility with any ``OmicsMatrix`` — use
    ``validate_design_against_omics()`` from ``compatibility.py``
    for cross-layer checks.
    """
    violations: list[str] = []
    violations.extend(_check_identifiers(design))
    violations.extend(_check_reference_integrity(design))
    violations.extend(_check_value_types(design))
    violations.extend(_check_contrasts(design))
    violations.extend(_check_treatments(design))
    violations.extend(_check_covariates(design))
    violations.extend(_check_replicates(design))
    violations.extend(_check_pairing(design))
    return violations


# ---------------------------------------------------------------------------
# Identifier uniqueness  (I-01 … I-12)
# ---------------------------------------------------------------------------


def _check_identifiers(design: ExperimentDesign) -> list[str]:
    """Ensure all identifiers are non-empty, non-whitespace, and unique."""
    violations: list[str] = []

    # I-12: experiment_id
    if not design.experiment_id.strip():
        violations.append("Empty experiment_id")

    # I-03: duplicate sample_id
    violations.extend(_check_unique_ids(
        [s.sample_id for s in design.samples],
        "sample_id",
    ))
    # I-01: empty sample_id
    for i, sample in enumerate(design.samples):
        if sample.sample_id == "":
            violations.append(f"Empty sample_id at index {i}")
    # I-02: whitespace-only sample_id
    for sample in design.samples:
        if sample.sample_id != "" and not sample.sample_id.strip():
            violations.append(
                f"Whitespace-only sample_id: {sample.sample_id!r}"
            )

    # I-04: duplicate group_id
    violations.extend(_check_unique_ids(
        [g.group_id for g in design.groups],
        "group_id",
    ))
    # I-08: empty group_id
    for g in design.groups:
        if not g.group_id.strip():
            violations.append("Empty group_id")

    # I-05: duplicate factor_id
    violations.extend(_check_unique_ids(
        [f.factor_id for f in design.factors],
        "factor_id",
    ))
    # I-09: empty factor_id
    for f in design.factors:
        if not f.factor_id.strip():
            violations.append("Empty factor_id")

    # I-06: duplicate covariate_id
    violations.extend(_check_unique_ids(
        [c.covariate_id for c in design.covariates],
        "covariate_id",
    ))
    # I-10: empty covariate_id
    for c in design.covariates:
        if not c.covariate_id.strip():
            violations.append("Empty covariate_id")

    # I-07: duplicate contrast_id
    violations.extend(_check_unique_ids(
        [ct.contrast_id for ct in design.contrasts],
        "contrast_id",
    ))
    # I-11: empty contrast_id
    for ct in design.contrasts:
        if not ct.contrast_id.strip():
            violations.append("Empty contrast_id")

    return violations


# ---------------------------------------------------------------------------
# Reference integrity  (R-01 … R-05)
# ---------------------------------------------------------------------------


def _check_reference_integrity(design: ExperimentDesign) -> list[str]:
    """Ensure all foreign-key references resolve to defined objects."""
    violations: list[str] = []

    group_ids = {g.group_id for g in design.groups}
    factor_ids = {f.factor_id for f in design.factors}
    covariate_ids = {c.covariate_id for c in design.covariates}

    # R-01: sample → group
    for sample in design.samples:
        if sample.group_id not in group_ids:
            violations.append(
                f"Sample {sample.sample_id!r} references unknown "
                f"group_id {sample.group_id!r}"
            )

    # R-02, R-03: contrast → groups
    for ct in design.contrasts:
        if ct.comparison_group_id not in group_ids:
            violations.append(
                f"Contrast {ct.contrast_id!r} references unknown "
                f"comparison_group_id {ct.comparison_group_id!r}"
            )
        if ct.reference_group_id not in group_ids:
            violations.append(
                f"Contrast {ct.contrast_id!r} references unknown "
                f"reference_group_id {ct.reference_group_id!r}"
            )

    # R-04: sample.factor_values keys → factor_ids
    for sample in design.samples:
        for key in sample.factor_values:
            if key not in factor_ids:
                violations.append(
                    f"Sample {sample.sample_id!r} has factor_values key "
                    f"{key!r} not in defined factors"
                )

    # R-05: sample.covariate_values keys → covariate_ids
    for sample in design.samples:
        for key in sample.covariate_values:
            if key not in covariate_ids:
                violations.append(
                    f"Sample {sample.sample_id!r} has covariate_values key "
                    f"{key!r} not in defined covariates"
                )

    return violations


# ---------------------------------------------------------------------------
# Value type consistency  (V-01 … V-05)
# ---------------------------------------------------------------------------


def _check_value_types(design: ExperimentDesign) -> list[str]:
    """Ensure sample values match their declared factor/covariate types."""
    violations: list[str] = []

    factor_map = {f.factor_id: f for f in design.factors}
    covariate_map = {c.covariate_id: c for c in design.covariates}

    # V-05: categorical factor must have levels defined
    for f in design.factors:
        if f.factor_type == FactorType.CATEGORICAL and f.levels is None:
            violations.append(
                f"Factor {f.factor_id!r} is categorical but has no "
                f"levels defined"
            )

    # Per-sample checks
    for sample in design.samples:
        for key, value in sample.factor_values.items():
            factor = factor_map.get(key)
            if factor is None:
                continue  # R-04 already catches unknown keys
            violations.extend(
                _check_value_against_factor(key, value, factor, sample.sample_id)
            )

        for key, value in sample.covariate_values.items():
            covariate = covariate_map.get(key)
            if covariate is None:
                continue  # R-05 already catches unknown keys
            violations.extend(
                _check_value_against_covariate(
                    key, value, covariate, sample.sample_id
                )
            )

    return violations


def _check_value_against_factor(
    factor_id: str,
    value: object,
    factor: ExperimentalFactor,
    sample_id: str,
) -> list[str]:
    """Check a single factor value against its declared type."""
    violations: list[str] = []

    if factor.factor_type == FactorType.CATEGORICAL:
        if not isinstance(value, str):
            violations.append(
                f"Sample {sample_id!r} factor {factor_id!r} value "
                f"{value!r} is not a string for type categorical"
            )
        elif factor.levels is not None and value not in factor.levels:
            violations.append(
                f"Sample {sample_id!r} factor {factor_id!r} value "
                f"{value!r} not in levels {factor.levels}"
            )

    elif factor.factor_type == FactorType.CONTINUOUS:
        if isinstance(value, bool) or not isinstance(value, (int, float)):
            violations.append(
                f"Sample {sample_id!r} factor {factor_id!r} value "
                f"{value!r} is not numeric for type continuous"
            )

    elif factor.factor_type == FactorType.ORDINAL:
        if not isinstance(value, str):
            violations.append(
                f"Sample {sample_id!r} factor {factor_id!r} value "
                f"{value!r} is not a string for type ordinal"
            )

    return violations


def _check_value_against_covariate(
    covariate_id: str,
    value: object,
    covariate: CovariateDefinition,
    sample_id: str,
) -> list[str]:
    """Check a single covariate value against its declared type."""
    violations: list[str] = []

    if covariate.value_type == CovariateValueType.CATEGORICAL:
        if not isinstance(value, str):
            violations.append(
                f"Sample {sample_id!r} covariate {covariate_id!r} value "
                f"{value!r} is not categorical"
            )

    elif covariate.value_type == CovariateValueType.CONTINUOUS:
        if isinstance(value, bool) or not isinstance(value, (int, float)):
            violations.append(
                f"Sample {sample_id!r} covariate {covariate_id!r} value "
                f"{value!r} is not numeric for type continuous"
            )

    elif covariate.value_type == CovariateValueType.ORDINAL:
        if not isinstance(value, str):
            violations.append(
                f"Sample {sample_id!r} covariate {covariate_id!r} value "
                f"{value!r} is not a string for type ordinal"
            )

    return violations


# ---------------------------------------------------------------------------
# Contrast semantics  (C-01 … C-03)
# ---------------------------------------------------------------------------


def _check_contrasts(design: ExperimentDesign) -> list[str]:
    """Validate contrast internal consistency."""
    violations: list[str] = []

    group_sample_count = _count_samples_per_group(design.samples)

    # C-01: comparison_group == reference_group
    for ct in design.contrasts:
        if ct.comparison_group_id == ct.reference_group_id:
            violations.append(
                f"Contrast {ct.contrast_id!r} has same group on both "
                f"sides: {ct.comparison_group_id!r}"
            )

    # C-02: referenced group has no samples
    for ct in design.contrasts:
        comp_count = group_sample_count.get(ct.comparison_group_id, 0)
        if comp_count == 0:
            violations.append(
                f"Contrast {ct.contrast_id!r} references comparison group "
                f"{ct.comparison_group_id!r} with no samples"
            )
        ref_count = group_sample_count.get(ct.reference_group_id, 0)
        if ref_count == 0:
            violations.append(
                f"Contrast {ct.contrast_id!r} references reference group "
                f"{ct.reference_group_id!r} with no samples"
            )

    # C-03: duplicate contrast with same (comparison, reference) pair
    seen: dict[tuple[str, str], str] = {}
    for ct in design.contrasts:
        pair = (ct.comparison_group_id, ct.reference_group_id)
        if pair in seen:
            violations.append(
                f"Duplicate contrast with same comparison "
                f"{ct.comparison_group_id!r} vs reference "
                f"{ct.reference_group_id!r}: {seen[pair]!r} and "
                f"{ct.contrast_id!r}"
            )
        else:
            seen[pair] = ct.contrast_id

    return violations


# ---------------------------------------------------------------------------
# Treatment validity  (T-01 … T-07)
# ---------------------------------------------------------------------------


def _check_treatments(design: ExperimentDesign) -> list[str]:
    """Validate treatment fields on samples."""
    violations: list[str] = []

    for sample in design.samples:
        if sample.treatment is None:
            continue
        t = sample.treatment
        sid = sample.sample_id

        # T-01: empty compound
        if not t.compound.strip():
            violations.append(
                f"Sample {sid!r} treatment has empty compound"
            )

        # T-02 … T-07: dose checks
        if t.dose is not None:
            val = t.dose.value
            unit = t.dose.unit

            if math.isnan(val):
                violations.append(
                    f"Sample {sid!r} treatment dose is NaN"
                )
            elif math.isinf(val):
                violations.append(
                    f"Sample {sid!r} treatment dose is infinity"
                )
            else:
                if val < 0:
                    violations.append(
                        f"Sample {sid!r} treatment dose value {val} "
                        f"is negative"
                    )
                if val == 0:
                    violations.append(
                        f"Sample {sid!r} treatment dose value is zero"
                    )

            if unit == "":
                violations.append(
                    f"Sample {sid!r} treatment dose has empty unit"
                )
            elif not unit.strip():
                violations.append(
                    f"Sample {sid!r} treatment dose unit is whitespace-only"
                )

    return violations


# ---------------------------------------------------------------------------
# Covariate special values  (CV-01 … CV-02)
# ---------------------------------------------------------------------------


def _check_covariates(design: ExperimentDesign) -> list[str]:
    """Check for NaN / infinity in continuous covariate values."""
    violations: list[str] = []

    covariate_map = {c.covariate_id: c for c in design.covariates}

    for sample in design.samples:
        for key, value in sample.covariate_values.items():
            covariate = covariate_map.get(key)
            if covariate is None:
                continue
            if covariate.value_type != CovariateValueType.CONTINUOUS:
                continue
            if isinstance(value, bool) or not isinstance(value, (int, float)):
                continue  # V-04 already catches type mismatch
            if math.isnan(value):
                violations.append(
                    f"Sample {sample.sample_id!r} covariate {key!r} "
                    f"value is NaN for continuous type"
                )
            elif math.isinf(value):
                violations.append(
                    f"Sample {sample.sample_id!r} covariate {key!r} "
                    f"value is infinity for continuous type"
                )

    return violations


# ---------------------------------------------------------------------------
# Replicate constraints  (RP-01 … RP-03)
# ---------------------------------------------------------------------------


def _check_replicates(design: ExperimentDesign) -> list[str]:
    """Replicate validation — no enforcement at schema level.

    RP-01: empty replicate IDs are allowed.
    RP-02: replicate IDs may repeat across groups.
    RP-03: no minimum replicate count is enforced.
    """
    return []


# ---------------------------------------------------------------------------
# Pairing structure  (P-01 … P-03)
# ---------------------------------------------------------------------------


def _check_pairing(design: ExperimentDesign) -> list[str]:
    """Validate pairing consistency."""
    violations: list[str] = []

    # P-01: pairing defined but no samples have pair_id
    if design.pairing is not None:
        paired_samples = [s for s in design.samples if s.pair_id is not None]
        if not paired_samples:
            violations.append(
                "Pairing is defined but no samples have a pair_id"
            )

        # P-02: each pair_id must have at least 2 samples
        pair_counts = _count_by_pair_id(design.samples)
        for pair_id, count in sorted(pair_counts.items()):
            if count < 2:
                violations.append(
                    f"Pair {pair_id!r} has only {count} sample"
                )

    return violations


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------


def _check_unique_ids(ids: list[str], label: str) -> list[str]:
    """Return violation strings for duplicate values in *ids*."""
    dups = _find_duplicates(ids)
    if not dups:
        return []
    examples = sorted(dups)[:5]
    return [f"Duplicate {label}: {examples}"]


def _find_duplicates(ids: list[str]) -> set[str]:
    """Return the set of duplicated values in *ids*."""
    seen: set[str] = set()
    dups: set[str] = set()
    for item in ids:
        if item in seen:
            dups.add(item)
        seen.add(item)
    return dups


def _count_samples_per_group(
    samples: list[DesignSample],
) -> dict[str, int]:
    """Count how many samples belong to each group_id."""
    counts: dict[str, int] = {}
    for s in samples:
        counts[s.group_id] = counts.get(s.group_id, 0) + 1
    return counts


def _count_by_pair_id(
    samples: list[DesignSample],
) -> dict[str, int]:
    """Count how many samples share each pair_id (excluding None)."""
    counts: dict[str, int] = {}
    for s in samples:
        if s.pair_id is not None:
            counts[s.pair_id] = counts.get(s.pair_id, 0) + 1
    return counts
