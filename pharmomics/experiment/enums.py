"""PharmOmics experiment domain — controlled vocabularies.

Defines the stable enum types used by the experiment design schema.
These enums describe *what* groups, factors, and covariates represent
in an experiment, independent of any analysis backend.

See ADR 0004 (contrast semantics) and ADR 0005 (batch as covariate)
for the rationale behind these categories.
"""

from __future__ import annotations

from enum import StrEnum

# ---------------------------------------------------------------------------
# GroupRole — semantic role of an experimental group
# ---------------------------------------------------------------------------


class GroupRole(StrEnum):
    """Semantic role of an experimental group within an experiment.

    These roles describe the *experimental fact* of what a group
    represents (e.g. a treated population, a vehicle control), **not**
    the analytical direction of any comparison.  Comparison intent is
    expressed by ``Contrast.comparison_group_id`` and
    ``Contrast.reference_group_id`` (ADR 0004).
    """

    TREATMENT = "treatment"
    """Group receiving an active intervention (drug, compound, …)."""

    CONTROL = "control"
    """Group serving as a baseline (vehicle, sham, untreated, …)."""

    OBSERVATIONAL = "observational"
    """Group not receiving any intervention; observed as-is."""

    OTHER = "other"
    """Group that does not fit the above categories."""


# ---------------------------------------------------------------------------
# FactorType — what kind of values a factor takes
# ---------------------------------------------------------------------------


class FactorType(StrEnum):
    """What kind of values an experimental factor takes."""

    CATEGORICAL = "categorical"
    """Finite set of named levels (e.g. drug identity, cell line)."""

    CONTINUOUS = "continuous"
    """Numeric values on a continuum (e.g. dose amount, time)."""

    ORDINAL = "ordinal"
    """Ordered categories (e.g. disease stage I/II/III)."""


# ---------------------------------------------------------------------------
# CovariateValueType — what kind of values a covariate holds
# ---------------------------------------------------------------------------


class CovariateValueType(StrEnum):
    """What kind of values a covariate holds per sample."""

    CATEGORICAL = "categorical"
    """Named categories (e.g. batch, sex)."""

    CONTINUOUS = "continuous"
    """Numeric values (e.g. age, weight)."""

    ORDINAL = "ordinal"
    """Ordered categories (e.g. ECOG score)."""


# ---------------------------------------------------------------------------
# CovariateRole — why a covariate is tracked
# ---------------------------------------------------------------------------


class CovariateRole(StrEnum):
    """Why a covariate is tracked in the experiment.

    See ADR 0005 for the rationale: batch is a covariate role, not a
    separate ``BatchDesign`` type.
    """

    BATCH = "batch"
    """Technical batch identifier (sequencing run, library prep date, …)."""

    CLINICAL = "clinical"
    """Subject-level clinical variable (age, sex, diagnosis, …)."""

    TECHNICAL = "technical"
    """Technical measurement not tied to biology (RIN, read depth, …)."""

    OTHER = "other"
    """Covariate that does not fit the above categories."""
