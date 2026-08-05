"""Pydantic schemas for the PharmOmics experiment domain.

Defines the declarative domain objects used to describe an experimental
design: samples, groups, factors, treatments, covariates, pairing, and
contrasts.

These schemas are **purely declarative** — they carry no validation
state, workflow state, or analysis configuration.  Validation is
performed by a separate ``validate()`` function (Phase 2B.2).

See the Phase 2B.0 schema contract and ADR 0003–0006 for the frozen
field definitions and semantic rules.
"""

from __future__ import annotations

from pydantic import BaseModel, ConfigDict, Field

from pharmomics.experiment.enums import (
    CovariateRole,
    CovariateValueType,
    FactorType,
    GroupRole,
)

# ---------------------------------------------------------------------------
# Type aliases
# ---------------------------------------------------------------------------

ScalarValue = str | int | float | bool
"""JSON-compatible scalar used in ``factor_values`` and
``covariate_values`` on ``DesignSample``."""

JSONScalar = str | int | float | bool | None
"""Any JSON-compatible value used in extension ``metadata`` dicts."""


# ---------------------------------------------------------------------------
# Quantity — value object for dose / concentration
# ---------------------------------------------------------------------------


class Quantity(BaseModel):
    """A measured quantity with a numeric value and a unit.

    Used for treatment dose and other values where magnitude and unit
    both carry meaning.  No unit conversion or equivalence checking is
    performed at the schema level.
    """

    model_config = ConfigDict(frozen=True)

    value: float
    """Numeric magnitude (e.g. ``500.0``)."""

    unit: str
    """Unit string (e.g. ``"nM"``, ``"mg/kg"``).  Stored as-is; no
    canonicalisation or conversion is applied."""


# ---------------------------------------------------------------------------
# Treatment
# ---------------------------------------------------------------------------


class Treatment(BaseModel):
    """Treatment applied to a sample.

    Records the compound and dose/duration that a sample was exposed to.
    This is a value object owned by ``DesignSample``.
    """

    model_config = ConfigDict(frozen=True)

    compound: str
    """Name of the drug, compound, or intervention."""

    dose: Quantity | None = None
    """Dose amount with units (e.g. ``Quantity(value=500, unit="nM")``)."""

    duration: str | None = None
    """Exposure duration with units (e.g. ``"24h"``, ``"72h"``)."""

    description: str | None = None
    """Human-readable description of the treatment protocol."""


# ---------------------------------------------------------------------------
# ExperimentalGroup
# ---------------------------------------------------------------------------


class ExperimentalGroup(BaseModel):
    """A named experimental condition group.

    Groups are referenced by ``DesignSample.group_id`` and by
    ``Contrast.comparison_group_id`` / ``Contrast.reference_group_id``.
    """

    model_config = ConfigDict(frozen=True)

    group_id: str
    """Unique identifier within the experiment."""

    label: str
    """Short human-readable label."""

    description: str | None = None
    """Longer description of the group."""

    role: GroupRole
    """The semantic role of this group (treatment, control, …).
    This expresses the experimental fact of what the group represents,
    **not** the analytical direction of a comparison (ADR 0004)."""


# ---------------------------------------------------------------------------
# ExperimentalFactor
# ---------------------------------------------------------------------------


class ExperimentalFactor(BaseModel):
    """A manipulated variable in the experiment.

    Examples: drug identity, time point, dose level.
    """

    model_config = ConfigDict(frozen=True)

    factor_id: str
    """Unique identifier within the experiment."""

    factor_type: FactorType
    """Whether the factor is categorical, continuous, or ordinal."""

    description: str | None = None

    levels: list[str] | None = None
    """Allowed level names for categorical / ordinal factors.
    ``None`` for continuous factors."""


# ---------------------------------------------------------------------------
# CovariateDefinition
# ---------------------------------------------------------------------------


class CovariateDefinition(BaseModel):
    """A measured variable that may confound or modulate results.

    See ADR 0005: batch is represented as a covariate with
    ``role=CovariateRole.BATCH``, not as a separate ``BatchDesign``.
    """

    model_config = ConfigDict(frozen=True)

    covariate_id: str
    """Unique identifier within the experiment."""

    role: CovariateRole
    """Why this covariate is tracked (batch, clinical, technical, …)."""

    value_type: CovariateValueType
    """Whether values are categorical, continuous, or ordinal."""

    description: str | None = None


# ---------------------------------------------------------------------------
# DesignSample
# ---------------------------------------------------------------------------


class DesignSample(BaseModel):
    """One experimental sample.

    Each sample belongs to exactly one group, has values for each
    experimental factor, and optionally carries treatment, replicate,
    covariate, and pairing information.
    """

    model_config = ConfigDict(frozen=True)

    sample_id: str
    """Unique identifier within the experiment."""

    group_id: str
    """References ``ExperimentalGroup.group_id``."""

    factor_values: dict[str, ScalarValue] = Field(default_factory=dict)
    """Values keyed by ``ExperimentalFactor.factor_id``."""

    treatment: Treatment | None = None

    biological_replicate: str | None = None
    """Replicate group identifier (string, not auto-inferred)."""

    technical_replicate: str | None = None

    covariate_values: dict[str, ScalarValue] = Field(default_factory=dict)
    """Values keyed by ``CovariateDefinition.covariate_id``.
    Batch is included here when ``role=CovariateRole.BATCH``."""

    pair_id: str | None = None
    """Pairing group identifier.  ``None`` means the sample is unpaired."""


# ---------------------------------------------------------------------------
# Contrast
# ---------------------------------------------------------------------------


class Contrast(BaseModel):
    """A comparison between two experimental groups.

    See ADR 0004: comparison/reference semantics.

    A **positive effect** means:
    ``comparison_group`` expression > ``reference_group`` expression.

    A vs B and B vs A are two distinct ``Contrast`` objects.
    """

    model_config = ConfigDict(frozen=True)

    contrast_id: str
    """Unique identifier within the experiment."""

    comparison_group_id: str
    """The 'A' in A vs B.  Must reference an ``ExperimentalGroup.group_id``."""

    reference_group_id: str
    """The 'B' in A vs B.  Must reference an ``ExperimentalGroup.group_id``."""

    description: str | None = None
    """Human-readable description of the comparison intent."""


# ---------------------------------------------------------------------------
# PairingDefinition
# ---------------------------------------------------------------------------


class PairingDefinition(BaseModel):
    """Describes the paired-sample structure of an experiment.

    Examples: before/after measurements on the same subject, matched
    tumour/normal pairs.
    """

    model_config = ConfigDict(frozen=True)

    pairing_type: str
    """Type of pairing, e.g. ``"before_after"``, ``"matched"``."""

    description: str | None = None


# ---------------------------------------------------------------------------
# ExperimentDesign
# ---------------------------------------------------------------------------


class ExperimentDesign(BaseModel):
    """Aggregate root: a complete declarative description of an experiment.

    Owns samples, groups, factors, contrasts, covariates, and pairing.
    Does **not** hold or reference any ``OmicsMatrix`` (ADR 0003).

    This model is **purely declarative** — it carries no validation
    state, workflow state, or analysis configuration (ADR 0006).
    """

    model_config = ConfigDict(frozen=True)

    schema_version: str = "1.0.0"
    """Schema version, e.g. ``"1.0.0"``."""

    experiment_id: str
    """Unique experiment identifier."""

    description: str | None = None
    """Human-readable description of the experiment."""

    samples: list[DesignSample] = Field(default_factory=list)
    """All samples in the experiment."""

    groups: list[ExperimentalGroup] = Field(default_factory=list)
    """Group definitions."""

    factors: list[ExperimentalFactor] = Field(default_factory=list)
    """Experimental factors."""

    contrasts: list[Contrast] = Field(default_factory=list)
    """Comparison intents.  May be empty."""

    covariates: list[CovariateDefinition] = Field(default_factory=list)
    """Covariate definitions (including batch)."""

    pairing: PairingDefinition | None = None
    """Paired-sample structure, if applicable."""
