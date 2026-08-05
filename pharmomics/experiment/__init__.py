"""PharmOmics experiment design domain.

Phase 2B.1 provides the declarative schema types for describing
experimental designs: samples, groups, factors, treatments, covariates,
pairing, and contrasts.

These types are independent of the omics data layer (ADR 0003) and
carry no validation or workflow state (ADR 0006).
"""

from pharmomics.experiment.enums import (
    CovariateRole,
    CovariateValueType,
    FactorType,
    GroupRole,
)
from pharmomics.experiment.schemas import (
    Contrast,
    CovariateDefinition,
    DesignSample,
    ExperimentalFactor,
    ExperimentalGroup,
    ExperimentDesign,
    JSONScalar,
    PairingDefinition,
    Quantity,
    ScalarValue,
    Treatment,
)

__all__ = [
    # Enums
    "CovariateRole",
    "CovariateValueType",
    "FactorType",
    "GroupRole",
    # Schema types
    "Contrast",
    "CovariateDefinition",
    "DesignSample",
    "ExperimentDesign",
    "ExperimentalFactor",
    "ExperimentalGroup",
    "JSONScalar",
    "PairingDefinition",
    "Quantity",
    "ScalarValue",
    "Treatment",
]
