# ExperimentDesign Schema Contract (Phase 2B.0 — Frozen)

**Status:** Frozen
**Date:** 2026-08-05

This document defines the canonical schema for `ExperimentDesign` and
its nested types. These contracts are **frozen** — changes require a
new ADR and explicit architect approval.

Related ADRs:
- [ADR 0003](adr/0003-design-matrix-separation.md): Design and matrix are independent.
- [ADR 0004](adr/0004-contrast-semantics.md): Contrast uses comparison/reference semantics.
- [ADR 0005](adr/0005-batch-as-covariate.md): Batch is a covariate role.
- [ADR 0006](adr/0006-experiment-design-pure.md): ExperimentDesign is purely declarative.

---

## Aggregate Root: ExperimentDesign

ExperimentDesign is the aggregate root. It owns samples, groups,
factors, contrasts, covariates, and pairing definitions. It does **not**
own or reference OmicsMatrix.

```python
class ExperimentDesign(BaseModel):
    model_config = ConfigDict(frozen=True)

    schema_version: str                       # e.g. "1.0.0"
    experiment_id: str                        # unique identifier
    description: str | None                   # human-readable description

    samples: list[DesignSample]               # all samples in the experiment
    groups: list[ExperimentalGroup]           # group definitions
    factors: list[ExperimentalFactor]         # experimental factors
    contrasts: list[Contrast]                 # comparison intents (may be empty)
    covariates: list[CovariateDefinition]     # covariates including batch
    pairing: PairingDefinition | None         # paired-sample structure
    metadata: dict[str, JSONScalar] | None    # non-core extensions only
```

### Constraints

- **Aggregate root.** ExperimentDesign owns all nested value objects.
- **No OmicsMatrix.** Does not hold, reference, or embed numeric data.
  Cross-layer alignment with OmicsMatrix is checked via a future
  `validate_design_against_omics()` function using `sample_id`.
- **Empty contrasts allowed.** `contrasts` may be `[]` when no
  comparison has been defined yet.
- **Metadata is for extensions only.** Core semantic fields must not
  hide inside `metadata`.
- **Immutability.** `frozen=True` — no mutation after construction.

---

## DesignSample

One row per experimental sample.

```python
class DesignSample(BaseModel):
    model_config = ConfigDict(frozen=True)

    sample_id: str                            # unique within this experiment
    group_id: str                             # references ExperimentalGroup.group_id
    factor_values: dict[str, ScalarValue]     # values for each ExperimentalFactor
    treatment: Treatment | None               # treatment details (if applicable)
    biological_replicate: str | None          # replicate group identifier
    technical_replicate: str | None           # technical replicate identifier
    covariate_values: dict[str, ScalarValue]  # all covariates including batch
    pair_id: str | None                       # pairing group identifier
```

### Constraints

- **One primary group_id.** A sample belongs to exactly one group in its
  current version. Multi-group membership is not supported in v1.
- **Replicate IDs are strings.** No automatic inference from names.
- **pair_id is a grouping key.** It identifies which pairing unit a
  sample belongs to, not the other members of the pair.
- **No treatment_sample_id / control_sample_id.** The PairGroup pattern
  of labeling samples as treatment vs control is not used.
- **No auto-inference.** Group, pair, and treatment are explicit — never
  derived from sample name patterns.
- **covariate_values is the single covariate channel.** Batch and other
  covariates share this dict (see ADR 0005).

---

## ExperimentalGroup

Defines a named experimental condition group.

```python
class ExperimentalGroup(BaseModel):
    model_config = ConfigDict(frozen=True)

    group_id: str                             # unique within this experiment
    label: str                                # short human-readable label
    description: str | None                   # longer description
```

---

## ExperimentalFactor

A manipulated variable in the experiment (e.g. drug, time, dose).

```python
class ExperimentalFactor(BaseModel):
    model_config = ConfigDict(frozen=True)

    factor_id: str                            # unique within this experiment
    factor_type: FactorType                   # categorical | continuous | ordinal
    description: str | None
    levels: list[str] | None                  # for categorical factors
```

---

## Contrast

Defines a comparison between two groups. See ADR 0004 for full rationale.

```python
class Contrast(BaseModel):
    model_config = ConfigDict(frozen=True)

    contrast_id: str                          # unique within this experiment
    comparison_group_id: str                  # the "A" in A vs B
    reference_group_id: str                   # the "B" in A vs B
    description: str | None
```

### Semantic rule

A **positive effect** means:

> `comparison_group` expression > `reference_group` expression

A vs B and B vs A are two different Contrast objects.

---

## CovariateDefinition

A measured variable that may confound or modulate results.
See ADR 0005 for batch-as-covariate rationale.

```python
class CovariateDefinition(BaseModel):
    model_config = ConfigDict(frozen=True)

    covariate_id: str                         # unique within this experiment
    role: CovariateRole                       # batch | clinical | technical | other
    value_type: CovariateValueType            # categorical | continuous | ordinal
    description: str | None
```

---

## PairingDefinition

Describes the paired-sample structure (e.g. before/after, matched
tumour/normal).

```python
class PairingDefinition(BaseModel):
    model_config = ConfigDict(frozen=True)

    pairing_type: str                         # e.g. "before_after", "matched"
    description: str | None
```

---

## Treatment

Treatment applied to a sample.

```python
class Treatment(BaseModel):
    model_config = ConfigDict(frozen=True)

    compound: str                             # drug or compound name
    dose: str | None                          # dose with units, e.g. "500nM"
    duration: str | None                      # e.g. "24h", "72h"
    description: str | None
```

---

## Value Types

```python
# ScalarValue: any JSON-compatible scalar used in factor_values and
# covariate_values.
ScalarValue = str | int | float | bool

# JSONScalar: any JSON-compatible value used in metadata extensions.
JSONScalar = str | int | float | bool | None
```

---

## Validation Contract

```python
def validate(design: ExperimentDesign) -> list[str]:
    """Validate an ExperimentDesign.

    Parameters
    ----------
    design : ExperimentDesign
        The design to validate.

    Returns
    -------
    list[str]
        Empty list if valid; otherwise violation descriptions.

    Notes
    -----
    Pure in-memory validation. Does not modify the input.
    Does not write results back to the schema.
    """
```

Validation rules (non-exhaustive, to be refined in Phase 2B.1):

- All `group_id` references in samples exist in `groups`.
- All `comparison_group_id` and `reference_group_id` in contrasts exist
  in `groups`.
- All `factor_values` keys reference valid `factor_id` values.
- `sample_id` values are unique.
- If `pairing` is defined, at least two samples share a `pair_id`.
- `covariate_values` keys reference valid `covariate_id` values.

---

## Type Reference Table

| Type | Owns | Referenced by |
|---|---|---|
| ExperimentDesign | samples, groups, factors, contrasts, covariates, pairing | — (aggregate root) |
| DesignSample | treatment (optional) | ExperimentDesign.samples |
| ExperimentalGroup | — | ExperimentDesign.groups, DesignSample.group_id |
| ExperimentalFactor | — | ExperimentDesign.factors, DesignSample.factor_values keys |
| Contrast | — | ExperimentDesign.contrasts |
| CovariateDefinition | — | ExperimentDesign.covariates, DesignSample.covariate_values keys |
| PairingDefinition | — | ExperimentDesign.pairing, DesignSample.pair_id |
| Treatment | — | DesignSample.treatment |

---

## Alignment with OmicsMatrix (future)

```python
def validate_design_against_omics(
    design: ExperimentDesign,
    omics: OmicsMatrix,
) -> list[str]:
    """Check that the design's sample_ids are compatible with the
    matrix's sample_ids. Returns violation strings.

    This function does NOT modify either input. It is not a method
    on ExperimentDesign or OmicsMatrix — it is a standalone validator
    that reads both objects.

    Typical checks:
    - Every sample_id in design.samples appears in omics.sample_ids.
    - No duplicate sample_ids across the two objects.
    """
```

See [ADR 0003](adr/0003-design-matrix-separation.md) for rationale.

---

## Module Placement

### Source package

The `experiment` package is a **top-level sibling** to `omics/`, `ingestion/`,
`agents/`, `run_store/`, and `config.py`. It is **not** a submodule of any
existing package.

```
pharmomics/
├── __init__.py
├── main.py
├── config.py
├── omics/                  # Phase 2A — existing
│   ├── __init__.py
│   ├── adapter.py
│   ├── enums.py
│   ├── schemas.py
│   └── validation.py
├── ingestion/              # Phase 1 — existing
│   ├── __init__.py
│   ├── loader.py
│   └── schemas.py
├── agents/                 # Phase 1 — existing
│   └── ...
├── run_store/              # Phase 1 — existing
│   └── ...
└── experiment/             # Phase 2B — NEW top-level package
    ├── __init__.py
    ├── enums.py            # stable domain enums
    ├── schemas.py          # pure domain schema only
    ├── validation.py       # ExperimentDesign internal consistency
    └── compatibility.py    # ExperimentDesign ↔ OmicsMatrix cross-layer
```

### Explicit placement rules

| Rule | Decision |
|---|---|
| `experiment` as a top-level package | **Required** — sibling to `omics/` |
| Place in `pharmomics/omics/` | **Forbidden** — omics and experiment are separate domains |
| Place in `pharmomics/analysis/` | **Forbidden** — analysis is a separate downstream consumer |
| Create `selectors.py` | **Forbidden** — no selector/query language in this phase |
| `schemas.py` contents | **Pure domain schema only** — no validation logic, no adapters |
| `enums.py` contents | **Stable domain enums** — FactorType, CovariateRole, CovariateValueType |
| `validation.py` contents | **ExperimentDesign internal consistency** — identifier uniqueness, reference integrity, constraint checks |
| `compatibility.py` contents | **Cross-layer validation** — `validate_design_against_omics(design, matrix) -> list[str]` |
| Create Python files this phase | **Not allowed** — position frozen, implementation deferred to Phase 2B.1 |

### Test package

Tests live under `tests/experiment/`, parallel to the source package.

```
tests/
├── ... existing Phase 1/2A tests
└── experiment/
    ├── conftest.py
    ├── test_schemas.py
    ├── test_validation_identifiers.py
    ├── test_validation_references.py
    ├── test_validation_groups.py
    ├── test_validation_contrasts.py
    ├── test_validation_treatments.py
    ├── test_validation_covariates.py
    ├── test_validation_replicates.py
    ├── test_validation_pairing.py
    └── test_omics_compatibility.py
```

### File-to-responsibility mapping

| File | Responsibility |
|---|---|
| `enums.py` | `FactorType`, `CovariateRole`, `CovariateValueType` |
| `schemas.py` | `ExperimentDesign`, `DesignSample`, `ExperimentalGroup`, `ExperimentalFactor`, `Contrast`, `CovariateDefinition`, `PairingDefinition`, `Treatment` |
| `validation.py` | `validate(ExperimentDesign) -> list[str]` and internal `_check_*` helpers |
| `compatibility.py` | `validate_design_against_omics(ExperimentDesign, OmicsMatrix) -> list[str]` |
| `conftest.py` | Shared fixtures for all experiment tests |

---

## Testing Strategy

### Principles

1. **Each validation rule gets at least one valid case and one invalid case.**
2. **No monolithic fixture.** Each test receives a minimal, composable fixture.
3. **Error strings must be stable and predictable.** Tests assert on specific
   substring patterns, not just "non-empty list".
4. **Validation is pure.** Tests verify that the input object is not mutated
   after calling `validate()`.
5. **Determinism.** Same input always produces the same ordered error list.

### Fixture strategy

Fixtures are provided in `tests/experiment/conftest.py`. Each fixture is
small, focused, and composable.

| Fixture | Purpose |
|---|---|
| `minimal_valid_design` | Smallest valid ExperimentDesign: 1 sample, 1 group, no contrasts |
| `simple_drug_control_design` | 2 groups (drug, vehicle), 3 samples each, 1 contrast |
| `batch_adjusted_design` | 2 batches as covariates, samples distributed across batches |
| `paired_before_after_design` | 3 pairs, each with before/after samples |
| `multi_factor_design` | 2 factors (drug, timepoint), factorial group structure |
| `invalid_reference_design` | Design with broken references (unknown group_id, etc.) |

**Rules:**
- No single fixture that covers every scenario.
- Tests compose by modifying a base fixture, not by reaching for a giant one.
- Property-based testing (e.g. Hypothesis) is **optional** and not required for
  Phase 2B initial implementation.

### Schema tests

| Test | File |
|---|---|
| minimal valid construction | `test_schemas.py` |
| full valid construction | `test_schemas.py` |
| optional fields (None allowed) | `test_schemas.py` |
| default values | `test_schemas.py` |
| enum validation (invalid enum values rejected) | `test_schemas.py` |
| invalid field types (wrong type rejected) | `test_schemas.py` |
| mutable default isolation (mutating one instance does not affect another) | `test_schemas.py` |
| serialization round-trip (`model_dump()` → rebuild → equality) | `test_schemas.py` |

### Identifier validation

| Test | File |
|---|---|
| empty sample_id rejected | `test_validation_identifiers.py` |
| whitespace-only sample_id rejected | `test_validation_identifiers.py` |
| duplicate sample_id rejected | `test_validation_identifiers.py` |
| duplicate group_id rejected | `test_validation_identifiers.py` |
| duplicate factor_id rejected | `test_validation_identifiers.py` |
| duplicate covariate_id rejected | `test_validation_identifiers.py` |
| duplicate contrast_id rejected | `test_validation_identifiers.py` |
| IDs are case-sensitive ("A" ≠ "a") | `test_validation_identifiers.py` |
| no automatic trimming of whitespace | `test_validation_identifiers.py` |

### Reference integrity

| Test | File |
|---|---|
| sample references unknown group_id | `test_validation_references.py` |
| sample factor_values key not in factors | `test_validation_references.py` |
| factor_values value not in factor.levels | `test_validation_references.py` |
| sample covariate_values key not in covariates | `test_validation_references.py` |
| contrast references unknown comparison_group_id | `test_validation_references.py` |
| contrast references unknown reference_group_id | `test_validation_references.py` |
| pairing references condition_factor not a valid factor_id | `test_validation_references.py` |
| factor_id and covariate_id collision (same string in both) | `test_validation_references.py` |

### Contrast tests

| Test | File |
|---|---|
| valid comparison/reference contrast | `test_validation_contrasts.py` |
| same group on both sides (comparison == reference) rejected | `test_validation_contrasts.py` |
| comparison_group_id points to unknown group | `test_validation_contrasts.py` |
| reference_group_id points to unknown group | `test_validation_contrasts.py` |
| referenced group has no samples (empty group) | `test_validation_contrasts.py` |
| duplicate contrast_id with same direction rejected | `test_validation_contrasts.py` |
| A vs B and B vs A both allowed as separate contrasts | `test_validation_contrasts.py` |
| direction semantics remain stable across validation runs | `test_validation_contrasts.py` |

### Treatment and Quantity tests

| Test | File |
|---|---|
| treatment with agent only (no dose, no duration) | `test_validation_treatments.py` |
| treatment with dose and duration | `test_validation_treatments.py` |
| vehicle control (agent="vehicle" or empty) | `test_validation_treatments.py` |
| empty agent string rejected | `test_validation_treatments.py` |
| negative quantity rejected | `test_validation_treatments.py` |
| zero quantity rejected | `test_validation_treatments.py` |
| empty unit string rejected | `test_validation_treatments.py` |
| NaN rejected | `test_validation_treatments.py` |
| infinity rejected | `test_validation_treatments.py` |
| no unit conversion performed (units stored as-is) | `test_validation_treatments.py` |

### Covariate tests

| Test | File |
|---|---|
| categorical covariate accepts string values | `test_validation_covariates.py` |
| continuous covariate accepts numeric values | `test_validation_covariates.py` |
| integer accepted for continuous covariate | `test_validation_covariates.py` |
| boolean accepted for categorical covariate | `test_validation_covariates.py` |
| string accepted for categorical covariate | `test_validation_covariates.py` |
| allowed_values constraint enforced | `test_validation_covariates.py` |
| bool rejected as integer-type covariate value | `test_validation_covariates.py` |
| NaN rejected for continuous covariate | `test_validation_covariates.py` |
| batch is represented as `CovariateDefinition(role="batch")` | `test_validation_covariates.py` |

### Replicate tests

| Test | File |
|---|---|
| biological_replicate ID string accepted | `test_validation_replicates.py` |
| technical_replicate ID string accepted | `test_validation_replicates.py` |
| empty replicate ID string allowed (optional field) | `test_validation_replicates.py` |
| same biological_replicate ID appears across different groups | `test_validation_replicates.py` |
| paired samples can reuse biological_replicate ID | `test_validation_replicates.py` |
| technical replicate not mapped to multiple biological replicates | `test_validation_replicates.py` |
| no minimum replicate count enforced at schema level | `test_validation_replicates.py` |

### Pairing tests

| Test | File |
|---|---|
| valid before/after pair structure | `test_validation_pairing.py` |
| sample with missing pair_id is valid (unpaired sample) | `test_validation_pairing.py` |
| pair with only one sample flagged | `test_validation_pairing.py` |
| multiple independent pairs in one design | `test_validation_pairing.py` |
| pairing references condition_factor that is not a valid factor_id | `test_validation_pairing.py` |
| duplicated condition value within a pair flagged | `test_validation_pairing.py` |
| incomplete pair (one sample missing) flagged | `test_validation_pairing.py` |
| require_complete_pairs flag behavior | `test_validation_pairing.py` |
| no automatic pair inference from sample names | `test_validation_pairing.py` |

### Compatibility tests

| Test | File |
|---|---|
| exact sample match between design and matrix | `test_omics_compatibility.py` |
| same samples in different order (order-independent matching) | `test_omics_compatibility.py` |
| design is a subset of matrix samples | `test_omics_compatibility.py` |
| matrix is a subset of design samples | `test_omics_compatibility.py` |
| partial overlap (some samples missing) → violations reported | `test_omics_compatibility.py` |
| no overlap at all → violations reported | `test_omics_compatibility.py` |
| duplicate IDs in either object flagged | `test_omics_compatibility.py` |
| one design applied to multiple modalities (transcriptomics + proteomics) | `test_omics_compatibility.py` |
| compatibility validator does not mutate the OmicsMatrix | `test_omics_compatibility.py` |
| compatibility validator does not mutate the ExperimentDesign | `test_omics_compatibility.py` |

### Purity and determinism tests

| Test | File |
|---|---|
| `validate()` does not mutate the input ExperimentDesign | `test_validation_*.py` (each file) |
| `validate()` returns identical error list on repeated calls | `test_validation_*.py` (each file) |
| error ordering is stable across runs | `test_validation_*.py` (each file) |
| validation performs no file I/O | `test_validation_*.py` (each file) |

---

## Explicitly Deferred Topics

The following topics are **not rejected** — they are deferred to a future ADR
or phase decision. Phase 2B.0 freezes the minimum viable contract only.
Future extensions **must not break the semantics of any field defined in this
document**.

### Statistical analysis (deferred)

| Topic | Reason for deferral |
|---|---|
| differential analysis | Requires external analysis backend |
| statistical formulas | Belongs in analysis layer, not schema |
| design matrix generation | Analysis-layer concern |
| contrast matrix generation | Analysis-layer concern |
| coefficient selection | Depends on analysis backend |
| interaction terms | Statistical model concern |
| random effects | Mixed-model concern |
| mixed models | Analysis-layer concern |
| longitudinal statistical models | Out of scope for v1 |
| repeated-measures analysis | Analysis-layer concern |

### Quality and power (deferred)

| Topic | Reason for deferral |
|---|---|
| power analysis | Requires external power calculation tool |
| minimum replicate recommendation | Statistical guideline, not schema rule |
| replicate adequacy scoring | Analysis-layer concern |
| effect-size thresholds | Analysis configuration |
| p-value thresholds | Analysis configuration |
| multiple-testing configuration | Analysis configuration |

### Workflow state (deferred)

| Topic | Reason for deferral |
|---|---|
| workflow status | Orchestration-layer concern |
| approval state | Orchestration-layer concern |
| validation results stored in schema | Violates ADR 0006 purity constraint |

### Data manipulation (deferred)

| Topic | Reason for deferral |
|---|---|
| omics matrix mutation | Phase 2A immutability contract applies |
| sample filtering | Analysis-layer concern |
| sample reordering | Analysis-layer concern |
| selector/query language | Explicitly rejected for this phase |

### Automatic inference (deferred / explicitly rejected for this phase)

| Topic | Reason for deferral |
|---|---|
| ConditionSelector | Rejected — violates declarative principle (ADR 0006) |
| matrix_index | Rejected — dynamic resolution, not declarative |
| automatic group inference from sample names | Rejected — violates explicit-data principle |
| automatic control inference | Rejected — violates explicit-data principle |
| automatic pair inference from sample names | Rejected — violates explicit-data principle |
| automatic contrast generation from group labels | Rejected — analysis-layer concern |
| sample-name parsing | Requires domain-specific heuristics |

### Domain knowledge (deferred)

| Topic | Reason for deferral |
|---|---|
| unit conversion | Requires unit library / equivalence database |
| unit equivalence | Requires domain knowledge base |
| drug knowledge-base integration | External service dependency |
| target and MOA annotation | External knowledge-base dependency |
| pathway analysis | Separate analysis module |
| ontology mapping | External ontology service dependency |
| combination therapy model | Complex schema extension |
| pharmacokinetic fields | Out of scope for transcriptomics v1 |

### Infrastructure (deferred)

| Topic | Reason for deferral |
|---|---|
| schema migration tooling | Needed only after v1 ships |
| backward-compatibility migration | Needed only after schema evolution |
| JSON serialization format standardization beyond basic `model_dump()` | Deferred until cross-system interoperability required |
| enum extension policy beyond the frozen initial set | Deferred until new enum values are needed |
| structured error/warning/info issue model | Deferred until downstream consumers require structured issues |

### Critical guarantees for future extensions

When any deferred topic is addressed in a future phase, the following
constraints apply:

1. **No field semantics change.** The meaning of any field frozen in this
   document must remain stable. New fields may be added, but existing field
   meanings cannot be redefined.
2. **Backward compatibility.** Existing valid ExperimentDesign objects must
   remain valid under any schema extension. Adding optional fields with
   defaults is permitted; removing or changing required fields is not.
3. **Purity preserved.** Validation and compatibility checking remain pure
   functions — they do not mutate inputs or carry side effects.
4. **Separation preserved.** ExperimentDesign and OmicsMatrix remain
   independent domain objects (ADR 0003).

---

## Phase 2B.0 Acceptance Criteria

| # | Criterion | Status |
|---|---|---|
| 1 | Schema types and fields frozen | Done |
| 2 | ADR 0003: Design-Matrix separation documented | Done |
| 3 | ADR 0004: Contrast semantics documented | Done |
| 4 | ADR 0005: Batch-as-covariate documented | Done |
| 5 | ADR 0006: Pure declarative design documented | Done |
| 6 | Validation contract defined | Done |
| 7 | Module placement frozen per module | Done |
| 8 | Future test file layout recorded | Done |
| 9 | Fixture strategy recorded | Done |
| 10 | Deferred topics fully enumerated | Done |
| 11 | Schema migration, serialization, enum extension explicitly deferred | Done |
| 12 | Zero Python, test, or config file modifications | Done |
