# ADR 0006: ExperimentDesign Is a Pure Declarative Domain Object

**Status:** Accepted (Phase 2B.0 — Architecture Freeze)
**Date:** 2026-08-05

## Context

ExperimentDesign captures the facts of an experiment: what samples
exist, how they are grouped, what factors were manipulated, and what
comparisons are intended. It is easy to conflate this descriptive role
with validation state, analysis configuration, or workflow state.

## Decision

**ExperimentDesign describes only experimental facts.** It carries no
validation state, no workflow state, no analysis parameters, and no
statistical configuration.

### What ExperimentDesign describes

- Experimental identity (experiment_id, description)
- Sample list and group membership
- Experimental factors and their values
- Contrast definitions (comparison intent)
- Covariate definitions
- Pairing structure
- Treatment information
- Replicate identifiers
- Non-core extensions via a typed metadata dict

### What ExperimentDesign does NOT contain

The following are explicitly forbidden in the schema:

| Category | Forbidden fields |
|---|---|
| Validation state | `validation_errors`, `status`, `validated`, `approval_state` |
| Workflow state | `workflow_state`, `phase`, `step` |
| Statistical thresholds | `effect_size_threshold`, `p_value_threshold`, `fdr_threshold` |
| Outcome definitions | `primary_outcome`, `secondary_outcomes` |
| Statistical formulas | `design_matrix`, `contrast_matrix`, `coefficient`, `model_type` |
| Analysis backend config | `analysis_backend`, `deseq2_config`, `limma_config`, `edgeR_config` |
| Quality gates | `minimum_per_group`, `replicate_adequacy`, `power_calculation` |

### Validation contract

Validation is a separate pure function:

```python
def validate(design: ExperimentDesign) -> list[str]:
    """Return violation strings. Empty list = valid.
    Does not modify the input design."""
    ...
```

- Validation is pure in-memory: no file I/O, no network calls.
- Validation does not mutate its input.
- Validation results are not written back into the schema.
- The design object is immutable (`frozen=True` in Pydantic).

## Consequences

### Positive

- ExperimentDesign is a stable, serialisable description of
  experimental intent. It can be saved, compared, versioned, and
  reviewed without carrying ephemeral state.
- Validation can evolve independently (new rules, stricter checks)
  without changing the schema.
- Analysis backends can be swapped without modifying the design.
- The schema remains small, reviewable, and understandable.
- Clear boundary: the schema says "what we want to know"; validation
  says "is this well-formed"; analysis says "how we compute it."

### Trade-offs

- The validation results must be carried alongside the design in
  memory or in a wrapper object, not inside the design itself.
- Workflow state (e.g. "pending human review") must be managed by a
  separate orchestration layer, not embedded in the domain object.
