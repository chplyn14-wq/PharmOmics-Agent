# ADR 0005: Batch Is a Covariate Role

**Status:** Accepted (Phase 2B.0 — Architecture Freeze)
**Date:** 2026-08-05

## Context

Batch effects are a common confounder in omics experiments. Some
systems model batch as a first-class concept separate from other
covariates (e.g. age, sex, library prep date), leading to parallel
structures:

```
BatchDesign(batch_id, batch_values)
CovariateDefinition(covariate_id, covariate_values)
```

This duplication forces downstream code to merge two parallel sources
of sample-level annotations and creates ambiguity about which field a
given annotation belongs to.

## Decision

**Batch is modeled as a CovariateDefinition with `role="batch"`.**
There is no separate BatchDesign type.

### Formal modeling

```python
CovariateDefinition(
    covariate_id="batch",
    role="batch",
    value_type="categorical",
)
```

Samples use the same field for all covariates, including batch:

```python
DesignSample(
    sample_id="...",
    covariate_values={"batch": "B1", "age": 45},
)
```

### Anti-patterns (explicitly forbidden)

- A standalone `BatchDesign` or `BatchDefinition` model
- `batch_values` as a separate field alongside `covariate_values` on
  DesignSample
- Special-case batch handling in the schema layer (batch adjustment is
  an analysis concern, not a schema concern)
- Automatic batch detection from sample naming conventions

## Consequences

### Positive

- Single data structure for all sample-level annotations:
  `covariate_values`. Downstream code iterates one dict, not two.
- Batch is treated like any other covariate at the schema level,
  reflecting the scientific reality that batch is just another variable
  that may confound results.
- Analysis-layer batch correction (e.g. ComBat, limma removeBatchEffect)
  can still special-case `role="batch"` without requiring schema
  duplication.

### Trade-offs

- Analysis code must check `role` to identify batch covariates, rather
  than relying on a dedicated type. This is a minor runtime cost for
  significant schema simplicity.
- The convention that `covariate_id="batch"` signals a batch covariate
  relies on string equality, not type safety. A future enum could
  tighten this without changing the schema shape.
