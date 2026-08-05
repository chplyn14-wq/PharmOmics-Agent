# ADR 0003: ExperimentDesign and OmicsMatrix Are Independent Domain Objects

**Status:** Accepted (Phase 2B.0 — Architecture Freeze)
**Date:** 2026-08-05

## Context

ExperimentDesign describes the intent and structure of an experiment
(samples, groups, factors, contrasts). OmicsMatrix holds the numeric
data and its provenance. During Phase 2A we built OmicsMatrix as a
frozen Pydantic domain object with its own validation pipeline.

A tempting coupling is to embed OmicsMatrix instances inside
ExperimentDesign, or vice versa. That creates a single aggregate that
conflates experimental intent with measured data.

## Problem

If ExperimentDesign owns OmicsMatrix:

- One experiment design can only ever reference one matrix modality, but
  real experiments may produce transcriptomics, proteomics, and
  metabolomics for the same sample set.
- The design schema grows to carry large numeric dataframes it does not
  own or produce.
- Validation logic must reason about cross-cutting concerns inside one
  model, blurring the boundary between the design layer and the data
  layer.
- A design cannot be authored before data is available, breaking the
  human-in-the-loop workflow where a researcher specifies intent first.

## Decision

**ExperimentDesign and OmicsMatrix are independent domain objects.**
They are aligned only through stable `sample_id` strings.

### Formal rules

1. **No ownership.** ExperimentDesign does not hold, reference, or
   modify OmicsMatrix instances. OmicsMatrix does not reference
   ExperimentDesign.

2. **Alignment by sample_id.** The only link between the two objects
   is that sample_ids appearing in ExperimentDesign.samples may overlap
   with sample_ids in OmicsMatrix.sample_ids. Alignment is checked, not
   assumed.

3. **One design, many matrices.** A single ExperimentDesign may be
   compatible with multiple OmicsMatrix objects of different modalities
   (e.g. transcriptomics and proteomics for the same experiment).

4. **Compatibility validator.** A future standalone function will
   check alignment:

   ```python
   def validate_design_against_omics(
       design: ExperimentDesign,
       omics: OmicsMatrix,
   ) -> list[str]:
       """Return violation strings where design and matrix are
       incompatible (e.g. missing sample_ids). Does not modify
       either input."""
       ...
   ```

   The validator is pure: it reads both objects and returns a list
   of violation descriptions. It does not mutate its inputs.

### Anti-patterns (explicitly forbidden)

- `ExperimentDesign.omics_matrices: list[OmicsMatrix]`
- `OmicsMatrix.experiment_design: ExperimentDesign`
- Embedding numeric data inside the design schema
- Embedding design metadata inside the matrix schema
- Auto-generating an ExperimentDesign from an OmicsMatrix
- Auto-generating an OmicsMatrix from an ExperimentDesign

## Consequences

### Positive

- ExperimentDesign can be authored before any data exists, supporting
  pre-registration of experimental intent.
- Multiple modalities can share the same design without duplication.
- Clear separation of concerns: the design layer describes intent; the
  matrix layer holds measured data.
- The compatibility validator can be swapped or extended independently
  (e.g. adding tolerance for partial sample overlap) without changing
  either domain object.

### Trade-offs

- Cross-referencing requires an explicit validation step rather than
  being guaranteed by type composition.
- Downstream analysis code must receive both objects separately and
  perform alignment checks itself, or rely on the validator.

### Neutral

- The sample_id string is the single contract between layers. If
  sample_id conventions change, both objects must be updated
  consistently — this is a future concern for the compatibility
  validator, not a structural coupling.
