i # ADR 0007: Analysis Specification Contract

## Status

Proposed

## Context

ExperimentDesign defines experimental structure.

OmicsMatrix defines measured data.

Compatibility layer verifies that ExperimentDesign and OmicsMatrix can interact.

However, there is no domain object representing the analytical question.

## Decision

Introduce AnalysisSpecification as an independent declarative domain object.

ExperimentDesign describes what happened in the experiment.

AnalysisSpecification describes what analysis question should be asked.

OmicsMatrix remains the source of measured data.

## Ownership

Contrast belongs to AnalysisSpecification.

A single experiment may support multiple analytical questions, therefore contrasts are analysis-specific.

## Non-goals

This ADR does not define:

- statistical models
- design matrix generation
- execution engine
- automatic inference
- batch correction strategy

## Sample Identity

AnalysisSpecification MUST NOT contain sample identifiers.

OmicsMatrix.sample_ids remains the single source of truth for sample identity.
