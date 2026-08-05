i # AnalysisSpecification Schema Contract

## Purpose

AnalysisSpecification represents analytical intent.

It does not represent experimental structure or measured data.

## Relationship

ExperimentDesign
    |
    + AnalysisSpecification
    |
OmicsMatrix

## Fields

AnalysisSpecification MAY contain:

- analysis intent
- factor references
- contrast references
- declarative parameters
## Immutability

AnalysisSpecification is a declarative object.

Validation MUST NOT modify AnalysisSpecification.
## Ownership

Contrast belongs to AnalysisSpecification.

ExperimentDesign does not contain analysis-specific contrasts.

## Constraints

AnalysisSpecification MUST NOT:

- contain sample identifiers
- generate design matrix
- select statistical models
- infer treatment/control/batch/pairing

## Open Questions

- Should analysis intent be enum or free-form?
- Should factor references use string ids or typed references?
- Should parameters be extensible?
