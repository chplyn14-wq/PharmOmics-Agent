# ADR 0004: Contrast Uses comparison/reference Group Semantics

**Status:** Accepted (Phase 2B.0 — Architecture Freeze)
**Date:** 2026-08-05

## Context

A Contrast describes which experimental groups are compared in a
differential expression analysis. Multiple naming conventions exist:

- `treatment_group` / `control_group` — assumes a treatment/control
  paradigm that does not always apply (e.g. time-series, dose-response,
  multi-factor designs).
- `numerator` / `denominator` — borrows from fold-change math but
  conflates representation with direction.
- `case` / `control` — specific to case-control study designs.
- `group_a` / `group_b` — ambiguous without external documentation.

## Decision

**Contrast uses `comparison_group_id` and `reference_group_id` as its
two group identifiers.**

### Formal fields

```python
class Contrast(BaseModel):
    contrast_id: str
    comparison_group_id: str   # the "A" in A vs B
    reference_group_id: str    # the "B" in A vs B
```

### Fixed semantic rule

A **positive effect** (e.g. positive log2 fold-change) means:

> **comparison group expression > reference group expression**

This is the single, invariant interpretation rule. All downstream
analysis must respect it.

### Direction

"A vs B" and "B vs A" are two distinct Contrast objects. They are not
interchangeable; reversing direction flips the sign of every effect.

### Anti-patterns (explicitly forbidden)

- `treatment_group_id` / `control_group_id` in Contrast
- `numerator` / `denominator`
- `ConditionSelector`, `matrix_index`, or `annotation_filter` fields
  that dynamically resolve groups at runtime
- Automatic sample query based on contrast — group membership comes
  from ExperimentDesign.samples, not from the contrast
- Automatic contrast generation from group labels
- Storing a list of sample IDs inside Contrast
- Embedding statistical test parameters (p-value, effect size) in
  Contrast

## Consequences

### Positive

- Neutral terminology works for any experimental paradigm — treatment,
  time-course, dose-response, factorial designs — without semantic
  awkwardness.
- The direction rule is unambiguous and can be stated in plain English
  in any output table or report.
- Contrast remains a simple value object: two group IDs and a label.
  No runtime resolution logic, no hidden queries.

### Trade-offs

- Researchers accustomed to "treatment vs control" must map their
  intent to comparison/reference explicitly. This is a feature, not a
  bug: it forces clarity about which direction is being tested.
- When groups are not naturally ordered (e.g. three drug arms), the
  choice of comparison vs reference is an analytical decision, not a
  property of the experiment itself.
