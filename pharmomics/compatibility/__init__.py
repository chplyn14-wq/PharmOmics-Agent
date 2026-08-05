"""PharmOmics compatibility layer.

Provides cross-domain validation between independent domain objects,
ensuring that an ExperimentDesign is compatible with an OmicsMatrix
without coupling the two schemas.

See ADR 0003 for the rationale: ExperimentDesign and OmicsMatrix are
independent domain objects aligned only through sample_id strings.
"""
