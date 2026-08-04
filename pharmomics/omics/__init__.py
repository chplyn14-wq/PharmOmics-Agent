"""PharmOmics omics core — domain types and schemas.

Phase 2A provides the multi-omics compatible ``OmicsMatrix`` data object
and its lightweight descriptor counterpart.  The enums defined here are
the controlled vocabularies shared across all omics modalities.

This module is deliberately independent of the ingestion layer — it
defines the *domain* shapes; an adapter (added later) bridges the gap
between ``ExpressionLoadResult`` / ``IngestionResult`` and these
schemas.
"""
