"""Agent/LLM-domain Pydantic schemas for PharmOmics.

Provides schemas for LLM call provenance and agent run metadata.
These are separate from run provenance (run_store) and ingestion
(ingestion.schemas).
"""

from __future__ import annotations

from pydantic import BaseModel, ConfigDict


class LLMCallRecord(BaseModel):
    """Log entry for a single LLM API call."""

    model_config = ConfigDict(frozen=True)

    run_id: str
    call_id: str
    model: str
    prompt_hash: str  # SHA-256 of the prompt text
    timestamp: str  # ISO-8601 datetime string
    token_count_input: int
    token_count_output: int
    cost_usd: float


class AgentRunMetadata(BaseModel):
    """Metadata for an agent task within a run."""

    model_config = ConfigDict(frozen=True)

    run_id: str
    agent_name: str
    task: str
    started_at: str  # ISO-8601 datetime string
    completed_at: str | None = None  # ISO-8601 datetime string
    status: str  # e.g. "running", "completed", "failed"
    input_file: str | None = None
    output_file: str | None = None
