"""PharmOmics run store — utilities and provenance schemas.

Milestone 1 provides:
- Run ID generation
- Run directory creation
- Relative path resolution (with escape protection)
- SHA-256 file hashing
- Atomic JSON I/O
- Pydantic schemas for run provenance and ingestion manifests.

Ingestion-domain schemas (gene mapping, ingested dataset metadata) live in
``pharmomics.ingestion.schemas``.
Agent/LLM provenance schemas live in ``pharmomics.agents.schemas``.
"""

from __future__ import annotations

import hashlib
import json
import secrets
import tempfile
from datetime import UTC, datetime
from enum import StrEnum
from pathlib import Path
from typing import Any

from pydantic import BaseModel, ConfigDict, Field

# ---------------------------------------------------------------------------
# Run ID generation
# ---------------------------------------------------------------------------


def generate_run_id() -> str:
    """Return a unique run ID string.

    Format: ``run-YYYYMMDD-HHMMSS-<8-hex-chars>``.
    """
    now = datetime.now(UTC)
    timestamp = now.strftime("%Y%m%d-%H%M%S")
    suffix = secrets.token_hex(4)
    return f"run-{timestamp}-{suffix}"


# ---------------------------------------------------------------------------
# Run directory management
# ---------------------------------------------------------------------------


class RunDirectoryExistsError(RuntimeError):
    """Raised when a run directory already exists."""


def create_run_directory(store_dir: Path, run_id: str) -> Path:
    """Create a new run directory inside *store_dir*.

    Parameters
    ----------
    store_dir : Path
        The parent directory that holds all runs.
    run_id : str
        The unique identifier for this run (typically from ``generate_run_id``).

    Returns
    -------
    Path
        The absolute path of the newly created run directory.

    Raises
    ------
    RunDirectoryExistsError
        If ``store_dir / run_id`` already exists.
    """
    run_dir = store_dir / run_id
    if run_dir.exists():
        raise RunDirectoryExistsError(f"Run directory already exists: {run_dir}")
    run_dir.mkdir(parents=True, exist_ok=False)
    return run_dir.resolve()


# ---------------------------------------------------------------------------
# Path resolution
# ---------------------------------------------------------------------------


class PathEscapeError(ValueError):
    """Raised when a relative path resolves outside the allowed base directory."""


def resolve_relative_path(base: Path, rel: str) -> Path:
    """Resolve *rel* against *base*, ensuring the result stays inside *base*.

    Parameters
    ----------
    base : Path
        The anchor directory (must be absolute).
    rel : str
        A relative path string.

    Returns
    -------
    Path
        The resolved absolute path.

    Raises
    ------
    PathEscapeError
        If the resolved path escapes the base directory.
    ValueError
        If *rel* is an absolute path or *base* is not absolute.
    """
    if not base.is_absolute():
        raise ValueError(f"base must be an absolute path, got: {base}")
    if rel.startswith("/") or (len(rel) > 1 and rel[1] == ":"):
        raise ValueError(f"rel must be a relative path, got: {rel}")

    resolved = (base / rel).resolve()
    # Use string comparison for containment check (handles Windows drive letters)
    try:
        resolved.relative_to(base)
    except ValueError:
        raise PathEscapeError(
            f"Path '{rel}' resolves outside of base directory '{base}'"
        ) from None
    return resolved


# ---------------------------------------------------------------------------
# File hashing
# ---------------------------------------------------------------------------


def hash_file_sha256(path: Path) -> str:
    """Return the lowercase hexadecimal SHA-256 digest of a file.

    Parameters
    ----------
    path : Path
        Path to the file to hash.

    Returns
    -------
    str
        64-character lowercase hex digest.
    """
    h = hashlib.sha256()
    with open(path, "rb") as f:
        while True:
            chunk = f.read(65536)
            if not chunk:
                break
            h.update(chunk)
    return h.hexdigest()


# ---------------------------------------------------------------------------
# JSON I/O
# ---------------------------------------------------------------------------


def write_json(path: Path, data: Any) -> None:
    """Write *data* as JSON to *path* atomically (temp file + rename).

    Parameters
    ----------
    path : Path
        Destination file path.
    data : Any
        A JSON-serialisable object.
    """
    path = path.resolve()
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp = tempfile.mkstemp(dir=path.parent, suffix=".tmp")
    try:
        with open(fd, "w", encoding="utf-8") as f:
            json.dump(data, f, indent=2, default=str)
            f.write("\n")
        Path(tmp).replace(path)
    except Exception:
        Path(tmp).unlink(missing_ok=True)
        raise


def read_json(path: Path) -> Any:
    """Read and parse a JSON file.

    Parameters
    ----------
    path : Path
        Path to the JSON file.

    Returns
    -------
    Any
        The parsed JSON content.
    """
    with open(path, encoding="utf-8") as f:
        return json.load(f)


# ---------------------------------------------------------------------------
# Enums
# ---------------------------------------------------------------------------


class RunStatus(StrEnum):
    """Possible states of a run."""

    CREATED = "created"
    RUNNING = "running"
    COMPLETED = "completed"
    FAILED = "failed"


# ---------------------------------------------------------------------------
# Pydantic schemas
# ---------------------------------------------------------------------------


class IngestionManifest(BaseModel):
    """Record of an ingested input file.

    Tracks the file-level provenance (hash, size, source accession) that
    is attached to a :class:`RunProvenance`.  This lives in run_store
    because it is part of the run's provenance chain, not the ingestion
    domain (which concerns gene mappings and dataset-level metadata).
    """

    model_config = ConfigDict(frozen=True)

    file_path: str
    sha256: str
    file_size: int
    uploaded_at: str  # ISO-8601 datetime string
    gse_accession: str


class RunProvenance(BaseModel):
    """Provenance record for a single run."""

    model_config = ConfigDict(frozen=True)

    run_id: str
    created_at: str  # ISO-8601 datetime string
    settings_snapshot: dict[str, Any] = Field(default_factory=dict)
    input_manifest: IngestionManifest | None = None
    status: RunStatus = RunStatus.CREATED
