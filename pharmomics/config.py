"""Configuration management for PharmOmics.

Provides a Pydantic Settings model that reads defaults from code and can be
overridden via environment variables prefixed with PHARMOMICS_.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

from pydantic import Field
from pydantic_settings import BaseSettings, SettingsConfigDict

# Project root is the parent of the pharmomics/ package.
_PROJECT_ROOT = Path(__file__).resolve().parent.parent


class Settings(BaseSettings):
    """Application configuration.

    All fields can be overridden by environment variables with the prefix
    ``PHARMOMICS_`` (e.g. ``PHARMOMICS_DATA_DIR=/tmp/data``).
    """

    model_config = SettingsConfigDict(
        env_prefix="PHARMOMICS_",
        env_file=".env",
        env_file_encoding="utf-8",
        extra="ignore",
    )

    data_dir: Path = Field(
        default=Path("data"),
        description="Directory for input data files.",
    )
    run_store_dir: Path = Field(
        default=Path("runs"),
        description="Directory for versioned run outputs.",
    )
    llm_api_key: str | None = Field(
        default=None,
        description="API key for the LLM provider.",
    )
    demo_gse_accession: str = Field(
        default="GSE193258",
        description="GEO accession for the demo dataset.",
    )
    demo_expression_file: str = Field(
        default="GSE193258_RNAseq_estimated_counts.tsv.gz",
        description="Filename of the demo expression matrix.",
    )

    def resolved_data_dir(self) -> Path:
        """Return data_dir resolved against the project root."""
        return _resolve(self.data_dir)

    def resolved_run_store_dir(self) -> Path:
        """Return run_store_dir resolved against the project root."""
        return _resolve(self.run_store_dir)

    def snapshot(self) -> dict[str, Any]:
        """Return a JSON-serialisable snapshot of current settings."""
        return {
            "data_dir": str(self.data_dir),
            "run_store_dir": str(self.run_store_dir),
            "llm_api_key": "***" if self.llm_api_key else None,
            "demo_gse_accession": self.demo_gse_accession,
            "demo_expression_file": self.demo_expression_file,
        }


def _resolve(path: Path) -> Path:
    """Resolve a path, treating relative paths as relative to the project root."""
    if path.is_absolute():
        return path
    return (_PROJECT_ROOT / path).resolve()


def load_settings() -> Settings:
    """Load settings from defaults, .env file, and environment variables."""
    return Settings()
