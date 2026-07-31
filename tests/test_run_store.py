"""Tests for pharmomics.run_store — run ID generation and directory creation."""

from __future__ import annotations

import re
from pathlib import Path

import pytest

from pharmomics.run_store import (
    RunDirectoryExistsError,
    create_run_directory,
    generate_run_id,
)


class TestGenerateRunId:
    """Verify run ID generation."""

    _PATTERN = re.compile(
        r"^run-\d{8}-\d{6}-[0-9a-f]{8}$"
    )

    def test_format(self) -> None:
        run_id = generate_run_id()
        assert self._PATTERN.match(run_id), f"Unexpected format: {run_id}"

    def test_unique(self) -> None:
        id1 = generate_run_id()
        id2 = generate_run_id()
        assert id1 != id2

    def test_starts_with_run_prefix(self) -> None:
        assert generate_run_id().startswith("run-")


class TestCreateRunDirectory:
    """Verify run directory creation."""

    def test_creates_directory(self, tmp_path: Path) -> None:
        store = tmp_path / "runs"
        store.mkdir()
        run_dir = create_run_directory(store, "run-20250101-120000-abcdef01")
        assert run_dir.exists()
        assert run_dir.is_dir()
        assert run_dir.name == "run-20250101-120000-abcdef01"

    def test_creates_parent_dirs(self, tmp_path: Path) -> None:
        store = tmp_path / "deeply" / "nested" / "runs"
        run_dir = create_run_directory(store, "run-20250101-120000-abcdef01")
        assert run_dir.exists()
        assert run_dir.is_dir()

    def test_raises_on_existing(self, tmp_path: Path) -> None:
        store = tmp_path / "runs"
        store.mkdir()
        run_dir = store / "run-20250101-120000-abcdef01"
        run_dir.mkdir()
        with pytest.raises(RunDirectoryExistsError):
            create_run_directory(store, "run-20250101-120000-abcdef01")

    def test_returns_resolved_path(self, tmp_path: Path) -> None:
        store = tmp_path / "runs"
        store.mkdir()
        run_dir = create_run_directory(store, "run-20250101-120000-abcdef01")
        assert run_dir.is_absolute()
