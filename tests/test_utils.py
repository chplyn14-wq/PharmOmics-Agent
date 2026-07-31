"""Tests for pharmomics.run_store utilities — hashing, path resolution, JSON I/O."""

from __future__ import annotations

from pathlib import Path

import pytest

from pharmomics.run_store import (
    PathEscapeError,
    hash_file_sha256,
    read_json,
    resolve_relative_path,
    write_json,
)


class TestHashFileSha256:
    """Verify SHA-256 file hashing."""

    def test_known_content(self, tmp_path: Path) -> None:
        f = tmp_path / "test.txt"
        f.write_bytes(b"hello\n")
        digest = hash_file_sha256(f)
        expected = (
            "5891b5b522d5df086d0ff0b110fbd9d2"
            "1bb4fc7163af34d08286a2e846f6be03"
        )
        assert digest == expected

    def test_empty_file(self, tmp_path: Path) -> None:
        f = tmp_path / "empty.txt"
        f.write_bytes(b"")
        digest = hash_file_sha256(f)
        expected = (
            "e3b0c44298fc1c149afbf4c8996fb924"
            "27ae41e4649b934ca495991b7852b855"
        )
        assert digest == expected

    def test_lowercase_hex(self, tmp_path: Path) -> None:
        f = tmp_path / "data.bin"
        f.write_bytes(b"\xff\x00")
        digest = hash_file_sha256(f)
        assert digest == digest.lower()

    def test_length(self, tmp_path: Path) -> None:
        f = tmp_path / "x.txt"
        f.write_text("x", encoding="utf-8")
        assert len(hash_file_sha256(f)) == 64


class TestResolveRelativePath:
    """Verify relative path resolution with escape protection."""

    def test_valid_relative(self, tmp_path: Path) -> None:
        base = tmp_path / "data"
        base.mkdir()
        resolved = resolve_relative_path(base, "sub/file.tsv")
        assert resolved == base / "sub" / "file.tsv"

    def test_dot_relative(self, tmp_path: Path) -> None:
        base = tmp_path / "data"
        base.mkdir()
        resolved = resolve_relative_path(base, "file.tsv")
        assert resolved == base / "file.tsv"

    def test_escape_raises(self, tmp_path: Path) -> None:
        base = tmp_path / "data"
        base.mkdir()
        with pytest.raises(PathEscapeError):
            resolve_relative_path(base, "../escape.txt")

    def test_deep_escape_raises(self, tmp_path: Path) -> None:
        base = tmp_path / "data"
        base.mkdir()
        with pytest.raises(PathEscapeError):
            resolve_relative_path(base, "a/../../escape.txt")

    def test_absolute_rel_raises(self, tmp_path: Path) -> None:
        base = tmp_path / "data"
        base.mkdir()
        with pytest.raises(ValueError, match="relative path"):
            resolve_relative_path(base, "/absolute/path.txt")

    def test_base_not_absolute_raises(self, tmp_path: Path) -> None:
        with pytest.raises(ValueError, match="absolute path"):
            resolve_relative_path(Path("relative"), "file.txt")


class TestJsonRoundTrip:
    """Verify JSON write/read round-trip."""

    def test_dict_roundtrip(self, tmp_path: Path) -> None:
        data = {"key": "value", "number": 42, "nested": {"a": [1, 2, 3]}}
        path = tmp_path / "test.json"
        write_json(path, data)
        result = read_json(path)
        assert result == data

    def test_list_roundtrip(self, tmp_path: Path) -> None:
        data = [1, "two", {"three": 3}]
        path = tmp_path / "test.json"
        write_json(path, data)
        result = read_json(path)
        assert result == data

    def test_creates_parent_dirs(self, tmp_path: Path) -> None:
        data = {"x": 1}
        path = tmp_path / "a" / "b" / "c" / "test.json"
        write_json(path, data)
        assert path.exists()
        assert read_json(path) == data
