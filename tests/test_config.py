"""Tests for pharmomics.config."""

from __future__ import annotations

from pathlib import Path

from pharmomics.config import _PROJECT_ROOT, Settings, load_settings


class TestSettingsDefaults:
    """Verify default field values."""

    def test_data_dir_default(self) -> None:
        s = Settings()
        assert s.data_dir == Path("data")

    def test_run_store_dir_default(self) -> None:
        s = Settings()
        assert s.run_store_dir == Path("runs")

    def test_llm_api_key_default(self) -> None:
        s = Settings()
        assert s.llm_api_key is None

    def test_demo_gse_accession_default(self) -> None:
        s = Settings()
        assert s.demo_gse_accession == "GSE193258"

    def test_demo_expression_file_default(self) -> None:
        s = Settings()
        assert s.demo_expression_file == "GSE193258_RNAseq_estimated_counts.tsv.gz"


class TestSettingsResolvedPaths:
    """Verify path resolution against project root."""

    def test_resolved_data_dir(self) -> None:
        s = Settings()
        resolved = s.resolved_data_dir()
        assert resolved.is_absolute()
        assert resolved == (_PROJECT_ROOT / "data").resolve()

    def test_resolved_run_store_dir(self) -> None:
        s = Settings()
        resolved = s.resolved_run_store_dir()
        assert resolved.is_absolute()
        assert resolved == (_PROJECT_ROOT / "runs").resolve()


class TestSettingsSnapshot:
    """Verify the snapshot method."""

    def test_snapshot_keys(self) -> None:
        s = Settings()
        snap = s.snapshot()
        assert set(snap.keys()) == {
            "data_dir",
            "run_store_dir",
            "llm_api_key",
            "demo_gse_accession",
            "demo_expression_file",
        }

    def test_snapshot_hides_api_key(self) -> None:
        s = Settings(llm_api_key="secret-key-123")
        snap = s.snapshot()
        assert snap["llm_api_key"] == "***"

    def test_snapshot_null_api_key(self) -> None:
        s = Settings()
        snap = s.snapshot()
        assert snap["llm_api_key"] is None


class TestSettingsEnvOverride:
    """Verify environment-variable override behavior."""

    def test_env_override_data_dir(self, monkeypatch: object) -> None:
        monkeypatch.setenv("PHARMOMICS_DATA_DIR", "/tmp/custom_data")
        s = Settings()
        assert s.data_dir == Path("/tmp/custom_data")

    def test_env_override_gse_accession(self, monkeypatch: object) -> None:
        monkeypatch.setenv("PHARMOMICS_DEMO_GSE_ACCESSION", "GSE99999")
        s = Settings()
        assert s.demo_gse_accession == "GSE99999"

    def test_env_override_llm_api_key(self, monkeypatch: object) -> None:
        monkeypatch.setenv("PHARMOMICS_LLM_API_KEY", "env-key-456")
        s = Settings()
        assert s.llm_api_key == "env-key-456"


class TestLoadSettings:
    """Verify the load_settings convenience function."""

    def test_returns_settings(self) -> None:
        s = load_settings()
        assert isinstance(s, Settings)
