"""Tests for PharmOmics Pydantic schemas across all domains."""

from __future__ import annotations

import pytest

from pharmomics.agents.schemas import AgentRunMetadata, LLMCallRecord
from pharmomics.ingestion.schemas import (
    GeneMappingRecord,
    GeneMappingSet,
    IngestedDataset,
)
from pharmomics.run_store import IngestionManifest, RunProvenance, RunStatus

# ---------------------------------------------------------------------------
# Ingestion schemas
# ---------------------------------------------------------------------------


class TestGeneMappingRecord:
    """Verify GeneMappingRecord schema."""

    def test_valid_full(self) -> None:
        g = GeneMappingRecord(
            gene_symbol="EGFR",
            ensembl_id="ENSG00000146648",
            entrez_id="1956",
            organism="Homo sapiens",
            source="HGNC",
        )
        assert g.gene_symbol == "EGFR"
        assert g.ensembl_id == "ENSG00000146648"

    def test_valid_minimal(self) -> None:
        g = GeneMappingRecord(gene_symbol="TP53")
        assert g.gene_symbol == "TP53"
        assert g.ensembl_id is None
        assert g.entrez_id is None
        assert g.organism == "Homo sapiens"
        assert g.source == "HGNC"

    def test_frozen(self) -> None:
        g = GeneMappingRecord(gene_symbol="TP53")
        with pytest.raises(Exception):
            g.gene_symbol = "BRCA1"

    def test_json_roundtrip(self) -> None:
        g = GeneMappingRecord(
            gene_symbol="EGFR",
            ensembl_id="ENSG00000146648",
            entrez_id="1956",
        )
        data = g.model_dump()
        roundtrip = GeneMappingRecord(**data)
        assert roundtrip == g


class TestGeneMappingSet:
    """Verify GeneMappingSet schema."""

    def test_valid(self) -> None:
        gs = GeneMappingSet(
            mappings=[
                GeneMappingRecord(gene_symbol="EGFR", ensembl_id="ENSG00000146648"),
                GeneMappingRecord(gene_symbol="TP53"),
            ]
        )
        assert len(gs.mappings) == 2

    def test_empty(self) -> None:
        gs = GeneMappingSet()
        assert gs.mappings == []

    def test_mutable_default_isolation(self) -> None:
        """Two instances must not share the same default list."""
        gs1 = GeneMappingSet()
        gs2 = GeneMappingSet()
        assert gs1.mappings == []
        assert gs2.mappings == []
        assert gs1.mappings is not gs2.mappings

    def test_json_roundtrip(self) -> None:
        gs = GeneMappingSet(
            mappings=[
                GeneMappingRecord(gene_symbol="MYC", entrez_id="4609"),
            ]
        )
        data = gs.model_dump()
        roundtrip = GeneMappingSet(**data)
        assert roundtrip == gs


class TestIngestedDataset:
    """Verify IngestedDataset schema."""

    def test_valid(self) -> None:
        ds = IngestedDataset(
            gse_accession="GSE193258",
            expression_file="GSE193258_RNAseq_estimated_counts.tsv.gz",
            n_genes=19712,
            n_samples=60,
            gene_id_type="HGNC gene symbols",
            value_type="estimated counts",
            ingested_at="2025-01-01T00:00:00Z",
            sha256="abc123",
        )
        assert ds.gse_accession == "GSE193258"
        assert ds.n_genes == 19712
        assert ds.n_samples == 60

    def test_frozen(self) -> None:
        ds = IngestedDataset(
            gse_accession="GSE00000",
            expression_file="x.tsv",
            n_genes=1,
            n_samples=1,
            gene_id_type="HGNC",
            value_type="counts",
            ingested_at="2025-01-01T00:00:00Z",
            sha256="abc",
        )
        with pytest.raises(Exception):
            ds.n_genes = 999

    def test_json_roundtrip(self) -> None:
        ds = IngestedDataset(
            gse_accession="GSE_SYNTHETIC",
            expression_file="synthetic.tsv",
            n_genes=5,
            n_samples=6,
            gene_id_type="HGNC",
            value_type="synthetic",
            ingested_at="2025-01-01T00:00:00Z",
            sha256="def456",
        )
        data = ds.model_dump()
        roundtrip = IngestedDataset(**data)
        assert roundtrip == ds


# ---------------------------------------------------------------------------
# Agent schemas
# ---------------------------------------------------------------------------


class TestLLMCallRecord:
    """Verify LLMCallRecord schema."""

    def test_valid(self) -> None:
        log = LLMCallRecord(
            run_id="run-20250101-120000-abcdef01",
            call_id="call-001",
            model="gpt-4",
            prompt_hash="def456",
            timestamp="2025-01-01T12:05:00Z",
            token_count_input=500,
            token_count_output=200,
            cost_usd=0.015,
        )
        assert log.model == "gpt-4"
        assert log.token_count_input == 500
        assert log.cost_usd == 0.015

    def test_frozen(self) -> None:
        log = LLMCallRecord(
            run_id="run-20250101-120000-abcdef01",
            call_id="call-001",
            model="gpt-4",
            prompt_hash="def456",
            timestamp="2025-01-01T12:05:00Z",
            token_count_input=500,
            token_count_output=200,
            cost_usd=0.015,
        )
        with pytest.raises(Exception):
            log.cost_usd = 0.99

    def test_json_roundtrip(self) -> None:
        log = LLMCallRecord(
            run_id="run-20250101-120000-abcdef01",
            call_id="call-001",
            model="gpt-4",
            prompt_hash="abc",
            timestamp="2025-01-01T12:05:00Z",
            token_count_input=100,
            token_count_output=50,
            cost_usd=0.001,
        )
        data = log.model_dump()
        roundtrip = LLMCallRecord(**data)
        assert roundtrip == log


class TestAgentRunMetadata:
    """Verify AgentRunMetadata schema."""

    def test_valid(self) -> None:
        m = AgentRunMetadata(
            run_id="run-20250101-120000-abcdef01",
            agent_name="qc_pipeline",
            task="validate expression matrix",
            started_at="2025-01-01T12:00:00Z",
            completed_at="2025-01-01T12:01:00Z",
            status="completed",
            input_file="/data/matrix.tsv",
            output_file="/runs/run-xxx/qc_report.json",
        )
        assert m.agent_name == "qc_pipeline"
        assert m.status == "completed"

    def test_valid_running(self) -> None:
        m = AgentRunMetadata(
            run_id="run-20250101-120000-abcdef01",
            agent_name="hypothesis_gen",
            task="generate hypotheses",
            started_at="2025-01-01T12:00:00Z",
            status="running",
        )
        assert m.completed_at is None

    def test_frozen(self) -> None:
        m = AgentRunMetadata(
            run_id="run-20250101-120000-abcdef01",
            agent_name="qc",
            task="validate",
            started_at="2025-01-01T12:00:00Z",
            status="running",
        )
        with pytest.raises(Exception):
            m.status = "completed"

    def test_json_roundtrip(self) -> None:
        m = AgentRunMetadata(
            run_id="run-20250101-120000-abcdef01",
            agent_name="qc",
            task="validate",
            started_at="2025-01-01T12:00:00Z",
            status="completed",
        )
        data = m.model_dump()
        roundtrip = AgentRunMetadata(**data)
        assert roundtrip == m


# ---------------------------------------------------------------------------
# Run-store schemas
# ---------------------------------------------------------------------------


class TestIngestionManifest:
    """Verify IngestionManifest schema."""

    def test_valid(self) -> None:
        m = IngestionManifest(
            file_path="/data/GSE193258_counts.tsv.gz",
            sha256="abc123",
            file_size=8413632,
            uploaded_at="2025-01-01T00:00:00Z",
            gse_accession="GSE193258",
        )
        assert m.file_path == "/data/GSE193258_counts.tsv.gz"
        assert m.sha256 == "abc123"
        assert m.file_size == 8413632
        assert m.gse_accession == "GSE193258"

    def test_frozen(self) -> None:
        m = IngestionManifest(
            file_path="/data/x.tsv",
            sha256="abc",
            file_size=100,
            uploaded_at="2025-01-01T00:00:00Z",
            gse_accession="GSE00000",
        )
        with pytest.raises(Exception):
            m.file_path = "/data/y.tsv"

    def test_json_roundtrip(self) -> None:
        m = IngestionManifest(
            file_path="/data/x.tsv",
            sha256="abc",
            file_size=100,
            uploaded_at="2025-01-01T00:00:00Z",
            gse_accession="GSE00000",
        )
        data = m.model_dump()
        roundtrip = IngestionManifest(**data)
        assert roundtrip == m


class TestRunProvenance:
    """Verify RunProvenance schema."""

    def test_valid_minimal(self) -> None:
        p = RunProvenance(
            run_id="run-20250101-120000-abcdef01",
            created_at="2025-01-01T12:00:00Z",
        )
        assert p.run_id == "run-20250101-120000-abcdef01"
        assert p.status == RunStatus.CREATED
        assert p.input_manifest is None
        assert p.settings_snapshot == {}

    def test_valid_full(self) -> None:
        manifest = IngestionManifest(
            file_path="/data/x.tsv",
            sha256="abc",
            file_size=100,
            uploaded_at="2025-01-01T00:00:00Z",
            gse_accession="GSE00000",
        )
        p = RunProvenance(
            run_id="run-20250101-120000-abcdef01",
            created_at="2025-01-01T12:00:00Z",
            settings_snapshot={"data_dir": "/data"},
            input_manifest=manifest,
            status=RunStatus.RUNNING,
        )
        assert p.input_manifest is not None
        assert p.input_manifest.gse_accession == "GSE00000"
        assert p.status == RunStatus.RUNNING

    def test_frozen(self) -> None:
        p = RunProvenance(
            run_id="run-20250101-120000-abcdef01",
            created_at="2025-01-01T12:00:00Z",
        )
        with pytest.raises(Exception):
            p.status = RunStatus.FAILED

    def test_json_roundtrip(self) -> None:
        p = RunProvenance(
            run_id="run-20250101-120000-abcdef01",
            created_at="2025-01-01T12:00:00Z",
            settings_snapshot={"key": "value"},
        )
        data = p.model_dump()
        roundtrip = RunProvenance(**data)
        assert roundtrip == p


class TestRunStatus:
    """Verify RunStatus enum."""

    def test_values(self) -> None:
        assert RunStatus.CREATED.value == "created"
        assert RunStatus.RUNNING.value == "running"
        assert RunStatus.COMPLETED.value == "completed"
        assert RunStatus.FAILED.value == "failed"

    def test_string_comparison(self) -> None:
        assert RunStatus.CREATED == "created"


# ---------------------------------------------------------------------------
# Module-boundary enforcement
# ---------------------------------------------------------------------------


class TestModuleBoundaries:
    """Verify schemas live in the correct modules."""

    def test_run_store_has_no_ingestion_models(self) -> None:
        import pharmomics.run_store as rs

        assert not hasattr(rs, "GeneMappingRecord")
        assert not hasattr(rs, "GeneMappingSet")
        assert not hasattr(rs, "IngestedDataset")

    def test_run_store_has_no_agent_models(self) -> None:
        import pharmomics.run_store as rs

        assert not hasattr(rs, "LLMCallRecord")
        assert not hasattr(rs, "AgentRunMetadata")

    def test_ingestion_has_gene_mapping(self) -> None:
        from pharmomics.ingestion import schemas as ing

        assert hasattr(ing, "GeneMappingRecord")
        assert hasattr(ing, "GeneMappingSet")

    def test_agents_has_llm_call_record(self) -> None:
        from pharmomics.agents import schemas as ag

        assert hasattr(ag, "LLMCallRecord")
        assert hasattr(ag, "AgentRunMetadata")
