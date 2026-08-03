# PharmOmics Agent

> **AI-native Scientific Research Copilot for Transcriptomics-driven Drug Discovery**

PharmOmics Agent is an open-source platform that helps researchers transform transcriptomics data into evidence-backed scientific hypotheses and experiment recommendations.

Unlike traditional RNA-seq pipelines that primarily generate statistical outputs, PharmOmics Agent is designed to bridge computational biology, biomedical knowledge retrieval, and AI-assisted scientific reasoning into a unified research workflow.

> **AI assists scientific reasoning. Scientists make the decisions.**

---

## Why PharmOmics Agent?

High-throughput transcriptomics experiments routinely identify thousands of differentially expressed genes.

The real challenge begins **after** the analysis:

- Which genes deserve further investigation?
- Which pathways are most biologically relevant?
- What does the current literature say?
- Are there existing drugs targeting these genes?
- Have related clinical trials already been conducted?
- Which experiments should be prioritized next?

Answering these questions typically requires researchers to manually combine information from multiple independent resources.

PharmOmics Agent aims to streamline this process by integrating transcriptomics, biomedical knowledge, retrieval-augmented generation (RAG), and AI-assisted reasoning into a single evidence-first workflow.

---

## Long-term Workflow

```text
Transcriptomics
        │
        ▼
Knowledge Retrieval
        │
        ▼
Evidence Integration
        │
        ▼
Retrieval-Augmented Generation (RAG)
        │
        ▼
AI Scientific Reasoning
        │
        ▼
Hypothesis Generation
        │
        ▼
Experiment Recommendation
        │
        ▼
Human Review
        │
        ▼
Scientific Decision Support
```

---

## Project Philosophy

PharmOmics Agent follows six guiding principles:

- 🧠 **AI-native**
- 📚 **Evidence-first**
- 👩‍🔬 **Human-in-the-loop**
- 🔁 **Reproducible**
- 🧩 **Modular**
- 🌍 **Open Source**

AI augments scientific reasoning—it does not replace scientific expertise.

---

## Development Status

**Current milestone:** Foundation (Milestone 1)

### ✅ Implemented

- Project architecture
- Python package structure
- CLI framework
- Metadata schema
- Dataset schema
- Pydantic validation
- Testing framework
- Ruff linting
- GitHub Actions CI
- Documentation skeleton

### 🚧 Planned

The following capabilities are planned but **not yet implemented**:

- GEO downloader
- Differential expression analysis (DEG)
- Gene Set Enrichment Analysis (GSEA)
- Pathway enrichment
- PubMed retrieval
- DrugBank integration
- ClinicalTrials retrieval
- Retrieval-Augmented Generation (RAG)
- LLM-based scientific reasoning
- Multi-agent workflow

---

## Roadmap

| Milestone | Status |
|------------|--------|
| Milestone 1 — Foundation | ✅ Completed |
| Milestone 2 — Transcriptomics Analysis | 🚧 Planned |
| Milestone 3 — Biomedical Knowledge Retrieval | 🚧 Planned |
| Milestone 4 — Evidence Integration | 🚧 Planned |
| Milestone 5 — Retrieval-Augmented Generation | 🚧 Planned |
| Milestone 6 — AI Scientific Reasoning | 🚧 Planned |
| Milestone 7 — Multi-Agent Scientific Workflow | 🚧 Planned |

See **ROADMAP.md** for detailed planning.

---

## Repository Structure

```text
pharmomics-agent/
│
├── pharmomics/
├── tests/
├── docs/
│   ├── architecture.md
│   ├── design_principles.md
│   ├── scientific_scope.md
│   ├── project_philosophy.md
│   └── development_notes.md
│
├── README.md
├── ROADMAP.md
├── CONTRIBUTING.md
├── CHANGELOG.md
└── pyproject.toml
```

---

## Documentation

| Document | Description |
|----------|-------------|
| README.md | Project overview |
| ROADMAP.md | Development roadmap |
| CONTRIBUTING.md | Contribution guide |
| architecture.md | System architecture |
| design_principles.md | Design philosophy |
| scientific_scope.md | Scientific scope |
| project_philosophy.md | Core philosophy |

---

## Local Data Ingestion

Ingest and validate a local expression matrix and sample metadata:

```bash
uv run pharmomics ingest \
  --expression-file data/matrix.tsv.gz \
  --metadata-file data/metadata.json \
  --source-id GSE193258 \
  --run-dir runs/my-run
```

### Expected Expression File Format

- **Format:** TSV or CSV, optionally gzip-compressed (`.gz` extension)
- **First column:** Gene identifiers (Ensembl IDs, HGNC symbols, Entrez IDs)
- **Remaining columns:** Numeric expression values, one per sample
- **Header row:** Required — first row contains column names
- **Sample names:** Must be unique, no duplicates allowed

Example TSV:

```
gene	Sample_A_1	Sample_A_2	Sample_B_1	Sample_B_2
EGFR	1000	1050	200	210
TP53	300	310	305	295
```

### Expected Metadata File Format

- **Format:** JSON, TSV, or CSV (gzip-compressed variants supported)
- **Required fields:** `sample_id`, `condition`
- **Optional fields:** `cell_line`, `treatment`, `replicate`, `batch`
- Every expression sample must have exactly one metadata row
- No duplicate sample IDs allowed
- No unexplained extra metadata samples

Example JSON:

```json
{
  "gse_accession": "GSE193258",
  "samples": {
    "PC9_DMSO_1": {"condition": "DMSO", "cell_line": "PC9", "replicate": 1},
    "PC9_DMSO_2": {"condition": "DMSO", "cell_line": "PC9", "replicate": 2},
    "PC9_osi_DTP_1": {"condition": "osi_DTP", "cell_line": "PC9", "replicate": 1}
  }
}
```

### Optional Arguments

| Flag | Description |
|---|---|
| `--value-type` | Override value classification (`raw_integer_counts`, `non_integer_estimated_counts`, `normalized_nonnegative_values`, `transformed_values`, `unknown`) |
| `--gene-id-type` | Override gene ID classification (`ensembl_ids`, `hgnc_symbols`, `entrez_ids`, `mixed`, `unknown`) |
| `--contrast-control` | Control condition name for contrast validation |
| `--contrast-treatment` | Treatment condition name for contrast validation |

### Important Notes

- The ingestion classification step **does not select or approve** a differential-expression method. Method selection requires independent verification of the quantification source and normalization state.
- Ensembl version suffixes (e.g., `.10`) are stripped only in a separate normalized field; original identifiers are always preserved.
- Gene symbols are treated as annotations and are not silently collapsed or merged.
- Missing batch information is recorded as unknown, not inferred.
- All persisted paths are relative to the run directory.

---

## Quick Start

Clone the repository:

```bash
git clone https://github.com/<your-org>/pharmomics-agent.git
cd pharmomics-agent
```

Create a virtual environment:

```bash
uv sync
```

Run the CLI:

```bash
uv run pharmomics --help
```

> **Note**
>
> Transcriptomics analysis and AI reasoning modules are currently under active development.
> The current release focuses on establishing the project's foundational architecture.

---

## Contributing

Community contributions are welcome.

Whether your background is in:

- Bioinformatics
- Computational Biology
- Artificial Intelligence
- Software Engineering
- Biomedical Knowledge Graphs

we encourage you to contribute.

Please read **CONTRIBUTING.md** before submitting a pull request.

---

## License

This project will be released under an open-source license.

---

## Disclaimer

PharmOmics Agent is intended for scientific research support only.

It does **not** provide medical advice, diagnosis, or treatment recommendations.

All AI-generated outputs should be reviewed and validated by qualified researchers before being used in scientific or clinical decision-making.