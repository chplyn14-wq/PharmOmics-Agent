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