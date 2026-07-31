# Architecture v0.1 — Milestone 0

## Overview

PharmOmics Agent v0.1 is a human-in-the-loop research assistant for drug-response
transcriptomics. It performs deterministic QC and statistical analysis, retrieves
PubMed evidence, generates evidence-grounded hypotheses, and proposes validation
experiments.

## Demo dataset

- **GEO accession:** GSE193258
- **Publication:** Criscione et al., *NPJ Precision Oncology* 6:95 (2022).
  DOI: 10.1038/s41698-022-00337-w
- **Supplementary file:** `GSE193258_RNAseq_estimated_counts.tsv.gz`
- **SHA-256:** `51C92720B7BF4A0D7D29F6B7F33E304BD9F858BA8182572260DDF7F1C4DE6D33`
- **File size:** 8,413,632 bytes

## Verified data properties

See `docs/demo_data_report.md` for the full data report. Summary:

| Property | Value |
|---|---|
| Format | Tab-separated, gzip-compressed |
| Rows | 19,712 (gene rows) |
| Columns | 60 sample columns + 1 gene column = 61 |
| Gene identifiers | HGNC gene symbols (e.g., A1BG, ZYX) |
| Value type | Non-integer estimated counts (83.9% non-integer, 16.1% zero) |
| Value range | 0 to 1,142,959.08 |
| Median | 277.83 |
| Unique genes | 19,712 (0 duplicates) |

## Sample design

**4 cell lines** × **5 conditions** × **3 biological replicates** = 60 samples

| Cell line | Conditions |
|---|---|
| HCC2935 | DMSO, osi_DTP, osi_acute, short_wash, long_wash |
| HCC827 | DMSO, osi_DTP, osi_acute, short_wash, long_wash |
| H1975 | DMSO, osi_DTP, osi_acute, short_wash, long_wash |
| PC9 | DMSO, osi_DTP, osi_acute, short_wash, long_wash |

Each condition has exactly **3 biological replicates** (suffixes _1, _2, _3).

## Primary demo comparison

**PC9 DMSO vs PC9 osi_DTP**

- Control: `PC9_DMSO_1`, `PC9_DMSO_2`, `PC9_DMSO_3`
- Treatment: `PC9_osi_DTP_1`, `PC9_osi_DTP_2`, `PC9_osi_DTP_3`

This is the smallest self-contained comparison that demonstrates the
full agent pipeline (QC → DE → hypothesis → experiment).

## DE method selection

### Recommendation: limma-voom

**Rationale:**
- The file contains estimated (non-integer) counts, not raw integer counts.
- DESeq2 and edgeR assume integer raw counts from feature-counting pipelines.
- limma-voom accepts non-integer count-like data and applies a
  mean-variance trend to compute observation-level precision weights.
- The original publication (Criscione et al., 2022) used **Limma** for
  differential expression analysis, providing precedent for this choice
  on the same dataset.

### Alternative considerations

- **DESeq2** with rounding: possible but rounding estimated counts is
  not standard practice and introduces discretization artifacts.
- **DESeq2 on raw counts** would require access to the original FASTQ
  files and a full re-quantification pipeline, which is outside scope
  for Milestone 0.

### Uncertainty

- The exact quantification tool (RSEM, Salmon, Kallisto) used by the
  original authors is **unknown** — the GEO series matrix file and
  supplementary metadata do not specify it. The term "estimated counts"
  is conventionally associated with RSEM output, but this has not been
  confirmed from the original source.
- Batch information (sequencing batch, library prep date) is **not
  present** in the expression file column names and is **unknown**
  whether it was recorded in the original study metadata.

## Component architecture

```
app.py (Streamlit)
  ├── qc/             — file validation, integrity checks
  ├── analysis/       — DE analysis (limma-voom), result tables
  ├── evidence/       — PubMed/NCBI retrieval, citation linking
  ├── hypothesis/     — LLM-based hypothesis generation with constraints
  └── experiment/     — validation experiment proposals (human approval required)
```

## Human-in-the-loop checkpoints

1. LLMs may interpret results but must never modify calculated values.
2. Numerical results come from deterministic analysis code only.
3. Every hypothesis includes limitations and falsification criteria.
4. No outputs described as clinically validated.
5. Experiment recommendations require explicit human approval.
