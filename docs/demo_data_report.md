# GSE193258 Demo Data Report — Milestone 0

## Data provenance

| Field | Value |
|---|---|
| GEO accession | GSE193258 |
| Supplementary file | GSE193258_RNAseq_estimated_counts.tsv.gz |
| Download URL | NCBI GEO (GSE193nnn/GSE193258/matrix/) |
| SHA-256 | 51C92720B7BF4A0D7D29F6B7F33E304BD9F858BA8182572260DDF7F1C4DE6D33 |
| File size | 8,413,632 bytes |
| Associated publication | Criscione et al., NPJ Precision Oncology 6:95 (2022) |

All properties below were verified by direct file inspection.

## File structure

| Property | Value |
|---|---|
| Format | TSV, gzip-compressed |
| Encoding | UTF-8 |
| Total rows (incl. header) | 19,713 |
| Data rows (genes) | 19,712 |
| Total columns | 61 (1 gene + 60 samples) |
| Delimiter | Tab |
| Gene identifier type | HGNC gene symbols (verified: A1BG, A2M, ZYX, ZZEF1, etc.) |
| Duplicate gene rows | 0 |

## Expression values

| Property | Value |
|---|---|
| Value type | Non-integer estimated counts |
| Non-integer fraction | 83.9% (992,546 / 1,182,720) |
| Exactly zero | 16.1% (190,174 / 1,182,720) |
| Integer-valued (non-zero) | ~0% (all integer-valued entries are zero) |
| Minimum | 0.0 |
| Maximum | 1,142,959.076 |
| Mean | 1,596.67 |
| Median | 277.83 |
| Q1 | 5.62 |
| Q3 | 1,338.33 |

Sample of first 5 genes (first 4 sample values each):

```
A1BG:    0.568, 5.155, 2.477, 9.383
A1CF:    1.690, 5.044, 0, 0.853
A2M:     4.592, 0, 0, 0
A2ML1:  10.993, 23.050, 26.685, 22.353
A3GALT2: 1.024, 0, 0, 0
```

Values are clearly non-integer (e.g., 0.568264844844155), confirming these
are estimated counts from a probabilistic quantification method, not raw
integer counts from a feature-counting pipeline.

The original quantification tool, normalization state, and intended
downstream statistical workflow remain unverified. Until these are
confirmed, no differential-expression method (limma-voom, DESeq2,
edgeR, or other) is approved for this file.

## Sample metadata

### Full sample list (60 samples)

Naming convention: `{cell_line}_{condition}_{replicate}`

**4 cell lines:** HCC2935, HCC827, H1975, PC9

**5 conditions per cell line:**

| Condition | Description (inferred) |
|---|---|
| DMSO | Vehicle control |
| osi_DTP | Osimertinib drug-tolerant persister |
| osi_acute | Acute osimertinib treatment |
| short_wash | Short washout after treatment |
| long_wash | Long washout after treatment |

**3 biological replicates per condition** (suffixes _1, _2, _3).

### PC9 samples (primary demo subset)

| Condition | Sample names |
|---|---|
| DMSO | PC9_DMSO_1, PC9_DMSO_2, PC9_DMSO_3 |
| osi_DTP | PC9_osi_DTP_1, PC9_osi_DTP_2, PC9_osi_DTP_3 |
| osi_acute | PC9_osi_acute_1, PC9_osi_acute_2, PC9_osi_acute_3 |
| short_wash | PC9_short_wash_1, PC9_short_wash_2, PC9_short_wash_3 |
| long_wash | PC9_long_wash_1, PC9_long_wash_2, PC9_long_wash_3 |

### Replicate structure

Every cell line × condition combination has exactly 3 biological replicates.
The replicate numbering (_1, _2, _3) is consistent across all conditions
and cell lines, suggesting these are matched biological replicates from the
same experimental batch.

### Batch metadata

- Batch variables (sequencing lane, library prep date, technician,
  etc.) are **not present** in the expression matrix column names.
- Whether batch information exists in the GEO sample metadata or the
  publication supplementary tables remains **unverified**.
- Treatment/batch confounding risk **cannot yet be excluded**.

## Original study context

**Source:** Criscione et al. (2022), *NPJ Precision Oncology* 6:95.
DOI: 10.1038/s41698-022-00337-w

**Study description (from GEO):** Profiled osimertinib drug-tolerant
persister (DTP) cells using RNA-seq, ChIP-seq, and ATAC-seq across four
EGFR-mutant NSCLC cell lines to identify common gene regulatory changes
and therapeutic vulnerabilities.

**Analysis method reported in publication:** The authors used Gene Set
Variation Analysis (GSVA) and the Limma R package for differential
pathway analysis.

**Quantification method:** Exact quantification tool: UNKNOWN.
The supplementary filename alone is insufficient to identify the
quantification software. No methods-section verification has been
completed.

## Suitability for demo

| Criterion | Status |
|---|---|
| File downloadable and intact | PASS (SHA-256 verified) |
| ≥3 replicates per condition | PASS (3 per condition) |
| Both conditions present (DMSO, osi_DTP) | PASS |
| Gene identifiers resolvable | PASS (HGNC symbols) |
| Values suitable for DE analysis | PROVISIONAL — matrix contains
non-integer estimated counts, but quantification tool, normalization
state, and intended downstream statistical workflow remain unverified |
| Batch metadata available in file | UNKNOWN (not present in file) |
| Exact quantification tool known | UNKNOWN — supplementary filename alone is
insufficient to identify the quantification software |

## Blockers

**Blocker M0-B1:** The differential-expression method cannot be finalized
until the quantification method and normalization state are verified from
an authoritative source or a better-characterized input matrix is selected.

## Milestone status

- File structure, integrity, and sample dimensions: **PASS**
- Expression values and gene identifiers: **PASS**
- Sample naming and replicate structure: **PASS**
- DE method selection: **BLOCKED (M0-B1)**
- Batch metadata: **UNKNOWN**
- Quantification tool: **UNKNOWN**

**Milestone 0: PROVISIONAL PASS WITH METHOD BLOCKER**

Milestone 1 (repository structure, ingestion, provenance, and validation)
may proceed. Milestone 2 (DE implementation) must not begin until M0-B1
is resolved.
