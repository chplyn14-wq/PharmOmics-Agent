# Demo Data — GSE193258

## Source

- **GEO Accession:** [GSE193258](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE193258)
- **Publication:** Criscione et al., *NPJ Precision Oncology* 6:95 (2022)
- **DOI:** 10.1038/s41698-022-00337-w
- **Supplementary file:** `GSE193258_RNAseq_estimated_counts.tsv.gz`

## File

The expression file is **not committed to the repository**. Download it from
the NCBI GEO series matrix directory and place it here, or keep it in a
temporary location and reference it by path.

- **SHA-256:** `51C92720B7BF4A0D7D29F6B7F33E304BD9F858BA8182572260DDF7F1C4DE6D33`
- **Size:** 8,413,632 bytes

## Content

- 19,712 genes (HGNC symbols)
- 60 sample columns (4 cell lines × 5 conditions × 3 replicates)
- Non-integer estimated counts

## Primary demo comparison

**PC9 DMSO vs PC9 osi_DTP**

| Group | Samples |
|---|---|
| PC9 DMSO | PC9_DMSO_1, PC9_DMSO_2, PC9_DMSO_3 |
| PC9 osi_DTP | PC9_osi_DTP_1, PC9_osi_DTP_2, PC9_osi_DTP_3 |

## All samples

```
HCC2935: DMSO_1-3, osi_DTP_1-3, osi_acute_1-3, short_wash_1-3, long_wash_1-3
HCC827:  DMSO_1-3, osi_DTP_1-3, osi_acute_1-3, short_wash_1-3, long_wash_1-3
H1975:   DMSO_1-3, osi_DTP_1-3, osi_acute_1-3, short_wash_1-3, long_wash_1-3
PC9:     DMSO_1-3, osi_DTP_1-3, osi_acute_1-3, short_wash_1-3, long_wash_1-3
```

## Notes

- Values are estimated counts (non-integer), not raw counts.
- Recommended DE method: **limma-voom** (compatible with non-integer
  count-like data).
- Full data report: `docs/demo_data_report.md`
- Architecture: `docs/architecture_v0.1.md`
