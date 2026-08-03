"""PharmOmics ingestion — expression matrix and sample metadata loading.

Milestone 1B provides:
- Expression matrix loading from TSV/CSV (optionally gzip-compressed)
- Sample metadata loading and validation
- Conservative expression-value classification
- Gene identifier inspection
- Ingestion manifest generation

This module does NOT select a differential-expression method.
"""

from __future__ import annotations

import csv
import gzip
import hashlib
import io
import re
import shutil
from dataclasses import dataclass, field
from datetime import UTC, datetime
from enum import StrEnum
from pathlib import Path
from typing import Any

import pandas as pd

import pharmomics
from pharmomics.run_store import (
    hash_file_sha256,
    write_json,
)

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

ENSEMBL_PATTERN = r"^ENS[GTP][0-9]+(\.[0-9]+)?$"
ENTREZ_PATTERN = r"^[0-9]+$"
HGNC_PATTERN = r"^[A-Z][A-Z0-9]*(-?[A-Z0-9]+)*$"

# ---------------------------------------------------------------------------
# Enums
# ---------------------------------------------------------------------------


class ValueType(StrEnum):
    """Conservative expression-value classifications."""

    RAW_INTEGER_COUNTS = "raw_integer_counts"
    NON_INTEGER_ESTIMATED_COUNTS = "non_integer_estimated_counts"
    NORMALIZED_NONNEGATIVE_VALUES = "normalized_nonnegative_values"
    TRANSFORMED_VALUES = "transformed_values"
    UNKNOWN = "unknown"


class GeneIdType(StrEnum):
    """Gene identifier classifications."""

    ENSEMBL = "ensembl_ids"
    HGNC_SYMBOLS = "hgnc_symbols"
    ENTREZ_IDS = "entrez_ids"
    MIXED = "mixed"
    UNKNOWN = "unknown"


class CompressionType(StrEnum):
    """File compression types."""

    GZIP = "gzip"
    NONE = "none"


# ---------------------------------------------------------------------------
# Data classes for ingestion results
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class ExpressionLoadResult:
    """Result of loading an expression matrix."""

    dataframe: pd.DataFrame
    n_genes: int
    n_samples: int
    sample_ids: list[str]
    gene_ids: list[str]
    delimiter: str
    compression: CompressionType
    original_filename: str


@dataclass(frozen=True)
class MetadataLoadResult:
    """Result of loading sample metadata."""

    dataframe: pd.DataFrame
    sample_ids: set[str]
    conditions: dict[str, str]  # sample_id -> condition
    cell_lines: dict[str, str | None]  # sample_id -> cell_line
    replicates: dict[str, int | None]  # sample_id -> replicate
    batch_values: dict[str, str | None]  # sample_id -> batch
    original_filename: str


@dataclass(frozen=True)
class GeneIdInspectionResult:
    """Result of gene identifier inspection."""

    id_type: GeneIdType
    original_ids: list[str]
    normalized_ids: list[str]  # Ensembl version stripped if applicable
    duplicate_ids: list[str]
    missing_ids: list[str]
    ensembl_count: int
    hgnc_count: int
    entrez_count: int
    unknown_count: int


@dataclass(frozen=True)
class ValueClassificationResult:
    """Result of expression-value classification."""

    value_type: ValueType
    total_values: int
    integer_count: int
    non_integer_count: int
    zero_count: int
    negative_count: int
    min_value: float
    max_value: float
    warnings: list[str] = field(default_factory=list)


@dataclass(frozen=True)
class IngestionResult:
    """Complete ingestion result ready for manifest serialization."""

    source_id: str
    original_expression_file: str
    original_metadata_file: str
    expression_path: str  # relative to run_dir
    metadata_path: str  # relative to run_dir
    expression_sha256: str
    metadata_sha256: str
    expression_file_size: int
    metadata_file_size: int
    delimiter: str
    compression: CompressionType
    n_genes: int
    n_samples: int
    sample_ids: list[str]
    conditions: dict[str, str]
    replicate_counts: dict[str, int]  # condition -> count
    batch_values: dict[str, str | None]
    value_type: ValueType
    gene_id_type: GeneIdType
    gene_id_inspection: GeneIdInspectionResult
    value_classification: ValueClassificationResult
    warnings: list[str] = field(default_factory=list)
    ingested_at: str = ""
    configuration_hash: str = ""
    software_version: str = ""

    def __post_init__(self) -> None:
        if not self.ingested_at:
            object.__setattr__(
                self, "ingested_at", datetime.now(UTC).strftime("%Y-%m-%dT%H:%M:%SZ")
            )
        if not self.software_version:
            object.__setattr__(
                self, "software_version",
                getattr(pharmomics, "__version__", "0.1.0"),
            )

    def to_dict(self) -> dict[str, Any]:
        """Return a JSON-serialisable dictionary."""
        return {
            "source_id": self.source_id,
            "original_expression_file": self.original_expression_file,
            "original_metadata_file": self.original_metadata_file,
            "expression_path": self.expression_path,
            "metadata_path": self.metadata_path,
            "expression_sha256": self.expression_sha256,
            "metadata_sha256": self.metadata_sha256,
            "expression_file_size": self.expression_file_size,
            "metadata_file_size": self.metadata_file_size,
            "delimiter": self.delimiter,
            "compression": self.compression.value,
            "n_genes": self.n_genes,
            "n_samples": self.n_samples,
            "sample_ids": self.sample_ids,
            "conditions": self.conditions,
            "replicate_counts": self.replicate_counts,
            "batch_values": {k: v for k, v in self.batch_values.items()},
            "value_type": self.value_type.value,
            "gene_id_type": self.gene_id_type.value,
            "gene_id_inspection": {
                "id_type": self.gene_id_inspection.id_type.value,
                "original_ids": self.gene_id_inspection.original_ids,
                "normalized_ids": self.gene_id_inspection.normalized_ids,
                "duplicate_ids": self.gene_id_inspection.duplicate_ids,
                "missing_ids": self.gene_id_inspection.missing_ids,
                "ensembl_count": self.gene_id_inspection.ensembl_count,
                "hgnc_count": self.gene_id_inspection.hgnc_count,
                "entrez_count": self.gene_id_inspection.entrez_count,
                "unknown_count": self.gene_id_inspection.unknown_count,
            },
            "value_classification": {
                "value_type": self.value_classification.value_type.value,
                "total_values": self.value_classification.total_values,
                "integer_count": self.value_classification.integer_count,
                "non_integer_count": self.value_classification.non_integer_count,
                "zero_count": self.value_classification.zero_count,
                "negative_count": self.value_classification.negative_count,
                "min_value": self.value_classification.min_value,
                "max_value": self.value_classification.max_value,
                "warnings": self.value_classification.warnings,
            },
            "warnings": self.warnings,
            "ingested_at": self.ingested_at,
            "configuration_hash": self.configuration_hash,
            "software_version": self.software_version,
        }


# ---------------------------------------------------------------------------
# Expression matrix loading
# ---------------------------------------------------------------------------


def _detect_compression(path: Path) -> CompressionType:
    """Detect compression from filename extension."""
    if path.name.endswith(".gz"):
        return CompressionType.GZIP
    return CompressionType.NONE


def _detect_delimiter(path: Path) -> str:
    """Detect whether a file is tab- or comma-delimited.

    Reads the first non-comment line and checks for tab vs comma.
    """
    opener = gzip.open if _detect_compression(path) == CompressionType.GZIP else open
    with opener(path, "rt", encoding="utf-8") as f:
        first_line = f.readline()

    tab_count = first_line.count("\t")
    comma_count = first_line.count(",")

    if tab_count >= comma_count and tab_count > 0:
        return "\t"
    if comma_count > 0:
        return ","
    # Default to tab if neither found
    return "\t"


def _open_file(path: Path):
    """Return an appropriate file handle based on compression."""
    if _detect_compression(path) == CompressionType.GZIP:
        return gzip.open(path, "rt", encoding="utf-8")
    return open(path, encoding="utf-8")


class ExpressionFileError(ValueError):
    """Raised when an expression file fails validation."""


def load_expression_matrix(
    path: Path,
    *,
    value_type: ValueType | None = None,
) -> ExpressionLoadResult:
    """Load an expression matrix from TSV or CSV (optionally gzip-compressed).

    Parameters
    ----------
    path : Path
        Path to the expression file.
    value_type : ValueType, optional
        Override for value classification.  If None, classification is
        determined programmatically.

    Returns
    -------
    ExpressionLoadResult
        Loaded dataframe and metadata.

    Raises
    ------
    ExpressionFileError
        If the file is empty, malformed, or contains duplicate sample names.
    FileNotFoundError
        If the file does not exist.
    """
    if not path.exists():
        raise FileNotFoundError(f"Expression file not found: {path}")

    compression = _detect_compression(path)
    delimiter = _detect_delimiter(path)

    opener = gzip.open if compression == CompressionType.GZIP else open
    with opener(path, "rt", encoding="utf-8") as f:
        content = f.read()

    reader = csv.reader(io.StringIO(content), delimiter=delimiter)
    rows = list(reader)

    if not rows:
        raise ExpressionFileError("Expression file is empty")

    header = rows[0]
    if len(header) < 2:
        raise ExpressionFileError(
            "Expression file must have at least 2 columns"
            f" (gene + 1 sample), got {len(header)}"
        )

    # First column is the gene identifier column
    sample_names = [h.strip() for h in header[1:]]

    # Reject duplicate sample names
    seen: set[str] = set()
    for name in sample_names:
        if name in seen:
            raise ExpressionFileError(f"Duplicate sample column name: {name}")
        seen.add(name)

    # Reject empty sample names
    for name in sample_names:
        if not name:
            raise ExpressionFileError(
                "Expression file contains empty sample column name(s)"
            )

    # Build dataframe
    data_rows = []
    for i, row in enumerate(rows[1:], start=2):
        if len(row) != len(header):
            raise ExpressionFileError(
                f"Line {i}: expected {len(header)} columns, got {len(row)}"
            )
        data_rows.append(row)

    if not data_rows:
        raise ExpressionFileError("Expression file contains no data rows")

    # Build dataframe with proper dtypes from the start
    # First column as strings (gene IDs), rest as float
    gene_ids_str = [str(g).strip() for g in [row[0] for row in data_rows]]

    # Convert numeric columns properly - create new dataframe to avoid dtype issues
    df = pd.DataFrame(data_rows, columns=header)
    numeric_data = {}
    for col_name in header[1:]:
        try:
            numeric_data[col_name] = pd.to_numeric(df[col_name], errors="raise")
        except (ValueError, TypeError) as exc:
            raise ExpressionFileError(
                f"Non-numeric value found in column '{col_name}': {exc}"
            ) from exc
    numeric_df = pd.DataFrame(numeric_data)
    # Insert gene column
    numeric_df.insert(0, header[0], gene_ids_str)
    df = numeric_df

    # Reject missing gene identifiers
    missing = [g for g in gene_ids_str if not g]
    if missing:
        raise ExpressionFileError(f"Missing gene identifiers at rows: {missing[:10]}")

    result = ExpressionLoadResult(
        dataframe=df,
        n_genes=len(data_rows),
        n_samples=len(sample_names),
        sample_ids=sample_names,
        gene_ids=gene_ids_str,
        delimiter="\t" if delimiter == "\t" else ",",
        compression=compression,
        original_filename=path.name,
    )

    return result


# ---------------------------------------------------------------------------
# Sample metadata loading
# ---------------------------------------------------------------------------


class MetadataFileError(ValueError):
    """Raised when a metadata file fails validation."""


def _load_metadata_from_json(path: Path) -> pd.DataFrame:
    """Load sample metadata from a JSON file."""
    import json

    with open(path, encoding="utf-8") as f:
        data = json.load(f)

    if not isinstance(data, dict):
        raise MetadataFileError("JSON metadata must be an object")

    if "samples" not in data:
        raise MetadataFileError("JSON metadata must contain a 'samples' key")

    samples = data["samples"]
    if not isinstance(samples, dict):
        raise MetadataFileError(
            "'samples' must map sample_id to metadata"
        )

    records = []
    for sample_id, meta in samples.items():
        record = {"sample_id": sample_id}
        record.update(meta)
        records.append(record)

    return pd.DataFrame(records)


def _load_metadata_from_tsv(path: Path) -> pd.DataFrame:
    """Load sample metadata from a TSV or CSV file."""
    compression = _detect_compression(path)
    delimiter = _detect_delimiter(path)

    opener = gzip.open if compression == CompressionType.GZIP else open
    with opener(path, "rt", encoding="utf-8") as f:
        df = pd.read_csv(f, delimiter=delimiter)

    return df


def load_sample_metadata(
    path: Path,
    *,
    expression_sample_ids: list[str] | None = None,
) -> MetadataLoadResult:
    """Load and validate sample metadata.

    Parameters
    ----------
    path : Path
        Path to the metadata file (JSON, TSV, or CSV).
    expression_sample_ids : list of str, optional
        Sample IDs from the expression matrix for cross-validation.

    Returns
    -------
    MetadataLoadResult
        Loaded metadata and validation results.

    Raises
    ------
    MetadataFileError
        If the metadata is invalid.
    FileNotFoundError
        If the file does not exist.
    """
    if not path.exists():
        raise FileNotFoundError(f"Metadata file not found: {path}")

    if path.suffix == ".json":
        df = _load_metadata_from_json(path)
    else:
        df = _load_metadata_from_tsv(path)

    if df.empty:
        raise MetadataFileError("Metadata file contains no data rows")

    # Ensure sample_id column exists
    if "sample_id" not in df.columns:
        # Try using the first column as sample_id
        first_col = df.columns[0]
        df = df.rename(columns={first_col: "sample_id"})
        df["sample_id"] = df["sample_id"].astype(str).str.strip()

    # Check for duplicate sample IDs
    dupes = df[df.duplicated(subset=["sample_id"], keep=False)]
    if not dupes.empty:
        dup_ids = dupes["sample_id"].unique().tolist()
        raise MetadataFileError(f"Duplicate sample IDs in metadata: {dup_ids[:10]}")

    # Ensure condition column exists
    if "condition" not in df.columns:
        raise MetadataFileError("Metadata must contain a 'condition' column")

    sample_ids = set(df["sample_id"].tolist())
    conditions = dict(zip(df["sample_id"], df["condition"]))

    # Optional fields
    cell_lines: dict[str, str | None] = {}
    replicates: dict[str, int | None] = {}
    batch_values: dict[str, str | None] = {}

    for _, row in df.iterrows():
        sid = str(row["sample_id"])
        if "cell_line" in df.columns and pd.notna(row.get("cell_line")):
            cell_lines[sid] = str(row["cell_line"])
        else:
            cell_lines[sid] = None
        if "replicate" in df.columns and pd.notna(row.get("replicate")):
            replicates[sid] = int(row["replicate"])
        else:
            replicates[sid] = None
        if "batch" in df.columns and pd.notna(row.get("batch")):
            batch_values[sid] = str(row["batch"])
        else:
            batch_values[sid] = None

    # Cross-validation with expression samples
    if expression_sample_ids is not None:
        expr_set = set(expression_sample_ids)

        # Every expression sample must have a metadata row
        missing_meta = expr_set - sample_ids
        if missing_meta:
            raise MetadataFileError(
                f"Expression samples missing in metadata: {sorted(missing_meta)[:10]}"
            )

        # Check for unexplained extra metadata samples
        extra_meta = sample_ids - expr_set
        if extra_meta:
            extra = sorted(extra_meta)[:10]
            raise MetadataFileError(
                f"Metadata contains extra samples not in expression: {extra}"
            )

    return MetadataLoadResult(
        dataframe=df,
        sample_ids=sample_ids,
        conditions=conditions,
        cell_lines=cell_lines,
        replicates=replicates,
        batch_values=batch_values,
        original_filename=path.name,
    )


# ---------------------------------------------------------------------------
# Expression-value classification
# ---------------------------------------------------------------------------


def classify_expression_values(
    df: pd.DataFrame,
    *,
    value_type_override: ValueType | None = None,
) -> ValueClassificationResult:
    """Classify expression values conservatively.

    Parameters
    ----------
    df : pd.DataFrame
        Expression dataframe (first column is gene IDs, rest are numeric).
    value_type_override : ValueType, optional
        If provided, skip programmatic classification and use this value.

    Returns
    -------
    ValueClassificationResult
        Classification result with statistics.
    """
    numeric_cols = df.iloc[:, 1:]
    values = numeric_cols.values.flatten().astype(float)

    total = len(values)
    zeros = int((values == 0).sum())
    negatives = int((values < 0).sum())

    # Check for integer vs non-integer
    non_zero_mask = values != 0
    non_zero_values = values[non_zero_mask]

    if len(non_zero_values) == 0:
        integer_count = 0
        non_integer_count = 0
    else:
        integer_count = int((non_zero_values == non_zero_values.astype(int)).sum())
        non_integer_count = len(non_zero_values) - integer_count

    min_val = float(values.min())
    max_val = float(values.max())

    warnings: list[str] = []

    if value_type_override is not None:
        return ValueClassificationResult(
            value_type=value_type_override,
            total_values=total,
            integer_count=integer_count,
            non_integer_count=non_integer_count,
            zero_count=zeros,
            negative_count=negatives,
            min_value=min_val,
            max_value=max_val,
            warnings=warnings,
        )

    # Programmatic classification
    if negatives > 0:
        # Has negative values — likely transformed (log, z-score, etc.)
        value_type = ValueType.TRANSFORMED_VALUES
        if negatives < total * 0.01:
            warnings.append(
                f"Found {negatives} negative values ({negatives/total*100:.1f}%); "
                "classified as transformed_values"
            )
    elif integer_count == len(non_zero_values) and len(non_zero_values) > 0:
        value_type = ValueType.RAW_INTEGER_COUNTS
    elif non_integer_count > 0 and min_val >= 0:
        # Non-integer, non-negative
        nz_count = len(non_zero_values)
        non_integer_fraction = non_integer_count / nz_count if nz_count > 0 else 0
        if non_integer_fraction > 0.5:
            value_type = ValueType.NON_INTEGER_ESTIMATED_COUNTS
            warnings.append(
                f"{non_integer_fraction*100:.1f}% of non-zero values are non-integer; "
                "classified as non_integer_estimated_counts. "
                "DE method selection requires verification of quantification source."
            )
        else:
            # Mostly integer but some non-integer — could be normalized
            value_type = ValueType.NORMALIZED_NONNEGATIVE_VALUES
            warnings.append(
                "Mixture of integer and non-integer values; "
                "classified as normalized_nonnegative_values"
            )
    elif min_val >= 0 and integer_count == 0 and len(non_zero_values) > 0:
        value_type = ValueType.NORMALIZED_NONNEGATIVE_VALUES
    else:
        value_type = ValueType.UNKNOWN
        warnings.append(
            "Could not determine value type; classified as unknown. "
            "Provide --value-type to override."
        )

    return ValueClassificationResult(
        value_type=value_type,
        total_values=total,
        integer_count=integer_count,
        non_integer_count=non_integer_count,
        zero_count=zeros,
        negative_count=negatives,
        min_value=min_val,
        max_value=max_val,
        warnings=warnings,
    )


# ---------------------------------------------------------------------------
# Gene identifier inspection
# ---------------------------------------------------------------------------


def inspect_gene_identifiers(
    gene_ids: list[str],
    *,
    gene_id_type_override: GeneIdType | None = None,
) -> GeneIdInspectionResult:
    """Inspect and classify gene identifiers.

    Parameters
    ----------
    gene_ids : list of str
        List of gene identifiers from the first column of the expression matrix.
    gene_id_type_override : GeneIdType, optional
        If provided, skip detection and use this classification.

    Returns
    -------
    GeneIdInspectionResult
        Inspection result with classifications.
    """
    ensembl_re = re.compile(ENSEMBL_PATTERN)
    entrez_re = re.compile(ENTREZ_PATTERN)
    hgnc_re = re.compile(HGNC_PATTERN)

    ensembl_ids = []
    hgnc_symbols = []
    entrez_ids = []
    unknown_ids = []
    normalized_ids = []

    for gid in gene_ids:
        gid = str(gid).strip()
        if not gid:
            continue

        if ensembl_re.match(gid):
            ensembl_ids.append(gid)
            # Strip version suffix for normalized field
            normalized = re.sub(r"\.[0-9]+$", "", gid)
            normalized_ids.append(normalized)
        elif entrez_re.match(gid) and len(gid) <= 10:
            entrez_ids.append(gid)
            normalized_ids.append(gid)
        elif hgnc_re.match(gid):
            hgnc_symbols.append(gid)
            normalized_ids.append(gid)
        else:
            unknown_ids.append(gid)
            normalized_ids.append(gid)

    # Determine overall type
    categories_with_hits = sum([
        len(ensembl_ids) > 0,
        len(hgnc_symbols) > 0,
        len(entrez_ids) > 0,
        len(unknown_ids) > 0,
    ])

    if gene_id_type_override is not None:
        id_type = gene_id_type_override
    elif categories_with_hits > 1:
        id_type = GeneIdType.MIXED
    elif len(ensembl_ids) > 0:
        id_type = GeneIdType.ENSEMBL
    elif len(hgnc_symbols) > 0:
        id_type = GeneIdType.HGNC_SYMBOLS
    elif len(entrez_ids) > 0:
        id_type = GeneIdType.ENTREZ_IDS
    else:
        id_type = GeneIdType.UNKNOWN

    # Detect duplicates
    seen: set[str] = set()
    duplicates: set[str] = set()
    for gid in gene_ids:
        gid = str(gid).strip()
        if gid in seen:
            duplicates.add(gid)
        seen.add(gid)

    # Detect missing
    missing = [gid for gid in gene_ids if not str(gid).strip()]

    return GeneIdInspectionResult(
        id_type=id_type,
        original_ids=gene_ids,
        normalized_ids=normalized_ids,
        duplicate_ids=sorted(duplicates),
        missing_ids=missing,
        ensembl_count=len(ensembl_ids),
        hgnc_count=len(hgnc_symbols),
        entrez_count=len(entrez_ids),
        unknown_count=len(unknown_ids),
    )


# ---------------------------------------------------------------------------
# Replicate counting
# ---------------------------------------------------------------------------


def count_replicates_per_condition(conditions: dict[str, str]) -> dict[str, int]:
    """Count samples per condition.

    Parameters
    ----------
    conditions : dict
        Mapping of sample_id -> condition.

    Returns
    -------
    dict
        Mapping of condition -> sample count.
    """
    counts: dict[str, int] = {}
    for cond in conditions.values():
        counts[cond] = counts.get(cond, 0) + 1
    return counts


# ---------------------------------------------------------------------------
# Contrast validation
# ---------------------------------------------------------------------------


class ContrastError(ValueError):
    """Raised when a contrast cannot be validated."""


def validate_contrast(
    conditions: dict[str, str],
    control: str,
    treatment: str,
) -> dict[str, Any]:
    """Validate that a contrast has samples in both groups.

    Parameters
    ----------
    conditions : dict
        Mapping of sample_id -> condition.
    control : str
        Control condition name.
    treatment : str
        Treatment condition name.

    Returns
    -------
    dict
        Summary with control_count, treatment_count, and valid bool.

    Raises
    ------
    ContrastError
        If either group has no samples.
    """
    control_samples = [s for s, c in conditions.items() if c == control]
    treatment_samples = [s for s, c in conditions.items() if c == treatment]

    if not control_samples:
        raise ContrastError(f"No samples found for control condition: {control}")
    if not treatment_samples:
        raise ContrastError(f"No samples found for treatment condition: {treatment}")

    return {
        "control": control,
        "treatment": treatment,
        "control_count": len(control_samples),
        "treatment_count": len(treatment_samples),
        "control_samples": sorted(control_samples),
        "treatment_samples": sorted(treatment_samples),
        "valid": True,
    }


# ---------------------------------------------------------------------------
# Full ingestion pipeline
# ---------------------------------------------------------------------------


class IngestionError(RuntimeError):
    """Raised when ingestion fails."""


def ingest(
    expression_path: Path,
    metadata_path: Path,
    source_id: str,
    run_dir: Path,
    *,
    value_type_override: ValueType | None = None,
    gene_id_type_override: GeneIdType | None = None,
    contrast_control: str | None = None,
    contrast_treatment: str | None = None,
) -> IngestionResult:
    """Run the full ingestion pipeline.

    Parameters
    ----------
    expression_path : Path
        Path to the expression matrix file.
    metadata_path : Path
        Path to the sample metadata file.
    source_id : str
        Source identifier (e.g., GSE accession).
    run_dir : Path
        Run directory for storing artifacts.
    value_type_override : ValueType, optional
        Override value classification.
    gene_id_type_override : GeneIdType, optional
        Override gene ID classification.
    contrast_control : str, optional
        Control condition for contrast validation.
    contrast_treatment : str, optional
        Treatment condition for contrast validation.

    Returns
    -------
    IngestionResult
        Complete ingestion result.
    """
    run_dir.mkdir(parents=True, exist_ok=True)

    # Store artifacts
    expr_dest = run_dir / "input" / expression_path.name
    expr_dest.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(expression_path, expr_dest)

    meta_dest = run_dir / "input" / metadata_path.name
    shutil.copy2(metadata_path, meta_dest)

    # Compute hashes on original files
    expr_hash = hash_file_sha256(expression_path)
    meta_hash = hash_file_sha256(metadata_path)
    expr_size = expression_path.stat().st_size
    meta_size = metadata_path.stat().st_size

    # Load expression matrix
    expr_result = load_expression_matrix(expression_path)

    # Load and validate metadata
    meta_result = load_sample_metadata(
        metadata_path,
        expression_sample_ids=expr_result.sample_ids,
    )

    # Classify values
    value_class = classify_expression_values(
        expr_result.dataframe,
        value_type_override=value_type_override,
    )

    # Inspect gene IDs
    gene_inspection = inspect_gene_identifiers(
        expr_result.gene_ids,
        gene_id_type_override=gene_id_type_override,
    )

    # Replicate counts
    replicate_counts = count_replicates_per_condition(meta_result.conditions)

    # Contrast validation (optional)
    if contrast_control and contrast_treatment:
        validate_contrast(
            meta_result.conditions, contrast_control, contrast_treatment
        )

    warnings = list(value_class.warnings)
    if gene_inspection.duplicate_ids:
        warnings.append(
            "Duplicate gene identifiers found:"
            f" {gene_inspection.duplicate_ids[:10]}"
        )
    if gene_inspection.id_type == GeneIdType.MIXED:
        warnings.append(
            "Mixed gene identifier types detected;"
            " identifier normalization may be incomplete"
        )

    # Configuration hash
    config_str = f"{source_id}:{value_type_override}:{gene_id_type_override}"
    config_hash = hashlib.sha256(config_str.encode()).hexdigest()[:16]

    result = IngestionResult(
        source_id=source_id,
        original_expression_file=expression_path.name,
        original_metadata_file=metadata_path.name,
        expression_path=str(expr_dest.relative_to(run_dir)).replace("\\", "/"),
        metadata_path=str(meta_dest.relative_to(run_dir)).replace("\\", "/"),
        expression_sha256=expr_hash,
        metadata_sha256=meta_hash,
        expression_file_size=expr_size,
        metadata_file_size=meta_size,
        delimiter=expr_result.delimiter,
        compression=expr_result.compression,
        n_genes=expr_result.n_genes,
        n_samples=expr_result.n_samples,
        sample_ids=expr_result.sample_ids,
        conditions=meta_result.conditions,
        replicate_counts=replicate_counts,
        batch_values=meta_result.batch_values,
        value_type=value_class.value_type,
        gene_id_type=gene_inspection.id_type,
        gene_id_inspection=gene_inspection,
        value_classification=value_class,
        warnings=warnings,
        configuration_hash=config_hash,
    )

    return result


# ---------------------------------------------------------------------------
# Manifest writing
# ---------------------------------------------------------------------------


def write_ingestion_manifest(
    run_dir: Path,
    result: IngestionResult,
) -> Path:
    """Write the ingestion manifest JSON to the run directory.

    Parameters
    ----------
    run_dir : Path
        The run directory.
    result : IngestionResult
        The ingestion result to serialize.

    Returns
    -------
    Path
        Path to the written manifest file.
    """
    manifest_path = run_dir / "ingestion_manifest.json"
    write_json(manifest_path, result.to_dict())
    return manifest_path
