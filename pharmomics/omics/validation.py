"""Validation rules for ``OmicsMatrix`` domain objects.

All validators operate purely on in-memory objects — no file I/O, no
network calls.  Each function returns a list of human-readable violation
strings (empty list = valid).
"""

from __future__ import annotations

from pharmomics.omics.schemas import OmicsMatrix

KNOWN_SCHEMA_VERSIONS = frozenset({"1.0.0"})


def validate(matrix: OmicsMatrix) -> list[str]:
    """Run all domain validations and return any violations.

    Parameters
    ----------
    matrix : OmicsMatrix
        The omics matrix to validate.

    Returns
    -------
    list[str]
        Empty list if the matrix is valid; otherwise a list of
        violation descriptions.
    """
    violations: list[str] = []
    violations.extend(_check_schema_version(matrix))
    violations.extend(_check_counts_consistency(matrix))
    violations.extend(_check_feature_metadata_coverage(matrix))
    violations.extend(_check_sample_metadata_coverage(matrix))
    violations.extend(_check_dataframe_shape(matrix))
    violations.extend(_check_unique_ids(matrix))
    violations.extend(_check_nonempty_ids(matrix))
    return violations


# ---------------------------------------------------------------------------
# Individual validators (private)
# ---------------------------------------------------------------------------


def _check_schema_version(matrix: OmicsMatrix) -> list[str]:
    """Reject unknown schema versions."""
    if matrix.schema_version not in KNOWN_SCHEMA_VERSIONS:
        return [
            f"Unknown schema version: {matrix.schema_version!r} "
            f"(known: {sorted(KNOWN_SCHEMA_VERSIONS)})"
        ]
    return []


def _check_counts_consistency(matrix: OmicsMatrix) -> list[str]:
    """Ensure n_features / n_samples match feature_ids / sample_ids lengths."""
    violations: list[str] = []
    if matrix.n_features != len(matrix.feature_ids):
        violations.append(
            f"n_features ({matrix.n_features}) does not match "
            f"feature_ids length ({len(matrix.feature_ids)})"
        )
    if matrix.n_samples != len(matrix.sample_ids):
        violations.append(
            f"n_samples ({matrix.n_samples}) does not match "
            f"sample_ids length ({len(matrix.sample_ids)})"
        )
    return violations


def _check_feature_metadata_coverage(matrix: OmicsMatrix) -> list[str]:
    """Every feature_id must have a corresponding feature_metadata entry."""
    missing = set(matrix.feature_ids) - set(matrix.feature_metadata.keys())
    if missing:
        example = sorted(missing)[:5]
        return [f"{len(missing)} feature(s) missing from feature_metadata: {example}"]
    return []


def _check_sample_metadata_coverage(matrix: OmicsMatrix) -> list[str]:
    """Every sample_id must have a corresponding sample_metadata entry."""
    missing = set(matrix.sample_ids) - set(matrix.sample_metadata.keys())
    if missing:
        example = sorted(missing)[:5]
        return [f"{len(missing)} sample(s) missing from sample_metadata: {example}"]
    return []


def _check_dataframe_shape(matrix: OmicsMatrix) -> list[str]:
    """DataFrame must have n_features rows and n_samples+1 columns."""
    df = matrix.dataframe
    violations: list[str] = []

    expected_rows = matrix.n_features
    expected_cols = matrix.n_samples + 1  # feature_id column + samples

    if len(df) != expected_rows:
        violations.append(
            f"DataFrame has {len(df)} rows but n_features={expected_rows}"
        )
    if len(df.columns) != expected_cols:
        violations.append(
            f"DataFrame has {len(df.columns)} columns but expected "
            f"{expected_cols} (n_samples={matrix.n_samples} + 1)"
        )
    return violations


def _check_unique_ids(matrix: OmicsMatrix) -> list[str]:
    """feature_ids and sample_ids must not contain duplicates."""
    violations: list[str] = []

    dup_features = _find_duplicates(matrix.feature_ids)
    if dup_features:
        violations.append(f"Duplicate feature_ids: {sorted(dup_features)[:5]}")

    dup_samples = _find_duplicates(matrix.sample_ids)
    if dup_samples:
        violations.append(f"Duplicate sample_ids: {sorted(dup_samples)[:5]}")

    return violations


def _check_nonempty_ids(matrix: OmicsMatrix) -> list[str]:
    """No feature_id or sample_id may be an empty string."""
    violations: list[str] = []

    empty_features = [i for i, f in enumerate(matrix.feature_ids) if not f.strip()]
    if empty_features:
        violations.append(f"Empty feature_id at index/indices: {empty_features[:5]}")

    empty_samples = [i for i, s in enumerate(matrix.sample_ids) if not s.strip()]
    if empty_samples:
        violations.append(f"Empty sample_id at index/indices: {empty_samples[:5]}")

    return violations


# ---------------------------------------------------------------------------
# Utility
# ---------------------------------------------------------------------------


def _find_duplicates(ids: list[str]) -> set[str]:
    """Return the set of duplicated values in *ids*."""
    seen: set[str] = set()
    dups: set[str] = set()
    for item in ids:
        if item in seen:
            dups.add(item)
        seen.add(item)
    return dups
