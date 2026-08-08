"""CLI ``analyze`` command — run differential analysis on real data files.

Orchestration only: builds domain objects from raw files via existing
ingest/parser modules, then delegates to ``run_analysis()`` and returns
an ``AnalysisResult`` for downstream rendering.

No statistical or scientific logic is implemented here.
"""

from __future__ import annotations

from pathlib import Path

from pharmomics.analysis.results import AnalysisResult
from pharmomics.analysis.run import run_analysis
from pharmomics.experiment.enums import FactorType, GroupRole
from pharmomics.experiment.schemas import (
    Contrast,
    DesignSample,
    ExperimentalFactor,
    ExperimentalGroup,
    ExperimentDesign,
)
from pharmomics.ingestion.loader import (
    classify_expression_values,
    inspect_gene_identifiers,
    load_expression_matrix,
    load_sample_metadata,
)
from pharmomics.omics.adapter import from_load_results

# ---------------------------------------------------------------------------
# Domain-object builders
# ---------------------------------------------------------------------------


def _slug(value: str) -> str:
    """Derive a safe group/contrast identifier from a condition string.

    Replaces non-alphanumeric characters with underscores and lowercases.
    """
    return "".join(c if c.isalnum() else "_" for c in value).lower().rstrip("_")


def build_experiment_design(
    sample_ids: list[str],
    conditions: dict[str, str],
    *,
    experiment_id: str = "exp-local",
    description: str | None = None,
    contrast_control: str = "",
    contrast_treatment: str = "",
) -> ExperimentDesign:
    """Construct an ``ExperimentDesign`` from flat metadata.

    Parameters
    ----------
    sample_ids : list[str]
        Ordered sample IDs (from the expression matrix header).
    conditions : dict[str, str]
        Mapping of sample_id → condition string.
    experiment_id : str, optional
        Experiment identifier.
    description : str, optional
        Human-readable description.
    contrast_control : str
        Condition name to use as the control (reference) group.
    contrast_treatment : str
        Condition name to use as the treatment (comparison) group.

    Returns
    -------
    ExperimentDesign
        A fully populated design with one ``ExperimentalGroup`` per
        unique condition and a single ``Contrast``.

    Raises
    ------
    ValueError
        If a sample_id has no matching condition, or if the control
        or treatment condition has no samples.
    """
    if not sample_ids:
        raise ValueError("sample_ids must not be empty")

    # --- Validate all samples have conditions ---
    for sid in sample_ids:
        if sid not in conditions:
            raise ValueError(f"Sample '{sid}' has no condition assigned")

    # --- Collect unique conditions and build groups ---
    unique_conditions: list[str] = list(dict.fromkeys(conditions.values()))

    groups: list[ExperimentalGroup] = []
    for cond in unique_conditions:
        group_id = _slug(cond)
        role = GroupRole.CONTROL if cond == contrast_control else GroupRole.TREATMENT
        if cond not in (contrast_control, contrast_treatment):
            role = GroupRole.OBSERVATIONAL
        groups.append(
            ExperimentalGroup(
                group_id=group_id,
                label=cond,
                role=role,
            )
        )

    # --- Validate control and treatment conditions have samples ---
    observed_conditions = set(conditions.get(sid) for sid in sample_ids)
    if contrast_control and contrast_control not in observed_conditions:
        raise ValueError(f"Control condition '{contrast_control}' has no samples")
    if contrast_treatment and contrast_treatment not in observed_conditions:
        raise ValueError(f"Treatment condition '{contrast_treatment}' has no samples")

    # --- Build samples ---
    samples: list[DesignSample] = []
    for sid in sample_ids:
        cond = conditions[sid]
        group_id = _slug(cond)
        samples.append(
            DesignSample(
                sample_id=sid,
                group_id=group_id,
                factor_values={"condition": cond},
            )
        )

    # --- Build factor ---
    factor = ExperimentalFactor(
        factor_id="condition",
        factor_type=FactorType.CATEGORICAL,
        description="Treatment condition",
        levels=unique_conditions,
    )

    # --- Build contrast ---
    control_group_id = _slug(contrast_control)
    treatment_group_id = _slug(contrast_treatment)
    contrast_id = f"{treatment_group_id}_vs_{control_group_id}"

    contrasts = [
        Contrast(
            contrast_id=contrast_id,
            comparison_group_id=treatment_group_id,
            reference_group_id=control_group_id,
            description=f"{contrast_treatment} vs {contrast_control}",
        )
    ]

    return ExperimentDesign(
        experiment_id=experiment_id,
        description=description,
        samples=samples,
        groups=groups,
        factors=[factor],
        contrasts=contrasts,
    )


def build_analysis_spec(
    *,
    analysis_type: str = "differential_analysis",
    contrast_id: str = "",
    fdr_threshold: float = 0.05,
) -> object:
    """Construct an ``AnalysisSpecification``.

    Parameters
    ----------
    analysis_type : str
        Type of analysis to run (default: ``"differential_analysis"``).
    contrast_id : str
        The contrast ID that was created in ``build_experiment_design``.
    fdr_threshold : float
        BH-FDR significance threshold (default: 0.05).

    Returns
    -------
    AnalysisSpecification
        The analysis intent for ``run_analysis()``.
    """
    from pharmomics.analysis.schemas import AnalysisSpecification

    return AnalysisSpecification(
        analysis_type=analysis_type,
        factor_references=["condition"],
        contrast_references=[contrast_id],
        parameters={"fdr_threshold": fdr_threshold},
    )


# ---------------------------------------------------------------------------
# CLI orchestration
# ---------------------------------------------------------------------------


def analyze(
    expression_file: Path,
    metadata_file: Path,
    contrast_control: str,
    contrast_treatment: str,
    output: Path,
    *,
    source_id: str = "local",
) -> AnalysisResult:
    """Run the full analysis pipeline on real data files.

    Parameters
    ----------
    expression_file : Path
        Path to the expression matrix file (TSV/CSV/gzip).
    metadata_file : Path
        Path to the sample metadata file (JSON/TSV/CSV).
    contrast_control : str
        Condition name for the reference group.
    contrast_treatment : str
        Condition name for the comparison group.
    output : Path
        Output path for the Markdown report (unused for the return
        value; kept for API compatibility).
    source_id : str, optional
        Provenance source identifier.

    Returns
    -------
    AnalysisResult
        The completed analysis result.  All values are computed by the
        analysis engine; this function does no statistical computation.

    Raises
    ------
    FileNotFoundError
        If input files do not exist.
    ExpressionFileError
        If the expression matrix is malformed.
    MetadataFileError
        If the metadata is malformed or inconsistent.
    AnalysisValidationError
        If the analysis inputs fail cross-domain validation.
    ValueError
        If control or treatment condition has no samples.
    """
    # Load expression matrix
    expr_result = load_expression_matrix(expression_file)

    # Load and cross-validate metadata
    meta_result = load_sample_metadata(
        metadata_file,
        expression_sample_ids=expr_result.sample_ids,
    )

    # Classify expression values
    value_class = classify_expression_values(expr_result.dataframe)

    # Inspect gene identifiers
    gene_inspection = inspect_gene_identifiers(expr_result.gene_ids)

    # Build OmicsMatrix via existing adapter
    omics = from_load_results(
        expr_result,
        meta_result,
        value_class=value_class,
        gene_inspection=gene_inspection,
        source_id=source_id,
    )

    # Build ExperimentDesign
    design = build_experiment_design(
        sample_ids=expr_result.sample_ids,
        conditions=meta_result.conditions,
        contrast_control=contrast_control,
        contrast_treatment=contrast_treatment,
    )

    # Build AnalysisSpecification
    contrast_id = design.contrasts[0].contrast_id
    spec = build_analysis_spec(contrast_id=contrast_id)

    # Run analysis (includes full validation chain)
    result = run_analysis(spec, design, omics)

    return result
