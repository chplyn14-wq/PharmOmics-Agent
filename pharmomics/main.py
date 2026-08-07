"""PharmOmics CLI entry point.

Provides the ``pharmomics`` command via Typer.  Milestone 1 implements
the ``init`` and ``status`` subcommands plus the ``ingest`` command
for local data ingestion and validation (Milestone 1B).
"""

from __future__ import annotations

from pathlib import Path

import typer

from pharmomics.analysis.example_data import make_demo_inputs
from pharmomics.analysis.render import render_markdown_report
from pharmomics.analysis.run import run_analysis
from pharmomics.analysis.runner import AnalysisValidationError
from pharmomics.cli.analyze import analyze as run_analyze
from pharmomics.config import load_settings
from pharmomics.ingestion.loader import (
    GeneIdType,
    IngestionError,
    ValueType,
    write_ingestion_manifest,
)
from pharmomics.ingestion.loader import ingest as run_ingest
from pharmomics.run_store import create_run_directory, generate_run_id

app = typer.Typer(
    name="pharmomics",
    help="Human-in-the-loop research assistant for drug-response transcriptomics.",
    add_completion=False,
)


@app.command()
def init() -> None:
    """Create a new run directory and print its path."""
    settings = load_settings()
    run_id = generate_run_id()
    run_dir = create_run_directory(settings.resolved_run_store_dir(), run_id)
    typer.echo(f"Created run directory: {run_dir}")
    typer.echo(f"Run ID: {run_id}")


@app.command()
def status() -> None:
    """Print current configuration summary."""
    settings = load_settings()
    typer.echo("PharmOmics configuration:")
    for key, value in settings.snapshot().items():
        typer.echo(f"  {key}: {value}")


@app.command()
def ingest(
    expression_file: Path = typer.Option(
        ...,
        "--expression-file",
        "-e",
        help="Path to the expression matrix file"
        " (TSV/CSV, optionally gzip-compressed).",
    ),
    metadata_file: Path = typer.Option(
        ...,
        "--metadata-file",
        "-m",
        help="Path to the sample metadata file (JSON, TSV, or CSV).",
    ),
    source_id: str = typer.Option(
        ...,
        "--source-id",
        "-s",
        help="Source identifier (e.g., GEO accession).",
    ),
    run_dir: Path = typer.Option(
        ...,
        "--run-dir",
        "-r",
        help="Path to the run directory for storing artifacts.",
    ),
    value_type: ValueType | None = typer.Option(
        None,
        "--value-type",
        help="Override expression-value classification.",
    ),
    gene_id_type: GeneIdType | None = typer.Option(
        None,
        "--gene-id-type",
        help="Override gene identifier classification.",
    ),
    contrast_control: str | None = typer.Option(
        None,
        "--contrast-control",
        help="Control condition name for contrast validation.",
    ),
    contrast_treatment: str | None = typer.Option(
        None,
        "--contrast-treatment",
        help="Treatment condition name for contrast validation.",
    ),
) -> None:
    """Ingest and validate a local expression matrix and sample metadata.

    Loads the expression file, validates sample metadata, classifies
    expression values, inspects gene identifiers, and writes an ingestion
    manifest to the run directory.

    Returns a non-zero exit code on validation failure.
    """
    # Validate input files exist
    if not expression_file.exists():
        typer.echo(f"Error: Expression file not found: {expression_file}", err=True)
        raise typer.Exit(1)

    if not metadata_file.exists():
        typer.echo(f"Error: Metadata file not found: {metadata_file}", err=True)
        raise typer.Exit(1)

    # Create run directory
    run_dir.mkdir(parents=True, exist_ok=True)

    try:
        result = run_ingest(
            expression_path=expression_file,
            metadata_path=metadata_file,
            source_id=source_id,
            run_dir=run_dir,
            value_type_override=value_type,
            gene_id_type_override=gene_id_type,
            contrast_control=contrast_control,
            contrast_treatment=contrast_treatment,
        )
    except (ValueError, IngestionError) as exc:
        typer.echo(f"Ingestion failed: {exc}", err=True)
        raise typer.Exit(1) from exc

    # Write manifest
    manifest_path = write_ingestion_manifest(run_dir, result)

    # Print summary
    typer.echo("=" * 60)
    typer.echo("PharmOmics Ingestion Summary")
    typer.echo("=" * 60)
    typer.echo(f"  Source ID:              {result.source_id}")
    typer.echo(f"  Expression file:        {result.original_expression_file}")
    typer.echo(f"  Metadata file:          {result.original_metadata_file}")
    typer.echo(f"  Genes:                  {result.n_genes}")
    typer.echo(f"  Samples:                {result.n_samples}")
    delim = "tab" if result.delimiter == "\t" else result.delimiter
    typer.echo(f"  Delimiter:              {delim}")
    typer.echo(f"  Compression:            {result.compression.value}")
    typer.echo(f"  Value type:             {result.value_type.value}")
    typer.echo(f"  Gene ID type:           {result.gene_id_type.value}")
    conds = ", ".join(sorted(result.replicate_counts.keys()))
    typer.echo(f"  Conditions:             {conds}")
    typer.echo(f"  Replicates:             {dict(result.replicate_counts)}")
    if result.warnings:
        typer.echo(f"  Warnings ({len(result.warnings)}):")
        for w in result.warnings:
            typer.echo(f"    - {w}")
    typer.echo(f"  Run directory:          {run_dir}")
    typer.echo(f"  Manifest:               {manifest_path}")
    typer.echo("=" * 60)
    typer.echo("Ingestion complete.")


@app.command("analyze-demo")
def analyze_demo(
    output: Path = typer.Option(
        Path("report.md"),
        "--output",
        "-o",
        help="Output path for the Markdown report.",
    ),
) -> None:
    """Run differential analysis on built-in demo data and write a report.

    Uses deterministic demo data (6 genes, 6 samples) to exercise the
    full analysis pipeline: validation, differential analysis, Markdown
    rendering, and report persistence.

    Returns a non-zero exit code on failure.
    """
    try:
        omics, design, spec = make_demo_inputs()
        result = run_analysis(spec, design, omics)
    except AnalysisValidationError as exc:
        typer.echo(f"Analysis validation error: {exc}", err=True)
        raise typer.Exit(1) from exc

    markdown = render_markdown_report(result)

    output_parent = output.parent
    if not output_parent.exists():
        typer.echo(
            f"Error: Output directory does not exist: {output_parent}",
            err=True,
        )
        raise typer.Exit(1)

    output.write_text(markdown, encoding="utf-8")
    typer.echo(f"Report written to {output.resolve()}")


@app.command("analyze")
def analyze_cmd(
    expression_file: Path = typer.Option(
        ...,
        "--expression-file",
        "-e",
        help="Path to the expression matrix file"
        " (TSV/CSV, optionally gzip-compressed).",
    ),
    metadata_file: Path = typer.Option(
        ...,
        "--metadata-file",
        "-m",
        help="Path to the sample metadata file (JSON, TSV, or CSV).",
    ),
    contrast_control: str = typer.Option(
        ...,
        "--contrast-control",
        help="Control condition name (reference group).",
    ),
    contrast_treatment: str = typer.Option(
        ...,
        "--contrast-treatment",
        help="Treatment condition name (comparison group).",
    ),
    output: Path = typer.Option(
        Path("report.md"),
        "--output",
        "-o",
        help="Output path for the Markdown report.",
    ),
    source_id: str = typer.Option(
        "local",
        "--source-id",
        "-s",
        help="Source identifier for provenance (default: local).",
    ),
) -> None:
    """Run differential analysis on real data files and write a report.

    Loads an expression matrix and sample metadata, constructs the
    experimental design and analysis specification, runs the full
    analysis pipeline, and writes a Markdown report.

    Returns a non-zero exit code on failure.
    """
    if not expression_file.exists():
        typer.echo(
            f"Error: Expression file not found: {expression_file}",
            err=True,
        )
        raise typer.Exit(1)

    if not metadata_file.exists():
        typer.echo(
            f"Error: Metadata file not found: {metadata_file}",
            err=True,
        )
        raise typer.Exit(1)

    try:
        markdown = run_analyze(
            expression_file=expression_file,
            metadata_file=metadata_file,
            contrast_control=contrast_control,
            contrast_treatment=contrast_treatment,
            output=output,
            source_id=source_id,
        )
    except (
        FileNotFoundError,
        AnalysisValidationError,
        ValueError,
    ) as exc:
        typer.echo(f"Analysis failed: {exc}", err=True)
        raise typer.Exit(1) from exc
    except Exception as exc:
        typer.echo(f"Unexpected error: {exc}", err=True)
        raise typer.Exit(1) from exc

    output_parent = output.parent
    if not output_parent.exists():
        typer.echo(
            f"Error: Output directory does not exist: {output_parent}",
            err=True,
        )
        raise typer.Exit(1)

    output.write_text(markdown, encoding="utf-8")
    typer.echo(f"Report written to {output.resolve()}")


if __name__ == "__main__":
    app()
