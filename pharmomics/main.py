"""PharmOmics CLI entry point.

Provides the ``pharmomics`` command via Typer.  Milestone 1 only implements
the ``init`` subcommand which creates a new run directory.
"""

from __future__ import annotations

import typer

from pharmomics.config import load_settings
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


if __name__ == "__main__":
    app()
