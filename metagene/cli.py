#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright © 2023 Ye Chang yech1990@gmail.com
# Distributed under terms of the GNU license.
#
# Command-line interface for metagene analysis

import logging
import sys

import click
import polars as pl
from rich.console import Console
from rich.progress import (
    BarColumn,
    Progress,
    SpinnerColumn,
    TextColumn,
    TimeElapsedColumn,
)

from .annotation import map_to_transcripts, normalize_positions
from .config import BUILTIN_REFERENCES
from .download import download_references
from .download import list_references as show_references
from .gtf import load_gtf
from .io import load_reference, load_sites
from .plotting import plot_profile
from .utils import NewlineRichHandler

# Set up rich console
console = Console()


# Set up logging with custom Rich handler
logging.basicConfig(
    level=logging.INFO,
    format="%(message)s",
    datefmt="[%X]",
    handlers=[NewlineRichHandler(console=console)],
)

# Configure root logger
root_logger = logging.getLogger()
root_logger.handlers = []  # Remove any existing handlers
root_logger.addHandler(NewlineRichHandler())

# Configure metagene logger
logger = logging.getLogger("metagene")
logger.handlers = []  # Remove any existing handlers
logger.addHandler(NewlineRichHandler())
logger.propagate = False  # Prevent propagation to root logger


def parse_comma_separated_ints(ctx, param, value):
    """Parse comma-separated integers."""
    if not value:
        return []
    try:
        return [int(x.strip()) for x in value.split(",")]
    except ValueError:
        raise click.BadParameter("Must be comma-separated integers")


def parse_comma_separated_strings(ctx, param, value):
    """Parse comma-separated strings."""
    if not value:
        return []
    return [x.strip() for x in value.split(",")]


def update_progress_description(progress, task, description: str) -> None:
    """Update progress bar description."""
    progress.update(task, description=f"[cyan]{description}...")


@click.command(
    help="Run metagene analysis on genomic sites.",
    context_settings=dict(help_option_names=["-h", "--help"]),
)
@click.version_option()
@click.option(
    "--input",
    "-i",
    "input_file",
    type=click.Path(exists=True),
    help="Input file path (BED, GTF, TSV or CSV, etc.)",
)
@click.option(
    "--output",
    "-o",
    "output_file",
    type=click.Path(),
    help="Output file path (TSV, CSV)",
)
@click.option(
    "--output-score",
    "-s",
    type=click.Path(),
    help="Output file for binned score statistics",
)
@click.option(
    "--output-figure",
    "-p",
    type=click.Path(),
    help="Output file for metagene plot",
)
@click.option(
    "--reference",
    "-r",
    type=str,
    help="Built-in reference genome to use (e.g., GRCh38, GRCm39)",
)
@click.option(
    "--gtf",
    "-g",
    type=click.Path(exists=True),
    help="GTF/GFF file path for custom reference",
)
@click.option(
    "--region",
    type=click.Choice(["all", "5utr", "cds", "3utr"]),
    default="all",
    help="Region to analyze (default: all)",
)
@click.option(
    "--bins",
    "-b",
    type=int,
    default=100,
    help="Number of bins for analysis (default: 100)",
)
@click.option(
    "--with-header",
    "-H",
    is_flag=True,
    help="Input file has header line",
)
@click.option(
    "--separator",
    "-S",
    type=str,
    default="\t",
    help="Separator for input file (default: tab)",
)
@click.option(
    "--meta-columns",
    "-m",
    type=str,
    default="1,2,3,6",
    callback=parse_comma_separated_ints,
    help="Input column indices (1-based) for genomic coordinates. The columns should contain Chromosome,Start,End,Strand or Chromosome,Site,Strand",
)
@click.option(
    "--weight-columns",
    "-w",
    type=str,
    default="",
    callback=parse_comma_separated_ints,
    help="Input column indices (1-based) for weight/score values",
)
@click.option(
    "--weight-names",
    "-n",
    type=str,
    default="",
    callback=parse_comma_separated_strings,
    help="Names for weight columns",
)
@click.option(
    "--score-transform",
    type=click.Choice(["none", "log2", "log10"]),
    default="none",
    help="Transform to apply to scores (default: none)",
)
@click.option(
    "--normalize",
    is_flag=True,
    help="Normalize scores by transcript length",
)
@click.option(
    "--list",
    "list_references_flag",
    is_flag=True,
    help="List all available built-in references and exit",
)
@click.option(
    "--download",
    type=str,
    help="Download a specific reference (e.g., GRCh38) or 'all' for all references",
)
def cli(
    input_file: str,
    output_file: str,
    output_score: str | None,
    output_figure: str | None,
    reference: str | None,
    gtf: str | None,
    region: str,
    bins: int,
    with_header: bool,
    separator: str,
    meta_columns: list[int],
    weight_columns: list[int],
    weight_names: list[str],
    score_transform: str,
    normalize: bool,
    list_references_flag: bool,
    download: str | None,
) -> None:
    """Run metagene analysis on genomic sites."""
    # Check if no arguments provided (except defaults) and show help
    ctx = click.get_current_context()
    if (
        not list_references_flag
        and not download
        and not input_file
        and not output_file
        and not reference
        and not gtf
    ):
        # Check if user just ran 'metagene' with no arguments
        if len(sys.argv) == 1:
            console.print(ctx.get_help())
            return

    # Handle list references option
    if list_references_flag:
        show_references(console)
        return

    # Handle download option
    if download:
        try:
            console.print(f"[bold cyan]Downloading {download}...[/bold cyan]")
            download_references(download)
            console.print(f"[green]✓[/green] Downloaded {download}")
        except Exception as e:
            console.print(f"[red]✗[/red] Failed to download {download}: {e}")
        return

    # Validate required options for analysis
    if not input_file:
        console.print(
            "[red]✗[/red] Input file is required for analysis (use -i/--input)"
        )
        console.print(
            "Use --list to see available references or --download to download references"
        )
        sys.exit(1)

    if not output_file:
        console.print(
            "[red]✗[/red] Output file is required for analysis (use -o/--output)"
        )
        sys.exit(1)

    try:
        # Initialize variable with proper type
        exon_ref: pl.DataFrame

        # Validate that exactly one of reference or gtf is provided
        if reference and gtf:
            console.print(
                "[red]✗[/red] Cannot specify both --reference and --gtf options"
            )
            sys.exit(1)
        elif not reference and not gtf:
            console.print(
                "[red]✗[/red] Must specify either --reference or --gtf option"
            )
            sys.exit(1)

        # Pre-load reference data BEFORE progress bar to handle any download prompts
        exon_ref: pl.DataFrame
        if reference:
            # Use built-in reference
            if reference in BUILTIN_REFERENCES:
                console.print(f"[cyan]Checking reference '{reference}'...")
                exon_ref_result = load_reference(reference)
                if exon_ref_result is None or not isinstance(
                    exon_ref_result, pl.DataFrame
                ):
                    console.print(
                        f"[red]✗[/red] Failed to load reference '{reference}'"
                    )
                    sys.exit(1)
                exon_ref = exon_ref_result
                console.print(f"[green]✓[/green] Reference '{reference}' ready")
            else:
                console.print(f"[red]✗[/red] Unknown built-in reference: {reference}")
                console.print(
                    f"Available references: {list(BUILTIN_REFERENCES.keys())}"
                )
                sys.exit(1)
        else:
            # Use custom GTF file
            console.print(f"[cyan]Loading GTF file '{gtf}'...")
            if gtf is None:
                console.print("[red]✗[/red] GTF file path is required")
                sys.exit(1)
            exon_ref = load_gtf(gtf)
            console.print("[green]✓[/green] GTF file loaded")

        # Create progress display
        with Progress(
            SpinnerColumn(),
            TextColumn("[progress.description]{task.description:<40}"),
            BarColumn(bar_width=40),
            TextColumn("[progress.percentage]{task.percentage:>3.0f}%"),
            TimeElapsedColumn(),
            console=console,
            transient=False,
        ) as progress:
            # Add main analysis task
            task = progress.add_task("[cyan]Running metagene analysis...", total=100)

            # Step 1: Load input sites
            update_progress_description(progress, task, "Loading input sites")
            progress.update(task, completed=20)

            # Convert 1-based to 0-based indices for meta_columns
            meta_col_index = [col - 1 for col in meta_columns]
            weight_col_index = [col - 1 for col in weight_columns]

            input_df = load_sites(
                input_file,
                with_header=with_header,
                meta_col_index=meta_col_index,
                separator=separator,
            )
            progress.console.log(f"[green]✓[/green] Loaded {len(input_df)} input sites")

            # Step 2: Reference data already loaded
            update_progress_description(progress, task, "Using reference data")
            progress.update(task, completed=40)
            progress.console.log(
                f"[green]✓[/green] Using reference with {len(exon_ref)} records"
            )

            # Step 3: Annotate with transcripts
            update_progress_description(progress, task, "Annotating transcripts")
            progress.update(task, completed=50)

            annotated_df = map_to_transcripts(input_df, exon_ref)
            progress.console.log("[green]✓[/green] Annotated transcripts")

            # Step 4: Normalize positions
            update_progress_description(progress, task, "Normalizing positions")
            progress.update(task, completed=70)

            if output_score or output_figure:
                gene_bins, gene_stats, gene_splits = normalize_positions(
                    annotated_df,
                    split_strategy="median",
                    bin_number=bins,
                    weight_col_index=weight_col_index,
                )
                progress.console.log(
                    f"[green]✓[/green] Normalized {gene_stats['5UTR'] + gene_stats['CDS'] + gene_stats['3UTR']} positions"
                )
                progress.console.log(
                    f"Gene splits - 5'UTR: {gene_splits[0]:.3f}, CDS: {gene_splits[1]:.3f}, 3'UTR: {gene_splits[2]:.3f}"
                )

            # Step 5: Save results
            update_progress_description(progress, task, "Saving results")
            progress.update(task, completed=90)

            # Save annotated data
            if output_file:
                annotated_df.write_csv(output_file, separator=separator)
                progress.console.log(
                    f"[green]✓[/green] Saved annotated intervals to: {output_file}"
                )

            # Save score statistics (if requested)
            if output_score:
                # insert (insert_column) a new feature_type column after the first column
                gene_bins.insert_column(
                    0,
                    pl.when(pl.col("feature_midpoint") < gene_splits[0])
                    .then(pl.lit("5UTR"))
                    .when(
                        (pl.col("feature_midpoint") > gene_splits[0] + gene_splits[1])
                    )
                    .then(pl.lit("3UTR"))
                    .otherwise(pl.lit("CDS"))
                    .alias("feature_type"),
                ).write_csv(output_score, separator=separator)
                progress.console.log(
                    f"[green]✓[/green] Saved binned statistics to: {output_score}"
                )

            # Step 6: Generate plot
            if output_figure:
                update_progress_description(progress, task, "Generating plot")
                plot_profile(gene_bins, gene_splits, output_figure)
                progress.console.log(f"[green]✓[/green] Saved plot to: {output_figure}")

            # Update progress to complete
            update_progress_description(progress, task, "Analysis complete")
            progress.update(task, completed=100)

            # Print summary
            progress.console.log("\n[bold cyan]Analysis Summary:[/bold cyan]")
            progress.console.log(f"[blue]Input file:[/blue] {input_file}")
            progress.console.log(f"[blue]Output file:[/blue] {output_file}")
            if output_score:
                progress.console.log(f"[blue]Score output file:[/blue] {output_score}")
            if output_figure:
                progress.console.log(
                    f"[blue]Figure output file:[/blue] {output_figure}"
                )
            if reference:
                progress.console.log(f"[blue]Reference:[/blue] {reference}")
            if gtf:
                progress.console.log(f"[blue]GTF file:[/blue] {gtf}")
            progress.console.log(f"[blue]Region:[/blue] {region}")
            progress.console.log(f"[blue]Number of bins:[/blue] {bins}")
            progress.console.log(f"[blue]Has header:[/blue] {with_header}")
            progress.console.log(f"[blue]Meta columns:[/blue] {meta_columns}")
            progress.console.log(f"[blue]Weight columns:[/blue] {weight_columns}")
            if weight_names:
                progress.console.log(f"[blue]Weight names:[/blue] {weight_names}")
            progress.console.log(f"[blue]Score transform:[/blue] {score_transform}")
            progress.console.log(f"[blue]Normalize:[/blue] {normalize}")
            if output_figure:
                progress.console.log("\n[bold cyan]Plot Settings:[/bold cyan]")

    except Exception as e:
        logger.error(f"Error: {e}")
        sys.exit(1)


if __name__ == "__main__":
    cli()
