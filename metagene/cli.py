#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright © 2023 Ye Chang yech1990@gmail.com
# Distributed under terms of the GNU license.
#
# Command-line interface for metagene analysis

import click
import os
import sys
from pathlib import Path
from typing import Optional, List, Dict, Any
from rich.progress import Progress, SpinnerColumn, TextColumn, BarColumn, TimeElapsedColumn
from rich.console import Console
from rich.logging import RichHandler
import logging



from .io import load_sites, load_reference
from .annotation import map_to_transcripts, normalize_positions, show_summary_stats
from .gtf import load_gtf
from .download import download_references, list_references as show_references, download_reference_cli
from .config import BUILTIN_REFERENCES
from .plotting import simple_metagene_plot

# Set up rich console
console = Console()

class NewlineRichHandler(RichHandler):
    def __init__(
        self,
        console=console,
        rich_tracebacks=True,
        markup=True,
        show_time=True,
        show_path=False,
        show_level=True,
        enable_link_path=False,
        **kwargs
    ):
        super().__init__(
            console=console,
            rich_tracebacks=rich_tracebacks,
            markup=markup,
            show_time=show_time,
            show_path=show_path,
            show_level=show_level,
            enable_link_path=enable_link_path,
            **kwargs
        )

    def emit(self, record):
        try:
            message = self.format(record)
            self.console.print("\n" + message + "\n")
        except Exception:
            self.handleError(record)

# Set up logging with custom Rich handler
logging.basicConfig(
    level=logging.INFO,
    format="%(message)s",
    datefmt="[%X]",
    handlers=[NewlineRichHandler()]
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
        return [int(x.strip()) for x in value.split(',')]
    except ValueError:
        raise click.BadParameter('Must be comma-separated integers')


def parse_comma_separated_strings(ctx, param, value):
    """Parse comma-separated strings."""
    if not value:
        return []
    return [x.strip() for x in value.split(',')]


def parse_comma_separated_floats(ctx, param, value):
    """Parse comma-separated floats."""
    if not value:
        return None
    try:
        return [float(x.strip()) for x in value.split(',')]
    except ValueError:
        raise click.BadParameter('Must be comma-separated floats')


def check_disk_space(path, required_mb):
    """Check if there's enough disk space available."""
    import shutil
    total, used, free = shutil.disk_usage(path)
    free_mb = free // (1024 * 1024)
    return free_mb >= required_mb


def print_analysis_summary(params):
    """Print a summary of analysis parameters."""
    logger.info("Analysis Parameters:")
    logger.info("------------------")
    logger.info(f"Input file: {params['input_file']}")
    logger.info(f"Reference: {params['builtin_features'] or params['features']}")
    logger.info(f"Meta columns: {params['meta_columns']}")
    logger.info(f"Weight columns: {params['weight_columns']}")
    logger.info(f"Bin number: {params['bin_number']}")
    logger.info(f"Output files:")
    if params['output_file']:
        logger.info(f"  - Annotated intervals: {params['output_file']}")
    if params['output_score']:
        logger.info(f"  - Binned statistics: {params['output_score']}")
    if params['output_figure']:
        logger.info(f"  - Plot: {params['output_figure']}")
    logger.info("------------------")


def update_progress_description(progress: Any, task: Any, description: str) -> None:
    """Update progress bar description."""
    progress.update(task, description=f"[cyan]{description}...")


@click.command(
    help="Run metagene analysis on genomic sites.",
    context_settings=dict(help_option_names=['-h', '--help'])
)
@click.version_option()
@click.option(
    "--input",
    "-i",
    "input_file",
    type=click.Path(exists=True),
    help="Input file path (BED, GTF, or TSV format)",
)
@click.option(
    "--output",
    "-o",
    "output_file",
    type=click.Path(),
    help="Output file path (TSV format)",
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
    is_flag=True,
    help="Input file has header line",
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
    default="5",
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
    "--plot",
    is_flag=True,
    help="Generate plots",
)
@click.option(
    "--plot-dir",
    type=click.Path(),
    help="Directory to save plots (default: same as output file)",
)
@click.option(
    "--plot-format",
    type=click.Choice(["png", "pdf", "svg"]),
    default="png",
    help="Plot file format (default: png)",
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
    output_score: Optional[str],
    output_figure: Optional[str],
    reference: Optional[str],
    gtf: Optional[str],
    region: str,
    bins: int,
    with_header: bool,
    meta_columns: List[int],
    weight_columns: List[int],
    weight_names: List[str],
    score_transform: str,
    normalize: bool,
    plot: bool,
    plot_dir: Optional[str],
    plot_format: str,
    list_references_flag: bool,
    download: Optional[str],
) -> None:
    """Run metagene analysis on genomic sites."""
    # Check if no arguments provided (except defaults) and show help
    ctx = click.get_current_context()
    if (not list_references_flag and not download and not input_file and not output_file and 
        not reference and not gtf):
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
        download_reference_cli(download, console)
        return

    # Validate required options for analysis
    if not input_file:
        console.print("[red]✗[/red] Input file is required for analysis (use -i/--input)")
        console.print("Use --list to see available references or --download to download references")
        sys.exit(1)
    
    if not output_file:
        console.print("[red]✗[/red] Output file is required for analysis (use -o/--output)")
        sys.exit(1)

    try:
        # Validate that exactly one of reference or gtf is provided
        if reference and gtf:
            console.print("[red]✗[/red] Cannot specify both --reference and --gtf options")
            sys.exit(1)
        elif not reference and not gtf:
            console.print("[red]✗[/red] Must specify either --reference or --gtf option")
            sys.exit(1)
        
        # Pre-load reference data BEFORE progress bar to handle any download prompts
        if reference:
            # Use built-in reference
            if reference in BUILTIN_REFERENCES:
                console.print(f"[cyan]Checking reference '{reference}'...")
                exon_ref = load_reference(reference)
                if exon_ref is None:
                    console.print(f"[red]✗[/red] Failed to load reference '{reference}'")
                    sys.exit(1)
                console.print(f"[green]✓[/green] Reference '{reference}' ready")
            else:
                console.print(f"[red]✗[/red] Unknown built-in reference: {reference}")
                console.print(f"Available references: {list(BUILTIN_REFERENCES.keys())}")
                sys.exit(1)
        else:
            # Use custom GTF file
            console.print(f"[cyan]Loading GTF file '{gtf}'...")
            if gtf is None:
                console.print("[red]✗[/red] GTF file path is required")
                sys.exit(1)
            exon_ref = load_gtf(gtf)
            console.print(f"[green]✓[/green] GTF file loaded")
        
        # Create progress display
        with Progress(
            SpinnerColumn(),
            TextColumn("[progress.description]{task.description:<40}"),
            BarColumn(bar_width=40),
            TextColumn("[progress.percentage]{task.percentage:>3.0f}%"),
            TimeElapsedColumn(),
            console=console,
            transient=False
        ) as progress:
            # Add main analysis task
            task = progress.add_task("[cyan]Running metagene analysis...", total=100)

            # Step 1: Load input sites
            update_progress_description(progress, task, "Loading input sites")
            progress.update(task, completed=20)
            
            # Convert 1-based to 0-based indices for meta_columns
            meta_col_index = [col - 1 for col in meta_columns]
            
            input_df = load_sites(
                input_file, 
                with_header=with_header, 
                meta_col_index=meta_col_index
            )
            progress.console.log(f"[green]✓[/green] Loaded {len(input_df)} input sites")

            # Step 2: Reference data already loaded
            update_progress_description(progress, task, "Using reference data")
            progress.update(task, completed=40)
            progress.console.log(f"[green]✓[/green] Using reference with {len(exon_ref)} records")

            # Step 3: Annotate with transcripts
            update_progress_description(progress, task, "Annotating transcripts")
            progress.update(task, completed=50)
            
            annotated_df = map_to_transcripts(input_df, exon_ref)
            progress.console.log(f"[green]✓[/green] Annotated transcripts")

            # Step 4: Normalize positions
            update_progress_description(progress, task, "Normalizing positions")
            progress.update(task, completed=70)
            
            final_df, gene_splits = normalize_positions(annotated_df, strategy="median")
            progress.console.log(f"[green]✓[/green] Normalized {len(final_df)} positions")
            progress.console.log(f"Gene splits - 5'UTR: {gene_splits[0]:.3f}, CDS: {gene_splits[1]:.3f}, 3'UTR: {gene_splits[2]:.3f}")

            # Step 5: Save results
            update_progress_description(progress, task, "Saving results")
            progress.update(task, completed=90)
            
            # Save annotated data
            if output_file:
                final_df.write_csv(output_file)
                progress.console.log(f"[green]✓[/green] Saved annotated intervals to: {output_file}")
            
            # Save score statistics (if requested)
            if output_score:
                # Create simple binned statistics
                import polars as pl
                import numpy as np
                
                # Create bins
                bin_edges = np.linspace(0, 1, bins + 1)
                bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
                
                # Count sites in each bin
                counts, _ = np.histogram(final_df["feature_pos"], bins=bin_edges)
                
                # Create score DataFrame
                score_df = pl.DataFrame({
                    "bin_center": bin_centers,
                    "count": counts
                })
                score_df.write_csv(output_score)
                progress.console.log(f"[green]✓[/green] Saved binned statistics to: {output_score}")

            # Step 6: Generate plot
            if output_figure:
                update_progress_description(progress, task, "Generating plot")
                simple_metagene_plot(final_df, gene_splits, output_figure)
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
                progress.console.log(f"[blue]Figure output file:[/blue] {output_figure}")
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
                if plot_dir:
                    progress.console.log(f"[blue]Plot directory:[/blue] {plot_dir}")
                progress.console.log(f"[blue]Plot format:[/blue] {plot_format}")

    except Exception as e:
        logger.error(f"Error: {e}")
        sys.exit(1)


def print_analysis_summary(
    input_file: str,
    output_file: str,
    output_score: Optional[str],
    output_figure: Optional[str],
    reference: Optional[str],
    gtf: Optional[str],
    region: str,
    bins: int,
    with_header: bool,
    meta_columns: List[int],
    weight_columns: List[int],
    weight_names: List[str],
    score_transform: str,
    normalize: bool,
    plot: bool,
    plot_dir: Optional[str],
    plot_format: str,
) -> None:
    """Print summary of analysis parameters."""
    logger.info("\nAnalysis Summary:")
    logger.info(f"Input file: {input_file}")
    logger.info(f"Output file: {output_file}")
    if output_score:
        logger.info(f"Score output file: {output_score}")
    if output_figure:
        logger.info(f"Figure output file: {output_figure}")
    if reference:
        logger.info(f"Reference: {reference}")
    if gtf:
        logger.info(f"GTF file: {gtf}")
    logger.info(f"Region: {region}")
    logger.info(f"Number of bins: {bins}")
    logger.info(f"Has header: {with_header}")
    logger.info(f"Meta columns: {meta_columns}")
    logger.info(f"Weight columns: {weight_columns}")
    if weight_names:
        logger.info(f"Weight names: {weight_names}")
    logger.info(f"Score transform: {score_transform}")
    logger.info(f"Normalize: {normalize}")
    if plot:
        logger.info("\nPlot Settings:")
        if plot_dir:
            logger.info(f"Plot directory: {plot_dir}")
        logger.info(f"Plot format: {plot_format}")


if __name__ == "__main__":
    cli()
