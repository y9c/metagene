#!/usr/bin/env python3
"""
Test script for demonstrating the use of built-in GRCh38 reference.
This script shows how to:
1. Load the built-in GRCh38 reference
2. Analyze sample data with the built-in reference
3. Generate and visualize metagene profiles
"""

from pathlib import Path
import logging
import polars as pl
from metagene.io import load_sites, load_reference
from metagene.annotation import (
    map_to_transcripts,
    normalize_positions,
    show_summary_stats,
)
from metagene.plotting import simple_metagene_plot

# Force reconfigure logging with timestamps after imports
logging.getLogger().handlers.clear()  # Clear existing handlers
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
    force=True,  # Force reconfiguration
)


def test_buildin_grch38():
    logging.info("=== Starting Metagene Analysis Demo with Built-in GRCh38 ===")

    # Define file paths
    current_dir = Path(__file__).parent
    gtf_file = current_dir / "example.gtf.gz"
    sites_file = current_dir / "sites.tsv.gz"
    meta_col_index = [0, 1, 2]
    weight_col_index = [4]

    # Step 1: Load built-in reference (GRCh38)
    logging.info("Step 1: Loading built-in reference (GRCh38)...")
    exon_ref = load_reference("GRCh38")

    # Step 2: Load input sites
    logging.info("Step 2: Loading input sites...")
    input_df = load_sites(
        str(sites_file), with_header=True, meta_col_index=meta_col_index
    )
    logging.info(f"Loaded {len(input_df)} input sites")

    # Step 3: Annotate genes using package function
    logging.info("Step 3: Annotating genes with transcript information...")
    annotated_df = map_to_transcripts(input_df, exon_ref)
    logging.info(
        f"Annotated {len(annotated_df.filter(pl.col('transcript_id').is_not_null()))} sites"
    )

    # Step 4: Normalize feature positions using package function
    logging.info("Step 4: Normalizing feature positions...")
    gene_bins, gene_stats, gene_splits = normalize_positions(
        annotated_df, split_strategy="median", bin_number=200, weight_col_index=None
    )
    logging.info(
        f"Normalized {gene_stats['5UTR'] + gene_stats['CDS'] + gene_stats['3UTR']} positions"
    )
    logging.info(f"Gene splits:")
    logging.info(f"  5'UTR: {gene_splits[0]:.3f}")
    logging.info(f"  CDS: {gene_splits[1]:.3f}")
    logging.info(f"  3'UTR: {gene_splits[2]:.3f}")
    logging.info(f"Gene stats:")
    logging.info(f"  Unknown: {gene_stats['None']}")
    logging.info(f"  5'UTR: {gene_stats['5UTR']}")
    logging.info(f"  CDS: {gene_stats['CDS']}")
    logging.info(f"  3'UTR: {gene_stats['3UTR']}")

    # Step 5: Generate plot
    logging.info("Step 5: Generating metagene profile plot...")
    output_path = Path(__file__).parent / "metagene_demo_package.png"
    simple_metagene_plot(gene_bins, gene_splits, output_path)
    logging.info(f"Plot saved to: {output_path}")

    # # Step 6: Show summary statistics
    # logging.info("Step 6: Summary statistics...")
    # show_summary_stats(gene_stats)

    logging.info("Metagene Analysis Demo completed successfully!")


if __name__ == "__main__":
    test_buildin_grch38()
