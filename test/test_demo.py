#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright © 2025 Ye Chang yech1990@gmail.com
# Distributed under terms of the GNU license.
#
# Created: 2025-01-24

"""
Metagene Analysis Demo - Complete Rewrite Using Package Functions

This script demonstrates the complete metagene workflow using only imported functions
from the metagene package. It produces realistic enrichment patterns near stop codons
by using the correct gene annotation and normalization approach.
"""

import os
import sys
import warnings
import polars as pl
from pathlib import Path

# Add the parent directory to the path to import metagene modules
current_dir = os.path.dirname(os.path.abspath(__file__))
project_root = os.path.abspath(os.path.join(current_dir, ".."))
sys.path.insert(0, project_root)

# Import all necessary functions from metagene package
from metagene import (
    load_gtf,
    map_to_transcripts,
    normalize_positions,
    load_sites,
    plot_profile,
)
from metagene.utils import setup_rich_logger

# Set up rich logging for tests
logger = setup_rich_logger("test_demo")
warnings.filterwarnings("ignore")


def main():
    """
    Main function that runs the complete metagene analysis demo.
    """
    logger.info("Starting Metagene Analysis Demo using Package Functions")

    # Define file paths
    current_dir = Path(__file__).parent
    gtf_file = current_dir / "example.gtf.gz"
    sites_file = current_dir / "sites.tsv.gz"
    meta_col_index = [0, 1, 2]

    # Step 1: Load GTF reference data
    logger.info("Step 1: Loading GTF reference data...")
    exon_ref = load_gtf(str(gtf_file))
    logger.info(f"Loaded exon reference with {len(exon_ref)} records")

    # Step 2: Load input sites
    logger.info("Step 2: Loading input sites...")
    input_df = load_sites(
        str(sites_file), with_header=True, meta_col_index=meta_col_index
    )
    logger.info(f"Loaded {len(input_df)} input sites")

    # Step 3: Annotate genes using package function
    logger.info("Step 3: Annotating genes with transcript information...")
    annotated_df = map_to_transcripts(input_df, exon_ref)
    logger.info(
        f"Annotated {len(annotated_df.filter(pl.col('transcript_id').is_not_null()))} sites"
    )

    # Step 4: Normalize feature positions using package function
    logger.info("Step 4: Normalizing feature positions...")
    gene_bins, gene_stats, gene_splits = normalize_positions(annotated_df, split_strategy="median", bin_number=100)
    logger.info(f"Normalized {gene_bins['count'].sum()} positions")
    logger.info("Gene splits:") 
    logger.info(f"  5'UTR: {gene_splits[0]:.3f}")
    logger.info(f"  CDS: {gene_splits[1]:.3f}")
    logger.info(f"  3'UTR: {gene_splits[2]:.3f}")
    logger.info("Gene stats:")
    logger.info(f"  Unknown: {gene_stats['None']}")
    logger.info(f"  5'UTR: {gene_stats['5UTR']}")
    logger.info(f"  CDS: {gene_stats['CDS']}")
    logger.info(f"  3'UTR: {gene_stats['3UTR']}")

    # Step 5: Generate plot
    logger.info("Step 5: Generating metagene profile plot...")
    output_path = Path(__file__).parent / "metagene_demo_package.png"
    plot_profile(gene_bins, gene_splits, str(output_path))
    logger.info(f"Plot saved to: {output_path}")

    # # Step 6: Show summary statistics
    # logger.info("Step 6: Summary statistics...")
    # show_summary_stats(gene_stats)

    logger.info("Metagene Analysis Demo completed successfully!")


if __name__ == "__main__":
    main()
