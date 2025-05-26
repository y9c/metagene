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
import logging
import warnings
import polars as pl
import matplotlib.pyplot as plt
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
    simple_metagene_plot,
    show_summary_stats,
)

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
)
warnings.filterwarnings("ignore")


def main():
    """
    Main function that runs the complete metagene analysis demo.
    """
    logging.info("Starting Metagene Analysis Demo using Package Functions")

    # Define file paths
    current_dir = Path(__file__).parent
    gtf_file = current_dir / "example.gtf.gz"
    sites_file = current_dir / "sites.tsv.gz"
    meta_col_index = [0, 1, 2]

    # Step 1: Load GTF reference data
    logging.info("Step 1: Loading GTF reference data...")
    exon_ref = load_gtf(str(gtf_file))
    logging.info(f"Loaded exon reference with {len(exon_ref)} records")

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
    gene_bins, gene_stats, gene_splits = normalize_positions(annotated_df, split_strategy="median", bin_number=100)
    logging.info(f"Normalized {gene_bins['count'].sum()} positions")
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
    main()
