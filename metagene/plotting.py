#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright © 2023 Ye Chang yech1990@gmail.com
# Distributed under terms of the GNU license.
#
# Plotting Module - Visualization functions for metagene analysis

import matplotlib.pyplot as plt
import polars as pl
from .utils import setup_logger

# Set up logger
logger = setup_logger(__name__)


def plot_profile(
    gene_bins: pl.DataFrame,
    gene_splits: tuple[float, float, float],
    output_path: str,
    figsize: tuple[int, int] = (10, 5),
):
    """
    Create the metagene profile plot using matplotlib.

    Args:
        gene_bins: DataFrame with normalized positions
        gene_splits: Tuple of gene region ratios
    """
    # Create the plot
    import matplotlib

    matplotlib.use("Agg")

    fig, ax = plt.subplots(figsize=figsize)

    # Plot histogram of feature positions
    for col in gene_bins.columns:
        if col.startswith("count"):
            ax.plot(
                gene_bins["feature_midpoint"],
                gene_bins[col],
                alpha=1,
                linewidth=0.5,
                label=col.removeprefix("count_"),
            )
            if len(gene_bins.columns) == 2:
                ax.fill_between(
                    gene_bins["feature_midpoint"], 0, gene_bins[col], alpha=0.5
                )

    # Add vertical lines to mark gene regions
    ax.axvline(gene_splits[0], color="red", linestyle="dashed", linewidth=2)
    ax.axvline(
        gene_splits[0] + gene_splits[1], color="green", linestyle="dashed", linewidth=2
    )

    # Add region labels
    ax.text(
        gene_splits[0] / 2,
        ax.get_ylim()[1] * 0.9,
        "5'UTR",
        ha="center",
        va="top",
        fontsize=12,
        fontweight="bold",
    )
    ax.text(
        gene_splits[0] + gene_splits[1] / 2,
        ax.get_ylim()[1] * 0.9,
        "CDS",
        ha="center",
        va="top",
        fontsize=12,
        fontweight="bold",
    )
    ax.text(
        gene_splits[0] + gene_splits[1] + gene_splits[2] / 2,
        ax.get_ylim()[1] * 0.9,
        "3'UTR",
        ha="center",
        va="top",
        fontsize=12,
        fontweight="bold",
    )

    # Set x-axis limits to 0-1
    ax.set_xlim(0, 1)

    # Set labels and title
    ax.set_xlabel("Normalized Gene Position", fontsize=12)
    ax.set_ylabel("Density", fontsize=12)
    ax.set_title("Metagene Profile", fontsize=14, fontweight="bold")
    ax.legend()
    ax.grid(True, alpha=0.3)

    plt.tight_layout()

    # Save plot
    plt.savefig(output_path, dpi=300, bbox_inches="tight")

    # Close the plot instead of showing it
    plt.close()
