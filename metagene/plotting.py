#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright © 2023 Ye Chang yech1990@gmail.com
# Distributed under terms of the GNU license.
#
# Plotting Module - Visualization functions for metagene analysis

import matplotlib.pyplot as plt
import polars as pl
import numpy as np
import pandas as pd
import os
from typing import Tuple, Optional, Union, List
from matplotlib.figure import Figure
from matplotlib.axes import Axes
from .utils import setup_logger

# Set up logger
logger = setup_logger(__name__)


def plot_profile(
    feature_pos_series: Union[pl.Series, pd.Series, np.ndarray],
    gene_splits: Tuple[float, float, float],  # (prop_5utr, prop_cds, prop_3utr)
    title: str = "Distribution of Sites Across Normalized Transcript Features",
    output_plot_path: Optional[str] = None,
    figsize: Tuple[int, int] = (10, 5),
    bins: int = 100,
) -> None:
    """
    Create a metagene profile plot showing distribution across transcript features.

    Args:
        feature_pos_series: Series of normalized feature positions (0-1 scale)
        gene_splits: Tuple of (5'UTR proportion, CDS proportion, 3'UTR proportion)
        title: Plot title
        output_plot_path: Path to save plot (optional)
        figsize: Figure size (width, height)
        bins: Number of histogram bins
    """
    # Convert to numpy array if needed
    if hasattr(feature_pos_series, "to_numpy"):
        if hasattr(feature_pos_series, "drop_nulls"):  # Polars Series
            data = feature_pos_series.drop_nulls().to_numpy()
        else:  # Pandas Series
            data = feature_pos_series.dropna().to_numpy()
    else:  # Already numpy array
        data = feature_pos_series[~np.isnan(feature_pos_series)]

    # Create figure and axis
    fig, ax = plt.subplots(figsize=figsize)

    # Create histogram
    n, bins_edges, patches = ax.hist(
        data,
        bins=bins,
        color="blue",
        alpha=0.7,
        edgecolor="black",
        linewidth=0.5,
        density=True,
    )

    # Add vertical lines for feature boundaries
    line_5utr_cds = gene_splits[0]
    line_cds_3utr = gene_splits[0] + gene_splits[1]

    ax.axvline(
        line_5utr_cds,
        color="red",
        linestyle="dashed",
        linewidth=2,
        label=f"5'UTR-CDS boundary ({line_5utr_cds:.2f})",
    )
    ax.axvline(
        line_cds_3utr,
        color="green",
        linestyle="dashed",
        linewidth=2,
        label=f"CDS-3'UTR boundary ({line_cds_3utr:.2f})",
    )

    # Customize plot
    ax.set_xlabel("Normalized Position Across Transcript")
    ax.set_ylabel("Density")
    ax.set_title(title)
    ax.legend()
    ax.grid(True, alpha=0.3)

    # Add feature region labels
    ax.text(
        line_5utr_cds / 2,
        ax.get_ylim()[1] * 0.9,
        "5'UTR",
        ha="center",
        va="center",
        fontweight="bold",
    )
    ax.text(
        (line_5utr_cds + line_cds_3utr) / 2,
        ax.get_ylim()[1] * 0.9,
        "CDS",
        ha="center",
        va="center",
        fontweight="bold",
    )
    ax.text(
        (line_cds_3utr + 1) / 2,
        ax.get_ylim()[1] * 0.9,
        "3'UTR",
        ha="center",
        va="center",
        fontweight="bold",
    )

    # Save or show plot
    if output_plot_path:
        plt.tight_layout()
        plt.savefig(output_plot_path, dpi=300, bbox_inches="tight")
        logger.info(f"Plot saved to: {output_plot_path}")
    else:
        plt.show()

    plt.close()


def plot_multiple_metagene_profiles(
    data_dict: dict,
    gene_splits: Tuple[float, float, float],
    title: str = "Metagene Profiles Comparison",
    output_plot_path: Optional[str] = None,
    figsize: Tuple[int, int] = (12, 8),
    bins: int = 100,
) -> None:
    """
    Create multiple metagene profile plots for comparison.

    Args:
        data_dict: Dictionary mapping sample names to position series
        gene_splits: Tuple of (5'UTR proportion, CDS proportion, 3'UTR proportion)
        title: Plot title
        output_plot_path: Path to save plot (optional)
        figsize: Figure size (width, height)
        bins: Number of histogram bins
    """
    # Create figure and axis
    fig, ax = plt.subplots(figsize=figsize)

    colors = plt.cm.Set1(np.linspace(0, 1, len(data_dict)))

    for i, (sample_name, feature_pos_series) in enumerate(data_dict.items()):
        # Convert to numpy array if needed
        if hasattr(feature_pos_series, "to_numpy"):
            if hasattr(feature_pos_series, "drop_nulls"):  # Polars Series
                data = feature_pos_series.drop_nulls().to_numpy()
            else:  # Pandas Series
                data = feature_pos_series.dropna().to_numpy()
        else:  # Already numpy array
            data = feature_pos_series[~np.isnan(feature_pos_series)]

        # Create histogram
        ax.hist(
            data,
            bins=bins,
            alpha=0.6,
            color=colors[i],
            label=sample_name,
            density=True,
            histtype="step",
            linewidth=2,
        )

    # Add vertical lines for feature boundaries
    line_5utr_cds = gene_splits[0]
    line_cds_3utr = gene_splits[0] + gene_splits[1]

    ax.axvline(
        line_5utr_cds,
        color="red",
        linestyle="dashed",
        linewidth=2,
        label=f"5'UTR-CDS boundary ({line_5utr_cds:.2f})",
    )
    ax.axvline(
        line_cds_3utr,
        color="green",
        linestyle="dashed",
        linewidth=2,
        label=f"CDS-3'UTR boundary ({line_cds_3utr:.2f})",
    )

    # Customize plot
    ax.set_xlabel("Normalized Position Across Transcript")
    ax.set_ylabel("Density")
    ax.set_title(title)
    ax.legend(bbox_to_anchor=(1.05, 1), loc="upper left")
    ax.grid(True, alpha=0.3)

    # Add feature region labels
    ax.text(
        line_5utr_cds / 2,
        ax.get_ylim()[1] * 0.9,
        "5'UTR",
        ha="center",
        va="center",
        fontweight="bold",
    )
    ax.text(
        (line_5utr_cds + line_cds_3utr) / 2,
        ax.get_ylim()[1] * 0.9,
        "CDS",
        ha="center",
        va="center",
        fontweight="bold",
    )
    ax.text(
        (line_cds_3utr + 1) / 2,
        ax.get_ylim()[1] * 0.9,
        "3'UTR",
        ha="center",
        va="center",
        fontweight="bold",
    )

    # Save or show plot
    if output_plot_path:
        plt.tight_layout()
        plt.savefig(output_plot_path, dpi=300, bbox_inches="tight")
        logger.info(f"Plot saved to: {output_plot_path}")
    else:
        plt.show()

    plt.close()


def plot_binned_statistics(
    df_score: pd.DataFrame,
    gene_splits: Tuple[float, float, float],
    title: str = "Metagene Profile - Binned Statistics",
    output_plot_path: Optional[str] = None,
    figsize: Tuple[int, int] = (10, 6),
) -> None:
    """
    Plot binned statistics from metagene analysis.

    Args:
        df_score: DataFrame with binned statistics
        gene_splits: Tuple of (5'UTR proportion, CDS proportion, 3'UTR proportion)
        title: Plot title
        output_plot_path: Path to save plot (optional)
        figsize: Figure size (width, height)
    """
    fig, ax = plt.subplots(figsize=figsize)

    # Plot mean values for each weight column
    mean_cols = [
        col for col in df_score.columns if col.endswith("_mean") or col == "mean"
    ]

    for col in mean_cols:
        label = col.replace("_mean", "") if col != "mean" else "default"
        ax.plot(df_score.index, df_score[col], label=label, linewidth=2)

    # Add vertical lines for feature boundaries
    line_5utr_cds = gene_splits[0]
    line_cds_3utr = gene_splits[0] + gene_splits[1]

    ax.axvline(
        line_5utr_cds,
        color="red",
        linestyle="dashed",
        linewidth=2,
        alpha=0.7,
        label=f"5'UTR-CDS boundary ({line_5utr_cds:.2f})",
    )
    ax.axvline(
        line_cds_3utr,
        color="green",
        linestyle="dashed",
        linewidth=2,
        alpha=0.7,
        label=f"CDS-3'UTR boundary ({line_cds_3utr:.2f})",
    )

    # Customize plot
    ax.set_xlabel("Normalized Position Across Transcript")
    ax.set_ylabel("Mean Score")
    ax.set_title(title)
    ax.legend()
    ax.grid(True, alpha=0.3)

    # Add feature region labels
    ax.text(
        line_5utr_cds / 2,
        ax.get_ylim()[1] * 0.9,
        "5'UTR",
        ha="center",
        va="center",
        fontweight="bold",
    )
    ax.text(
        (line_5utr_cds + line_cds_3utr) / 2,
        ax.get_ylim()[1] * 0.9,
        "CDS",
        ha="center",
        va="center",
        fontweight="bold",
    )
    ax.text(
        (line_cds_3utr + 1) / 2,
        ax.get_ylim()[1] * 0.9,
        "3'UTR",
        ha="center",
        va="center",
        fontweight="bold",
    )

    # Save or show plot
    if output_plot_path:
        plt.tight_layout()
        plt.savefig(output_plot_path, dpi=300, bbox_inches="tight")
        logger.info(f"Plot saved to: {output_plot_path}")
    else:
        plt.show()

    plt.close()


def simple_metagene_plot(gene_bins: pl.DataFrame, gene_splits: tuple, output_path: str):
    """
    Create the metagene profile plot using matplotlib.

    Args:
        gene_bins: DataFrame with normalized positions
        gene_splits: Tuple of gene region ratios
    """
    # Create the plot
    fig, ax = plt.subplots(figsize=(12, 6))

    # Plot histogram of feature positions
    for col in gene_bins.columns:
        if col.startswith("count"):
            ax.plot(gene_bins["feature_midpoint"], gene_bins[col], alpha=1, linewidth=0.5, label=col.removeprefix("count_"))
            if len(gene_bins.columns) == 2:
                ax.fill_between(gene_bins["feature_midpoint"], 0, gene_bins[col], alpha=0.5)

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
