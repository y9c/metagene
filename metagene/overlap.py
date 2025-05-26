#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright © 2023 Ye Chang yech1990@gmail.com
# Distributed under terms of the GNU license.
#
# Feature Overlap Module - Handles genomic overlap analysis with PyRanges v1

import sys
import numpy as np
import pandas as pd
import pyranges as pr


def calculate_bin_statistics(
    data: np.ndarray, 
    weights: np.ndarray, 
    num_bins: int = 100, 
    bin_range: tuple[float, float] = (0, 1), 
    suffix: str = ""
) -> pd.DataFrame:
    """
    Calculate binned statistics for metagene analysis.
    
    Args:
        data: Array of normalized positions
        weights: Array of weights for each position
        num_bins: Number of bins for analysis
        bin_range: Range for binning (default 0-1)
        suffix: Suffix for column names
    
    Returns:
        DataFrame with bin statistics
    """
    data = np.asarray(data)
    
    # Generate the bin edges
    bins = np.linspace(bin_range[0], bin_range[1], num_bins + 1)

    # Initialize arrays for sums and counts
    bin_sums = np.zeros(num_bins)
    bin_counts = np.zeros(num_bins)

    # Calculate which bin each data point falls into
    bin_indices = np.digitize(data, bins) - 1

    # Update the bin_sums and bin_counts arrays
    for w, b in zip(weights, bin_indices):
        if 0 <= b < num_bins:
            bin_sums[b] += w
            bin_counts[b] += 1

    # Calculate the mean for each bin
    bin_means = np.divide(
        bin_sums,
        bin_counts,
        out=np.zeros_like(bin_sums),
        where=bin_counts != 0,
    )

    df = pd.DataFrame(
        {"count": bin_counts, "sum": bin_sums, "mean": bin_means},
        index=(bins[1:] + bins[:-1]) / 2,
    )
    
    if suffix:
        df.columns = [c + suffix for c in df.columns]
    return df


def annotate_with_features(
    df_input: pd.DataFrame,
    df_feature: pd.DataFrame,
    bin_number: int = 100,
    type_ratios: list[float] | None = None,
    annot_name: bool = False,
    keep_all: bool = False,
    by_strand: bool = False,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Annotate input genomic intervals with feature information using PyRanges v1.
    
    Args:
        df_input: Input genomic intervals
        df_feature: Feature annotations 
        bin_number: Number of bins for analysis
        type_ratios: Custom ratios for feature types
        annot_name: Whether to include feature names
        keep_all: Whether to keep all columns
        by_strand: Whether to consider strand
    
    Returns:
        Tuple of (annotated dataframe, binned statistics)
    """
    # Perform overlap using PyRanges v1 join_ranges
    df_joined = (
        pr.PyRanges(df_input)
        .join_ranges(
            pr.PyRanges(df_feature),
            report_overlap_column="Overlap",
            strand_behavior="same" if by_strand else "ignore",
        )
    )

    # Get the best overlap for each input interval
    df = (
        df_joined.sort_values(by=["Overlap"], ascending=False)
        .groupby(["Chromosome", "Start", "End"], as_index=False, group_keys=False)
        .head(1)
        .assign(
            d=lambda x: np.where(
                x.Strand_b == "+",  # Fixed: use _b suffix for PyRanges v1
                (np.minimum((x.Start + x.End) // 2, x.End_b) - x.Start_b) 
                / x.len_of_feature
                + x.frac_of_feature,
                (x.End_b - np.maximum((x.Start + x.End) // 2, x.Start_b))
                / x.len_of_feature
                + x.frac_of_feature,
            )
        )
    )

    # Calculate feature type ratios
    if type_ratios is not None:
        type_ratios_array = np.array(type_ratios)
        type2ratio = dict(
            zip(["5UTR", "CDS", "3UTR"], type_ratios_array / np.sum(type_ratios_array))
        )
    else:
        # Calculate ratios from overlapping transcripts
        s = (
            df_feature[df_feature["Transcript"].isin(df["Transcript"].unique())]
            .groupby(["Transcript", "Type"])["len_of_window"]
            .sum()
            .reset_index()
            .groupby(["Type"])["len_of_window"]
            .median()
        )
        type2ratio = (s / s.sum()).to_dict()
        
        # Check that all required feature types are present
        required_types = {"5UTR", "CDS", "3UTR"}
        if not required_types.issubset(type2ratio.keys()):
            missing_types = required_types - set(type2ratio.keys())
            sys.exit(
                f"Error: Given ranges do not overlap all 3 features: {missing_types}! "
                "Please use the type_ratios argument instead."
            )

    # Normalize distances based on feature types
    df["d_norm"] = np.where(
        df["Type"] == "5UTR",
        df["d"] * type2ratio["5UTR"],
        np.where(
            df["Type"] == "CDS",
            df["d"] * type2ratio["CDS"] + type2ratio["5UTR"],
            df["d"] * type2ratio["3UTR"] + type2ratio["5UTR"] + type2ratio["CDS"],
        ),
    )
    
    if annot_name:
        df["Name"] = df["Name_b"]  # Fixed: use _b suffix

    # Calculate bin statistics for each weight column
    weight_cols = [c for c in df.columns if c.startswith("Weight_")]
    df_list = []
    
    if len(weight_cols) > 0:
        for c in weight_cols:
            df_list.append(
                calculate_bin_statistics(
                    df["d_norm"],
                    weights=df[c].fillna(0),
                    num_bins=bin_number,
                    bin_range=(0, 1),
                    suffix="_" + c.replace("Weight_", ""),
                )
            )
    else:
        df_list.append(
            calculate_bin_statistics(
                df["d_norm"],
                weights=np.ones(len(df["d_norm"])),
                num_bins=bin_number,
                bin_range=(0, 1),
            )
        )
    
    df_score = pd.concat(df_list, axis=1)

    # Select output columns
    if not keep_all:
        df = df.loc[:, ["Chromosome", "Start", "End", "Name", "d_norm", "Strand"]]
    
    # Store metadata in dataframe attributes
    df.attrs.update(type2ratio.items())
    df_score.attrs.update(type2ratio.items())
    
    return df, df_score
