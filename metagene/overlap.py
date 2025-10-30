#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright © 2023 Ye Chang yech1990@gmail.com
# Distributed under terms of the GNU license.
#
# Feature Overlap Module - Handles genomic overlap analysis with ruranges

import sys
import numpy as np
import polars as pl
import ruranges


def calculate_bin_statistics(
    data: np.ndarray, 
    weights: np.ndarray, 
    num_bins: int = 100, 
    bin_range: tuple[float, float] = (0, 1), 
    suffix: str = ""
) -> pl.DataFrame:
    """
    Calculate binned statistics for metagene analysis.
    
    Args:
        data: Array of normalized positions
        weights: Array of weights for each position
        num_bins: Number of bins for analysis
        bin_range: Range for binning (default 0-1)
        suffix: Suffix for column names
    
    Returns:
        Polars DataFrame with bin statistics
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

    # Create column names
    col_names = {
        "count": "count",
        "sum": "sum",
        "mean": "mean"
    }
    if suffix:
        col_names = {k: v + suffix for k, v in col_names.items()}

    df = pl.DataFrame({
        col_names["count"]: bin_counts,
        col_names["sum"]: bin_sums,
        col_names["mean"]: bin_means,
    }).with_columns(
        pl.Series("index", (bins[1:] + bins[:-1]) / 2)
    )
    
    return df


def annotate_with_features(
    df_input: pl.DataFrame,
    df_feature: pl.DataFrame,
    bin_number: int = 100,
    type_ratios: list[float] | None = None,
    annot_name: bool = False,
    keep_all: bool = False,
    by_strand: bool = False,
) -> tuple[pl.DataFrame, pl.DataFrame]:
    """
    Annotate input genomic intervals with feature information using ruranges.
    
    Args:
        df_input: Input genomic intervals (Polars DataFrame)
        df_feature: Feature annotations (Polars DataFrame)
        bin_number: Number of bins for analysis
        type_ratios: Custom ratios for feature types
        annot_name: Whether to include feature names
        keep_all: Whether to keep all columns
        by_strand: Whether to consider strand
    
    Returns:
        Tuple of (annotated polars dataframe, binned statistics polars dataframe)
    """
    # Prepare arrays for overlap detection
    input_starts = df_input["Start"].cast(pl.Int64).to_numpy()
    input_ends = df_input["End"].cast(pl.Int64).to_numpy()
    input_chroms = df_input["Chromosome"].to_numpy()
    
    feature_starts = df_feature["Start"].cast(pl.Int64).to_numpy()
    feature_ends = df_feature["End"].cast(pl.Int64).to_numpy()
    feature_chroms = df_feature["Chromosome"].to_numpy()
    
    # Create group IDs for chromosomes (and optionally strand)
    if by_strand:
        input_strands = df_input["Strand"].to_numpy()
        feature_strands = df_feature["Strand"].to_numpy()
        
        unique_chr_strand = np.unique(
            np.concatenate([
                np.char.add(input_chroms.astype(str), input_strands.astype(str)),
                np.char.add(feature_chroms.astype(str), feature_strands.astype(str))
            ])
        )
        chr_strand_to_id = {cs: i for i, cs in enumerate(unique_chr_strand)}
        
        input_groups = np.array([
            chr_strand_to_id[c + s] 
            for c, s in zip(input_chroms.astype(str), input_strands.astype(str))
        ], dtype=np.uint32)
        
        feature_groups = np.array([
            chr_strand_to_id[c + s] 
            for c, s in zip(feature_chroms.astype(str), feature_strands.astype(str))
        ], dtype=np.uint32)
    else:
        unique_chroms = np.unique(np.concatenate([input_chroms, feature_chroms]))
        chrom_to_id = {chrom: i for i, chrom in enumerate(unique_chroms)}
        input_groups = np.array([chrom_to_id[c] for c in input_chroms], dtype=np.uint32)
        feature_groups = np.array([chrom_to_id[c] for c in feature_chroms], dtype=np.uint32)
    
    # Find overlaps
    idx_input, idx_feature = ruranges.overlaps(
        starts=input_starts,
        ends=input_ends,
        starts2=feature_starts,
        ends2=feature_ends,
        groups=input_groups,
        groups2=feature_groups,
    )
    
    # Calculate overlap lengths
    overlap_starts = np.maximum(input_starts[idx_input], feature_starts[idx_feature])
    overlap_ends = np.minimum(input_ends[idx_input], feature_ends[idx_feature])
    overlap_lengths = overlap_ends - overlap_starts
    
    # Build overlapping pairs dataframe
    overlaps_df = pl.DataFrame({
        "input_idx": idx_input,
        "feature_idx": idx_feature,
        "Overlap": overlap_lengths,
    })
    
    # Join with original dataframes
    input_indexed = df_input.with_row_index("input_idx")
    feature_indexed = df_feature.with_row_index("feature_idx")
    
    df_joined = (
        overlaps_df
        .join(input_indexed, on="input_idx")
        .join(feature_indexed, on="feature_idx", suffix="_b")
    )
    
    # Get the best overlap for each input interval (sort and take first per group)
    df = (
        df_joined
        .sort("Overlap", descending=True)
        .group_by(["Chromosome", "Start", "End"], maintain_order=True)
        .first()
    )
    
    # Calculate position within feature
    df = df.with_columns(
        pl.when(pl.col("Strand_b") == "+")
        .then(
            (
                (pl.col("Start") + pl.col("End")) // 2
            ).clip(pl.col("Start_b"), pl.col("End_b"))
            - pl.col("Start_b")
        )
        .otherwise(
            pl.col("End_b")
            - (
                (pl.col("Start") + pl.col("End")) // 2
            ).clip(pl.col("Start_b"), pl.col("End_b"))
        )
        .alias("_midpoint_offset")
    ).with_columns(
        (
            pl.col("_midpoint_offset") / pl.col("len_of_feature")
            + pl.col("frac_of_feature")
        ).alias("d")
    )

    # Calculate feature type ratios
    if type_ratios is not None:
        type_ratios_array = np.array(type_ratios)
        type2ratio = dict(
            zip(["5UTR", "CDS", "3UTR"], type_ratios_array / np.sum(type_ratios_array))
        )
    else:
        # Calculate ratios from overlapping transcripts
        matching_features = df_feature.filter(
            pl.col("Transcript").is_in(df["Transcript"].unique())
        )
        
        s = (
            matching_features
            .group_by(["Transcript", "Type"])
            .agg(pl.col("len_of_window").sum())
            .group_by("Type")
            .agg(pl.col("len_of_window").median())
        )
        
        total = s["len_of_window"].sum()
        type2ratio = dict(zip(s["Type"], s["len_of_window"] / total))
        
        # Check that all required feature types are present
        required_types = {"5UTR", "CDS", "3UTR"}
        if not required_types.issubset(type2ratio.keys()):
            missing_types = required_types - set(type2ratio.keys())
            sys.exit(
                f"Error: Given ranges do not overlap all 3 features: {missing_types}! "
                "Please use the type_ratios argument instead."
            )

    # Normalize distances based on feature types
    df = df.with_columns(
        pl.when(pl.col("Type") == "5UTR")
        .then(pl.col("d") * type2ratio["5UTR"])
        .when(pl.col("Type") == "CDS")
        .then(pl.col("d") * type2ratio["CDS"] + type2ratio["5UTR"])
        .otherwise(
            pl.col("d") * type2ratio["3UTR"] + type2ratio["5UTR"] + type2ratio["CDS"]
        )
        .alias("d_norm")
    )
    
    if annot_name:
        df = df.with_columns(pl.col("Name_b").alias("Name"))

    # Calculate bin statistics for each weight column
    weight_cols = [c for c in df.columns if c.startswith("Weight_")]
    df_list = []
    
    if len(weight_cols) > 0:
        for c in weight_cols:
            df_list.append(
                calculate_bin_statistics(
                    df["d_norm"].to_numpy(),
                    weights=df[c].fill_null(0).to_numpy(),
                    num_bins=bin_number,
                    bin_range=(0, 1),
                    suffix="_" + c.replace("Weight_", ""),
                )
            )
    else:
        df_list.append(
            calculate_bin_statistics(
                df["d_norm"].to_numpy(),
                weights=np.ones(len(df)),
                num_bins=bin_number,
                bin_range=(0, 1),
            )
        )
    
    # Concatenate horizontally (column-wise)
    df_score = pl.concat(df_list, how="horizontal")

    # Select output columns
    if not keep_all:
        df = df.select(["Chromosome", "Start", "End", "Name", "d_norm", "Strand"])
    
    # Store metadata in dataframe attributes (Polars doesn't support attrs, but we can return it separately)
    # For backward compatibility, we could store it as metadata in the schema, but simpler to just document it
    
    return df, df_score
