#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Annotation and Metagene Analysis Functions

import pyranges as pr
import polars as pl
import numpy as np


def map_to_transcripts(
    input_sites: pl.DataFrame, exon_ref: pr.PyRanges
) -> pl.DataFrame:
    """
    Annotate input sites with transcript information using exon reference.
    Returns a Polars DataFrame with transcript mapping.
    """
    # Prepare query with row index for join
    qry = input_sites.loc[:, ["Chromosome", "Start", "End", "Strand"]].assign(
        _tmp_row_index=np.arange(0, input_sites.shape[0])
    )
    annot = exon_ref.join_ranges(qry, suffix="_qry", join_type="inner")

    annot = pl.DataFrame(annot)
    # Add reference columns
    annot = annot.with_columns(
        pl.col("Chromosome").alias("Chromosome_ref"),
        pl.col("Start").alias("Start_ref"),
        pl.col("End").alias("End_ref"),
        pl.col("Strand").alias("Strand_ref"),
    )

    # Calculate new Start/End based on strand
    annot = annot.with_columns(
        [
            pl.when(pl.col("Strand_ref") == "+")
            .then(
                (pl.col("Start_qry") - pl.col("Start_ref")).clip(
                    0, pl.col("End_exon") - pl.col("Start_exon")
                )
                + pl.col("Start_exon")
            )
            .otherwise(
                (pl.col("End_ref") - pl.col("End_qry")).clip(
                    0, pl.col("End_exon") - pl.col("Start_exon")
                )
                + pl.col("Start_exon")
            )
            .alias("transcript_start"),
            pl.when(pl.col("Strand_ref") == "+")
            .then(
                (pl.col("End_qry") - pl.col("Start_ref")).clip(
                    0, pl.col("End_exon") - pl.col("Start_exon")
                )
                + pl.col("Start_exon")
            )
            .otherwise(
                (pl.col("End_ref") - pl.col("End_qry")).clip(
                    0, pl.col("End_exon") - pl.col("Start_exon")
                )
                + pl.col("Start_exon")
            )
            .alias("transcript_end"),
        ]
    )

    # Use Polars groupby.apply to pick best transcript per gene
    def pick_best_transcript(df: pl.DataFrame) -> pl.DataFrame:
        min_level = df["transcript_level"].min()
        best = df.filter(pl.col("transcript_level") == min_level)
        max_length = best["transcript_length"].max()
        best = best.filter(pl.col("transcript_length") == max_length)
        best_transcript_id = best["transcript_id"][0]
        # If multiple, pick the first
        return best.filter(pl.col("transcript_id") == best_transcript_id)

    annot = annot.group_by("gene_id").map_groups(pick_best_transcript)

    annotation_cols = [
        "gene_id",
        "transcript_id",
        "transcript_start",
        "transcript_end",
        "transcript_length",
        "start_codon_pos",
        "stop_codon_pos",
        "exon_number",
        "Start_exon",
        "End_exon",
    ]
    annot = annot.select(["_tmp_row_index"] + annotation_cols)

    # Join annotation back to input_sites
    annotated_sites = (
        pl.DataFrame(input_sites)
        .with_row_index("_tmp_row_index")
        .drop(annotation_cols, strict=False)
        .join(annot, on="_tmp_row_index", how="left")
        .with_columns(
            pl.col("transcript_start").cast(pl.Int64),
            pl.col("transcript_end").cast(pl.Int64),
            pl.col("transcript_length").cast(pl.Int64),
            pl.col("start_codon_pos").cast(pl.Int64),
            pl.col("stop_codon_pos").cast(pl.Int64),
            pl.col("Start_exon").cast(pl.Int64),
            pl.col("End_exon").cast(pl.Int64),
        )
        .with_columns(record_id=pl.col("_tmp_row_index"))
        .drop("_tmp_row_index")
    )
    return annotated_sites


def calcualte_gene_splits(
    annotated_sites: pr.PyRanges, strategy: str = "mean"
) -> tuple:
    """
    Calculate gene region splits (5'UTR, CDS, 3'UTR) from annotated sites.
    """
    df = (
        annotated_sites.select(
            "transcript_id", "transcript_length", "start_codon_pos", "stop_codon_pos"
        )
        .drop_nulls()
        .unique()
    )
    if strategy == "mean":
        len_5utr = df["start_codon_pos"].mean()
        len_cds = df["stop_codon_pos"].mean() - df["start_codon_pos"].mean()
        len_3utr = df["transcript_length"].mean() - df["stop_codon_pos"].mean()
    elif strategy == "median":
        len_5utr = df["start_codon_pos"].median()
        len_cds = df["stop_codon_pos"].median() - df["start_codon_pos"].median()
        len_3utr = df["transcript_length"].median() - df["stop_codon_pos"].median()
    len_total = len_5utr + len_cds + len_3utr
    return len_5utr / len_total, len_cds / len_total, len_3utr / len_total


def normalize_positions(
    annotated_sites: pl.DataFrame, split_strategy: str = "median", bin_number: int = 100, weight_col_index: list[int] | None = None
) -> tuple[pl.DataFrame, dict, tuple]:
    """
    Normalize transcript positions to relative feature positions (0-1 scale).
    Returns the normalized DataFrame and the gene splits.
    """
    # check if the "transcript_id", "transcript_start" and  "transcript_end" in the dataframe columns
    # use the mid point of transcript_start and transcript_end as transcript_pos

    gene_splits = calcualte_gene_splits(annotated_sites, split_strategy)

    gene_stats = (
        annotated_sites.with_columns(
            transcript_pos=(pl.col("transcript_start") + pl.col("transcript_end")) // 2
        )
        .with_columns(
            feature_type=pl.when(pl.col("transcript_pos").is_null())
            .then(pl.lit("None"))
            .when(pl.col("transcript_pos") < pl.col("start_codon_pos"))
            .then(pl.lit("5UTR"))
            .when(pl.col("transcript_pos") > pl.col("stop_codon_pos"))
            .then(pl.lit("3UTR"))
            .otherwise(pl.lit("CDS"))
        )
        .with_columns(feature_weight=1 / pl.len().over("record_id"))
        .group_by("feature_type")
        .agg(count=pl.col("feature_weight").sum())
    )
    gene_stats = dict( zip(gene_stats["feature_type"], gene_stats["count"]))

    gene_bins = (
        annotated_sites.with_columns(
            transcript_pos=(pl.col("transcript_start") + pl.col("transcript_end")) // 2
        )
        .filter(pl.col("transcript_pos").is_not_null())
        .with_columns(
            feature_pos=pl.when(pl.col("transcript_pos") < pl.col("start_codon_pos"))
            .then(pl.col("transcript_pos") / pl.col("start_codon_pos") * gene_splits[0])
            .when(pl.col("transcript_pos") > pl.col("stop_codon_pos"))
            .then(
                gene_splits[0]
                + gene_splits[1]
                + (pl.col("transcript_pos") - pl.col("stop_codon_pos"))
                / (pl.col("transcript_length") - pl.col("stop_codon_pos"))
                * gene_splits[2]
            )
            .otherwise(
                gene_splits[0]
                + (pl.col("transcript_pos") - pl.col("start_codon_pos"))
                / (pl.col("stop_codon_pos") - pl.col("start_codon_pos"))
                * gene_splits[1]
            )
        )
        .with_columns(feature_weight=1 / pl.len().over("record_id"))
        .with_columns(
            feature_bin=pl.col("feature_pos").cut(
                breaks=np.linspace(0, 1, bin_number + 1)
            )
        )
    )
    n2c = {}
    if weight_col_index is None:
        bin_counts, _ = np.histogram(
            gene_bins["feature_pos"],
            weights=gene_bins["feature_weight"],
            bins=np.linspace(0, 1, bin_number + 1),
        )
        n2c["count"] = bin_counts
    else:
        for col_index in weight_col_index:
            col_name = annotated_sites.columns[col_index]
            bin_counts, _ = np.histogram(
                gene_bins["feature_pos"],
                weights=gene_bins["feature_weight"] * gene_bins[col_name],
                bins=np.linspace(0, 1, bin_number + 1),
            )
            n2c[f"count_{col_name}"] = bin_counts
    bin_midpoints = np.linspace(0, 1, bin_number) + 0.5 / bin_number
    gene_bins = pl.DataFrame({
        "feature_midpoint": bin_midpoints,
        **n2c,
    })
    return gene_bins, gene_stats, gene_splits


def show_summary_stats(df: pl.DataFrame) -> str:
    """
    Generate summary statistics of the analysis.
    
    Args:
        df_normalized: Final DataFrame with all annotations
        
    Returns:
        A formatted string containing the summary statistics
    """
    # filter record with feature_type is not null, and show the proportion passed the filter
    total_passed = df.height
    total_sites = df.height
    pass_percentage = (total_passed / total_sites * 100) if total_sites > 0 else 0

    # Count by feature type
    feature_counts = df.group_by("feature_type").len().sort("feature_type")
    
    # Build feature distribution string
    feature_dist = []
    for row in feature_counts.iter_rows():
        feature_type, count = row
        percentage = (count / total_passed * 100) if total_passed > 0 else 0
        feature_dist.append(f"{feature_type}: {count} sites ({percentage:.1f}%)")

    # Calculate position statistics
    feature_positions = df["feature_pos"]
    pos_stats = []
    if len(feature_positions) > 0:
        pos_stats = [
            f"Mean: {feature_positions.mean():.3f}",
            f"Median: {feature_positions.median():.3f}",
            f"Min: {feature_positions.min():.3f}",
            f"Max: {feature_positions.max():.3f}"
        ]
    else:
        pos_stats = ["No valid position statistics available (all values are null)"]

    # Combine all parts into a single string
    summary = (
        f"Total sites passed the filter: {total_passed} / {total_sites} ({pass_percentage:.1f}%)\n\n"
        f"Feature Distribution:\n  " + "\n  ".join(feature_dist) + "\n\n"
        f"Position Statistics:\n  " + "\n  ".join(pos_stats)
    )
    
    return summary
