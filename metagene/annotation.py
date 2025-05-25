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
            .alias("Start"),
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
            .alias("End"),
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

    # Add transcript_pos column
    annot = annot.with_columns([pl.col("End").alias("transcript_pos")])

    annotation_cols = [
        "gene_id",
        "transcript_id",
        "transcript_pos",
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
            [
                pl.col("transcript_pos").cast(pl.Int64),
                pl.col("transcript_length").cast(pl.Int64),
                pl.col("start_codon_pos").cast(pl.Int64),
                pl.col("stop_codon_pos").cast(pl.Int64),
                pl.col("Start_exon").cast(pl.Int64),
                pl.col("End_exon").cast(pl.Int64),
            ]
        )
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
        .unique()
        .drop_nulls()
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
    annotated_sites: pl.DataFrame, strategy: str = "median"
) -> tuple[pl.DataFrame, tuple]:
    """
    Normalize transcript positions to relative feature positions (0-1 scale).
    Returns the normalized DataFrame and the gene splits.
    """
    gene_splits = calcualte_gene_splits(annotated_sites, strategy)
    normalized_sites = annotated_sites.with_columns(
        feature_type=pl.when(pl.col("transcript_pos").is_null())
        .then(pl.lit(None))
        .when(pl.col("transcript_pos") < pl.col("start_codon_pos"))
        .then(pl.lit("5UTR"))
        .when(pl.col("transcript_pos") > pl.col("stop_codon_pos"))
        .then(pl.lit("3UTR"))
        .otherwise(pl.lit("CDS"))
    ).with_columns(
        feature_pos=pl.when(pl.col("transcript_pos").is_null())
        .then(pl.lit(None))
        .when(pl.col("transcript_pos") < pl.col("start_codon_pos"))
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
    return normalized_sites, gene_splits


def show_summary_stats(df_normalized: pl.DataFrame):
    import logging

    """
    Show summary statistics of the analysis.
    
    Args:
        df: Final DataFrame with all annotations
    """
    # filter record with feature_type is not null, and show the proportion passed the filter
    # of each feature_type (xx / xx (..%) passed)
    df = df_normalized.filter(pl.col("feature_pos").is_not_null())
    # logging.info(f"Total sites passed the filter: {df.height}")
    logging.info(
        f"Total sites passed the filter: {df.height} / {df_normalized.height} ({df.height / df_normalized.height * 100:.1f}%)"
    )

    # Count by feature type
    feature_counts = df.group_by("feature_type").len().sort("feature_type")
    total_sites = df.height

    logging.info("\nFeature Distribution:")
    for row in feature_counts.iter_rows():
        feature_type, count = row
        percentage = (count / total_sites) * 100
        logging.info(f"  {feature_type}: {count} sites ({percentage:.1f}%)")

    # Show position statistics - handle null values properly
    feature_positions = df["feature_pos"]
    if len(feature_positions) > 0:  # Only calculate stats if we have non-null values
        logging.info("Position Statistics:")
        logging.info(f"  Mean: {feature_positions.mean():.3f}")
        logging.info(f"  Median: {feature_positions.median():.3f}")
        logging.info(f"  Min: {feature_positions.min():.3f}")
        logging.info(f"  Max: {feature_positions.max():.3f}")
    else:
        logging.info("No valid position statistics available (all values are null)")
