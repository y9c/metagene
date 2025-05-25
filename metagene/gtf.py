#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright © 2023 Ye Chang yech1990@gmail.com
# Distributed under terms of the GNU license.
#
# Created: 2023-03-08 20:46

import numpy as np
import pyranges as pr
import os
import hashlib
import polars as pl
import pandas as pd
from .utils import get_cache_dir, setup_logger, get_file_hash, ensure_dir

# Set up logger
logger = setup_logger(__name__)


def prepare_exon_ref(gtf_file: str) -> pr.PyRanges:
    """
    Prepares a comprehensive exon reference from a GTF file.
    Uses Polars for initial GTF parsing for speed.
    """

    try:
        pl_df = pl.read_csv(
            gtf_file,
            separator="\t",
            has_header=False,
            comment_prefix="#",
            schema={
                "Chromosome": pl.Utf8,
                "source": pl.Utf8,
                "Feature": pl.Utf8,
                "Start": pl.Int64,
                "End": pl.Int64,
                "score": pl.Utf8,
                "Strand": pl.Utf8,
                "frame": pl.Utf8,
                "attribute": pl.Utf8,
            },
        )
    except Exception as e:
        logger.error(f"Error reading GTF with Polars: {e}. Consider checking GTF format.")
        raise

    # Filter out non-exon features
    pl_df = pl_df.filter(pl.col("Feature").is_in(["exon", "start_codon", "stop_codon"]))

    # Add row ID for join operations and split the attribute column
    # Only keep the relevant attributes, to avoid unnecessary data
    pl_df = (
        pl_df.with_row_index("_temp_row_id")
        .with_columns(
            pl.col("attribute")
            .str.extract_all(r'(\w+)\s+"([^"]+)"')
            .list.eval(
                pl.element().str.extract_groups(
                    r'(?<attribute_key>\w+)\s+"(?<attribute_value>[^"]+)"'
                )
            )
        )
        .explode("attribute")
        .unnest("attribute")
        .filter(
            pl.col("attribute_key").is_in(
                [
                    "gene_id",
                    "transcript_id",
                    "exon_number",
                    "transcript_support_level",
                    "tag",
                ]
            )
        )
        .pivot(
            on="attribute_key",
            values="attribute_value",
            aggregate_function="first",
        )
        .drop("_temp_row_id")
    )

    # GTF: 1-based, closed [Start, End]
    # PyRanges internal: 0-based, half-open [Start, End)
    pl_df = pl_df.with_columns(pl.col("Start") - 1)

    # now parse the exon records only
    pl_df_exon = pl_df.filter(pl.col("Feature") == "exon")

    # drop unnecessary columns
    pl_df_exon = pl_df_exon.drop(["source", "Feature", "score", "frame"])

    pl_df_exon = pl_df_exon.with_columns(
        ((pl.col("End") - pl.col("Start")).sum().alias("transcript_length")).over(
            "gene_id", "transcript_id"
        )
    )

    if "tag" not in pl_df_exon.columns:
        pl_df_exon = pl_df_exon.with_columns(pl.lit(None).cast(pl.Utf8).alias("tag"))
    if "transcript_support_level" not in pl_df_exon.columns:
        pl_df_exon = pl_df_exon.with_columns(
            pl.lit(None).cast(pl.Utf8).alias("transcript_support_level")
        )
    pl_df_exon = pl_df_exon.with_columns(
        pl.when(pl.col("tag") == "Ensembl_canonical")
        .then(0)
        .otherwise(
            pl.col("transcript_support_level")
            .str.split(" ")
            .list.get(0)
            .replace({"NA": "6", "NaN": "6"})
            .cast(pl.Int64, strict=False)
            .fill_null(10)
        )
        .alias("transcript_level")
    ).drop(["tag", "transcript_support_level"])
    # then convert to PyRanges
    pr_exon = pr.PyRanges(pl_df_exon.to_pandas())

    # generate cumulative exon lengths
    pr_exon = pr_exon.group_cumsum(
        match_by=["gene_id", "transcript_id"],
        use_strand=True,
        cumsum_start_column="Start_exon",
        cumsum_end_column="End_exon",
        sort=True,
    )

    # Extract both start and stop codons together
    pl_df_codons = (
        pl_df.filter(pl.col("Feature").is_in(["start_codon", "stop_codon"]))
        .with_columns(
            pl.when(pl.col("Strand") == "+")
            .then(pl.col("Start"))
            .otherwise(pl.col("End") - 1)
            .alias("CodonStart")
        )
        .group_by(["Chromosome", "gene_id", "transcript_id", "Feature"])
        .agg(
            pl.col("Strand").first(),
            pl.when((pl.col("Feature") == "start_codon") == (pl.col("Strand") == "+"))
            .then(pl.col("CodonStart").min())
            .otherwise(pl.col("CodonStart").max())
            .alias("Start"),
        )
        .with_columns(pl.col("Start").list.first())
        .with_columns((pl.col("Start") + 1).alias("End"))
        .select(["Chromosome", "Start", "End", "Strand", "transcript_id", "Feature"])
    )

    pr_codons = pr.PyRanges(pl_df_codons.to_pandas())
    # Join with exon reference and calculate positions
    df_codons = (
        pr_exon.join_ranges(pr_codons, join_type="inner", match_by="transcript_id")
        .assign(
            codon_pos=lambda df: np.where(
                df["Strand"] == "+",
                df["Start_b"] - df["Start"] + df["Start_exon"],
                df["End"] - df["End_b"] + df["Start_exon"],
            ),
        )
        .loc[:, ["transcript_id", "Feature", "codon_pos"]]
    )

    # Pivot to get start_codon_pos and stop_codon_pos columns
    pos_codon_df = df_codons.pivot_table(
        index="transcript_id", columns="Feature", values="codon_pos", aggfunc="first"
    ).reset_index()

    # Rename columns to match expected format
    pos_codon_df.columns.name = None
    if "start_codon" in pos_codon_df.columns:
        pos_codon_df = pos_codon_df.rename(columns={"start_codon": "start_codon_pos"})
    if "stop_codon" in pos_codon_df.columns:
        pos_codon_df = pos_codon_df.rename(columns={"stop_codon": "stop_codon_pos"})

    # Merge with exon reference
    pr_tx = pr_exon.merge(pos_codon_df, on="transcript_id", how="left")
    # Start_exon', 'End_exon', 'start_codon_pos', 'stop_codon_pos' to Int32 to reduce size
    pr_tx["Start_exon"] = pr_tx["Start_exon"].astype("Int32")
    pr_tx["End_exon"] = pr_tx["End_exon"].astype("Int32")
    pr_tx["start_codon_pos"] = pr_tx["start_codon_pos"].astype("Int32")
    pr_tx["stop_codon_pos"] = pr_tx["stop_codon_pos"].astype("Int32")
    return pr_tx


def load_gtf(gtf_file: str, use_cache: bool = True) -> pr.PyRanges:
    """
    Load and process a GTF file to a PyRanges object with caching support.
    
    Args:
        gtf_file: Path to the GTF file
        use_cache: Whether to use caching for faster loading
        
    Returns:
        PyRanges object with processed exon information
    """
    if not use_cache:
        logger.info(f"Cache disabled. Processing GTF file: {gtf_file}")
        return prepare_exon_ref(gtf_file)

    gtf_file_abs = os.path.abspath(gtf_file)
    gtf_dir = os.path.dirname(gtf_file_abs)
    gtf_basename = os.path.basename(gtf_file_abs)

    local_cache_filename = f".{gtf_basename}.parquet"
    local_cache_filepath = os.path.join(gtf_dir, local_cache_filename)

    if os.path.exists(local_cache_filepath) and os.access(
        local_cache_filepath, os.R_OK
    ):
        try:
            logger.info(f"Loading local cached reference from: {local_cache_filepath}")
            # Check if pyranges has read_parquet or use DataFrame
            if hasattr(pr, "read_parquet"):
                return pr.read_parquet(local_cache_filepath)
            else:
                df = pd.read_parquet(local_cache_filepath)
                return pr.PyRanges(df)
        except Exception as e:
            logger.warning(
                f"Error loading local cache file {local_cache_filepath}: {e}. Attempting default cache."
            )

    try:
        cache_dir = get_cache_dir()
    except OSError as e:
        logger.error(f"Could not create cache directory: {e}. Processing without caching.")
        return prepare_exon_ref(gtf_file)

    path_hash = get_file_hash(gtf_file_abs)
    default_cache_filename = f"{gtf_basename}_{path_hash}.parquet"
    default_cache_filepath = os.path.join(cache_dir, default_cache_filename)

    if os.path.exists(default_cache_filepath) and os.access(
        default_cache_filepath, os.R_OK
    ):
        try:
            logger.info(f"Loading default cached reference from: {default_cache_filepath}")
            # Check if pyranges has read_parquet or use DataFrame
            if hasattr(pr, "read_parquet"):
                return pr.read_parquet(default_cache_filepath)
            else:
                df = pd.read_parquet(default_cache_filepath)
                return pr.PyRanges(df)
        except Exception as e:
            logger.warning(
                f"Error loading default cache file {default_cache_filepath}: {e}. Processing without caching."
            )

    # If we get here, we need to process the GTF file
    logger.info(f"Processing GTF file: {gtf_file}")
    pr_tx = prepare_exon_ref(gtf_file)

    # Save to cache
    try:
        logger.info(f"Saving processed reference to cache: {default_cache_filepath}")
        if hasattr(pr_tx, "to_parquet"):
            pr_tx.to_parquet(default_cache_filepath)
        else:
            pr_tx.df.to_parquet(default_cache_filepath)
    except Exception as e:
        logger.warning(f"Could not save to cache: {e}")

    return pr_tx


if __name__ == "__main__":
    # Example usage:
    test_gtf_path = os.path.join(
        os.path.dirname(__file__), "..", "test", "example.gtf.gz"
    )  # Adjusted path

    # Create a dummy test.gtf if it doesn't exist for testing purposes
    if not os.path.exists(test_gtf_path):
        logger.info(
            f"Test GTF file not found at {test_gtf_path}. Creating a dummy file for demonstration."
        )
        dummy_gtf_content = """chr1\tunknown\texon\t1000\t2000\t.\t+\t.\tgene_id "test_gene"; transcript_id "test_transcript"; exon_number "1";
chr1\tunknown\tstart_codon\t1200\t1202\t.\t+\t.\ttranscript_id "test_transcript";
chr1\tunknown\tstop_codon\t1800\t1802\t.\t+\t.\ttranscript_id "test_transcript";
"""
        ensure_dir(os.path.dirname(test_gtf_path))

        with open(test_gtf_path, "wt") as f:  # Save as gzipped file
            f.write(dummy_gtf_content)
        logger.info(f"Dummy GTF file created at {test_gtf_path}")

    logger.info(f"Processing GTF file with caching: {test_gtf_path}")
    df_result = load_gtf(test_gtf_path)
    logger.info("\nProcessing complete. Resulting PyRanges object (first call):")
    if df_result is not None and not df_result.empty:
        logger.info(df_result.head())
    else:
        logger.warning("No data or empty dataframe returned.")

    logger.info(
        f"\nProcessing GTF file with caching again (should use cache): {test_gtf_path}"
    )
    df_result_cached = load_gtf(test_gtf_path)
    logger.info("\nResulting PyRanges object (second call):")
    if df_result_cached is not None and not df_result_cached.empty:
        logger.info(df_result_cached.head())
    else:
        logger.warning("No data or empty dataframe returned on second call.")

    # Test with cache disabled
    logger.info(f"\nProcessing GTF file with cache disabled: {test_gtf_path}")
    df_result_no_cache = load_gtf(test_gtf_path, use_cache=False)
    logger.info("\nResulting PyRanges object (cache disabled):")
    if df_result_no_cache is not None and not df_result_no_cache.empty:
        logger.info(df_result_no_cache.head())
    else:
        logger.warning("No data or empty dataframe returned with cache disabled.")
