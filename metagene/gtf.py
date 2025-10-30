#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright © 2023 Ye Chang yech1990@gmail.com
# Distributed under terms of the GNU license.
#
# Created: 2023-03-08 20:46

import logging
import os

import numpy as np
import polars as pl
import ruranges
from rich.console import Console

from .utils import NewlineRichHandler, ensure_dir, get_cache_dir, get_file_hash

# Set up rich console
console = Console()

# Set up logger with rich console style
logger = logging.getLogger(__name__)
logger.handlers = []  # Remove any existing handlers
logger.addHandler(NewlineRichHandler(console=console))
logger.setLevel(logging.INFO)
logger.propagate = False  # Prevent propagation to root logger


def prepare_exon_ref(gtf_file: str) -> pl.DataFrame:
    """
    Prepares a comprehensive exon reference from a GTF file.
    Uses Polars for initial GTF parsing for speed.
    Returns a Polars DataFrame with exon reference data.
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
        logger.error(
            f"[red]✗[/red] Error reading GTF with Polars: {e}. Consider checking GTF format."
        )
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
    
    # Sort by chromosome, transcript_id, and position for proper cumsum calculation
    pl_df_exon = pl_df_exon.sort(["Chromosome", "gene_id", "transcript_id", "Start"])
    
    # Generate cumulative exon lengths using ruranges
    # Group by gene_id and transcript_id for cumsum
    def apply_group_cumsum(group_df: pl.DataFrame) -> pl.DataFrame:
        starts = group_df["Start"].cast(pl.Int64).to_numpy()
        ends = group_df["End"].cast(pl.Int64).to_numpy()
        negative_strand = (group_df["Strand"] == "-").to_numpy()
        
        # Call ruranges.group_cumsum
        idx, cumsum_start, cumsum_end = ruranges.group_cumsum(
            starts=starts,
            ends=ends,
            negative_strand=negative_strand,
            groups=None,  # All in same group within this function
            sort=True,
        )
        
        # Reorder based on idx to get back to original order
        result_df = group_df.with_row_index("_orig_idx")
        sort_indices = np.argsort(idx)
        result_df = result_df.with_columns([
            pl.Series("Start_exon", cumsum_start[sort_indices], dtype=pl.Int64),
            pl.Series("End_exon", cumsum_end[sort_indices], dtype=pl.Int64),
        ])
        return result_df.drop("_orig_idx")
    
    pl_df_exon = pl_df_exon.group_by(["gene_id", "transcript_id"], maintain_order=True).map_groups(
        apply_group_cumsum
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

    # Join with exon reference and calculate positions using ruranges
    # Prepare arrays for overlap detection
    exon_starts = pl_df_exon["Start"].cast(pl.Int64).to_numpy()
    exon_ends = pl_df_exon["End"].cast(pl.Int64).to_numpy()
    exon_chroms = pl_df_exon["Chromosome"].to_numpy()
    
    codon_starts = pl_df_codons["Start"].cast(pl.Int64).to_numpy()
    codon_ends = pl_df_codons["End"].cast(pl.Int64).to_numpy()
    codon_chroms = pl_df_codons["Chromosome"].to_numpy()
    
    # Create group IDs for chromosomes using vectorized operations
    unique_chroms = np.unique(np.concatenate([exon_chroms, codon_chroms]))
    chrom_to_id = {chrom: i for i, chrom in enumerate(unique_chroms)}
    exon_groups = np.vectorize(chrom_to_id.get)(exon_chroms).astype(np.uint32)
    codon_groups = np.vectorize(chrom_to_id.get)(codon_chroms).astype(np.uint32)
    
    # Find overlaps
    idx_exon, idx_codon = ruranges.overlaps(
        starts=exon_starts,
        ends=exon_ends,
        starts2=codon_starts,
        ends2=codon_ends,
        groups=exon_groups,
        groups2=codon_groups,
    )
    
    # Build overlapping pairs dataframe
    overlaps_df = pl.DataFrame({
        "exon_idx": idx_exon,
        "codon_idx": idx_codon,
    })
    
    # Join with original dataframes to get transcript_id and match by transcript_id
    exon_indexed = pl_df_exon.with_row_index("exon_idx")
    codon_indexed = pl_df_codons.with_row_index("codon_idx")
    
    df_codons = (
        overlaps_df
        .join(exon_indexed, on="exon_idx")
        .join(codon_indexed, on="codon_idx", suffix="_codon")
        .filter(pl.col("transcript_id") == pl.col("transcript_id_codon"))
        .with_columns(
            pl.when(pl.col("Strand") == "+")
            .then(pl.col("Start_codon") - pl.col("Start") + pl.col("Start_exon"))
            .otherwise(pl.col("End") - pl.col("End_codon") + pl.col("Start_exon"))
            .alias("codon_pos")
        )
        .select(["transcript_id", "Feature", "codon_pos"])
    )
    
    # Pivot to get start_codon_pos and stop_codon_pos columns
    pos_codon_df = (
        df_codons
        .pivot(
            index="transcript_id",
            on="Feature",
            values="codon_pos",
            aggregate_function="first",
        )
    )
    
    # Rename columns to match expected format
    if "start_codon" in pos_codon_df.columns:
        pos_codon_df = pos_codon_df.rename({"start_codon": "start_codon_pos"})
    if "stop_codon" in pos_codon_df.columns:
        pos_codon_df = pos_codon_df.rename({"stop_codon": "stop_codon_pos"})
    
    # Merge with exon reference
    pl_tx = pl_df_exon.join(pos_codon_df, on="transcript_id", how="left")
    
    # Convert to Int32 to reduce size
    pl_tx = pl_tx.with_columns([
        pl.col("Start_exon").cast(pl.Int32),
        pl.col("End_exon").cast(pl.Int32),
        pl.col("start_codon_pos").cast(pl.Int32),
        pl.col("stop_codon_pos").cast(pl.Int32),
    ])
    
    return pl_tx


def load_gtf(gtf_file: str, use_cache: bool = True) -> pl.DataFrame:
    """
    Load and process a GTF file to a Polars DataFrame with caching support.

    Args:
        gtf_file: Path to the GTF file
        use_cache: Whether to use caching for faster loading

    Returns:
        Polars DataFrame with processed exon information
    """
    if not use_cache:
        logger.info(f"[cyan]Cache disabled. Processing GTF file: {gtf_file}[/cyan]")
        df_result = prepare_exon_ref(gtf_file)
        return df_result

    gtf_file_abs = os.path.abspath(gtf_file)
    gtf_dir = os.path.dirname(gtf_file_abs)
    gtf_basename = os.path.basename(gtf_file_abs)

    local_cache_filename = f".{gtf_basename}.parquet"
    local_cache_filepath = os.path.join(gtf_dir, local_cache_filename)

    if os.path.exists(local_cache_filepath) and os.access(
        local_cache_filepath, os.R_OK
    ):
        # Check if cache file is newer than the GTF file
        gtf_mtime = os.path.getmtime(gtf_file_abs)
        cache_mtime = os.path.getmtime(local_cache_filepath)

        if cache_mtime > gtf_mtime:
            try:
                logger.info(
                    f"[cyan]Loading local cached reference from: {local_cache_filepath}[/cyan]"
                )
                df = pl.read_parquet(local_cache_filepath)
                return df
            except Exception as e:
                logger.warning(
                    f"[yellow]⚠[/yellow] Error loading local cache file {local_cache_filepath}: {e}. Attempting default cache."
                )
        else:
            logger.info(
                f"[yellow]⚠[/yellow] Local cache file is older than GTF file. Checking default cache for: {gtf_file}"
            )
            # Remove old local cache file
            try:
                os.remove(local_cache_filepath)
            except Exception as e:
                logger.warning(
                    f"[yellow]⚠[/yellow] Could not remove old local cache file: {e}"
                )

    try:
        cache_dir = get_cache_dir()
    except OSError as e:
        logger.error(
            f"[red]✗[/red] Could not create cache directory: {e}. Processing without caching."
        )
        df_result = prepare_exon_ref(gtf_file)
        return df_result

    path_hash = get_file_hash(gtf_file_abs)
    default_cache_filename = f"{gtf_basename}_{path_hash}.parquet"
    default_cache_filepath = os.path.join(cache_dir, default_cache_filename)

    if os.path.exists(default_cache_filepath) and os.access(
        default_cache_filepath, os.R_OK
    ):
        # Check if cache file is newer than the GTF file
        gtf_mtime = os.path.getmtime(gtf_file_abs)
        cache_mtime = os.path.getmtime(default_cache_filepath)

        if cache_mtime > gtf_mtime:
            try:
                logger.info(
                    f"[cyan]Loading default cached reference from: {default_cache_filepath}[/cyan]"
                )
                df = pl.read_parquet(default_cache_filepath)
                return df
            except Exception as e:
                logger.warning(
                    f"[yellow]⚠[/yellow] Error loading default cache file {default_cache_filepath}: {e}. Processing without caching."
                )
        else:
            logger.info(
                f"[yellow]⚠[/yellow] Cache file is older than GTF file. Rebuilding cache for: {gtf_file}"
            )
            # Remove old cache file
            try:
                os.remove(default_cache_filepath)
            except Exception as e:
                logger.warning(
                    f"[yellow]⚠[/yellow] Could not remove old cache file: {e}"
                )

    # If we get here, we need to process the GTF file
    logger.info(f"[cyan]Processing GTF file: {gtf_file}[/cyan]")
    df_result = prepare_exon_ref(gtf_file)

    # Try to save to local cache first (same directory as GTF file)
    local_cache_saved = False
    try:
        logger.info(
            f"[cyan]Attempting to save processed reference to local cache: {local_cache_filepath}[/cyan]"
        )
        df_result.write_parquet(local_cache_filepath)
        local_cache_saved = True
        logger.info(
            f"[green]✓[/green] Successfully saved to local cache: {local_cache_filepath}"
        )
    except Exception as e:
        logger.warning(
            f"[yellow]⚠[/yellow] Could not save to local cache {local_cache_filepath}: {e}. Trying default cache directory."
        )

    # If local cache failed, try default cache directory
    if not local_cache_saved:
        try:
            logger.info(
                f"[cyan]Saving processed reference to default cache: {default_cache_filepath}[/cyan]"
            )
            df_result.write_parquet(default_cache_filepath)
            logger.info(
                f"[green]✓[/green] Successfully saved to default cache: {default_cache_filepath}"
            )
        except Exception as e:
            logger.warning(
                f"[yellow]⚠[/yellow] Could not save to default cache: {e}. Proceeding without caching."
            )

    return df_result


if __name__ == "__main__":
    # Example usage:
    test_gtf_path = os.path.join(
        os.path.dirname(__file__), "..", "test", "example.gtf.gz"
    )  # Adjusted path

    # Create a dummy test.gtf if it doesn't exist for testing purposes
    if not os.path.exists(test_gtf_path):
        logger.info(
            f"[yellow]Test GTF file not found at {test_gtf_path}. Creating a dummy file for demonstration.[/yellow]"
        )
        dummy_gtf_content = """chr1\tunknown\texon\t1000\t2000\t.\t+\t.\tgene_id "test_gene"; transcript_id "test_transcript"; exon_number "1";
chr1\tunknown\tstart_codon\t1200\t1202\t.\t+\t.\ttranscript_id "test_transcript";
chr1\tunknown\tstop_codon\t1800\t1802\t.\t+\t.\ttranscript_id "test_transcript";
"""
        ensure_dir(os.path.dirname(test_gtf_path))

        with open(test_gtf_path, "wt") as f:  # Save as gzipped file
            f.write(dummy_gtf_content)
        logger.info(f"[green]✓[/green] Dummy GTF file created at {test_gtf_path}")

    logger.info(f"[cyan]Processing GTF file with caching: {test_gtf_path}[/cyan]")
    df_result = load_gtf(test_gtf_path)
    logger.info(
        "[green]✓[/green] Processing complete. Resulting Polars DataFrame (first call):"
    )
    if df_result is not None and len(df_result) > 0:
        logger.info(str(df_result.head()))
    else:
        logger.warning("[yellow]⚠[/yellow] No data or empty dataframe returned.")

    logger.info(
        f"[cyan]Processing GTF file with caching again (should use cache): {test_gtf_path}[/cyan]"
    )
    df_result_cached = load_gtf(test_gtf_path)
    logger.info("[green]✓[/green] Resulting Polars DataFrame (second call):")
    if df_result_cached is not None and len(df_result_cached) > 0:
        logger.info(str(df_result_cached.head()))
    else:
        logger.warning(
            "[yellow]⚠[/yellow] No data or empty dataframe returned on second call."
        )

    # Test with cache disabled
    logger.info(
        f"[cyan]Processing GTF file with cache disabled: {test_gtf_path}[/cyan]"
    )
    df_result_no_cache = load_gtf(test_gtf_path, use_cache=False)
    logger.info("[green]✓[/green] Resulting Polars DataFrame (cache disabled):")
    if df_result_no_cache is not None and len(df_result_no_cache) > 0:
        logger.info(str(df_result_no_cache.head()))
    else:
        logger.warning(
            "[yellow]⚠[/yellow] No data or empty dataframe returned with cache disabled."
        )
