#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Map to Local Coordinates - Convert global genomic coordinates to local reference coordinates

This module implements map_to_local functionality using ruranges,
inspired by pyranges1's map_to_local but using Rust-based ruranges operations.
"""

import numpy as np
import polars as pl
import ruranges


def map_to_local(
    query: pl.DataFrame,
    reference: pl.DataFrame,
    ref_id_col: str = "transcript_id",
    match_by: list[str] | str | None = None,
    keep_global_chrom: bool = False,
    keep_global_loc: bool = False,
) -> pl.DataFrame:
    """
    Map genomic intervals from global coordinates to local reference coordinates.

    This function transforms query intervals to local coordinates within reference intervals,
    similar to pyranges1's map_to_local. It uses ruranges for efficient overlap detection.

    Args:
        query: Polars DataFrame with genomic intervals to map (must have Chromosome, Start, End, Strand)
        reference: Polars DataFrame with reference intervals (must have Chromosome, Start, End, Strand, and ref_id_col)
        ref_id_col: Column name in reference that identifies unique references (e.g., "transcript_id", "gene_id")
        match_by: Column name(s) to match between query and reference (optional)
        keep_global_chrom: If True, keep original global chromosome as "Chromosome_global"
        keep_global_loc: If True, keep original global Start/End as "Start_global"/"End_global"

    Returns:
        Polars DataFrame with local coordinates

    Example:
        >>> # Map SNPs to transcript coordinates
        >>> snps = pl.DataFrame({
        ...     "Chromosome": ["chr1", "chr1"],
        ...     "Start": [1000, 2000],
        ...     "End": [1001, 2001],
        ...     "Strand": ["+", "+"],
        ... })
        >>> transcripts = pl.DataFrame({
        ...     "Chromosome": ["chr1"],
        ...     "Start": [900],
        ...     "End": [2100],
        ...     "Strand": ["+"],
        ...     "transcript_id": ["TX1"],
        ... })
        >>> result = map_to_local(snps, transcripts, ref_id_col="transcript_id")
    """
    # Validate required columns
    required_query_cols = ["Chromosome", "Start", "End", "Strand"]
    required_ref_cols = ["Chromosome", "Start", "End", "Strand", ref_id_col]

    for col in required_query_cols:
        if col not in query.columns:
            raise ValueError(f"Query DataFrame must have '{col}' column")

    for col in required_ref_cols:
        if col not in reference.columns:
            raise ValueError(f"Reference DataFrame must have '{col}' column")

    # Store original query columns for reordering
    original_query_cols = query.columns

    # Handle match_by parameter
    if match_by is not None:
        if isinstance(match_by, str):
            match_by = [match_by]
    else:
        match_by = []

    # Calculate cumulative positions
    # For each reference ID, calculate cumulative lengths in 5'->3' order
    # Plus strand: left to right (ascending Start)
    # Minus strand: right to left (descending Start)

    ref_sorted_for_cumsum = reference.sort([ref_id_col, "Start"])

    def calculate_cumsum_positions(group_df: pl.DataFrame) -> pl.DataFrame:
        """Calculate cumulative positions for a group of exons"""
        # Sort by coordinate position
        if group_df["Strand"][0] == "+":
            # Plus strand: sort ascending by Start (5' to 3' is left to right)
            sorted_df = group_df.sort("Start")
        else:
            # Minus strand: sort descending by Start (5' to 3' is right to left)
            sorted_df = group_df.sort("Start", descending=True)

        # Calculate lengths
        lengths = (sorted_df["End"] - sorted_df["Start"]).to_numpy()

        # Calculate cumulative starts and ends
        cumsum_start = np.zeros(len(lengths), dtype=np.int64)
        cumsum_end = np.zeros(len(lengths), dtype=np.int64)

        cumsum_start[0] = 0
        cumsum_end[0] = lengths[0]
        for i in range(1, len(lengths)):
            cumsum_start[i] = cumsum_end[i - 1]
            cumsum_end[i] = cumsum_start[i] + lengths[i]

        # Add cumsum columns
        result_df = sorted_df.with_columns(
            [
                pl.Series("_cumsum_start", cumsum_start, dtype=pl.Int64),
                pl.Series("_cumsum_end", cumsum_end, dtype=pl.Int64),
            ]
        )

        return result_df

    ref_with_cumsum = ref_sorted_for_cumsum.group_by(
        ref_id_col, maintain_order=True
    ).map_groups(calculate_cumsum_positions)

    # Add row indices for joining
    query_indexed = query.with_row_index("_query_idx")
    ref_indexed = ref_with_cumsum.with_row_index("_ref_idx")

    # Prepare arrays for overlap detection
    query_starts = query_indexed["Start"].cast(pl.Int64).to_numpy()
    query_ends = query_indexed["End"].cast(pl.Int64).to_numpy()
    query_chroms = query_indexed["Chromosome"].to_numpy()
    # query_strands not needed for overlap grouping

    ref_starts = ref_indexed["Start"].cast(pl.Int64).to_numpy()
    ref_ends = ref_indexed["End"].cast(pl.Int64).to_numpy()
    ref_chroms = ref_indexed["Chromosome"].to_numpy()
    # ref_strands not needed for overlap grouping

    # Create group IDs - ignore strand for overlap detection (we'll handle strand later)
    unique_chroms = np.unique(np.concatenate([query_chroms, ref_chroms]))
    chrom_to_id = {chrom: i for i, chrom in enumerate(unique_chroms)}
    query_groups = np.vectorize(chrom_to_id.get)(query_chroms).astype(np.uint32)
    ref_groups = np.vectorize(chrom_to_id.get)(ref_chroms).astype(np.uint32)

    # Find overlaps
    idx_query, idx_ref = ruranges.overlaps(
        starts=query_starts,
        ends=query_ends,
        starts2=ref_starts,
        ends2=ref_ends,
        groups=query_groups,
        groups2=ref_groups,
    )

    if len(idx_query) == 0:
        # No overlaps found - return empty result with expected columns
        result_cols = ["Chromosome", "Start", "End", "Strand"]
        if keep_global_chrom:
            result_cols.append("Chromosome_global")
        if keep_global_loc:
            result_cols.extend(["Start_global", "End_global", "Strand_global"])
        return pl.DataFrame({col: [] for col in result_cols})

    # Build overlapping pairs dataframe
    overlaps_df = pl.DataFrame(
        {
            "_query_idx": idx_query,
            "_ref_idx": idx_ref,
        }
    )

    # Join with original dataframes
    result = overlaps_df.join(query_indexed, on="_query_idx").join(
        ref_indexed, on="_ref_idx", suffix="_ref"
    )

    # Filter by match_by columns if specified
    if match_by:
        for col in match_by:
            if col in result.columns and f"{col}_ref" in result.columns:
                result = result.filter(pl.col(col) == pl.col(f"{col}_ref"))

    # Calculate intersection of query and reference intervals
    result = result.with_columns(
        [
            pl.max_horizontal(pl.col("Start"), pl.col("Start_ref")).alias(
                "_intersect_start"
            ),
            pl.min_horizontal(pl.col("End"), pl.col("End_ref")).alias("_intersect_end"),
        ]
    )

    # Store global coordinates if requested
    if keep_global_chrom:
        result = result.with_columns(pl.col("Chromosome").alias("Chromosome_global"))

    if keep_global_loc:
        result = result.with_columns(
            [
                pl.col("Start").alias("Start_global"),
                pl.col("End").alias("End_global"),
                pl.col("Strand").alias("Strand_global"),
            ]
        )

    # Transform coordinates to local reference coordinates
    # Handle strand-aware transformation
    ref_is_minus = result["Strand_ref"] == "-"

    # For minus strand references, transform differently
    result = result.with_columns(
        [
            pl.when(ref_is_minus)
            .then(
                pl.col("End_ref") - pl.col("_intersect_end") + pl.col("_cumsum_start")
            )
            .otherwise(
                pl.col("_intersect_start")
                - pl.col("Start_ref")
                + pl.col("_cumsum_start")
            )
            .alias("_local_start"),
            pl.when(ref_is_minus)
            .then(
                pl.col("End_ref") - pl.col("_intersect_start") + pl.col("_cumsum_start")
            )
            .otherwise(
                pl.col("_intersect_end") - pl.col("Start_ref") + pl.col("_cumsum_start")
            )
            .alias("_local_end"),
        ]
    )

    # Transform strand: if query and ref have same strand, result is +, otherwise -
    result = result.with_columns(
        pl.when(pl.col("Strand") == pl.col("Strand_ref"))
        .then(pl.lit("+"))
        .otherwise(pl.lit("-"))
        .alias("_local_strand")
    )

    # Replace coordinates with local coordinates
    # ref_id_col should be in the reference dataframe and may or may not have _ref suffix
    # ref_id_col is present in joined frame as provided

    result = result.with_columns(
        [
            pl.col(ref_id_col).alias("Chromosome"),  # ref_id becomes new chromosome
            pl.col("_local_start").alias("Start"),
            pl.col("_local_end").alias("End"),
            pl.col("_local_strand").alias("Strand"),
        ]
    )

    # Select output columns
    output_cols = ["Chromosome", "Start", "End", "Strand"]

    # Add any additional query columns (excluding coordinate columns)
    for col in original_query_cols:
        if (
            col not in ["Chromosome", "Start", "End", "Strand"]
            and col in result.columns
        ):
            output_cols.append(col)

    if keep_global_chrom:
        output_cols.append("Chromosome_global")

    if keep_global_loc:
        output_cols.extend(["Start_global", "End_global", "Strand_global"])

    # Filter to available columns and return
    available_cols = [col for col in output_cols if col in result.columns]
    return result.select(available_cols)
