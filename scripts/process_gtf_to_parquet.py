#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
GTF to Parquet Converter for Metagene Built-in References

This script processes large GTF files and converts them to highly compressed
Parquet files for built-in species support in the metagene package.

This is the unified script that replaces the previous separate scripts:
- process_gtf_to_parquet.py (old)
- batch_process_gtf.py (old)
- test_gtf_processing.py

Usage:
    # Process single file
    python process_gtf_to_parquet.py --single input.gtf output_name

    # Process all predefined files (23 species/genome builds)
    python process_gtf_to_parquet.py --all

    # List all files to be processed
    python process_gtf_to_parquet.py --list

    # Create species mapping template
    python process_gtf_to_parquet.py --template

Examples:
    python process_gtf_to_parquet.py --single Homo_sapiens.GRCh38.109.gtf GRCh38
    python process_gtf_to_parquet.py --all
"""

import argparse
import os
import sys
import time
from pathlib import Path
from typing import Optional

import numpy as np
import polars as pl
import pyranges as pr


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
        print(f"Error reading GTF with Polars: {e}. Consider checking GTF format.")
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
        group_by=["gene_id", "transcript_id"],
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
    # make sure the codon positions are integers
    # Start_exon', 'End_exon', 'start_codon_pos', 'stop_codon_pos' to Int32 to reduce size
    pr_tx["Start_exon"] = pr_tx["Start_exon"].astype("Int32")
    pr_tx["End_exon"] = pr_tx["End_exon"].astype("Int32")
    pr_tx["start_codon_pos"] = pr_tx["start_codon_pos"].astype("Int32")
    pr_tx["stop_codon_pos"] = pr_tx["stop_codon_pos"].astype("Int32")
    return pr_tx


def compress_and_save_parquet(pr_data: pr.PyRanges, output_path: str) -> None:
    """
    Save PyRanges data as highly compressed Parquet file using Polars.

    Args:
        pr_data: PyRanges object to save
        output_path: Output file path
    """
    print("Converting to Polars DataFrame for compression...")

    # Convert PyRanges to pandas DataFrame first, then to Polars
    df_polars = pl.from_pandas(pr_data)  # type: ignore

    print(f"Original DataFrame shape: {df_polars.shape}")
    print(f"Columns: {list(df_polars.columns)}")  # type: ignore

    # Calculate memory usage
    memory_mb = df_polars.estimated_size("mb")  # type: ignore
    print(f"Memory usage: {memory_mb:.2f} MB")

    # Optimize data types for better compression
    # print("Optimizing data types...")

    # # Convert string columns to categorical where appropriate
    # for col in df_polars.columns:  # type: ignore
    #     if df_polars[col].dtype == pl.Utf8:  # type: ignore
    #         unique_count = df_polars[col].n_unique()  # type: ignore
    #         total_count = len(df_polars)
    #
    #         if unique_count < total_count * 0.1:  # Less than 10% unique values
    #             df_polars = df_polars.with_columns(  # type: ignore
    #                 pl.col(col).cast(pl.Categorical).alias(col)
    #             )
    #             print(f"  Converting {col} to Categorical ({unique_count} unique values)")

    # Calculate optimized memory usage
    optimized_memory_mb = df_polars.estimated_size("mb")  # type: ignore
    print(f"Optimized memory usage: {optimized_memory_mb:.2f} MB")

    # Save with maximum compression using Polars
    print(f"Saving to {output_path} with ZSTD compression (level 22)...")

    df_polars.write_parquet(  # type: ignore
        output_path, compression="zstd", compression_level=22
    )

    # Check file size
    file_size_mb = os.path.getsize(output_path) / (1024 * 1024)
    print(f"Compressed file size: {file_size_mb:.2f} MB")

    # Calculate compression ratio
    compression_ratio = optimized_memory_mb / file_size_mb if file_size_mb > 0 else 0
    print(f"Compression ratio: {compression_ratio:.1f}x")


def process_gtf_file(gtf_path: str, output_name: Optional[str] = None) -> str:
    """
    Process a GTF file and save as compressed Parquet.

    Args:
        gtf_path: Path to input GTF file
        output_name: Output name (without extension)

    Returns:
        Path to output Parquet file
    """
    gtf_path_obj = Path(gtf_path)

    if not gtf_path_obj.exists():
        raise FileNotFoundError(f"GTF file not found: {gtf_path}")

    if output_name is None:
        # Try to infer output name from filename
        output_name = gtf_path_obj.stem
        if output_name.endswith(".gtf"):
            output_name = output_name[:-4]

    output_path = f"{output_name}.parquet"

    print("=" * 60)
    print(f"PROCESSING GTF FILE: {gtf_path}")
    print(f"OUTPUT FILE: {output_path}")
    print("=" * 60)

    # Check input file size
    input_size_mb = gtf_path_obj.stat().st_size / (1024 * 1024)
    print(f"Input GTF file size: {input_size_mb:.2f} MB")

    # Process GTF file
    start_time = time.time()
    print("\nProcessing GTF file with prepare_exon_ref...")

    try:
        pr_result = prepare_exon_ref(gtf_path)
        processing_time = time.time() - start_time
        print(f"GTF processing completed in {processing_time:.2f} seconds")

        if pr_result.empty:
            raise ValueError("No data was extracted from the GTF file")

        print(f"Extracted {len(pr_result)} genomic features")

        # Save as compressed Parquet
        compress_start_time = time.time()
        compress_and_save_parquet(pr_result, output_path)
        compression_time = time.time() - compress_start_time

        print(f"Compression completed in {compression_time:.2f} seconds")

        # Calculate compression ratio
        output_size_mb = os.path.getsize(output_path) / (1024 * 1024)
        compression_ratio = input_size_mb / output_size_mb if output_size_mb > 0 else 0

        print("\n" + "=" * 60)
        print("COMPRESSION SUMMARY")
        print("=" * 60)
        print(f"Input size:        {input_size_mb:.2f} MB")
        print(f"Output size:       {output_size_mb:.2f} MB")
        print(f"Compression ratio: {compression_ratio:.1f}x")
        print(
            f"Space saved:       {((input_size_mb - output_size_mb) / input_size_mb * 100):.1f}%"
        )
        print(f"Total time:        {(processing_time + compression_time):.2f} seconds")

        return output_path

    except Exception as e:
        print(f"Error processing GTF file: {e}")
        raise


def create_species_mapping_template():
    """Create a template for species mapping configuration."""
    template = """# Species Mapping Configuration Template
# Add this to your metagene package after generating Parquet files

BUILTIN_SPECIES = {
    # Human
    'GRCh38': 'GRCh38.parquet',
    'hg38': 'hg38.parquet',
    'GRCh37': 'GRCh37.parquet', 
    'hg19': 'hg19.parquet',
    
    # Mouse
    'GRCm39': 'GRCm39.parquet',
    'mm39': 'mm39.parquet',
    'GRCm38': 'GRCm38.parquet',
    'mm10': 'mm10.parquet',
    'mm9': 'mm9.parquet',
    'NCBIM37': 'NCBIM37.parquet',
    
    # Model organisms
    'TAIR10': 'TAIR10.parquet',  # Arabidopsis
    'IRGSP-1.0': 'IRGSP-1.0.parquet',  # Rice
    'WBcel235': 'WBcel235.parquet',  # C. elegans
    'ce11': 'ce11.parquet',  # C. elegans
    'BDGP6.32': 'BDGP6.32.parquet',  # Drosophila
    'dm6': 'dm6.parquet',  # Drosophila
    'GRCz11': 'GRCz11.parquet',  # Zebrafish
    'GRCz10': 'GRCz10.parquet',  # Zebrafish
    'bGalGal1': 'bGalGal1.parquet',  # Chicken
    'Glycine_max_v2.1': 'Glycine_max_v2.1.parquet',  # Soybean
    'R64-1-1': 'R64-1-1.parquet',  # S. cerevisiae
    'sacCer3': 'sacCer3.parquet',  # S. cerevisiae
    'ASM294v2': 'ASM294v2.parquet',  # S. pombe
}
"""

    with open("species_mapping_template.py", "w") as f:
        f.write(template)

    print("Created species_mapping_template.py")


def main():
    """Main function to process all GTF files."""

    # Define all GTF files to process
    gtf_files = [
        # Human
        {
            "name": "GRCh38",
            "file": "Homo_sapiens/raw/GRCh38.release110.gtf.gz",
            "description": "Human genome GRCh38 (Ensembl release 110)",
        },
        {
            "name": "GRCh37",
            "file": "Homo_sapiens/raw/GRCh37.release75.gtf.gz",
            "description": "Human genome GRCh37 (Ensembl release 75)",
        },
        {
            "name": "hg38",
            "file": "Homo_sapiens/raw/hg38.20221028.gtf.gz",
            "description": "Human genome hg38 (UCSC 2022)",
        },
        {
            "name": "hg19",
            "file": "Homo_sapiens/raw/hg19.20210517.gtf.gz",
            "description": "Human genome hg19 (UCSC 2021)",
        },
        # Mouse
        {
            "name": "GRCm39",
            "file": "Mus_musculus/raw/GRCm39.release110.gtf.gz",
            "description": "Mouse genome GRCm39 (Ensembl release 110)",
        },
        {
            "name": "GRCm38",
            "file": "Mus_musculus/raw/GRCm38.release102.gtf.gz",
            "description": "Mouse genome GRCm38 (Ensembl release 102)",
        },
        {
            "name": "mm39",
            "file": "Mus_musculus/raw/mm39.20240214.gtf.gz",
            "description": "Mouse genome mm39 (UCSC 2024)",
        },
        {
            "name": "mm10",
            "file": "Mus_musculus/raw/mm10.20210423.gtf.gz",
            "description": "Mouse genome mm10 (UCSC 2021)",
        },
        {
            "name": "mm9",
            "file": "Mus_musculus/raw/mm9.20200110.gtf.gz",
            "description": "Mouse genome mm9 (UCSC 2020)",
        },
        {
            "name": "NCBIM37",
            "file": "Mus_musculus/raw/NCBIM37.release67.gtf.gz",
            "description": "Mouse genome NCBIM37 (Ensembl release 67)",
        },
        # Model organisms
        {
            "name": "TAIR10",
            "file": "Arabidopsis_thaliana/raw/TAIR10.release57.gtf.gz",
            "description": "Arabidopsis thaliana TAIR10",
        },
        {
            "name": "IRGSP-1.0",
            "file": "Oryza_sativa/raw/IRGSP-1.0.release57.gtf.gz",
            "description": "Rice IRGSP-1.0",
        },
        {
            "name": "WBcel235",
            "file": "Caenorhabditis_elegans/raw/WBcel235.release110.gtf.gz",
            "description": "C. elegans WBcel235",
        },
        {
            "name": "ce11",
            "file": "Caenorhabditis_elegans/raw/ce11.20200110.gtf.gz",
            "description": "C. elegans ce11 (UCSC)",
        },
        {
            "name": "BDGP6.32",
            "file": "Drosophila_melanogaster/raw/BDGP6.32.release110.gtf.gz",
            "description": "D. melanogaster BDGP6.32",
        },
        {
            "name": "dm6",
            "file": "Drosophila_melanogaster/raw/dm6.20210210.gtf.gz",
            "description": "D. melanogaster dm6 (UCSC)",
        },
        {
            "name": "GRCz11",
            "file": "Danio_rerio/raw/GRCz11.release110.gtf.gz",
            "description": "Zebrafish GRCz11",
        },
        {
            "name": "GRCz10",
            "file": "Danio_rerio/raw/GRCz10.release91.gtf.gz",
            "description": "Zebrafish GRCz10",
        },
        {
            "name": "bGalGal1",
            "file": "Gallus_gallus/raw/bGalGal1.mat.broiler.GRCg7b.release110.gtf.gz",
            "description": "Chicken bGalGal1.mat.broiler.GRCg7b",
        },
        {
            "name": "Glycine_max_v2.1",
            "file": "Glycine_max/raw/Glycine_max_v2.1.release57.gtf.gz",
            "description": "Soybean Glycine max v2.1",
        },
        {
            "name": "R64-1-1",
            "file": "Saccharomyces_cerevisiae/raw/R64-1-1.release57.gtf.gz",
            "description": "S. cerevisiae R64-1-1",
        },
        {
            "name": "sacCer3",
            "file": "Saccharomyces_cerevisiae/raw/sacCer3.20210210.gtf.gz",
            "description": "S. cerevisiae sacCer3 (UCSC)",
        },
        {
            "name": "ASM294v2",
            "file": "Schizosaccharomyces_pombe/raw/ASM294v2.release57.gtf.gz",
            "description": "S. pombe ASM294v2",
        },
    ]

    parser = argparse.ArgumentParser(
        description="Convert GTF files to compressed Parquet for metagene built-in references",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python process_gtf_to_parquet.py --single Homo_sapiens.GRCh38.109.gtf GRCh38
  python process_gtf_to_parquet.py --all  # Process all predefined GTF files
  python process_gtf_to_parquet.py --list  # List all files to be processed
        """,
    )

    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument(
        "--single",
        nargs=2,
        metavar=("GTF_FILE", "OUTPUT_NAME"),
        help="Process a single GTF file",
    )
    group.add_argument(
        "--all", action="store_true", help="Process all predefined GTF files"
    )
    group.add_argument(
        "--list", action="store_true", help="List all files to be processed"
    )
    group.add_argument(
        "--template", action="store_true", help="Create species mapping template file"
    )

    args = parser.parse_args()

    if args.template:
        create_species_mapping_template()
        return

    if args.list:
        print("Files to be processed:")
        print("=" * 80)
        for i, gtf_config in enumerate(gtf_files, 1):
            print(f"{i:2d}. {gtf_config['name']}")
            print(f"    File: {gtf_config['file']}")
            print(f"    Description: {gtf_config['description']}")
            print()
        return

    if args.single:
        # Process single file
        gtf_file, output_name = args.single
        try:
            output_path = process_gtf_file(gtf_file, output_name)
            print(f"\n✅ SUCCESS: Created {output_path}")
            print("\nNext steps:")
            print("1. Copy this Parquet file to your metagene package data directory")
            print("2. Update the BUILTIN_SPECIES mapping in your io.py module")
            print("3. Test with: from metagene.io import load_builtin_features")

        except Exception as e:
            print(f"\n❌ ERROR: {e}")
            sys.exit(1)

    elif args.all:
        # Process all predefined files
        print(f"Processing {len(gtf_files)} GTF files...")
        print("=" * 80)

        results = []

        for i, gtf_config in enumerate(gtf_files):
            name = gtf_config["name"]
            gtf_file = gtf_config["file"]
            description = gtf_config["description"]

            print(f"\n[{i + 1}/{len(gtf_files)}] Processing {name}: {description}")

            if not os.path.exists(gtf_file):
                print(f"  ❌ GTF file not found: {gtf_file}")
                results.append(
                    {"name": name, "status": "FAILED", "reason": "File not found"}
                )
                continue

            try:
                output_path = process_gtf_file(gtf_file, name)
                results.append(
                    {"name": name, "status": "SUCCESS", "output": output_path}
                )
                print(f"  ✅ Success: {output_path}")

            except Exception as e:
                print(f"  ❌ Failed: {e}")
                results.append({"name": name, "status": "FAILED", "reason": str(e)})

        # Print summary
        print("\n" + "=" * 80)
        print("BATCH PROCESSING SUMMARY")
        print("=" * 80)

        successful = [r for r in results if r["status"] == "SUCCESS"]
        failed = [r for r in results if r["status"] == "FAILED"]

        print(f"Total files:    {len(results)}")
        print(f"Successful:     {len(successful)}")
        print(f"Failed:         {len(failed)}")

        if successful:
            print("\n✅ Successful:")
            for result in successful:
                print(f"  {result['name']}: {result['output']}")

        if failed:
            print("\n❌ Failed:")
            for result in failed:
                print(f"  {result['name']}: {result['reason']}")

        # Create species mapping code
        if successful:
            print("\n" + "=" * 80)
            print("SPECIES MAPPING CODE")
            print("=" * 80)
            print("Add this to your metagene/io.py file:")
            print()
            print("BUILTIN_SPECIES = {")
            for result in successful:
                name = result["name"]
                filename = Path(result["output"]).name
                print(f"    '{name}': '{filename}',")
            print("}")

        if len(failed) == 0:
            print("\n🎉 All files processed successfully!")
            print("\nNext steps:")
            print(
                "1. Copy the generated .parquet files to your metagene/data/ directory"
            )
            print("2. Update the BUILTIN_SPECIES mapping in metagene/io.py")
            print("3. Test with: metagene.io.load_builtin_features('GRCh38')")
        else:
            print("\n⚠️  Some files failed to process. Check the summary above.")
            sys.exit(1)


if __name__ == "__main__":
    main()
