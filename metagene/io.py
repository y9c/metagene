#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright © 2023 Ye Chang yech1990@gmail.com
# Distributed under terms of the GNU license.
#
# Data I/O Module - Handles input/output operations

import polars as pl
import os
from pathlib import Path
from .utils import get_cache_dir
from .config import BUILTIN_REFERENCES


def load_sites(
    input_file_name: str,
    with_header: bool = False,
    meta_col_index: list[int] | None = None,
    separator: str = "\t",
) -> pl.DataFrame:
    """
    Load genomic sites from a file using Polars only.
    Returns:
        Two Polars DataFrame with processed site information
        df_input: Polars DataFrame with input site information
        df_meta: Polars DataFrame with meta information
    """
    df = pl.scan_csv(input_file_name, separator=separator, has_header=with_header)
    colnames = list(df.collect_schema())
    
    if meta_col_index is None:
        raise ValueError("meta_col_index must be provided")
    
    meta_col_names = [colnames[i] for i in meta_col_index]
    # check if the Chromosome, Start, End, Strand are in the colnames
    # if so add the _original_ prefix, and return a new list of colnames
    newnames = [
        "_original_" + col if col in ["Chromosome", "Start", "End", "Strand"] else col
        for col in colnames
    ]

    df = pl.scan_csv(
        input_file_name,
        separator=separator,
        has_header=with_header,
        new_columns=newnames,
        schema_overrides={meta_col_names[0]: pl.Utf8, meta_col_names[-1]: pl.Utf8},
    )

    if len(meta_col_names) == 4:
        df = df.with_columns(
            pl.col(meta_col_names[0]).alias("Chromosome"),
            pl.col(meta_col_names[1]).alias("Start"),
            pl.col(meta_col_names[2]).alias("End"),
            pl.col(meta_col_names[3]).alias("Strand"),
        )
    elif len(meta_col_names) == 3:
        df = df.with_columns(
            pl.col(meta_col_names[0]).alias("Chromosome"),
            (pl.col(meta_col_names[1]) - 1).alias("Start"),
            pl.col(meta_col_names[1]).alias("End"),
            pl.col(meta_col_names[2]).alias("Strand"),
        )
    else:
        raise ValueError("meta_col_index must specify either 3 or 4 column indices")

    return df.collect()


def parse_feature_file(feature_file_name: str) -> pl.DataFrame:
    """
    Parse feature file (BED or Parquet format) using Polars only.
    Returns:
        Polars DataFrame with processed feature information
    """
    df = pl.read_parquet(feature_file_name)
    return df


def load_reference(species: str | None = None) -> pl.DataFrame | dict:
    """
    Load built-in reference annotations for common species using Polars only.
    
    Args:
        species: Species name to load, or None to get available species
        
    Returns:
        Polars DataFrame with feature annotations, or dict of available species if species=None
    """
    if species is None:
        available = {}
        cache_dir = get_cache_dir()
        for species_name, info in BUILTIN_REFERENCES.items():
            cache_path = cache_dir / Path(info["parquet_file"]).name
            if cache_path.exists():
                file_size_mb = os.path.getsize(cache_path) / (1024 * 1024)
                available[species_name] = {
                    "file": info["parquet_file"],
                    "source": info["source_file"],
                    "description": info["description"],
                    "size_mb": round(file_size_mb, 2),
                    "location": "cache",
                }
        return available

    if species not in BUILTIN_REFERENCES:
        available_species = list(BUILTIN_REFERENCES.keys())
        raise ValueError(
            f"Species '{species}' not available. "
            f"Available species: {available_species}\n"
            f"Use load_reference() without arguments to see available species."
        )

    info = BUILTIN_REFERENCES[species]
    cache_dir = get_cache_dir()
    cache_path = cache_dir / Path(info["parquet_file"]).name
    
    if cache_path.exists():
        return parse_feature_file(str(cache_path))
    
    # File doesn't exist - prompt user to download
    try:
        import click
        from .download import download_references
        
        click.echo(f"Reference '{species}' not found locally.")
        click.echo(f"Description: {info['description']}")
        
        if click.confirm("Would you like to download it now?", default=True):
            try:
                download_references(species, silent=True)
                return parse_feature_file(str(cache_path))
            except Exception as e:
                raise RuntimeError(f"Failed to download {species}: {e}")
        else:
            raise ValueError(f"Reference '{species}' is required but not available locally.")
    except ImportError:
        # click not available, just raise an error
        raise ValueError(f"Reference '{species}' not found locally and cannot prompt for download.")
    
    raise ValueError(f"Reference '{species}' not found locally.")
