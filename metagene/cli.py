#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright © 2022 Ye Chang yech1990@gmail.com
# Distributed under terms of the GNU license.
#
# Created: 2022-07-09 01:18

"""metagene cli."""

import importlib.resources
import os
import sys

import click

from .data import __name__ as data_package
from .overlap import annotate_with_feature, parse_features, parse_input
from .read_gtf import gtf_to_bed


@click.command(
    help="metagene command line interface",
    no_args_is_help=True,
    context_settings=dict(help_option_names=["-h", "--help"]),
)
@click.option(
    "--input", "-i", type=click.Path(exists=True), help="Input file."
)
@click.option(
    "--columns",
    "-c",
    type=str,
    default="1,2,3,6",
    help="Input columns, Chromosome,Start,End,Strand or Chromosome,Site,Strand",
)
@click.option(
    "--with-header", "-H", is_flag=True, help="Input file with header."
)
@click.option(
    "--features", "-f", type=click.Path(exists=True), help="Freature file."
)
@click.option(
    "--threads", "-t", type=int, default=8, help="Number of threads."
)
@click.option(
    "--buildin-features", "-b", default="GRCh38", help="Buildin features."
)
def cli(input, with_header, columns, features, threads, buildin_features):
    if features is None:
        if buildin_features not in ["GRCh38", "GRCm39", "TAIR10", "IRGSP-1.0"]:
            sys.exit(
                f"Feature ({buildin_features}) has not been implemented yet..."
            )
        df_feature = parse_features(
            importlib.resources.files("metagene.data").joinpath(
                f"{buildin_features}.bed.parquet"
            )
        )
    else:
        if features.endswith(".gtf"):
            bed_file = features.rsplit(".", 1)[0] + ".bed.parquet"
            if not os.path.exists(bed_file):
                gtf_to_bed(features, bed_file)
            df_feature = parse_features(bed_file)
        else:
            df_feature = parse_features(features)
    df_input = parse_input(
        input, with_header, col_index=[int(x) - 1 for x in columns.split(",")]
    )
    annotate_with_feature(df_input, df_feature, nb_cpu=threads).to_csv(
        sys.stdout, sep="\t", index=False, header=False
    )
