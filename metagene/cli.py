#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright © 2022 Ye Chang yech1990@gmail.com
# Distributed under terms of the GNU license.
#
# Created: 2022-07-09 01:18

"""metagene cli."""

import importlib.resources
import sys

import click

from .data import __name__ as data_package
from .overlap import annotate_with_feature, parse_features, parse_input


@click.command(help="metagene command line interface", no_args_is_help=True)
@click.option(
    "--input", "-i", type=click.Path(exists=True), help="Input file."
)
@click.option(
    "--feature", "-f", type=click.Path(exists=True), help="Freature file."
)
@click.option(
    "--threads", "-t", type=int, default=8, help="Number of threads."
)
@click.option(
    "--buildin-features",
    "-b",
    type=click.Choice(["GRCh38", "GRCm38"]),
    help="Buildin features.",
)
def cli(input, feature, threads, buildin_features):
    if buildin_features is not None:
        df_feature, type2ratio = parse_features(
            importlib.resources.path(data_package, "GRCh38.bed.gz")
        )
    else:
        df_feature, type2ratio = parse_features(feature)
    df_input = parse_input(input)
    annotate_with_feature(
        df_input, df_feature, type2ratio, nb_cpu=threads
    ).to_csv(sys.stdout, sep="\t", index=False, header=False)
