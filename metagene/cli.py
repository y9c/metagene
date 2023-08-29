#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright © 2022 Ye Chang yech1990@gmail.com
# Distributed under terms of the GNU license.
#
# Created: 2022-07-09 01:18

"""metagene cli."""

import importlib.resources
import logging
import os
import sys

import rich_click as click
from rich.logging import RichHandler

from .data import __name__ as data_package
from .overlap import annotate_with_feature, parse_features, parse_input
from .read_gtf import gtf_to_bed

logger = logging.getLogger("metagene")
# logger = logging.getLogger(__name__)
logger.addHandler(RichHandler())
logger.propagate = False


@click.command(
    help="metagene command line interface",
    no_args_is_help=True,
    context_settings=dict(help_option_names=["-h", "--help"]),
)
@click.option(
    "--input", "-i", type=click.Path(exists=True), help="Input file."
)
@click.option("--output", "-o", default="-", help="Output file.")
@click.option(
    "--with-header", "-H", is_flag=True, help="Input file with header."
)
@click.option(
    "--columns",
    "-c",
    type=str,
    default="1,2,3,5,6",
    help="Input columns, [Chromosome,Start,End,Score,Strand], "
    "[Chromosome,Start,End,Strand] or [Chromosome,Site,Strand]",
)
@click.option(
    "--bin-number", "-b", type=int, default=None, help="Number of bins."
)
@click.option(
    "--features", "-f", type=click.Path(exists=True), help="Freature file."
)
@click.option(
    "--threads", "-t", type=int, default=8, help="Number of threads."
)
@click.option(
    "--buildin-features",
    "-F",
    default="GRCh38",
    type=click.Choice(
        [
            "GRCh38",
            "GRCm38",
            "GRCm39",
            "TAIR10",
            "IRGSP-1.0",
        ]
    ),
    help="Buildin features.",
)
def cli(
    input,
    output,
    with_header,
    columns,
    bin_number,
    features,
    threads,
    buildin_features,
):
    if features is None:
        logger.info("Parsing buildin features.")
        df_feature = parse_features(
            str(
                importlib.resources.files("metagene.data").joinpath(
                    f"{buildin_features}.bed.parquet"
                )
            )
        )
    else:
        if features.endswith(".gtf"):
            logger.info(f"Parsing features from gtf file ({features}).")
            bed_file = features.rsplit(".", 1)[0] + ".bed.parquet"
            if not os.path.exists(bed_file):
                logger.info("Generating cache file for gtf file.")
                gtf_to_bed(features, bed_file)
            df_feature = parse_features(bed_file)
        else:
            logger.info(f"Parsing features from file ({features}).")
            df_feature = parse_features(features)
    logger.info("Loading input data.")
    df_input = parse_input(
        input, with_header, col_index=[int(x) - 1 for x in columns.split(",")]
    )
    logger.info("Annotating input data using parsed feature data.")
    df_output = annotate_with_feature(
        df_input, df_feature, bin_number=bin_number, nb_cpu=threads
    )

    if output == "-":
        output = sys.stdout
    elif output.endswith(".gz"):
        import gzip

        output = gzip.open(output, "wt")
    else:
        output = open(output, "w")
    for k, v in df_output.attrs.items():
        if k in ["bin_y", "bin_x"]:
            logger.debug(f"{k}: {v}")
        else:
            logger.info(f"{k}: {v}")
    logger.info(
        ", and they are also written into the comment lines (#) of the output."
    )

    for k, v in df_output.attrs.items():
        print(f"# {k}: {v}", file=output)
    logger.info("Saving annotated output data.")
    df_output.to_csv(output, sep="\t", index=False, header=False)

    if "bin_y" in df_output.attrs:
        import asciichartpy

        logger.info("Plotting the distribution of the reuslts.")
        chart = asciichartpy.plot(
            df_output.attrs["bin_y"], {"height": 10, "format": "{:8.2f}"}
        )
        logger.info(chart)
