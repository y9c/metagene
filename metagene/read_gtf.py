#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright © 2023 Ye Chang yech1990@gmail.com
# Distributed under terms of the GNU license.
#
# Created: 2023-03-08 20:46


import polars as pl
from gtfparse import read_gtf


def tidy_record(group_df):
    name = (
        group_df.select(pl.col("transcript_id"))[0, 0]
        + "."
        + str(group_df.select(pl.col("transcript_version"))[0, 0])
    )
    col_chr = []
    col_start = []
    col_end = []
    col_type = []
    strand = group_df.select(pl.col("strand"))[0, 0]
    # r = []
    start_df = group_df.filter(pl.col("feature") == "start_codon")
    end_df = group_df.filter(pl.col("feature") == "start_codon")
    if start_df.shape[0] == 1 and end_df.shape[0] == 1:
        # start codon position
        cs = start_df.select(pl.col("start"))[0, 0]
        # stop codon position
        ce = start_df.select(pl.col("end"))[0, 0]

        if strand == "+":
            group_df = group_df.sort(pl.col("start"))
            # split a list of intervals into 5UTR, CDS, 3UTR
            for row in group_df.rows(named=True):
                if row["start"] <= cs <= row["end"]:
                    col_chr.append(row["seqname"])
                    col_start.append(row["start"])
                    col_end.append(cs - 1)
                    col_type.append("5UTR")
                    if row["end"] > ce:
                        # cds
                        col_chr.append(row["seqname"])
                        col_start.append(cs)
                        col_end.append(ce - 1)
                        col_type.append("CDS")
                        # 3 utr
                        col_chr.append(row["seqname"])
                        col_start.append(ce)
                        col_end.append(row["end"])
                        col_type.append("3UTR")
                    else:
                        col_chr.append(row["seqname"])
                        col_start.append(cs)
                        col_end.append(row["end"])
                        col_type.append("CDS")
                elif row["start"] <= ce <= row["end"]:
                    # cds
                    col_chr.append(row["seqname"])
                    col_start.append(row["start"])
                    col_end.append(ce - 1)
                    col_type.append("CDS")
                    # 3 utr
                    col_chr.append(row["seqname"])
                    col_start.append(ce)
                    col_end.append(row["end"])
                    col_type.append("3UTR")
                else:
                    col_chr.append(row["seqname"])
                    col_start.append(row["start"])
                    col_end.append(row["end"])
                    col_type.append("CDS")
        elif strand == "-":
            group_df = group_df.sort(pl.col("start"), descending=True)
            for row in group_df.rows(named=True):
                if row["start"] <= cs < row["end"]:
                    col_chr.append(row["seqname"])
                    col_start.append(row["start"])
                    col_end.append(cs - 1)
                    col_type.append("3UTR")
                    if row["end"] > ce:
                        # cds
                        col_chr.append(row["seqname"])
                        col_start.append(cs)
                        col_end.append(ce - 1)
                        col_type.append("CDS")
                        # 5 utr
                        col_chr.append(row["seqname"])
                        col_start.append(ce)
                        col_end.append(row["end"])
                        col_type.append("5UTR")
                    else:
                        col_chr.append(row["seqname"])
                        col_start.append(cs)
                        col_end.append(row["end"])
                        col_type.append("CDS")
                if row["end"] > ce:
                    # cds
                    col_chr.append(row["seqname"])
                    col_start.append(row["start"])
                    col_end.append(ce - 1)
                    col_type.append("CDS")
                    # 5 utr
                    col_chr.append(row["seqname"])
                    col_start.append(ce)
                    col_end.append(row["end"])
                    col_type.append("5UTR")
                else:
                    col_chr.append(row["seqname"])
                    col_start.append(row["start"])
                    col_end.append(row["end"])
                    col_type.append("CDS")
    return pl.DataFrame(
        [
            pl.Series("seqname", col_chr, dtype=pl.Utf8),
            pl.Series("start", col_start, dtype=pl.Int64),
            pl.Series("end", col_end, dtype=pl.Int64),
            pl.Series(
                "type", [name + ":" + i for i in col_type], dtype=pl.Utf8
            ),
            pl.Series("strand", [strand] * len(col_type), dtype=pl.Utf8),
        ]
    )


def add_order(group_df):
    feature = group_df.select(pl.col("feature"))[0, 0]
    strand = group_df.select(pl.col("strand"))[0, 0]
    if feature == "CDS":
        rank = group_df.select(
            pl.col("exon_number").alias("Index").cast(pl.Int64, strict=False)
        ).to_series(0)
    else:
        n = list(range(group_df.shape[0]))
        if strand == "+":
            rank = pl.Series("Index", n, dtype=pl.Int64)
        else:
            rank = pl.Series("Index", n[::-1], dtype=pl.Int64)

    # names=["Chromosome", "Start", "End", "Name", "Index", "Strand"],
    return pl.DataFrame(
        [
            group_df.select(pl.col("seqname").alias("Chromosome")).to_series(
                0
            ),
            group_df.select(
                pl.col("start").alias("Start").cast(pl.Int64, strict=False)
            ).to_series(0),
            group_df.select(
                pl.col("end").alias("End").cast(pl.Int64, strict=False)
            ).to_series(0),
            group_df.select(pl.col("name").alias("Name")).to_series(0),
            rank,
            group_df.select(pl.col("strand").alias("Strand")).to_series(0),
        ]
    )


def gtf_to_bed(gtf_file, bed_file):
    # 1. returns GTF with essential columns such as "feature", "seqname", "start", "end"
    #    alongside the names of any optional keys which appeared in the attribute column
    # 2. backup polars data frame to parquet file
    # df.write_parquet("hg38.parquet")
    # df = pl.read_parquet("hg38.parquet")
    # print(df.filter(pl.col("feature") == "gene"))

    df = (
        read_gtf(gtf_file, result_type="polars")
        .filter(
            pl.col("feature").cast(pl.Utf8)
            # .is_in(["start_codon", "CDS", "stop_codon", "exon"])
            .is_in(["five_prime_utr", "CDS", "three_prime_utr"])
        )
        # .with_columns(
        #    pl.concat_str(
        #        [
        #            pl.col("transcript_id"),
        #            pl.col("transcript_version"),
        #        ],
        #        separator=".",
        #    ).alias("name"),
        # )
        .with_columns(
            pl.col("transcript_id").alias("name"),
        )
        .with_columns(
            pl.concat_str(
                [
                    pl.col("name"),
                    pl.col("feature").map_dict(
                        {
                            "five_prime_utr": "5UTR",
                            "CDS": "CDS",
                            "three_prime_utr": "3UTR",
                        },
                        default="Unspecified",
                    ),
                ],
                separator=":",
            ).alias("name"),
        )
        .groupby(["name"], maintain_order=True)
        .apply(add_order)
    )
    df.write_parquet(bed_file)
    return df


if __name__ == "__main__":
    # df = pl.read_parquet("hg38.parquet")
    df = gtf_to_bed("./data/test.gtf", "./data/test.bed.parquet")
    with pl.Config() as cfg:
        cfg.set_tbl_rows(100)
        print(df)
    # d = (
    #     df.head(1000)
    #     .groupby(["transcript_id", "transcript_version", "strand"])
    #     .apply(tidy_record)
    # )
    # with pl.Config() as cfg:
    #     cfg.set_tbl_rows(1000)
