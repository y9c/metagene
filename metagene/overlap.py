#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright © 2022 Ye Chang yech1990@gmail.com
# Distributed under terms of the GNU license.
#
# Created: 2022-07-18 22:42

import os
import sys

import numpy as np
import pandas as pd
import pyranges as pr


def parse_features(feature_file_name: str) -> tuple[pd.DataFrame, dict]:
    df = (
        pd.read_csv(
            feature_file_name,
            sep="\t",
            names=["Chromosome", "Start", "End", "Name", "Index", "Strand"],
            comment="#",
        )
        .sort_values(["Name", "Index"], ascending=True)
        .assign(len_of_windown=lambda x: x["End"] - x["Start"])
    )
    df["len_of_feature"] = df.groupby("Name")["len_of_windown"].transform(
        "sum"
    )
    df["frac_of_feature"] = (
        df.groupby("Name")["len_of_windown"].transform("cumsum")
        - df["len_of_windown"]
    ) / df["len_of_feature"]
    df["Type"] = df["Name"].str.split(":").str[-1]
    s = df.groupby("Type")["len_of_windown"].sum()
    return df, (s / s.sum()).to_dict()


def parse_input(input_file_name: str) -> pd.DataFrame:
    df = pd.read_csv(
        input_file_name,
        sep="\t",
        usecols=[0, 1, 2, 5],
        names=["Chromosome", "Start", "End", "Strand"],
        comment="#",
    )

    df["Chromosome"] = (
        df["Chromosome"].str.replace("chrM", "MT").str.replace("chr", "")
    )
    return df


def annotate_with_feature(
    df_input: pd.DataFrame,
    df_feature: pd.DataFrame,
    type2ratio: dict,
    nb_cpu=8,
) -> pd.DataFrame:
    df = (
        pr.PyRanges(df_input)
        .join(
            pr.PyRanges(df_feature),
            suffix="_ref",
            report_overlap=True,
            nb_cpu=nb_cpu,
        )
        .df.groupby(
            ["Chromosome", "Start", "End"], as_index=False, group_keys=False
        )
        .apply(lambda x: x.nlargest(1, "Overlap"))
        .assign(
            d=lambda x: np.where(
                x.Strand_ref == "+",
                (np.minimum((x.Start + x.End) // 2, x.End_ref) - x.Start_ref)
                / x.len_of_feature
                + x.frac_of_feature,
                (x.End_ref - np.maximum((x.Start + x.End) // 2, x.Start_ref))
                / x.len_of_feature
                + x.frac_of_feature,
            )
        )
    )
    df["d_norm"] = np.where(
        df["Type"] == "5UTR",
        df["d"] * type2ratio["5UTR"],
        np.where(
            df["Type"] == "CDS",
            df["d"] * type2ratio["CDS"] + type2ratio["5UTR"],
            df["d"] * type2ratio["3UTR"]
            + type2ratio["5UTR"]
            + type2ratio["CDS"],
        ),
    )
    df = df.loc[:, ["Chromosome", "Start", "End", "Name", "d_norm", "Strand"]]
    # Use attrs property to store metadata in dataframe
    # DataFrame.attrs is an experimental feature, use be used with pandas >= 1.0
    print(type2ratio, file=sys.stderr)
    df.attrs.update(type2ratio)
    return df


if __name__ == "__main__":
    INPUT_FILE = "../data/input.bed.gz"
    FEATURE_FILE = "../data/features.bed.gz"
    df_feature, type2ratio = parse_features(FEATURE_FILE)
    df_input = parse_input(INPUT_FILE)
    annotate_with_feature(df_input, df_feature, type2ratio).to_csv(
        sys.stdout, sep="\t", index=False, header=False
    )
