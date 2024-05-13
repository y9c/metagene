#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright © 2022 Ye Chang yech1990@gmail.com
# Distributed under terms of the GNU license.
#
# Created: 2022-07-18 22:42

import sys

import numpy as np
import pandas as pd
import polars as pl
import pyranges as pr


def parse_features(feature_file_name: str) -> pd.DataFrame:
    if feature_file_name.endswith(".bed") or feature_file_name.endswith(".bed.gz"):
        df = pd.read_csv(
            feature_file_name,
            sep="\t",
            names=["Chromosome", "Start", "End", "Name", "Index", "Strand"],
            comment="#",
        )
    else:
        df = pl.read_parquet(feature_file_name).to_pandas()
    df = df.sort_values(["Name", "Index"], ascending=True).assign(
        len_of_window=lambda x: x["End"] - x["Start"]
    )
    df["len_of_feature"] = df.groupby("Name")["len_of_window"].transform("sum")
    df["frac_of_feature"] = (
        df.groupby("Name")["len_of_window"].transform("cumsum") - df["len_of_window"]
    ) / df["len_of_feature"]
    df[["Transcript", "Type"]] = df["Name"].str.split(":", expand=True)
    return df


def parse_input(
    input_file_name: str,
    with_header: bool = False,
    meta_col_index: list = [0, 1, 2, 5],
    weight_col_index: list = [4],
    weight_col_name: list = [],
) -> pd.DataFrame:
    if len(meta_col_index) == 4:
        col_name_dict = dict(
            zip(meta_col_index, ["Chromosome", "Start", "End", "Strand"])
        )
        col_type_dict = dict(zip(meta_col_index, [str, int, int, str]))
        for i, w in enumerate(weight_col_index):
            col_name_dict[w] = "Weight_" + (
                weight_col_name[i] if weight_col_name else str(i + 1)
            )
            col_type_dict[w] = np.float64
        col_name_dict = dict(sorted(col_name_dict.items()))

        df = pd.read_csv(
            input_file_name,
            sep="\t",
            usecols=col_name_dict.keys(),
            names=col_name_dict.values(),
            dtype={v: col_type_dict[k] for k, v in col_name_dict.items()},
            comment="#",
            skiprows=1 if with_header else 0,
        )
    elif len(meta_col_index) == 3:
        col_name_dict = dict(zip(meta_col_index, ["Chromosome", "End", "Strand"]))
        col_type_dict = dict(zip(meta_col_index, [str, int, str]))
        for i, w in enumerate(weight_col_index):
            col_name_dict[w] = "Weight_" + (
                weight_col_name[i] if weight_col_name else str(i + 1)
            )
            col_type_dict[w] = np.float64
        col_name_dict = dict(sorted(col_name_dict.items()))
        df = pd.read_csv(
            input_file_name,
            sep="\t",
            usecols=col_name_dict.keys(),
            names=col_name_dict.values(),
            dtype={v: col_type_dict[k] for k, v in col_name_dict.items()},
            comment="#",
            skiprows=1 if with_header else 0,
        ).assign(Start=lambda x: x["End"] - 1)
    df["Chromosome"] = df["Chromosome"].str.replace("chrM", "MT").str.replace("chr", "")
    return df


def cal_bin_means(data, weights, num_bins=100, bin_range=(0, 1), suffix=""):
    data = np.asarray(data)
    # Generate the bin edges
    bins = np.linspace(bin_range[0], bin_range[1], num_bins + 1)

    # Initialize an array to hold the sum of the values in each bin
    bin_sums = np.zeros(num_bins)

    # Initialize an array to hold the count of the number of values in each bin
    bin_counts = np.zeros(num_bins)

    # Calculate which bin each data point falls into
    bin_indices = np.digitize(data, bins) - 1

    # Update the bin_sums and bin_counts arrays
    for w, b in zip(weights, bin_indices):
        if b >= 0 and b < num_bins:
            bin_sums[b] += w
            bin_counts[b] += 1

    # Calculate the mean for each bin
    bin_means = np.divide(
        bin_sums,
        bin_counts,
        out=np.zeros_like(bin_sums),
        where=bin_counts != 0,
    )

    df = pd.DataFrame(
        {"count": bin_counts, "sum": bin_sums, "mean": bin_means},
        index=(bins[1:] + bins[:-1]) / 2,
    )
    # df["count"] = df["count"].astype("Int64")
    if suffix:
        df.columns = [c + suffix for c in df.columns]
    return df


def annotate_with_feature(
    df_input: pd.DataFrame,
    df_feature: pd.DataFrame,
    bin_number=100,
    nb_cpu=8,
    type_ratios=None,
    annot_name=False,
    keep_all=False,
    by_strand=False,
) -> (pd.DataFrame, pd.DataFrame):
    df = (
        pr.PyRanges(df_input)
        .join(
            pr.PyRanges(df_feature),
            suffix="_ref",
            report_overlap=True,
            nb_cpu=nb_cpu,
            strandedness="same" if by_strand else None,
        )
        .df
    )

    # df.groupby(
    #     ["Chromosome", "Start", "End"], as_index=False, group_keys=False
    # ).apply(lambda x: x.nlargest(1, "Overlap"))
    df = (
        df.sort_values(by=["Overlap"], ascending=False)
        .groupby(["Chromosome", "Start", "End"], as_index=False, group_keys=False)
        .head(1)
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

    if type_ratios is not None:
        type_ratios = np.array(type_ratios)
        type2ratio = dict(
            zip(["5UTR", "CDS", "3UTR"], type_ratios / np.sum(type_ratios))
        )
    else:
        # type to ratio is differ for different input bins
        s = (
            df_feature[df_feature["Transcript"].isin(df["Transcript"].unique())]
            .groupby(["Transcript", "Type"])["len_of_window"]
            .sum()
            .reset_index()
            .groupby(["Type"])["len_of_window"]
            .median()
        )
        type2ratio = (s / s.sum()).to_dict()
        if (
            "5UTR" not in type2ratio
            or "CDS" not in type2ratio
            or "3UTR" not in type2ratio
        ):
            sys.exit(
                "Error: Given ranges do not overlap all 3 features: "
                "5'UTR, CDS and 3'UTR !"
                "Please use the type_ratios argument instead."
            )
    df["d_norm"] = np.where(
        df["Type"] == "5UTR",
        df["d"] * type2ratio["5UTR"],
        np.where(
            df["Type"] == "CDS",
            df["d"] * type2ratio["CDS"] + type2ratio["5UTR"],
            df["d"] * type2ratio["3UTR"] + type2ratio["5UTR"] + type2ratio["CDS"],
        ),
    )
    if annot_name:
        df["Name"] = df["Name_ref"]

    weight_col = [c for c in df.columns if c.startswith("Weight_")]

    df_list = []
    if len(weight_col) > 0:
        for c in weight_col:
            df_list.append(
                cal_bin_means(
                    df["d_norm"],
                    weights=df[c].fillna(0),
                    # weights=np.ones(len(df["d_norm"])),
                    num_bins=bin_number,
                    bin_range=(0, 1),
                    suffix="_" + c.replace("Weight_", ""),
                )
            )
    else:
        df_list.append(
            cal_bin_means(
                df["d_norm"],
                weights=np.ones(len(df["d_norm"])),
                num_bins=bin_number,
                bin_range=(0, 1),
            )
        )
    df_score = pd.concat(df_list, axis=1)

    if not keep_all:
        df = df.loc[:, ["Chromosome", "Start", "End", "Name", "d_norm", "Strand"]]
    # Use attrs property to store metadata in dataframe
    # DataFrame.attrs is an experimental feature, use be used with pandas >= 1.0
    df.attrs.update(type2ratio)
    df_score.attrs.update(type2ratio)
    return df, df_score


if __name__ == "__main__":
    INPUT_FILE = "../data/input.bed.gz"
    FEATURE_FILE = "../data/features.bed.gz"
    df_feature = parse_features(FEATURE_FILE)
    df_input = parse_input(INPUT_FILE)
    annotate_with_feature(df_input, df_feature).to_csv(
        sys.stdout, sep="\t", index=False, header=False
    )
