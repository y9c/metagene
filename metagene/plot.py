#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright © 2023 Ye Chang yech1990@gmail.com
# Distributed under terms of the GNU license.
#
# Created: 2023-08-29 21:43

import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd


def plot_metagene(df, info):
    df = df.set_index("x")
    fig, ax = plt.subplots(figsize=(3.5, 2), dpi=300)
    for c in df.columns:
        n = "All" if df.columns[0] == "y" else c.split("_", 1)[1]
        ax.plot(df.index, df[c], lw=0.8, label=n)
        ax.fill_between(df.index, df[c], alpha=0.1)

    ax.axvline(info["5UTR"], ls="--", color="silver", lw=0.5)
    ax.axvline(info["5UTR"] + info["CDS"], ls="--", color="silver", lw=0.5)
    ax.set_xticks([0, info["5UTR"], info["5UTR"] + info["CDS"], 1])
    ax.set_xticklabels([], size="xx-small")
    ax.set_xlim(0, 1)
    ax.add_patch(
        mpl.patches.Rectangle(
            (0, 0),
            1,
            0.02,
            transform=ax.transAxes,
            facecolor="k",
            fill=True,
            lw=3,
            clip_on=False,
        )
    )
    ax.add_patch(
        mpl.patches.Rectangle(
            (info["5UTR"], -0.02),
            info["CDS"],
            0.06,
            transform=ax.transAxes,
            facecolor="k",
            fill=True,
            lw=3,
            clip_on=False,
        )
    )
    ax.yaxis.set_major_locator(mpl.ticker.MaxNLocator(4))
    ax.set_xlabel("Relative position on coding genes")
    ax.set_ylabel("Cumulative level")
    ax.legend(
        loc="lower left",
        fontsize="xx-small",
        frameon=False,
        bbox_to_anchor=(0, 1),
        ncol=4,
    )
    fig.tight_layout()
    return fig


if __name__ == "__main__":
    df = pd.read_csv("./test2.tsv", index_col=0, sep="\t")
    info = {
        "3UTR": 0.4299628942486085,
        "5UTR": 0.048701298701298704,
        "CDS": 0.5213358070500927,
    }
    fig = plot_metagene(df, info)
    fig.savefig("test.pdf", bbox_inches="tight")
