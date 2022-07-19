import numpy as np
import pandas as pd

df = pd.read_csv("canonical_transcript_exon.bed", sep="\t", low_memory=False)
for (t, s), g in df.groupby(["Transcript stable ID version", "Strand"]):
    g = g.sort_values("Exon region start (bp)", ascending=s == "+")
    f = -1
    n = 0

    def print_pos(p1, p2, n):
        if f == -1:
            f2str = ":5UTR"
        elif f == 0:
            f2str = ":CDS"
        elif f == 1:
            f2str = ":3UTR"
        else:
            f2str = ""
        print(chrom, p1, p2, t + f2str, n, s, sep="\t")
        n += 10000

    if np.isnan(g.iloc[:, 3]).all():
        continue
    if s == "+":
        for _, r in g.iloc[:, [0, 1, 2, 3, 4]].iterrows():
            chrom, exon_start, exon_end, coding_start, coding_end = r.values
            coding_start = (
                np.nan if np.isnan(coding_start) else int(coding_start)
            )
            coding_end = np.nan if np.isnan(coding_end) else int(coding_end)
            if f == -1 and coding_start >= exon_start:
                if coding_start > exon_start:
                    print_pos(exon_start, coding_start, n)
                    # n += coding_start - exon_start
                    n += 1
                f += 1
                n = 0
                print_pos(coding_start, coding_end, n)
                # n += coding_end - coding_start
                n += 1
                if coding_end < exon_end:
                    f += 1
                    n = 0
                    print_pos(coding_end, exon_end, n)
            elif f == 0 and coding_end < exon_end:
                print_pos(coding_start, coding_end, n)
                # n += coding_end - coding_start
                n += 1
                f += 1
                n = 0
                print_pos(coding_end, exon_end, n)
                # n += exon_end - coding_end
                n += 1
            else:
                print_pos(exon_start, exon_end, n)
                # n += exon_end - exon_start
                n += 1
    else:
        for _, r in g.iloc[:, [0, 1, 2, 3, 4]].iterrows():
            chrom, exon_start, exon_end, coding_start, coding_end = r.values
            coding_start = (
                np.nan if np.isnan(coding_start) else int(coding_start)
            )
            coding_end = np.nan if np.isnan(coding_end) else int(coding_end)
            if f == -1 and coding_end <= exon_end:
                if coding_end < exon_end:
                    print_pos(coding_end, exon_end, n)
                    # n += exon_end - coding_end
                    n += 1
                f += 1
                n = 0
                print_pos(coding_start, coding_end, n)
                # n += coding_end - coding_start
                n += 1
                if coding_start > exon_start:
                    f += 1
                    n = 0
                    print_pos(exon_start, coding_start, n)
            elif f == 0 and coding_start > exon_start:
                print_pos(coding_start, exon_end, n)
                # n += exon_end - coding_start
                n += 1
                f += 1
                n = 0
                print_pos(exon_start, coding_start, n)
                # n += coding_start - exon_start
                n += 1
            else:
                print_pos(exon_start, exon_end, n)
                # n += exon_end - exon_start
                n += 1
