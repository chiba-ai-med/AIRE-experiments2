#!/usr/bin/env python3
"""Visualize NMF component W columns as contact maps."""
import argparse
import os

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def w_to_matrix(w_col, feat, chrom, bin_start, bin_end):
    mask = feat["chrom"] == chrom
    feat_chr = feat[mask]
    w_chr = w_col[mask.values]
    size = bin_end - bin_start
    mat = np.zeros((size, size))
    for idx in range(len(feat_chr)):
        bi = feat_chr.iloc[idx]["bin_i"]
        bj = feat_chr.iloc[idx]["bin_j"]
        if bi >= bin_start and bi < bin_end and bj >= bin_start and bj < bin_end:
            li = bi - bin_start
            lj = bj - bin_start
            mat[li, lj] = w_chr[idx]
            mat[lj, li] = w_chr[idx]
    return mat


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True)
    parser.add_argument("--k", type=int, default=8)
    parser.add_argument("--resolution", type=int, default=100000)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    W = np.loadtxt(os.path.join(args.input, f"W_nmf_k{args.k}.tsv"),
                   delimiter="\t")
    feat = pd.read_csv(os.path.join(args.input, "features.tsv"), sep="\t",
                       header=None, names=["chrom", "bin_i", "bin_j"])
    H = np.loadtxt(os.path.join(args.input, f"H_nmf_k{args.k}.tsv"),
                   delimiter="\t")

    cells_path = os.path.join(args.input, "cells.tsv")
    if os.path.exists(cells_path):
        cells = pd.read_csv(cells_path, sep="\t", header=None,
                            names=["sample", "label"])
        sample_names = cells["sample"].tolist()
    else:
        sample_names = [f"s{i}" for i in range(H.shape[1])]

    H_frac = H / H.sum(axis=0, keepdims=True)
    res = args.resolution

    regions = [
        ("chr10", 170, 270, "Aire locus (chr10:17-27Mb)"),
        ("chr1", 500, 700, "chr1:50-70Mb"),
        ("chr6", 500, 650, "chr6:50-65Mb"),
    ]

    fig, axes = plt.subplots(len(regions), args.k,
                             figsize=(3.5 * args.k, 3.5 * len(regions)))

    for ri, (chrom, bstart, bend, title) in enumerate(regions):
        for k in range(args.k):
            ax = axes[ri, k]
            mat = w_to_matrix(W[:, k], feat, chrom, bstart, bend)
            start_mb = bstart * res / 1e6
            end_mb = bend * res / 1e6
            extent = [start_mb, end_mb, end_mb, start_mb]
            vmax = np.percentile(mat[mat > 0], 98) if np.any(mat > 0) else 1
            ax.imshow(mat, cmap="YlOrRd", vmin=0, vmax=vmax,
                      extent=extent, interpolation="none")
            if ri == 0:
                top = sample_names[np.argmax(H_frac[k, :])]
                ax.set_title(f"comp{k+1}\n({top})", fontsize=8)
            if k == 0:
                ax.set_ylabel(f"{title}\n(Mb)", fontsize=8)
            ax.tick_params(labelsize=6)

    plt.suptitle(f"Multi-band NMF (K={args.k}, {res//1000}kb)",
                 fontsize=13, y=1.01)
    plt.tight_layout()
    os.makedirs(os.path.dirname(args.output), exist_ok=True)
    plt.savefig(args.output, dpi=150, bbox_inches="tight")
    print(f"[saved] {args.output}", flush=True)


if __name__ == "__main__":
    main()
