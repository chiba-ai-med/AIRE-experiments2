#!/usr/bin/env python3
"""Differential contact maps at Aire locus (10kb).

Run from repo root: python src/plot_aire_locus_10kb_differential.py
Output: plot/aire_locus_10kb_differential.png
"""
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm

res = 10000
chrom = "chr10"
start_bp, end_bp = 16_000_000, 20_000_000
start_bin = start_bp // res
end_bin = end_bp // res
n_bins = end_bin - start_bin

W = np.loadtxt("output/nmf_10000/W_nmf_k8.tsv")
feat = pd.read_csv("output/nmf_10000/features.tsv", sep="\t", header=None,
                    names=["chrom","bin_i","bin_j"])
H = np.loadtxt("output/nmf_10000/H_nmf_k8.tsv")
cells = pd.read_csv("output/nmf_10000/cells.tsv", sep="\t", header=None,
                     names=["sample","label"])
samples = cells["sample"].tolist()
H_n = H / H.sum(axis=0, keepdims=True)

mask = feat["chrom"] == chrom
bi_all = feat[mask]["bin_i"].values
bj_all = feat[mask]["bin_j"].values

def build_mat(k):
    mat = np.zeros((n_bins, n_bins))
    w_chr = W[:, k][mask.values]
    for fi in range(len(w_chr)):
        i, j = bi_all[fi], bj_all[fi]
        if start_bin <= i < end_bin and start_bin <= j < end_bin:
            li = i - start_bin
            lj = j - start_bin
            mat[li, lj] = w_chr[fi]
            mat[lj, li] = w_chr[fi]
    return mat

# Build mats for key components
mat_mtec = build_mat(0)    # comp1 = mTEC
mat_dp = build_mat(3)      # comp4 = DP
mat_dn = build_mat(4)      # comp5 = DN
mat_cd4 = build_mat(1)     # comp2 = CD4_SP

# Normalize each for comparison (row-sum normalize)
def norm_mat(m):
    s = m.sum()
    return m / s if s > 0 else m

mat_mtec_n = norm_mat(mat_mtec)
mat_dp_n = norm_mat(mat_dp)
mat_dn_n = norm_mat(mat_dn)
mat_cd4_n = norm_mat(mat_cd4)

extent = [start_bp/1e6, end_bp/1e6, end_bp/1e6, start_bp/1e6]
aire_mb = 17.96

fig, axes = plt.subplots(2, 4, figsize=(20, 10))

# Row 1: raw contact maps
for ax, mat, title in zip(axes[0], 
    [mat_mtec, mat_dp, mat_dn, mat_cd4],
    ["comp1 (mTEC)", "comp4 (DP)", "comp5 (DN)", "comp2 (CD4_SP)"]):
    vmax = np.percentile(mat[mat > 0], 98) if np.any(mat > 0) else 1
    ax.imshow(mat, cmap="YlOrRd", vmin=0, vmax=vmax, extent=extent, interpolation="none")
    ax.set_title(title, fontsize=11)
    ax.axvline(aire_mb, color="cyan", lw=1, ls="--", alpha=0.8)
    ax.axhline(aire_mb, color="cyan", lw=1, ls="--", alpha=0.8)
    ax.tick_params(labelsize=7)
axes[0][0].set_ylabel("Contact maps\n(Mb)", fontsize=10)

# Row 2: differential (mTEC - other, normalized)
for ax, mat_other_n, title in zip(axes[1],
    [mat_mtec_n, mat_mtec_n - mat_dp_n, mat_mtec_n - mat_dn_n, mat_mtec_n - mat_cd4_n],
    ["mTEC (normalized)", "mTEC - DP", "mTEC - DN", "mTEC - CD4_SP"]):
    if "normalized" in title:
        vmax = np.percentile(mat_other_n[mat_other_n > 0], 98) if np.any(mat_other_n > 0) else 1
        ax.imshow(mat_other_n, cmap="YlOrRd", vmin=0, vmax=vmax, extent=extent, interpolation="none")
    else:
        vabs = np.percentile(np.abs(mat_other_n[mat_other_n != 0]), 95) if np.any(mat_other_n != 0) else 1
        norm = TwoSlopeNorm(vmin=-vabs, vcenter=0, vmax=vabs)
        ax.imshow(mat_other_n, cmap="RdBu_r", norm=norm, extent=extent, interpolation="none")
    ax.set_title(title, fontsize=11)
    ax.axvline(aire_mb, color="cyan", lw=1, ls="--", alpha=0.8)
    ax.axhline(aire_mb, color="cyan", lw=1, ls="--", alpha=0.8)
    ax.tick_params(labelsize=7)
axes[1][0].set_ylabel("Differential\n(Mb)", fontsize=10)

plt.suptitle("Aire locus (chr10:16-20Mb) — 10kb component comparison", fontsize=14, y=1.01)
plt.tight_layout()
plt.savefig("plot/aire_locus_10kb_differential.png", dpi=200, bbox_inches="tight")
print("[saved] plot/aire_locus_10kb_differential.png")
