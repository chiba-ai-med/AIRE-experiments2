#!/usr/bin/env python3
"""Selected 5 components at Aire locus (mTEC + major thymocytes, 10kb).

Run from repo root: python src/plot_aire_locus_10kb_selected.py
Output: plot/aire_locus_10kb_selected.png
"""
import numpy as np
import pandas as pd
import cooler
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

res = 10000
chrom = "chr10"
start_bp, end_bp = 15_000_000, 22_000_000
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

# Key components: comp1=mTEC, comp4=DP, comp5=DN, comp8=Thy1.1
show = [
    (0, "mTEC_WT"),
    (3, "DP_thymocyte"),
    (4, "DN_thymocyte"),
    (7, "Mature_Treg_Thy1.1"),
]

# Reference
clr_ref = cooler.Cooler(f"data/mtec/processed/mtec_hic.mcool::resolutions/{res}")
mat_ref = clr_ref.matrix(balance=True).fetch(f"{chrom}:{start_bp}-{end_bp}")
mat_ref = np.nan_to_num(mat_ref, nan=0.0)

fig, axes = plt.subplots(1, 5, figsize=(25, 5))
extent = [start_bp/1e6, end_bp/1e6, end_bp/1e6, start_bp/1e6]

# Reference
ax = axes[0]
vmax = np.percentile(mat_ref[mat_ref > 0], 98) if np.any(mat_ref > 0) else 1
ax.imshow(mat_ref, cmap="YlOrRd", vmin=0, vmax=vmax, extent=extent, interpolation="none")
ax.set_title("Reference\n(mTEC bulk Hi-C)", fontsize=11)
ax.set_ylabel("Mb", fontsize=10)

for idx, (k, label) in enumerate(show):
    ax = axes[idx + 1]
    mat = np.zeros((n_bins, n_bins))
    w_chr = W[:, k][mask.values]
    bi = feat[mask]["bin_i"].values
    bj = feat[mask]["bin_j"].values
    for fi in range(len(w_chr)):
        i, j = bi[fi], bj[fi]
        if start_bin <= i < end_bin and start_bin <= j < end_bin:
            li = i - start_bin
            lj = j - start_bin
            mat[li, lj] = w_chr[fi]
            mat[lj, li] = w_chr[fi]
    vmax_k = np.percentile(mat[mat > 0], 98) if np.any(mat > 0) else 1
    ax.imshow(mat, cmap="YlOrRd", vmin=0, vmax=vmax_k, extent=extent, interpolation="none")
    loading = H_n[k, samples.index(label)]
    ax.set_title(f"comp{k+1} ({label})\nloading={loading:.3f}", fontsize=10)

for ax in axes:
    ax.axvline(17.96, color="cyan", linewidth=1, linestyle="--", alpha=0.8)
    ax.axhline(17.96, color="cyan", linewidth=1, linestyle="--", alpha=0.8)
    ax.tick_params(labelsize=8)

plt.suptitle("Aire locus (chr10:15-22Mb) — 10kb Multi-band NMF", fontsize=14, y=1.02)
plt.tight_layout()
plt.savefig("plot/aire_locus_10kb_selected.png", dpi=200, bbox_inches="tight")
print("[saved] plot/aire_locus_10kb_selected.png")
