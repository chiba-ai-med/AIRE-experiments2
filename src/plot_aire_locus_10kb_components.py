#!/usr/bin/env python3
"""All 8 NMF components at Aire locus (10kb).

Run from repo root: python src/plot_aire_locus_10kb_components.py
Output: plot/aire_locus_10kb_components.png
"""
import numpy as np
import pandas as pd
import cooler
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# Aire locus: chr10 ~17-27Mb at 10kb
res = 10000
chrom = "chr10"
start_bp, end_bp = 17_000_000, 27_000_000
start_bin = start_bp // res
end_bin = end_bp // res

# Load NMF results
W = np.loadtxt("output/nmf_10000/W_nmf_k8.tsv")
feat = pd.read_csv("output/nmf_10000/features.tsv", sep="\t", header=None,
                    names=["chrom","bin_i","bin_j"])
H = np.loadtxt("output/nmf_10000/H_nmf_k8.tsv")
cells = pd.read_csv("output/nmf_10000/cells.tsv", sep="\t", header=None,
                     names=["sample","label"])
samples = cells["sample"].tolist()

H_n = H / H.sum(axis=0, keepdims=True)

# Build contact maps for each component
mask = feat["chrom"] == chrom
feat_chr = feat[mask]
n_bins = end_bin - start_bin

comp_labels = []
for k in range(8):
    dom = samples[np.argmax(H_n[k])]
    comp_labels.append(f"comp{k+1}\n({dom})")

# Reference mTEC
clr_ref = cooler.Cooler(f"data/mtec/processed/mtec_hic.mcool::resolutions/{res}")
mat_ref = clr_ref.matrix(balance=True).fetch(f"{chrom}:{start_bp}-{end_bp}")
mat_ref = np.nan_to_num(mat_ref, nan=0.0)

fig, axes = plt.subplots(1, 9, figsize=(36, 3.5))

# Plot reference
ax = axes[0]
vmax = np.percentile(mat_ref[mat_ref > 0], 98) if np.any(mat_ref > 0) else 1
extent = [start_bp/1e6, end_bp/1e6, end_bp/1e6, start_bp/1e6]
ax.imshow(mat_ref, cmap="YlOrRd", vmin=0, vmax=vmax, extent=extent, interpolation="none")
ax.set_title("Reference\n(mTEC bulk)", fontsize=9)
ax.set_ylabel("Mb", fontsize=8)
ax.tick_params(labelsize=6)

# Plot each NMF component
for k in range(8):
    ax = axes[k + 1]
    mat = np.zeros((n_bins, n_bins))
    w_chr = W[:, k][mask.values]
    bi = feat_chr["bin_i"].values
    bj = feat_chr["bin_j"].values
    for idx in range(len(feat_chr)):
        i, j = bi[idx], bj[idx]
        if start_bin <= i < end_bin and start_bin <= j < end_bin:
            li = i - start_bin
            lj = j - start_bin
            mat[li, lj] = w_chr[idx]
            mat[lj, li] = w_chr[idx]
    
    vmax_k = np.percentile(mat[mat > 0], 98) if np.any(mat > 0) else 1
    ax.imshow(mat, cmap="YlOrRd", vmin=0, vmax=vmax_k, extent=extent, interpolation="none")
    ax.set_title(comp_labels[k], fontsize=8)
    ax.tick_params(labelsize=6)
    ax.set_yticks([])

# Aire gene position (~17.96Mb on mm10)
for ax in axes:
    ax.axvline(17.96, color="blue", linewidth=0.5, linestyle="--", alpha=0.7)
    ax.axhline(17.96, color="blue", linewidth=0.5, linestyle="--", alpha=0.7)

plt.suptitle("chr10:17-27Mb (Aire locus) — 10kb resolution", fontsize=13, y=1.05)
plt.tight_layout()
plt.savefig("plot/aire_locus_10kb_components.png", dpi=200, bbox_inches="tight")
print("[saved] plot/aire_locus_10kb_components.png")
