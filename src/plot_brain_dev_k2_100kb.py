#!/usr/bin/env python3
"""Brain ES vs differentiation K=2 deconvolution (100kb).

Run from repo root: python src/plot_brain_dev_k2_100kb.py
Output: plot/brain_dev_k2_100kb.png
"""
import numpy as np
import pandas as pd
import cooler
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

res = 100000
W = np.loadtxt("output/nmf_brain_100000_k2/W_nmf_k2.tsv")
feat = pd.read_csv("output/nmf_brain_100000_k2/features.tsv", sep="\t", header=None,
                    names=["chrom","bin_i","bin_j"])

comp_labels = ["comp1 (Differentiated: NPC+CN)", "comp2 (Stem: ES)"]

loci = [
    ("chr3", 34_000_000, 44_000_000, "Sox2 (chr3:34-44Mb)"),
    ("chr2", 105_000_000, 115_000_000, "Pax6 (chr2:105-115Mb)"),
    ("chr5", 85_000_000, 95_000_000, "Neurog2 (chr5:85-95Mb)"),
]

refs = [
    ("Bonev_ES", "data/brain/processed/brain_hic_es.mcool::resolutions/100000"),
    ("Bonev_NPC", "data/brain/processed/brain_hic_npc.mcool::resolutions/100000"),
    ("Bonev_CN", "data/brain/processed/brain_hic_cn.mcool::resolutions/100000"),
]

fig, axes = plt.subplots(len(loci), 5, figsize=(25, 5*len(loci)))

for ri, (chrom, start_bp, end_bp, title) in enumerate(loci):
    start_bin = start_bp // res
    end_bin = end_bp // res
    n_bins = end_bin - start_bin
    extent = [start_bp/1e6, end_bp/1e6, end_bp/1e6, start_bp/1e6]
    mask = feat["chrom"] == chrom
    bi_all = feat[mask]["bin_i"].values
    bj_all = feat[mask]["bin_j"].values

    for ci, (rname, uri) in enumerate(refs):
        ax = axes[ri, ci]
        clr = cooler.Cooler(uri)
        has_bal = "weight" in clr.bins().columns
        mat = clr.matrix(balance=has_bal).fetch(f"{chrom}:{start_bp}-{end_bp}")
        mat = np.nan_to_num(mat, nan=0.0)
        vmax = np.percentile(mat[mat > 0], 98) if np.any(mat > 0) else 1
        ax.imshow(mat, cmap="YlOrRd", vmin=0, vmax=vmax, extent=extent, interpolation="none")
        ax.set_title(f"Ref: {rname}", fontsize=10)
        if ci == 0:
            ax.set_ylabel(f"{title}\n(Mb)", fontsize=9)
        ax.tick_params(labelsize=7)

    for k in range(2):
        ax = axes[ri, k + 3]
        mat = np.zeros((n_bins, n_bins))
        w_chr = W[:, k][mask.values]
        for fi in range(len(w_chr)):
            i, j = bi_all[fi], bj_all[fi]
            if start_bin <= i < end_bin and start_bin <= j < end_bin:
                mat[i - start_bin, j - start_bin] = w_chr[fi]
                mat[j - start_bin, i - start_bin] = w_chr[fi]
        vmax = np.percentile(mat[mat > 0], 98) if np.any(mat > 0) else 1
        ax.imshow(mat, cmap="YlOrRd", vmin=0, vmax=vmax, extent=extent, interpolation="none")
        ax.set_title(comp_labels[k], fontsize=10)
        ax.tick_params(labelsize=7)

plt.suptitle("Brain development: Reference vs Multi-band NMF (100kb, K=2)\nH: comp1=[ES:0.33, NPC:0.98, CN:1.00]  comp2=[ES:0.68, NPC:0.02, CN:0.00]",
             fontsize=13, y=1.03)
plt.tight_layout()
plt.savefig("plot/brain_dev_k2_100kb.png", dpi=200, bbox_inches="tight")
print("[saved] plot/brain_dev_k2_100kb.png")
