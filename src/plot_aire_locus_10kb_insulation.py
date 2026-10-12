#!/usr/bin/env python3
"""Insulation profile comparison at Aire locus (10kb).

Run from repo root: python src/plot_aire_locus_10kb_insulation.py
Output: plot/aire_locus_10kb_insulation.png
"""
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

chrom = "chr10"
start_bp, end_bp = 15_000_000, 22_000_000
aire_mb = 17.96

cells = pd.read_csv("output/nmf_10000/cells.tsv", sep="\t", header=None,
                     names=["sample","label"])
samples = cells["sample"].tolist()

H = np.loadtxt("output/nmf_10000/H_nmf_k8.tsv")
H_n = H / H.sum(axis=0, keepdims=True)

comps = [
    (0, "mTEC_WT", "comp1"),
    (3, "DP_thymocyte", "comp4"),
    (4, "DN_thymocyte", "comp5"),
    (1, "CD4_SP_thymocyte", "comp2"),
    (7, "Mature_Treg_Thy1.1", "comp8"),
]
colors = ["red", "blue", "green", "orange", "purple"]

fig, axes = plt.subplots(2, 1, figsize=(16, 8), sharex=True)

for (k, label, cname), color in zip(comps, colors):
    ins = pd.read_csv(f"output/tad_10000_k8/{cname.replace('comp','comp')}_insulation.tsv", sep="\t")
    # fix: file is comp1, comp2, etc
    fname = f"output/tad_10000_k8/comp{k+1}_insulation.tsv"
    ins = pd.read_csv(fname, sep="\t")
    chr_ins = ins[(ins["chrom"] == chrom) & 
                  (ins["start"] >= start_bp) & 
                  (ins["start"] <= end_bp)]
    
    x = chr_ins["start"].values / 1e6
    y = chr_ins["log2_insulation_score_500000"].values
    
    axes[0].plot(x, y, color=color, alpha=0.7, linewidth=1.2, label=f"{cname} ({label})")
    
    # boundary positions
    bcol = "is_boundary_500000"
    if bcol in chr_ins.columns:
        bounds = chr_ins[chr_ins[bcol] == True]
        axes[1].scatter(bounds["start"].values / 1e6, [k]*len(bounds), 
                       color=color, s=15, alpha=0.8)

axes[0].axvline(aire_mb, color="cyan", lw=1.5, ls="--", alpha=0.8, label="Aire")
axes[0].set_ylabel("log2 insulation score", fontsize=11)
axes[0].legend(fontsize=8, loc="lower right")
axes[0].set_title(f"Insulation profiles — chr10:{start_bp//1e6:.0f}-{end_bp//1e6:.0f}Mb (10kb)", fontsize=13)

axes[1].axvline(aire_mb, color="cyan", lw=1.5, ls="--", alpha=0.8)
axes[1].set_xlabel("Position (Mb)", fontsize=11)
axes[1].set_ylabel("Component", fontsize=11)
axes[1].set_yticks(range(len(comps)))
axes[1].set_yticklabels([f"{c[2]} ({c[1]})" for c in comps], fontsize=9)
axes[1].set_title("TAD boundary positions", fontsize=13)

plt.tight_layout()
plt.savefig("plot/aire_locus_10kb_insulation.png", dpi=200, bbox_inches="tight")
print("[saved] plot/aire_locus_10kb_insulation.png")
