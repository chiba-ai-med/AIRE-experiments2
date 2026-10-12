#!/usr/bin/env python3
"""Aire locus deep view: 8 components + reference + insulation (10kb).

Run from repo root: python src/plot_aire_locus_10kb_deep.py
Output: plot/aire_locus_10kb_deep.png
"""
import numpy as np
import pandas as pd
import cooler
import cooltools
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

chrom = "chr10"
start_bp, end_bp = 17_000_000, 19_500_000
aire_mb = 17.96
res = 10000

cells = pd.read_csv("output/nmf_10000/cells.tsv", sep="\t", header=None,
                     names=["sample","label"])
samples = cells["sample"].tolist()

H = np.loadtxt("output/nmf_10000/H_nmf_k8.tsv")
H_n = H / H.sum(axis=0, keepdims=True)

comps = [
    (0, "mTEC_WT", "red"),
    (3, "DP_thymocyte", "blue"),
    (4, "DN_thymocyte", "green"),
    (1, "CD4_SP_thymocyte", "orange"),
    (7, "Mature_Treg_Thy1.1", "purple"),
]

# Reference mTEC insulation at 10kb
clr_ref = cooler.Cooler(f"data/mtec/processed/mtec_hic.mcool::resolutions/{res}")
ins_ref = cooltools.insulation(clr_ref, [500000], verbose=False, ignore_diags=2)
ref_chr = ins_ref[(ins_ref["chrom"] == chrom) &
                  (ins_ref["start"] >= start_bp) &
                  (ins_ref["start"] <= end_bp)]

# Contact maps
W = np.loadtxt("output/nmf_10000/W_nmf_k8.tsv")
feat = pd.read_csv("output/nmf_10000/features.tsv", sep="\t", header=None,
                    names=["chrom","bin_i","bin_j"])
mask = feat["chrom"] == chrom
bi_all = feat[mask]["bin_i"].values
bj_all = feat[mask]["bin_j"].values
start_bin = start_bp // res
end_bin = end_bp // res
n_bins = end_bin - start_bin

def build_mat(k):
    mat = np.zeros((n_bins, n_bins))
    w_chr = W[:, k][mask.values]
    for fi in range(len(w_chr)):
        i, j = bi_all[fi], bj_all[fi]
        if start_bin <= i < end_bin and start_bin <= j < end_bin:
            mat[i - start_bin, j - start_bin] = w_chr[fi]
            mat[j - start_bin, i - start_bin] = w_chr[fi]
    return mat

mat_ref = clr_ref.matrix(balance=True).fetch(f"{chrom}:{start_bp}-{end_bp}")
mat_ref = np.nan_to_num(mat_ref, nan=0.0)

extent = [start_bp/1e6, end_bp/1e6, end_bp/1e6, start_bp/1e6]

fig = plt.figure(figsize=(22, 14))
gs = fig.add_gridspec(3, 5, height_ratios=[2, 2, 1.2], hspace=0.3, wspace=0.3)

# Row 1: contact maps (ref + 4 comps)
titles_row1 = ["Reference (mTEC bulk)"] + [f"comp{k+1} ({label})" for k, label, _ in comps[:4]]
mats_row1 = [mat_ref] + [build_mat(k) for k, _, _ in comps[:4]]
cmaps_row1 = ["YlOrRd"] * 5

for col, (mat, title) in enumerate(zip(mats_row1, titles_row1)):
    ax = fig.add_subplot(gs[0, col])
    vmax = np.percentile(mat[mat > 0], 98) if np.any(mat > 0) else 1
    ax.imshow(mat, cmap="YlOrRd", vmin=0, vmax=vmax, extent=extent, interpolation="none")
    ax.set_title(title, fontsize=10)
    ax.axvline(aire_mb, color="cyan", lw=1, ls="--", alpha=0.8)
    ax.axhline(aire_mb, color="cyan", lw=1, ls="--", alpha=0.8)
    ax.tick_params(labelsize=7)
    if col == 0:
        ax.set_ylabel("Mb", fontsize=9)

# Row 2: differential (mTEC - others) + Treg
from matplotlib.colors import TwoSlopeNorm
mat_mtec = build_mat(0)
def norm_mat(m):
    s = m.sum()
    return m / s if s > 0 else m
mat_mtec_n = norm_mat(mat_mtec)

diff_list = [
    ("mTEC (normalized)", mat_mtec_n, "YlOrRd", False),
    ("mTEC - DP", mat_mtec_n - norm_mat(build_mat(3)), "RdBu_r", True),
    ("mTEC - DN", mat_mtec_n - norm_mat(build_mat(4)), "RdBu_r", True),
    ("mTEC - CD4_SP", mat_mtec_n - norm_mat(build_mat(1)), "RdBu_r", True),
    (f"comp8 (Treg_Thy1.1)", build_mat(7), "YlOrRd", False),
]

for col, (title, mat, cmap, is_diff) in enumerate(diff_list):
    ax = fig.add_subplot(gs[1, col])
    if is_diff:
        vabs = np.percentile(np.abs(mat[mat != 0]), 95) if np.any(mat != 0) else 1
        norm = TwoSlopeNorm(vmin=-vabs, vcenter=0, vmax=vabs)
        ax.imshow(mat, cmap=cmap, norm=norm, extent=extent, interpolation="none")
    else:
        vmax = np.percentile(mat[mat > 0], 98) if np.any(mat > 0) else 1
        ax.imshow(mat, cmap=cmap, vmin=0, vmax=vmax, extent=extent, interpolation="none")
    ax.set_title(title, fontsize=10)
    ax.axvline(aire_mb, color="cyan", lw=1, ls="--", alpha=0.8)
    ax.axhline(aire_mb, color="cyan", lw=1, ls="--", alpha=0.8)
    ax.tick_params(labelsize=7)
    if col == 0:
        ax.set_ylabel("Mb", fontsize=9)

# Row 3: insulation profiles (spanning all 5 cols)
ax_ins = fig.add_subplot(gs[2, :])

# Reference
x_ref = ref_chr["start"].values / 1e6
y_ref = ref_chr["log2_insulation_score_500000"].values
ax_ins.plot(x_ref, y_ref, color="black", linewidth=2, alpha=0.8, label="Reference (mTEC bulk)")

ref_bounds = ref_chr[ref_chr["is_boundary_500000"] == True]
ax_ins.scatter(ref_bounds["start"].values / 1e6, ref_bounds["log2_insulation_score_500000"].values,
               color="black", s=40, zorder=5, marker="v")

# NMF components
for k, label, color in comps:
    ins = pd.read_csv(f"output/tad_10000_k8/comp{k+1}_insulation.tsv", sep="\t")
    chr_ins = ins[(ins["chrom"] == chrom) &
                  (ins["start"] >= start_bp) &
                  (ins["start"] <= end_bp)]
    x = chr_ins["start"].values / 1e6
    y = chr_ins["log2_insulation_score_500000"].values
    ax_ins.plot(x, y, color=color, alpha=0.7, linewidth=1.2, label=f"comp{k+1} ({label})")
    
    bounds = chr_ins[chr_ins["is_boundary_500000"] == True]
    ax_ins.scatter(bounds["start"].values / 1e6, bounds["log2_insulation_score_500000"].values,
                   color=color, s=25, zorder=4, marker="v")

ax_ins.axvline(aire_mb, color="cyan", lw=1.5, ls="--", alpha=0.8, label="Aire gene")
ax_ins.set_xlabel("Position (Mb)", fontsize=11)
ax_ins.set_ylabel("log2 insulation", fontsize=11)
ax_ins.legend(fontsize=8, loc="upper right", ncol=2)
ax_ins.set_xlim(start_bp/1e6, end_bp/1e6)
ax_ins.tick_params(labelsize=8)

plt.suptitle("Aire locus (chr10:17-19.5Mb) — 10kb Multi-band NMF", fontsize=14, y=1.01)
plt.savefig("plot/aire_locus_10kb_deep.png", dpi=200, bbox_inches="tight")
print("[saved] plot/aire_locus_10kb_deep.png")
