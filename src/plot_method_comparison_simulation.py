#!/usr/bin/env python3
"""3-method comparison: RIPHAMMER vs Standard NMF vs PCA.
NOTE: This script re-runs NMF fitting (may take several minutes).

Run from repo root: python src/plot_method_comparison_simulation.py
Output: plot/method_comparison_simulation.png
"""
import numpy as np
from scipy.optimize import nnls
from sklearn.decomposition import NMF, PCA
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

X = np.load("output/simulation_100000/X_bulk.npy")
feat = np.genfromtxt("output/simulation_100000/features.tsv", delimiter="\t", dtype=str)
dist = np.abs(feat[:, 2].astype(int) - feat[:, 1].astype(int))
true_frac = np.loadtxt("output/simulation_100000/true_fractions.tsv", delimiter="\t", skiprows=1)
cell_names = ["CD4_SP", "DN", "DP", "Treg"]
mix_names = ["equal", "CD4_dom", "CD4_DN", "Treg_rare", "gradient", "extreme"]
K = 4

# RIPHAMMER
H_r = np.loadtxt("output/nmf_simulation_100000/H_nmf_k4.tsv")

# Standard NMF
model_nmf = NMF(n_components=K, init="random", random_state=42, max_iter=500)
W_std = model_nmf.fit_transform(X)
H_s = model_nmf.components_

# PCA
pca = PCA(n_components=K, random_state=42)
pca.fit_transform(X)
H_p = np.abs(pca.components_)

def get_pred(H):
    H_pure, H_mix = H[:, :4], H[:, 4:]
    preds = []
    for i in range(6):
        f, _ = nnls(H_pure, H_mix[:, i])
        s = f.sum()
        preds.append(f / s if s > 0 else f)
    return np.array(preds)

pred_r = get_pred(H_r)
pred_s = get_pred(H_s)
pred_p = get_pred(H_p)

# H matrix heatmaps + scatter plots
fig, axes = plt.subplots(2, 3, figsize=(18, 10))

methods = [
    ("RIPHAMMER\n(Multi-band NMF)", H_r, pred_r, 0.023),
    ("Standard NMF\n(no band norm)", H_s, pred_s, 0.011),
    ("PCA\n(abs loadings)", H_p, pred_p, 0.193),
]

sample_labels = cell_names + mix_names

for col, (name, H, pred, mae) in enumerate(methods):
    # Row 1: H matrix heatmap
    ax = axes[0][col]
    H_n = H / H.sum(axis=0, keepdims=True)
    im = ax.imshow(H_n, cmap="YlOrRd", vmin=0, vmax=1, aspect="auto")
    ax.set_xticks(range(10))
    ax.set_xticklabels(sample_labels, fontsize=7, rotation=45, ha="right")
    ax.set_yticks(range(K))
    ax.set_yticklabels([f"comp{k+1}" for k in range(K)], fontsize=8)
    ax.set_title(name, fontsize=11, fontweight="bold")
    if col == 0:
        ax.set_ylabel("H matrix (normalized)", fontsize=10)
    ax.axvline(3.5, color="white", lw=2)
    for i in range(K):
        for j in range(10):
            ax.text(j, i, f"{H_n[i,j]:.2f}", ha="center", va="center", fontsize=5.5,
                    color="white" if H_n[i,j] > 0.5 else "black")

    # Row 2: scatter true vs pred
    ax = axes[1][col]
    colors = ["red", "green", "blue", "orange"]
    for j, (cell, color) in enumerate(zip(cell_names, colors)):
        ax.scatter(true_frac[:, j], pred[:, j], c=color, s=60, label=cell,
                   edgecolors="black", linewidth=0.5)
    ax.plot([0, 1], [0, 1], "k--", alpha=0.5)
    ax.set_xlabel("True fraction", fontsize=10)
    if col == 0:
        ax.set_ylabel("Predicted fraction", fontsize=10)
    ax.set_title(f"MAE = {mae:.3f}", fontsize=11)
    ax.set_xlim(-0.05, 1.05)
    ax.set_ylim(-0.05, 1.05)
    if col == 2:
        ax.legend(fontsize=8, loc="upper left")

plt.suptitle("Simulation: Method Comparison (100kb, K=4)", fontsize=14, y=1.01)
plt.tight_layout()
plt.savefig("plot/method_comparison_simulation.png", dpi=200, bbox_inches="tight")
print("[saved] plot/method_comparison_simulation.png")
