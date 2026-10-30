#!/usr/bin/env python3
"""Plot per-dataset scHi-C validation results."""
import argparse
import json

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--evaluation", required=True)
    parser.add_argument("--H", required=True)
    parser.add_argument("--true-fractions", required=True)
    parser.add_argument("--dataset", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    eval_df = pd.read_csv(args.evaluation, sep="\t")
    H = np.loadtxt(args.H, delimiter="\t")
    true_frac = pd.read_csv(args.true_fractions, sep="\t", index_col=0)

    perm = eval(eval_df["permutation"].iloc[0])
    K = H.shape[0]
    M = H.shape[1]
    celltypes = list(true_frac.columns)

    h_perm = H.T[:, list(perm)]

    dataset_labels = {
        "gage_seq": "GAGE-seq (Nat Genet 2024)",
        "droplet_hic": "Droplet Hi-C (Nat Biotech 2024)",
    }
    title = dataset_labels.get(args.dataset, args.dataset)

    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    fig.suptitle(f"RIPHAMMER Validation — {title}", fontsize=14, fontweight="bold")

    # Panel 1: Scatter — recovered vs true fractions
    ax = axes[0]
    colors = plt.cm.Set2(np.linspace(0, 1, K))
    for k in range(K):
        ax.scatter(true_frac.values[:, k], h_perm[:, k],
                   c=[colors[k]], label=celltypes[k], s=40, alpha=0.7)
    ax.plot([0, 1], [0, 1], "k--", alpha=0.3, lw=1)
    ax.set_xlabel("True fraction")
    ax.set_ylabel("Recovered fraction")
    ax.set_title(f"Fraction recovery (MAE={eval_df['MAE'].iloc[0]:.3f})")
    ax.legend(fontsize=7, loc="upper left", ncol=1)
    ax.set_xlim(-0.05, 1.05)
    ax.set_ylim(-0.05, 1.05)
    ax.set_aspect("equal")

    # Panel 2: Per-cell-type correlation bar
    ax = axes[1]
    corrs = [eval_df[f"corr_{ct}"].iloc[0] for ct in celltypes]
    bars = ax.bar(range(K), corrs, color=colors[:K])
    ax.set_xticks(range(K))
    ax.set_xticklabels(celltypes, rotation=45, ha="right", fontsize=8)
    ax.set_ylabel("Pearson r")
    ax.set_title(f"Per-type correlation (mean={eval_df['mean_correlation'].iloc[0]:.3f})")
    ax.set_ylim(-0.1, 1.1)
    ax.axhline(1.0, color="gray", ls="--", lw=0.5)

    # Panel 3: Stacked bar — H matrix composition per mixture
    ax = axes[2]
    bottom = np.zeros(M)
    for k in range(K):
        ax.bar(range(M), h_perm[:, k], bottom=bottom,
               color=colors[k], label=celltypes[k], width=0.8)
        bottom += h_perm[:, k]
    ax.set_xlabel("Mixture")
    ax.set_ylabel("Fraction")
    ax.set_title("Recovered H (stacked)")
    ax.legend(fontsize=7, loc="upper right", ncol=1)

    plt.tight_layout()
    fig.savefig(args.output, dpi=150, bbox_inches="tight")
    print(f"Saved {args.output}", flush=True)
    plt.close()


if __name__ == "__main__":
    main()
