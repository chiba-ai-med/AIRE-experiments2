#!/usr/bin/env python3
"""Combined scHi-C validation plot across datasets."""
import argparse

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--evaluations", nargs="+", required=True)
    parser.add_argument("--datasets", nargs="+", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    dataset_labels = {
        "gage_seq": "GAGE-seq",
        "droplet_hic": "Droplet Hi-C",
        "droplet_hic_merged": "Droplet Hi-C\n(merged)",
    }

    all_evals = []
    for path, ds in zip(args.evaluations, args.datasets):
        df = pd.read_csv(path, sep="\t")
        df["dataset"] = dataset_labels.get(ds, ds)
        all_evals.append(df)
    combined = pd.concat(all_evals, ignore_index=True)

    fig, axes = plt.subplots(1, 3, figsize=(12, 4))
    fig.suptitle("RIPHAMMER scHi-C Pseudo-bulk Validation", fontsize=13,
                 fontweight="bold")

    colors = {"GAGE-seq": "#2980b9", "Droplet Hi-C": "#e67e22",
              "Droplet Hi-C\n(merged)": "#27ae60"}

    xlabels = [f"{d}\n(K={k})" for d, k in
               zip(combined["dataset"], combined["K"])]

    # MAE comparison
    ax = axes[0]
    x = range(len(combined))
    ax.bar(x, combined["MAE"],
           color=[colors.get(d, "gray") for d in combined["dataset"]])
    ax.set_xticks(x)
    ax.set_xticklabels(xlabels, fontsize=9)
    ax.set_ylabel("MAE")
    ax.set_title("Mean Absolute Error")

    # Correlation comparison
    ax = axes[1]
    ax.bar(x, combined["mean_correlation"],
           color=[colors.get(d, "gray") for d in combined["dataset"]])
    ax.set_xticks(x)
    ax.set_xticklabels(xlabels, fontsize=9)
    ax.set_ylabel("Mean Pearson r")
    ax.set_title("Mean Correlation")
    ax.set_ylim(0, 1.1)

    # Pure separation
    ax = axes[2]
    pure = [int(s.split("/")[0]) / int(s.split("/")[1])
            for s in combined["pure_separation"]]
    ax.bar(x, pure,
           color=[colors.get(d, "gray") for d in combined["dataset"]])
    ax.set_xticks(x)
    ax.set_xticklabels(xlabels, fontsize=9)
    ax.set_ylabel("Fraction")
    ax.set_title("Pure Separation")
    ax.set_ylim(0, 1.1)

    plt.tight_layout()
    fig.savefig(args.output, dpi=150, bbox_inches="tight")
    print(f"Saved {args.output}", flush=True)
    plt.close()


if __name__ == "__main__":
    main()
