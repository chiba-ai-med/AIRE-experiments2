#!/usr/bin/env python3
"""Vectorize brain-only mcool files using RIPHAMMER."""
import argparse
import os

import yaml
import numpy as np

import riphammer


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    parser.add_argument("--resolution", type=int, required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--max-dist", type=int, default=5000000)
    args = parser.parse_args()

    with open(args.config) as f:
        cfg = yaml.safe_load(f)

    cool_uris = []
    labels = []
    group_labels = []
    for s in cfg["samples"]["brain"]["entries"]:
        cool_uris.append(f"{s['mcool']}::resolutions/{args.resolution}")
        labels.append(s["name"])
        group_labels.append("brain")

    print(f"[vectorize_brain] {len(cool_uris)} samples, "
          f"resolution={args.resolution}", flush=True)

    X, feat_df, dist_arr, sample_names = riphammer.mcool_to_features(
        cool_uris,
        resolution=args.resolution,
        max_dist_bp=args.max_dist,
        autosomes_only=True,
        labels=labels,
    )

    os.makedirs(args.output, exist_ok=True)
    np.save(os.path.join(args.output, "X_bulk.npy"), X)
    feat_df.to_csv(os.path.join(args.output, "features.tsv"),
                   sep="\t", header=False, index=False)

    with open(os.path.join(args.output, "cells.tsv"), "w") as f:
        for name, group in zip(sample_names, group_labels):
            f.write(f"{name}\t{group}\n")

    print(f"[saved] {args.output}/ (X: {X.shape})", flush=True)


if __name__ == "__main__":
    main()
