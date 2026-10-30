#!/usr/bin/env python3
"""Vectorize synthetic mixture mcool files into a feature matrix for RIPHAMMER."""
import argparse
import os

import numpy as np
import riphammer


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--mixtures-mcool", required=True,
                        help="File listing mcool paths, one per line")
    parser.add_argument("--sample-names", required=True,
                        help="cells.tsv with name and group columns")
    parser.add_argument("--resolution", type=int, required=True)
    parser.add_argument("--max-dist", type=int, default=5000000)
    parser.add_argument("--outdir", required=True)
    args = parser.parse_args()

    with open(args.mixtures_mcool) as f:
        mcool_paths = [line.strip() for line in f if line.strip()]

    cool_uris = [f"{p}::resolutions/{args.resolution}" for p in mcool_paths]

    samples = []
    with open(args.sample_names) as f:
        for line in f:
            name = line.strip().split("\t")[0]
            samples.append(name)

    print(f"[vectorize] {len(cool_uris)} mixtures, "
          f"resolution={args.resolution}", flush=True)

    X, feat_df, dist_arr, sample_names = riphammer.mcool_to_features(
        cool_uris,
        resolution=args.resolution,
        max_dist_bp=args.max_dist,
        autosomes_only=True,
        labels=samples,
    )

    os.makedirs(args.outdir, exist_ok=True)
    np.save(os.path.join(args.outdir, "X_bulk.npy"), X)
    feat_df.to_csv(os.path.join(args.outdir, "features.tsv"),
                   sep="\t", header=False, index=False)

    with open(os.path.join(args.outdir, "cells.tsv"), "w") as f:
        for name in sample_names:
            f.write(f"{name}\tscHiC_mixture\n")

    print(f"[saved] {args.outdir}/ (X: {X.shape})", flush=True)


if __name__ == "__main__":
    main()
