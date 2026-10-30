#!/usr/bin/env python3
"""Run RIPHAMMER NMF on synthetic scHi-C mixtures."""
import argparse
import json
import os

import numpy as np
import pandas as pd
import riphammer


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--X", required=True)
    parser.add_argument("--features", required=True)
    parser.add_argument("--cells", required=True)
    parser.add_argument("--true-fractions", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--max-iter", type=int, default=500)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--mode", default="single")
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    X = np.load(args.X)
    feat_df = pd.read_csv(args.features, sep="\t", header=None)
    distances = feat_df.iloc[:, -1].values.astype(int)

    true_frac = pd.read_csv(args.true_fractions, sep="\t", index_col=0)
    K = true_frac.shape[1]

    print(f"X shape: {X.shape}", flush=True)
    print(f"K (from true fractions): {K}", flush=True)
    print(f"Mode: {args.mode}", flush=True)

    model = riphammer.MultibandNMF(
        n_components=K,
        max_iter=args.max_iter,
        random_state=args.seed,
        mode=args.mode,
        normalize_h=True,
        normalize_bands=True,
    )
    model.fit(X, distances, verbose=True)

    H = model.H_
    W = model.get_W_full()

    np.savetxt(os.path.join(args.outdir, "H.tsv"), H, delimiter="\t")
    np.savetxt(os.path.join(args.outdir, "W.tsv"), W, delimiter="\t")

    meta = {
        "K": K,
        "celltypes": list(true_frac.columns),
        "samples": list(true_frac.index),
        "mode": args.mode,
        "max_iter": args.max_iter,
        "seed": args.seed,
        "reconstruction_error": model.reconstruction_error_,
        "X_shape": list(X.shape),
    }
    with open(os.path.join(args.outdir, "meta.json"), "w") as f:
        json.dump(meta, f, indent=2)

    print(f"\nH shape: {H.shape}", flush=True)
    print(f"W shape: {W.shape}", flush=True)
    print(f"Reconstruction error: {model.reconstruction_error_:.6f}", flush=True)
    print(f"Saved to {args.outdir}/", flush=True)


if __name__ == "__main__":
    main()
