#!/usr/bin/env python3
"""Evaluate RIPHAMMER deconvolution accuracy against known fractions.

Finds the best permutation mapping between NMF components and cell types,
then computes MAE, correlation, and pure-separation metrics.
"""
import argparse
import json

import numpy as np
import pandas as pd
from scipy import stats
from scipy.optimize import linear_sum_assignment


def find_best_permutation(H, true_frac):
    """Find the component permutation that best matches true fractions.

    Uses Hungarian algorithm (O(K^3)) instead of exhaustive search (O(K!)).
    """
    K = H.shape[0]
    M = H.shape[1]
    assert true_frac.shape == (M, K), f"Shape mismatch: H={H.shape}, frac={true_frac.shape}"

    true_mat = true_frac.values  # (M, K)
    h_mat = H.T  # (M, K)

    cost = np.zeros((K, K))
    for i in range(K):
        for j in range(K):
            r, _ = stats.pearsonr(h_mat[:, i], true_mat[:, j])
            cost[i, j] = -r  # negate for minimization

    row_ind, col_ind = linear_sum_assignment(cost)
    perm = tuple(col_ind[np.argsort(row_ind)])
    best_corr = -cost[row_ind, col_ind].mean()

    return perm, best_corr


def evaluate(H, true_frac, perm):
    """Compute evaluation metrics."""
    K = H.shape[0]
    M = H.shape[1]
    true_mat = true_frac.values  # (M, K)
    inv_perm = list(np.argsort(perm))
    h_mat = H.T[:, inv_perm]  # (M, K) permuted: column k = NMF component for cell type k

    mae = np.mean(np.abs(h_mat - true_mat))
    per_type_mae = np.mean(np.abs(h_mat - true_mat), axis=0)
    per_type_corr = [
        stats.pearsonr(h_mat[:, k], true_mat[:, k])[0]
        for k in range(K)
    ]
    mean_corr = np.mean(per_type_corr)

    pure_sep = 0
    for k in range(K):
        dominant = np.argmax(h_mat, axis=1)
        true_dominant = np.argmax(true_mat, axis=1)
        if np.any(true_dominant == k):
            mask = true_dominant == k
            if np.all(dominant[mask] == k):
                pure_sep += 1
    pure_sep_frac = pure_sep / K

    results = {
        "MAE": mae,
        "mean_correlation": mean_corr,
        "pure_separation": f"{pure_sep}/{K}",
    }
    for k, ct in enumerate(true_frac.columns):
        results[f"MAE_{ct}"] = per_type_mae[k]
        results[f"corr_{ct}"] = per_type_corr[k]

    return results


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--H", required=True)
    parser.add_argument("--true-fractions", required=True)
    parser.add_argument("--meta", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    H = np.loadtxt(args.H, delimiter="\t")
    true_frac = pd.read_csv(args.true_fractions, sep="\t", index_col=0)

    with open(args.meta) as f:
        meta = json.load(f)

    print(f"H shape: {H.shape}", flush=True)
    print(f"True fractions shape: {true_frac.shape}", flush=True)

    perm, best_corr = find_best_permutation(H, true_frac)
    print(f"Best permutation: {perm} (mean corr={best_corr:.4f})", flush=True)

    results = evaluate(H, true_frac, perm)
    results["permutation"] = str(perm)
    results["dataset"] = meta.get("dataset", "unknown")
    results["K"] = meta["K"]

    print("\n=== Evaluation Results ===", flush=True)
    for k, v in results.items():
        print(f"  {k}: {v}", flush=True)

    results_df = pd.DataFrame([results])
    results_df.to_csv(args.output, sep="\t", index=False)
    print(f"\nSaved to {args.output}", flush=True)


if __name__ == "__main__":
    main()
