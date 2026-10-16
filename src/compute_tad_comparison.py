#!/usr/bin/env python3
"""Compute TAD comparison across 3 methods on simulation data.

For each chromosome × cell type:
  - Build contact matrix from W (RIPHAMMER, Standard NMF) and THUNDER W
  - Compute insulation score, find TAD boundaries
  - Compare with ground truth (pure sample)

Inputs:
  output/simulation_100000/X_bulk.npy
  output/simulation_100000/features.tsv
  output/nmf_simulation_100000/W_nmf_k4.tsv
  output/nmf_simulation_100000/H_nmf_k4.tsv
  output/simulation_100000/thunder_W_k4.tsv
  output/simulation_100000/thunder_H_k4_raw.tsv

Output:
  output/simulation_100000/tad_comparison.npz
"""
import sys
import numpy as np
import pandas as pd
from scipy.stats import pearsonr
from sklearn.decomposition import NMF

sim_dir = sys.argv[1]   # output/simulation_100000
nmf_dir = sys.argv[2]   # output/nmf_simulation_100000
out_path = sys.argv[3]  # output/simulation_100000/tad_comparison.npz

X = np.load(f"{sim_dir}/X_bulk.npy")
feat = pd.read_csv(f"{sim_dir}/features.tsv", sep="\t", header=None,
                   names=["chrom", "bin1", "bin2"])
distances = np.abs(feat["bin2"].values - feat["bin1"].values)
resolution = 100000

W_riph = np.loadtxt(f"{nmf_dir}/W_nmf_k4.tsv")
H_riph = np.loadtxt(f"{nmf_dir}/H_nmf_k4.tsv")

model_std = NMF(n_components=4, init="nndsvd", random_state=42, max_iter=500)
W_std = model_std.fit_transform(np.maximum(X, 0))
H_std = model_std.components_

thunder_w_df = pd.read_csv(f"{sim_dir}/thunder_W_k4.tsv", sep="\t")
W_thunder_raw = thunder_w_df.iloc[:, 1:].values
thunder_features = thunder_w_df["feature_name"].values
thunder_feat_to_idx = {f: i for i, f in enumerate(thunder_features)}
H_thunder_raw = np.loadtxt(f"{sim_dir}/thunder_H_k4_raw.tsv",
                           delimiter="\t", skiprows=1)


def map_components(H, n_pure=4):
    H_norm = H / H.sum(axis=0, keepdims=True)
    comp_map = {np.argmax(H_norm[:, i]): i for i in range(n_pure)}
    col_order, used = [], set()
    for ct in range(n_pure):
        for c, m in comp_map.items():
            if m == ct and c not in used:
                col_order.append(c)
                used.add(c)
                break
        else:
            cands = [(c, H_norm[c, ct]) for c in range(H.shape[0])
                     if c not in used]
            best = max(cands, key=lambda x: x[1])[0]
            col_order.append(best)
            used.add(best)
    return col_order


riph_order = map_components(H_riph)
std_order = map_components(H_std)
thunder_order = map_components(H_thunder_raw)


def build_mat(b1, b2, vals, n):
    mat = np.zeros((n, n))
    for i, j, v in zip(b1, b2, vals):
        mat[i, j] = v
        mat[j, i] = v
    return mat


def insulation_score(mat, window=5):
    n = mat.shape[0]
    s = np.full(n, np.nan)
    for i in range(window, n - window):
        s[i] = np.mean(mat[i - window:i, i:i + window])
    v = ~np.isnan(s)
    if v.sum() > 0 and np.std(s[v]) > 0:
        s[v] = (s[v] - np.mean(s[v])) / np.std(s[v])
    return s


def find_tad_boundaries(ins, min_prom=0.5):
    bounds = []
    for i in range(1, len(ins) - 1):
        if np.isnan(ins[i]):
            continue
        if ins[i] < ins[i - 1] and ins[i] < ins[i + 1]:
            lm = np.nanmax(ins[max(0, i - 5):i])
            rm = np.nanmax(ins[i + 1:min(len(ins), i + 6)])
            if min(lm - ins[i], rm - ins[i]) > min_prom:
                bounds.append(i)
    return np.array(bounds)


def boundary_metrics(true_b, pred_b, tol=2):
    if len(true_b) == 0 or len(pred_b) == 0:
        return np.nan, np.nan, np.nan
    tp = sum(any(abs(p - t) <= tol for t in true_b) for p in pred_b)
    prec = tp / len(pred_b) if len(pred_b) > 0 else 0
    rec = tp / len(true_b) if len(true_b) > 0 else 0
    f1 = 2 * prec * rec / (prec + rec) if (prec + rec) > 0 else 0
    return prec, rec, f1


results = {k: [] for k in [
    "RIPHAMMER_precision", "RIPHAMMER_recall", "RIPHAMMER_f1",
    "RIPHAMMER_insul_corr",
    "THUNDER_precision", "THUNDER_recall", "THUNDER_f1",
    "THUNDER_insul_corr",
    "Standard NMF_precision", "Standard NMF_recall", "Standard NMF_f1",
    "Standard NMF_insul_corr",
]}

chroms = [f"chr{i}" for i in range(1, 20)]

for chrom in chroms:
    chr_mask = feat["chrom"] == chrom
    cf = feat[chr_mask]
    if len(cf) == 0:
        continue

    bins = sorted(set(cf["bin1"]) | set(cf["bin2"]))
    n_bins = len(bins)
    b2i = {b: i for i, b in enumerate(bins)}
    bin1_idx = np.array([b2i[b] for b in cf["bin1"]])
    bin2_idx = np.array([b2i[b] for b in cf["bin2"]])
    global_idx = cf.index.values

    for ct in range(4):
        mat_true = build_mat(bin1_idx, bin2_idx,
                             X[global_idx, ct], n_bins)
        ins_true = insulation_score(mat_true)
        tad_true = find_tad_boundaries(ins_true)

        mat_riph = build_mat(bin1_idx, bin2_idx,
                             W_riph[global_idx, riph_order[ct]], n_bins)
        ins_riph = insulation_score(mat_riph)
        tad_riph = find_tad_boundaries(ins_riph)

        mat_std = build_mat(bin1_idx, bin2_idx,
                            W_std[global_idx, std_order[ct]], n_bins)
        ins_std = insulation_score(mat_std)
        tad_std = find_tad_boundaries(ins_std)

        thunder_vals = np.zeros(len(global_idx))
        for k, gidx in enumerate(global_idx):
            row = feat.iloc[gidx]
            fname = (f"{row['chrom']}_{int(row['bin1']) * resolution}"
                     f"_{row['chrom']}_{int(row['bin2']) * resolution}")
            if fname in thunder_feat_to_idx:
                thunder_vals[k] = W_thunder_raw[
                    thunder_feat_to_idx[fname], thunder_order[ct]]
        mat_thunder = build_mat(bin1_idx, bin2_idx, thunder_vals, n_bins)
        ins_thunder = insulation_score(mat_thunder)
        tad_thunder = find_tad_boundaries(ins_thunder)

        for method, ins_pred, tad_pred in [
            ("RIPHAMMER", ins_riph, tad_riph),
            ("THUNDER", ins_thunder, tad_thunder),
            ("Standard NMF", ins_std, tad_std),
        ]:
            p, r, f1 = boundary_metrics(tad_true, tad_pred)
            results[f"{method}_precision"].append(p)
            results[f"{method}_recall"].append(r)
            results[f"{method}_f1"].append(f1)

            valid = ~(np.isnan(ins_true) | np.isnan(ins_pred))
            if valid.sum() > 10:
                corr, _ = pearsonr(ins_true[valid], ins_pred[valid])
            else:
                corr = np.nan
            results[f"{method}_insul_corr"].append(corr)

    print(f"  {chrom}: done ({n_bins} bins)", flush=True)

np.savez(out_path, **{k: np.array(v) for k, v in results.items()})
print(f"[saved] {out_path}")
