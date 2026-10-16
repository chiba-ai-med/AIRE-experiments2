#!/usr/bin/env python3
"""Compute 3-way comparison: RIPHAMMER vs Standard NMF vs PCA.

Inputs:
  output/simulation_100000/X_bulk.npy
  output/simulation_100000/features.tsv
  output/simulation_100000/true_fractions.tsv
  output/nmf_simulation_100000/H_nmf_k4.tsv

Output:
  output/simulation_100000/comparison_results.npz
"""
import sys
import numpy as np
from scipy.optimize import nnls
from scipy.stats import pearsonr
from sklearn.decomposition import NMF, PCA

sim_dir = sys.argv[1]   # output/simulation_100000
nmf_dir = sys.argv[2]   # output/nmf_simulation_100000
out_path = sys.argv[3]  # output/simulation_100000/comparison_results.npz

X = np.load(f"{sim_dir}/X_bulk.npy")
feat = np.genfromtxt(f"{sim_dir}/features.tsv", delimiter="\t", dtype=str)
true_frac = np.loadtxt(f"{sim_dir}/true_fractions.tsv",
                       delimiter="\t", skiprows=1)
K = 4

H_r = np.loadtxt(f"{nmf_dir}/H_nmf_k{K}.tsv")

model_nmf = NMF(n_components=K, init="random", random_state=42, max_iter=500)
model_nmf.fit_transform(X)
H_s = model_nmf.components_

pca = PCA(n_components=K, random_state=42)
pca.fit_transform(X)
H_p = np.abs(pca.components_)


def get_pred(H):
    H_pure, H_mix = H[:, :4], H[:, 4:]
    preds = []
    for i in range(H_mix.shape[1]):
        f, _ = nnls(H_pure, H_mix[:, i])
        s = f.sum()
        preds.append(f / s if s > 0 else f)
    return np.array(preds)


pred_r = get_pred(H_r)
pred_s = get_pred(H_s)
pred_p = get_pred(H_p)


def compute_metrics(true, pred):
    mae = np.mean(np.abs(true - pred))
    rmse = np.sqrt(np.mean((true - pred) ** 2))
    rs = []
    for ct in range(true.shape[1]):
        if np.std(true[:, ct]) > 1e-8 and np.std(pred[:, ct]) > 1e-8:
            r, _ = pearsonr(true[:, ct], pred[:, ct])
            rs.append(r)
    return mae, rmse, np.mean(rs) if rs else float("nan")


methods = ["RIPHAMMER", "Standard NMF", "PCA"]
preds = [pred_r, pred_s, pred_p]
mae_arr, rmse_arr, mean_r_arr = [], [], []
for name, pred in zip(methods, preds):
    mae, rmse, mr = compute_metrics(true_frac, pred)
    mae_arr.append(mae)
    rmse_arr.append(rmse)
    mean_r_arr.append(mr)
    print(f"  {name}: MAE={mae:.4f}, RMSE={rmse:.4f}, mean_r={mr:.3f}")

np.savez(
    out_path,
    pred_riphammer=pred_r,
    pred_stdnmf=pred_s,
    pred_pca=pred_p,
    true_frac=true_frac,
    mae=np.array(mae_arr),
    rmse=np.array(rmse_arr),
    mean_r=np.array(mean_r_arr),
)
print(f"[saved] {out_path}")
