#!/usr/bin/env python3
"""Noise/downsampling robustness (initial version; see robustness_analysis.py for final).
NOTE: This script re-runs NMF fitting (may take several minutes).

Run from repo root: python src/plot_noise_robustness.py
Output: plot/noise_robustness.png
"""

import numpy as np
import pandas as pd
import sys
import time
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.stats import pearsonr

from riphammer.core import MultibandNMF
from sklearn.decomposition import NMF

# Load data
X_orig = np.load('output/simulation_100000/X_bulk.npy')  # (1232364, 10)
feat = pd.read_csv('output/simulation_100000/features.tsv', sep='\t', header=None)
feat.columns = ['chrom', 'bin1', 'bin2']
distances = np.abs(feat['bin2'].values - feat['bin1'].values)

true_frac = np.loadtxt('output/simulation_100000/true_fractions.tsv', delimiter='\t', skiprows=1)
# true_frac: (6, 4) — CD4_SP, DN, DP, Treg_GFP

sample_names_pure = ['CD4_SP', 'DN', 'DP', 'Treg_GFP']

def map_components(H, n_pure=4):
    H_norm = H / H.sum(axis=0, keepdims=True)
    comp_map = {}
    for i in range(n_pure):
        comp_map[np.argmax(H_norm[:, i])] = i
    col_order, used = [], set()
    for ct in range(n_pure):
        for c, m in comp_map.items():
            if m == ct and c not in used:
                col_order.append(c); used.add(c); break
        else:
            cands = [(c, H_norm[c, ct]) for c in range(H.shape[0]) if c not in used]
            col_order.append(max(cands, key=lambda x: x[1])[0]); used.add(col_order[-1])
    return col_order

def evaluate(X, distances, true_frac):
    """Run RIPHAMMER and StdNMF, return MAE and pure separation for each."""
    results = {}
    
    # RIPHAMMER
    model_r = MultibandNMF(n_components=4, max_iter=500, random_state=42)
    model_r.fit(X, distances, verbose=False)
    H_r = model_r.H_
    H_r_norm = H_r / H_r.sum(axis=0, keepdims=True)
    order_r = map_components(H_r)
    H_r_ordered = H_r_norm[order_r, :]
    pred_r = H_r_ordered[:, 4:].T  # (6, 4)
    mae_r = np.mean(np.abs(true_frac - pred_r))
    # Pure separation
    pure_dom = [np.argmax(H_r_norm[:, i]) for i in range(4)]
    sep_r = len(set(pure_dom))
    results['riphammer'] = {'mae': mae_r, 'sep': sep_r, 'W': model_r.get_W_full(), 'H': H_r, 'order': order_r}
    
    # Standard NMF
    model_s = NMF(n_components=4, init='nndsvd', random_state=42, max_iter=500)
    W_s = model_s.fit_transform(np.maximum(X, 0))
    H_s = model_s.components_
    H_s_norm = H_s / H_s.sum(axis=0, keepdims=True)
    order_s = map_components(H_s)
    H_s_ordered = H_s_norm[order_s, :]
    pred_s = H_s_ordered[:, 4:].T
    mae_s = np.mean(np.abs(true_frac - pred_s))
    pure_dom_s = [np.argmax(H_s_norm[:, i]) for i in range(4)]
    sep_s = len(set(pure_dom_s))
    results['stdnmf'] = {'mae': mae_s, 'sep': sep_s, 'W': W_s, 'H': H_s, 'order': order_s}
    
    return results

# ============================================
# 1. Downsampling robustness
# ============================================
print("=== Downsampling robustness ===")
ds_fractions = [0.0, 0.1, 0.2, 0.3, 0.5, 0.7]
ds_seeds = [42, 123, 456]

ds_results = {m: {'mae_mean': [], 'mae_std': [], 'sep_mean': []} for m in ['riphammer', 'stdnmf']}

nz_mask = X_orig > 0
nz_indices = np.argwhere(nz_mask)
n_nz = len(nz_indices)
print(f"Total nonzero entries: {n_nz}")

for frac in ds_fractions:
    print(f"  Downsampling {frac*100:.0f}%...", flush=True)
    maes = {'riphammer': [], 'stdnmf': []}
    seps = {'riphammer': [], 'stdnmf': []}
    
    for seed in ds_seeds:
        rng = np.random.RandomState(seed)
        X_ds = X_orig.copy()
        if frac > 0:
            n_drop = int(n_nz * frac)
            drop_idx = rng.choice(n_nz, n_drop, replace=False)
            for idx in drop_idx:
                X_ds[nz_indices[idx, 0], nz_indices[idx, 1]] = 0
        
        res = evaluate(X_ds, distances, true_frac)
        for m in ['riphammer', 'stdnmf']:
            maes[m].append(res[m]['mae'])
            seps[m].append(res[m]['sep'])
    
    for m in ['riphammer', 'stdnmf']:
        ds_results[m]['mae_mean'].append(np.mean(maes[m]))
        ds_results[m]['mae_std'].append(np.std(maes[m]))
        ds_results[m]['sep_mean'].append(np.mean(seps[m]))
    
    print(f"    RIPHAMMER MAE: {np.mean(maes['riphammer']):.4f} ± {np.std(maes['riphammer']):.4f}, sep: {np.mean(seps['riphammer']):.1f}")
    print(f"    StdNMF    MAE: {np.mean(maes['stdnmf']):.4f} ± {np.std(maes['stdnmf']):.4f}, sep: {np.mean(seps['stdnmf']):.1f}")

# ============================================
# 2. Gaussian noise robustness
# ============================================
print("\n=== Gaussian noise robustness ===")
noise_fractions = [0.0, 0.1, 0.2, 0.5, 1.0]
mean_nz = X_orig[X_orig > 0].mean()
print(f"Mean nonzero value: {mean_nz:.6f}")

gn_results = {m: {'mae_mean': [], 'mae_std': [], 'sep_mean': []} for m in ['riphammer', 'stdnmf']}

for nf in noise_fractions:
    print(f"  Noise σ = {nf} × mean...", flush=True)
    maes = {'riphammer': [], 'stdnmf': []}
    seps = {'riphammer': [], 'stdnmf': []}
    
    for seed in ds_seeds:
        rng = np.random.RandomState(seed)
        X_noisy = X_orig.copy()
        if nf > 0:
            sigma = nf * mean_nz
            noise = rng.normal(0, sigma, X_noisy.shape)
            X_noisy = np.maximum(X_noisy + noise, 0)
        
        res = evaluate(X_noisy, distances, true_frac)
        for m in ['riphammer', 'stdnmf']:
            maes[m].append(res[m]['mae'])
            seps[m].append(res[m]['sep'])
    
    for m in ['riphammer', 'stdnmf']:
        gn_results[m]['mae_mean'].append(np.mean(maes[m]))
        gn_results[m]['mae_std'].append(np.std(maes[m]))
        gn_results[m]['sep_mean'].append(np.mean(seps[m]))
    
    print(f"    RIPHAMMER MAE: {np.mean(maes['riphammer']):.4f} ± {np.std(maes['riphammer']):.4f}, sep: {np.mean(seps['riphammer']):.1f}")
    print(f"    StdNMF    MAE: {np.mean(maes['stdnmf']):.4f} ± {np.std(maes['stdnmf']):.4f}, sep: {np.mean(seps['stdnmf']):.1f}")

# ============================================
# 3. Differential contact under noise (Aire locus)
# ============================================
print("\n=== Differential contact at Aire locus under noise ===")

def build_mat(b1, b2, vals, n):
    mat = np.zeros((n, n))
    for i, j, v in zip(b1, b2, vals):
        mat[i, j] = v; mat[j, i] = v
    return mat

def differential_score(mats, ct_idx, n_ct=4):
    target = mats[ct_idx].copy()
    others = np.mean([mats[i] for i in range(n_ct) if i != ct_idx], axis=0)
    eps = np.percentile(target[target > 0], 5) if (target > 0).sum() > 10 else 1e-10
    return np.log2((target + eps) / (others + eps))

# Aire locus: chr10, bins 770-820
chr_mask = feat['chrom'] == 'chr10'
cf = feat[chr_mask]
bins = sorted(set(cf['bin1']) | set(cf['bin2']))
n_bins_chr = len(bins)
b2i = {b: i for i, b in enumerate(bins)}

# Region: bins 770-820 in local index
region_bins = [b for b in bins if 770 <= b <= 820]
n_reg = len(region_bins)
rb2i = {b: i for i, b in enumerate(region_bins)}
cf_region = cf[(cf['bin1'].isin(region_bins)) & (cf['bin2'].isin(region_bins))]
bin1_local = np.array([rb2i[b] for b in cf_region['bin1']])
bin2_local = np.array([rb2i[b] for b in cf_region['bin2']])
global_idx_region = cf_region.index.values

# True differential (from clean data, focus on Treg_GFP = ct 3)
ct_focus = 3
true_mats = [build_mat(bin1_local, bin2_local, X_orig[global_idx_region, ct], n_reg) for ct in range(4)]
diff_true = differential_score(true_mats, ct_focus)

diff_corr_results = {'clean': {}, 'ds50': {}}

for condition, frac in [('clean', 0.0), ('ds50', 0.5)]:
    rng = np.random.RandomState(42)
    X_test = X_orig.copy()
    if frac > 0:
        n_drop = int(n_nz * frac)
        drop_idx = rng.choice(n_nz, n_drop, replace=False)
        for idx in drop_idx:
            X_test[nz_indices[idx, 0], nz_indices[idx, 1]] = 0
    
    res = evaluate(X_test, distances, true_frac)
    
    for method in ['riphammer', 'stdnmf']:
        W = res[method]['W']
        order = res[method]['order']
        pred_mats = [build_mat(bin1_local, bin2_local, W[global_idx_region, order[ct]], n_reg) for ct in range(4)]
        diff_pred = differential_score(pred_mats, ct_focus)
        
        mask = np.isfinite(diff_true) & np.isfinite(diff_pred)
        if mask.sum() > 10:
            r, _ = pearsonr(diff_true[mask].flatten(), diff_pred[mask].flatten())
        else:
            r = 0.0
        diff_corr_results[condition][method] = r
        print(f"  {condition} {method}: diff corr = {r:.3f}")

# ============================================
# Save results
# ============================================
np.savez('output/simulation_100000/noise_robustness.npz',
         ds_fractions=ds_fractions,
         ds_riphammer_mae_mean=ds_results['riphammer']['mae_mean'],
         ds_riphammer_mae_std=ds_results['riphammer']['mae_std'],
         ds_riphammer_sep=ds_results['riphammer']['sep_mean'],
         ds_stdnmf_mae_mean=ds_results['stdnmf']['mae_mean'],
         ds_stdnmf_mae_std=ds_results['stdnmf']['mae_std'],
         ds_stdnmf_sep=ds_results['stdnmf']['sep_mean'],
         gn_fractions=noise_fractions,
         gn_riphammer_mae_mean=gn_results['riphammer']['mae_mean'],
         gn_riphammer_mae_std=gn_results['riphammer']['mae_std'],
         gn_riphammer_sep=gn_results['riphammer']['sep_mean'],
         gn_stdnmf_mae_mean=gn_results['stdnmf']['mae_mean'],
         gn_stdnmf_mae_std=gn_results['stdnmf']['mae_std'],
         gn_stdnmf_sep=gn_results['stdnmf']['sep_mean'],
         diff_corr_clean_riphammer=diff_corr_results['clean']['riphammer'],
         diff_corr_clean_stdnmf=diff_corr_results['clean']['stdnmf'],
         diff_corr_ds50_riphammer=diff_corr_results['ds50']['riphammer'],
         diff_corr_ds50_stdnmf=diff_corr_results['ds50']['stdnmf'],
)
print("\n[saved] noise_robustness.npz")

# ============================================
# Figure
# ============================================
fig, axes = plt.subplots(2, 2, figsize=(14, 11))

# (a) Downsampling vs MAE
ax = axes[0, 0]
x_ds = [f*100 for f in ds_fractions]
ax.errorbar(x_ds, ds_results['riphammer']['mae_mean'], yerr=ds_results['riphammer']['mae_std'],
            marker='o', color='#2196F3', lw=2, capsize=4, label='RIPHAMMER')
ax.errorbar(x_ds, ds_results['stdnmf']['mae_mean'], yerr=ds_results['stdnmf']['mae_std'],
            marker='s', color='#4CAF50', lw=2, capsize=4, label='Standard NMF')
ax.set_xlabel('Downsampling fraction (%)', fontsize=11)
ax.set_ylabel('Proportion MAE', fontsize=11)
ax.set_title('(a) Robustness to downsampling', fontweight='bold', fontsize=12)
ax.legend(fontsize=10)
ax.grid(True, alpha=0.2)

# (b) Gaussian noise vs MAE
ax = axes[0, 1]
ax.errorbar(noise_fractions, gn_results['riphammer']['mae_mean'], yerr=gn_results['riphammer']['mae_std'],
            marker='o', color='#2196F3', lw=2, capsize=4, label='RIPHAMMER')
ax.errorbar(noise_fractions, gn_results['stdnmf']['mae_mean'], yerr=gn_results['stdnmf']['mae_std'],
            marker='s', color='#4CAF50', lw=2, capsize=4, label='Standard NMF')
ax.set_xlabel('Noise level (σ / mean)', fontsize=11)
ax.set_ylabel('Proportion MAE', fontsize=11)
ax.set_title('(b) Robustness to Gaussian noise', fontweight='bold', fontsize=12)
ax.legend(fontsize=10)
ax.grid(True, alpha=0.2)

# (c) Pure sample separation vs downsampling
ax = axes[1, 0]
ax.plot(x_ds, ds_results['riphammer']['sep_mean'], 'o-', color='#2196F3', lw=2, markersize=8, label='RIPHAMMER')
ax.plot(x_ds, ds_results['stdnmf']['sep_mean'], 's-', color='#4CAF50', lw=2, markersize=8, label='Standard NMF')
ax.set_xlabel('Downsampling fraction (%)', fontsize=11)
ax.set_ylabel('Unique dominant cell types (out of 4)', fontsize=11)
ax.set_title('(c) Cell-type separation under sparsity', fontweight='bold', fontsize=12)
ax.set_ylim(0.5, 4.5)
ax.set_yticks([1, 2, 3, 4])
ax.legend(fontsize=10)
ax.grid(True, alpha=0.2)

# (d) Differential contact at Aire: clean vs noisy
ax = axes[1, 1]
x_pos = np.arange(2)
width = 0.3
r_vals = [diff_corr_results['clean']['riphammer'], diff_corr_results['ds50']['riphammer']]
s_vals = [diff_corr_results['clean']['stdnmf'], diff_corr_results['ds50']['stdnmf']]
bars1 = ax.bar(x_pos - width/2, r_vals, width, color='#2196F3', alpha=0.85, label='RIPHAMMER')
bars2 = ax.bar(x_pos + width/2, s_vals, width, color='#4CAF50', alpha=0.85, label='Standard NMF')
for bars in [bars1, bars2]:
    for bar in bars:
        h = bar.get_height()
        y = h + 0.02 if h > 0 else h - 0.06
        ax.text(bar.get_x() + bar.get_width()/2, y, f'{h:.3f}',
                ha='center', fontsize=10, fontweight='bold')
ax.set_xticks(x_pos)
ax.set_xticklabels(['Clean (0%)', '50% downsampled'], fontsize=11)
ax.set_ylabel('Diff. contact correlation (Aire locus)', fontsize=10)
ax.set_title('(d) Biological signal under sparsity', fontweight='bold', fontsize=12)
ax.legend(fontsize=10)
ax.axhline(0, color='gray', alpha=0.3, ls='--')
ax.set_ylim(min(min(r_vals), min(s_vals)) - 0.15, max(max(r_vals), max(s_vals)) + 0.15)

fig.suptitle('Noise and Downsampling Robustness Analysis', fontsize=14, fontweight='bold', y=1.0)
plt.tight_layout()
plt.savefig('plot/noise_robustness.png', dpi=150, bbox_inches='tight')
plt.close()
print("[saved] plot/noise_robustness.png")
