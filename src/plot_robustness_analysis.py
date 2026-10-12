#!/usr/bin/env python3
"""K selection + noise/downsampling robustness (6 panels).
NOTE: This script re-runs NMF fitting (may take several minutes).

Run from repo root: python src/plot_robustness_analysis.py
Output: plot/robustness_analysis.png
"""
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sys
from riphammer.core import MultibandNMF
from sklearn.decomposition import NMF
from scipy.stats import pearsonr

X = np.load('output/simulation_100000/X_bulk.npy')
feat = pd.read_csv('output/simulation_100000/features.tsv', sep='\t', header=None)
feat.columns = ['chrom', 'bin1', 'bin2']
distances = np.abs(feat['bin2'].values - feat['bin1'].values)
true_frac = np.loadtxt('output/simulation_100000/true_fractions.tsv', delimiter='\t', skiprows=1)

sample_names = ['CD4_SP', 'DN', 'DP', 'Treg_GFP',
                'equal', 'CD4_dom', 'CD4_DN_dom', 'Treg_rare', 'gradient', 'extreme']

# Use 50% subsample for speed
rng = np.random.RandomState(42)
n_half = len(X) // 2
idx_half = np.sort(rng.choice(len(X), n_half, replace=False))
X_sub = X[idx_half]
d_sub = distances[idx_half]

def map_and_eval(H, true_frac):
    """Map components to cell types and compute MAE on mixed samples."""
    K = H.shape[0]
    H_norm = H / H.sum(axis=0, keepdims=True)
    
    # Pure sample separation
    used = set()
    comp_map = {}
    for i in range(4):
        dom = np.argmax(H_norm[:, i])
        if dom not in used:
            comp_map[dom] = i
            used.add(dom)
    n_unique = len(set(comp_map.values()))
    
    if n_unique < 4 or K < 4:
        return n_unique, float('nan')
    
    col_order = []
    used2 = set()
    for ct in range(4):
        for c, m in comp_map.items():
            if m == ct and c not in used2:
                col_order.append(c); used2.add(c); break
        else:
            cands = [(c, H_norm[c, ct]) for c in range(K) if c not in used2]
            if cands:
                best = max(cands, key=lambda x: x[1])[0]
                col_order.append(best); used2.add(best)
            else:
                return n_unique, float('nan')
    
    pred = H_norm[col_order, 4:10].T  # (6, 4)
    mae = np.mean(np.abs(true_frac - pred))
    return n_unique, mae

# ============================================
# K robustness (on 50% data for speed)
# ============================================
print("=== K Robustness ===")
ks = [2, 3, 4, 5, 6, 8]
k_results = {'k': [], 'n_unique': [], 'mae': [], 'recon_err': [], 'stability': []}

for K in ks:
    print(f"  K={K}...", end=' ', flush=True)
    
    # Main run
    m = MultibandNMF(n_components=K, max_iter=500, random_state=42)
    m.fit(X_sub, d_sub, verbose=False)
    W = m.get_W_full()
    H = m.H_
    
    recon = W @ H
    err = np.linalg.norm(X_sub - recon) / np.linalg.norm(X_sub)
    n_uniq, mae = map_and_eval(H, true_frac)
    
    # Stability: 3 seeds
    Hs = []
    for seed in [42, 123, 456]:
        m2 = MultibandNMF(n_components=K, max_iter=300, random_state=seed)
        m2.fit(X_sub, d_sub, verbose=False)
        H2 = m2.H_ / m2.H_.sum(axis=0, keepdims=True)
        Hs.append(H2)
    
    # Pairwise H correlation (match columns by max correlation)
    corrs = []
    for i in range(len(Hs)):
        for j in range(i+1, len(Hs)):
            # For each column of Hs[i], find best matching column in Hs[j]
            best_corrs = []
            for ci in range(K):
                max_r = max(abs(pearsonr(Hs[i][ci], Hs[j][cj])[0]) for cj in range(K))
                best_corrs.append(max_r)
            corrs.append(np.mean(best_corrs))
    stability = np.mean(corrs)
    
    k_results['k'].append(K)
    k_results['n_unique'].append(n_uniq)
    k_results['mae'].append(mae)
    k_results['recon_err'].append(err)
    k_results['stability'].append(stability)
    
    mae_str = f'{mae:.4f}' if not np.isnan(mae) else 'N/A'
    print(f"err={err:.4f}, unique={n_uniq}/4, MAE={mae_str}, stability={stability:.3f}")

# ============================================
# Noise robustness (on 50% data)
# ============================================
print("\n=== Noise Robustness ===")

# Downsampling
ds_fracs = [0.0, 0.1, 0.2, 0.3, 0.5, 0.7]
ds_results = {'frac': [], 'mae_riph': [], 'mae_std': [], 'mae_riph_std': [], 'mae_std_std': [],
              'sep_riph': [], 'sep_std': []}

nz_mask = X_sub > 0

for frac in ds_fracs:
    print(f"  Downsampling {frac:.0%}...", flush=True)
    maes_r, maes_s, seps_r, seps_s = [], [], [], []
    
    for seed in [42, 123, 456]:
        rng2 = np.random.RandomState(seed)
        X_ds = X_sub.copy()
        if frac > 0:
            drop_mask = rng2.random(X_ds.shape) < frac
            X_ds[nz_mask & drop_mask] = 0
        
        # RIPHAMMER
        m = MultibandNMF(n_components=4, max_iter=500, random_state=42)
        m.fit(X_ds, d_sub, verbose=False)
        sep_r, mae_r = map_and_eval(m.H_, true_frac)
        maes_r.append(mae_r); seps_r.append(sep_r)
        
        # StdNMF
        m2 = NMF(n_components=4, max_iter=500, random_state=42)
        m2.fit_transform(np.maximum(X_ds, 0))
        sep_s, mae_s = map_and_eval(m2.components_, true_frac)
        maes_s.append(mae_s); seps_s.append(sep_s)
    
    ds_results['frac'].append(frac)
    ds_results['mae_riph'].append(np.nanmean(maes_r))
    ds_results['mae_std'].append(np.nanmean(maes_s))
    ds_results['mae_riph_std'].append(np.nanstd(maes_r))
    ds_results['mae_std_std'].append(np.nanstd(maes_s))
    ds_results['sep_riph'].append(np.mean(seps_r))
    ds_results['sep_std'].append(np.mean(seps_s))
    print(f"    RIPHAMMER: MAE={np.nanmean(maes_r):.4f}±{np.nanstd(maes_r):.4f}, sep={np.mean(seps_r):.1f}/4")
    print(f"    StdNMF:    MAE={np.nanmean(maes_s):.4f}±{np.nanstd(maes_s):.4f}, sep={np.mean(seps_s):.1f}/4")

# Gaussian noise
noise_levels = [0.0, 0.1, 0.2, 0.5, 1.0]
noise_results = {'level': [], 'mae_riph': [], 'mae_std': []}
mean_nz = X_sub[X_sub > 0].mean()

for nl in noise_levels:
    print(f"  Noise σ={nl:.1f}×mean...", flush=True)
    X_noisy = X_sub + np.random.RandomState(42).randn(*X_sub.shape) * nl * mean_nz
    X_noisy = np.maximum(X_noisy, 0)
    
    m = MultibandNMF(n_components=4, max_iter=500, random_state=42)
    m.fit(X_noisy, d_sub, verbose=False)
    _, mae_r = map_and_eval(m.H_, true_frac)
    
    m2 = NMF(n_components=4, max_iter=500, random_state=42)
    m2.fit_transform(X_noisy)
    _, mae_s = map_and_eval(m2.components_, true_frac)
    
    noise_results['level'].append(nl)
    noise_results['mae_riph'].append(mae_r)
    noise_results['mae_std'].append(mae_s)
    print(f"    RIPHAMMER MAE={mae_r:.4f}, StdNMF MAE={mae_s:.4f}")

# ============================================
# Combined Figure: K robustness + Noise robustness
# ============================================
fig, axes = plt.subplots(2, 3, figsize=(18, 11))

# Row 1: K robustness
ax = axes[0, 0]
ax.plot(k_results['k'], k_results['recon_err'], 'o-', color='#2196F3', lw=2, ms=8)
ax.set_xlabel('K (number of components)')
ax.set_ylabel('Relative reconstruction error')
ax.set_title('(a) Reconstruction error vs K', fontweight='bold')
ax.grid(True, alpha=0.2)
ax.axvline(4, color='red', ls='--', alpha=0.5, label='True K=4')
ax.legend()

ax = axes[0, 1]
ax.bar(k_results['k'], k_results['n_unique'], color='#2196F3', alpha=0.8)
ax.axhline(4, color='red', ls='--', alpha=0.5)
ax.set_xlabel('K')
ax.set_ylabel('Unique cell types identified (out of 4)')
ax.set_title('(b) Cell type separation vs K', fontweight='bold')
ax.set_ylim(0, 5)

ax = axes[0, 2]
ax.plot(k_results['k'], k_results['stability'], 'o-', color='#2196F3', lw=2, ms=8)
ax.set_xlabel('K')
ax.set_ylabel('Mean pairwise H correlation')
ax.set_title('(c) Stability across random seeds', fontweight='bold')
ax.set_ylim(0.5, 1.05)
ax.grid(True, alpha=0.2)
ax.axvline(4, color='red', ls='--', alpha=0.5, label='True K=4')
ax.legend()

# Row 2: Noise robustness
ax = axes[1, 0]
ax.errorbar(ds_results['frac'], ds_results['mae_riph'], yerr=ds_results['mae_riph_std'],
            fmt='o-', color='#2196F3', lw=2, ms=8, capsize=4, label='RIPHAMMER')
ax.errorbar(ds_results['frac'], ds_results['mae_std'], yerr=ds_results['mae_std_std'],
            fmt='s-', color='#4CAF50', lw=2, ms=8, capsize=4, label='Standard NMF')
ax.set_xlabel('Downsampling fraction')
ax.set_ylabel('MAE')
ax.set_title('(d) Robustness to downsampling', fontweight='bold')
ax.legend()
ax.grid(True, alpha=0.2)

ax = axes[1, 1]
ax.plot(noise_results['level'], noise_results['mae_riph'], 'o-', color='#2196F3', lw=2, ms=8, label='RIPHAMMER')
ax.plot(noise_results['level'], noise_results['mae_std'], 's-', color='#4CAF50', lw=2, ms=8, label='Standard NMF')
ax.set_xlabel('Noise level (σ / mean)')
ax.set_ylabel('MAE')
ax.set_title('(e) Robustness to Gaussian noise', fontweight='bold')
ax.legend()
ax.grid(True, alpha=0.2)

ax = axes[1, 2]
ax.plot(ds_results['frac'], ds_results['sep_riph'], 'o-', color='#2196F3', lw=2, ms=8, label='RIPHAMMER')
ax.plot(ds_results['frac'], ds_results['sep_std'], 's-', color='#4CAF50', lw=2, ms=8, label='Standard NMF')
ax.set_xlabel('Downsampling fraction')
ax.set_ylabel('Unique cell types (out of 4)')
ax.set_title('(f) Separation under downsampling', fontweight='bold')
ax.set_ylim(0, 5)
ax.legend()
ax.grid(True, alpha=0.2)

plt.suptitle('Robustness Analysis: K selection and noise tolerance',
             fontsize=14, fontweight='bold', y=1.01)
plt.tight_layout()
plt.savefig('plot/robustness_analysis.png', dpi=150, bbox_inches='tight')
plt.close()

np.savez('output/simulation_100000/robustness_results.npz',
         k_vals=k_results['k'], k_recon=k_results['recon_err'],
         k_unique=k_results['n_unique'], k_stability=k_results['stability'],
         k_mae=k_results['mae'],
         ds_frac=ds_results['frac'], ds_mae_riph=ds_results['mae_riph'],
         ds_mae_std=ds_results['mae_std'], ds_sep_riph=ds_results['sep_riph'],
         ds_sep_std=ds_results['sep_std'],
         noise_level=noise_results['level'], noise_mae_riph=noise_results['mae_riph'],
         noise_mae_std=noise_results['mae_std'])

print("\n[saved] plot/robustness_analysis.png")
print("[saved] output/simulation_100000/robustness_results.npz")
