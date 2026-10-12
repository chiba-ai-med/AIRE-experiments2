#!/usr/bin/env python3
"""Computational scalability: data size × K × samples.
NOTE: This script re-runs NMF fitting (may take several minutes).

Run from repo root: python src/plot_scalability_benchmark.py
Output: plot/scalability_benchmark.png
"""
import numpy as np
import pandas as pd
import time
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sys
from riphammer.core import MultibandNMF
from sklearn.decomposition import NMF

X = np.load('output/simulation_100000/X_bulk.npy')
feat = pd.read_csv('output/simulation_100000/features.tsv', sep='\t', header=None)
feat.columns = ['chrom', 'bin1', 'bin2']
distances = np.abs(feat['bin2'].values - feat['bin1'].values)

# ============================================
# 1. Data size scaling (subsample features)
# ============================================
fracs = [0.1, 0.25, 0.5, 0.75, 1.0]
size_results = {'frac': [], 'n_features': [], 'riphammer': [], 'stdnmf': []}

for frac in fracs:
    n = int(len(X) * frac)
    rng = np.random.RandomState(42)
    idx = rng.choice(len(X), n, replace=False)
    idx.sort()
    X_sub = X[idx]
    d_sub = distances[idx]
    
    # RIPHAMMER
    times_r = []
    for seed in [42, 123, 456]:
        t0 = time.time()
        m = MultibandNMF(n_components=4, max_iter=200, random_state=seed)
        m.fit(X_sub, d_sub, verbose=False)
        times_r.append(time.time() - t0)
    
    # Std NMF
    times_s = []
    for seed in [42, 123, 456]:
        t0 = time.time()
        m = NMF(n_components=4, max_iter=200, random_state=seed)
        m.fit_transform(np.maximum(X_sub, 0))
        times_s.append(time.time() - t0)
    
    size_results['frac'].append(frac)
    size_results['n_features'].append(n)
    size_results['riphammer'].append(np.median(times_r))
    size_results['stdnmf'].append(np.median(times_s))
    print(f"  frac={frac:.0%} ({n:,} features): RIPHAMMER {np.median(times_r):.1f}s, StdNMF {np.median(times_s):.1f}s")

# ============================================
# 2. K scaling
# ============================================
ks = [2, 4, 6, 8]
k_results = {'k': [], 'riphammer': [], 'stdnmf': []}

# Use 50% subsample for speed
n_half = len(X) // 2
rng = np.random.RandomState(42)
idx_half = rng.choice(len(X), n_half, replace=False)
idx_half.sort()
X_half = X[idx_half]
d_half = distances[idx_half]

for k in ks:
    t0 = time.time()
    MultibandNMF(n_components=k, max_iter=200, random_state=42).fit(X_half, d_half, verbose=False)
    tr = time.time() - t0
    
    t0 = time.time()
    NMF(n_components=k, max_iter=200, random_state=42).fit_transform(np.maximum(X_half, 0))
    ts = time.time() - t0
    
    k_results['k'].append(k)
    k_results['riphammer'].append(tr)
    k_results['stdnmf'].append(ts)
    print(f"  K={k}: RIPHAMMER {tr:.1f}s, StdNMF {ts:.1f}s")

# ============================================
# 3. Sample scaling
# ============================================
sample_counts = [4, 6, 8, 10]
sample_results = {'n_samples': [], 'riphammer': [], 'stdnmf': []}

for ns in sample_counts:
    X_s = X_half[:, :ns]
    
    t0 = time.time()
    MultibandNMF(n_components=4, max_iter=200, random_state=42).fit(X_s, d_half, verbose=False)
    tr = time.time() - t0
    
    t0 = time.time()
    NMF(n_components=4, max_iter=200, random_state=42).fit_transform(np.maximum(X_s, 0))
    ts = time.time() - t0
    
    sample_results['n_samples'].append(ns)
    sample_results['riphammer'].append(tr)
    sample_results['stdnmf'].append(ts)
    print(f"  N_samples={ns}: RIPHAMMER {tr:.1f}s, StdNMF {ts:.1f}s")

# THUNDER reference times (from actual runs)
thunder_nrun1_time = 24.4 * 60 / 10  # ~146s per single run at 1.17M features
thunder_full_time = 24.4 * 60  # 24.4 min for nrun=10

# ============================================
# Figure
# ============================================
fig, axes = plt.subplots(1, 3, figsize=(18, 5.5))

# (a) Data size
ax = axes[0]
ax.plot([r/1e6 for r in size_results['n_features']], size_results['riphammer'], 
        'o-', color='#2196F3', lw=2, ms=8, label='RIPHAMMER')
ax.plot([r/1e6 for r in size_results['n_features']], size_results['stdnmf'],
        's-', color='#4CAF50', lw=2, ms=8, label='Standard NMF')
# THUNDER reference point
ax.axhline(thunder_nrun1_time, color='#FF9800', ls='--', lw=1.5, alpha=0.7, label=f'THUNDER (nrun=1, est.)')
ax.set_xlabel('Number of features (millions)', fontsize=11)
ax.set_ylabel('Runtime (seconds)', fontsize=11)
ax.set_title('(a) Scaling with data size\n(K=4, max_iter=200)', fontweight='bold')
ax.legend(fontsize=9)
ax.grid(True, alpha=0.2)

# (b) K
ax = axes[1]
ax.plot(k_results['k'], k_results['riphammer'], 'o-', color='#2196F3', lw=2, ms=8, label='RIPHAMMER')
ax.plot(k_results['k'], k_results['stdnmf'], 's-', color='#4CAF50', lw=2, ms=8, label='Standard NMF')
ax.set_xlabel('Number of components (K)', fontsize=11)
ax.set_ylabel('Runtime (seconds)', fontsize=11)
ax.set_title('(b) Scaling with K\n(616K features, max_iter=200)', fontweight='bold')
ax.legend(fontsize=9)
ax.grid(True, alpha=0.2)

# (c) Samples
ax = axes[2]
ax.plot(sample_results['n_samples'], sample_results['riphammer'], 'o-', color='#2196F3', lw=2, ms=8, label='RIPHAMMER')
ax.plot(sample_results['n_samples'], sample_results['stdnmf'], 's-', color='#4CAF50', lw=2, ms=8, label='Standard NMF')
ax.set_xlabel('Number of samples', fontsize=11)
ax.set_ylabel('Runtime (seconds)', fontsize=11)
ax.set_title('(c) Scaling with samples\n(616K features, K=4, max_iter=200)', fontweight='bold')
ax.legend(fontsize=9)
ax.grid(True, alpha=0.2)

plt.suptitle('Computational Scalability', fontsize=14, fontweight='bold', y=1.02)
plt.tight_layout()
plt.savefig('plot/scalability_benchmark.png', dpi=150, bbox_inches='tight')
plt.close()

# Save
np.savez('output/simulation_100000/scalability_results.npz',
         size_frac=fracs, size_nfeat=size_results['n_features'],
         size_riphammer=size_results['riphammer'], size_stdnmf=size_results['stdnmf'],
         k_vals=k_results['k'], k_riphammer=k_results['riphammer'], k_stdnmf=k_results['stdnmf'],
         sample_n=sample_results['n_samples'], sample_riphammer=sample_results['riphammer'],
         sample_stdnmf=sample_results['stdnmf'])

print("\n[saved] plot/scalability_benchmark.png")
print("[saved] output/simulation_100000/scalability_results.npz")
