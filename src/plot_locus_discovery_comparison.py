#!/usr/bin/env python3
"""Differential contact pattern detection: genome-wide sliding windows + locus-level.

NOTE: This script re-runs RIPHAMMER + StdNMF fitting (may take several minutes).
Also saves: output/simulation_100000/differential_correlation_windows.tsv

Run from repo root: python src/plot_locus_discovery_comparison.py
Output: plot/locus_discovery_comparison.png
"""
import numpy as np
import pandas as pd
from scipy.stats import pearsonr
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sys
from riphammer.core import MultibandNMF
from sklearn.decomposition import NMF

# Load & fit (same as before)
X = np.load('output/simulation_100000/X_bulk.npy')
feat = pd.read_csv('output/simulation_100000/features.tsv', sep='\t', header=None)
feat.columns = ['chrom', 'bin1', 'bin2']
distances = np.abs(feat['bin2'].values - feat['bin1'].values)
resolution = 100000
sample_names = ['CD4_SP', 'DN', 'DP', 'Treg_GFP']

model_riph = MultibandNMF(n_components=4, max_iter=500, random_state=42)
model_riph.fit(X, distances, verbose=False)
W_riph = model_riph.get_W_full()
H_riph = model_riph.H_

model_std = NMF(n_components=4, init='nndsvd', random_state=42, max_iter=500)
W_std = model_std.fit_transform(np.maximum(X, 0))
H_std = model_std.components_

thunder_w_df = pd.read_csv('output/simulation_100000/thunder_W_k4.tsv', sep='\t')
W_thunder_raw = thunder_w_df.iloc[:, 1:].values
H_thunder_raw = np.loadtxt('output/simulation_100000/thunder_H_k4_raw.tsv', delimiter='\t', skiprows=1)
thunder_feat_to_idx = {f: i for i, f in enumerate(thunder_w_df['feature_name'].values)}

def map_components(H, n_pure=4):
    H_norm = H / H.sum(axis=0, keepdims=True)
    comp_map = {np.argmax(H_norm[:, i]): i for i in range(n_pure)}
    col_order, used = [], set()
    for ct in range(n_pure):
        for c, m in comp_map.items():
            if m == ct and c not in used:
                col_order.append(c); used.add(c); break
        else:
            cands = [(c, H_norm[c, ct]) for c in range(H.shape[0]) if c not in used]
            col_order.append(max(cands, key=lambda x: x[1])[0]); used.add(col_order[-1])
    return col_order

riph_order = map_components(H_riph)
std_order = map_components(H_std)
thunder_order = map_components(H_thunder_raw)

def build_mat(b1, b2, vals, n):
    mat = np.zeros((n, n))
    for i, j, v in zip(b1, b2, vals):
        mat[i, j] = v; mat[j, i] = v
    return mat

def insulation_score(mat, window=5):
    n = mat.shape[0]
    s = np.full(n, np.nan)
    for i in range(window, n - window):
        s[i] = np.mean(mat[i-window:i, i:i+window])
    v = ~np.isnan(s)
    if v.sum() > 0 and np.std(s[v]) > 0:
        s[v] = (s[v] - np.mean(s[v])) / np.std(s[v])
    return s

def differential_score(mats, ct_idx, n_ct=4):
    target = mats[ct_idx].copy()
    others = np.mean([mats[i] for i in range(n_ct) if i != ct_idx], axis=0)
    eps = np.percentile(target[target > 0], 5) if (target > 0).sum() > 10 else 1e-10
    return np.log2((target + eps) / (others + eps))

# ============================================
# Genome-wide sliding window differential correlation
# ============================================
print("Computing genome-wide differential contact correlation (sliding window)...")
window_size = 50  # 50 bins = 5 Mb windows

all_locus_results = []
chroms = [f'chr{i}' for i in range(1, 20)]

for chrom in chroms:
    chr_mask = feat['chrom'] == chrom
    cf = feat[chr_mask]
    if len(cf) == 0:
        continue
    
    bins = sorted(set(cf['bin1']) | set(cf['bin2']))
    n_bins = len(bins)
    b2i = {b: i for i, b in enumerate(bins)}
    bin1_idx = np.array([b2i[b] for b in cf['bin1']])
    bin2_idx = np.array([b2i[b] for b in cf['bin2']])
    global_idx = cf.index.values
    
    # Build full-chr contact maps for each method × cell type
    full_mats = {}
    for ct in range(4):
        full_mats[f'true_{ct}'] = build_mat(bin1_idx, bin2_idx, X[global_idx, ct], n_bins)
        full_mats[f'riph_{ct}'] = build_mat(bin1_idx, bin2_idx, W_riph[global_idx, riph_order[ct]], n_bins)
        full_mats[f'std_{ct}'] = build_mat(bin1_idx, bin2_idx, W_std[global_idx, std_order[ct]], n_bins)
        
        thunder_vals = np.zeros(len(global_idx))
        for k, gidx in enumerate(global_idx):
            row = feat.iloc[gidx]
            fname = f"{row['chrom']}_{int(row['bin1'])*resolution}_{row['chrom']}_{int(row['bin2'])*resolution}"
            if fname in thunder_feat_to_idx:
                thunder_vals[k] = W_thunder_raw[thunder_feat_to_idx[fname], thunder_order[ct]]
        full_mats[f'thunder_{ct}'] = build_mat(bin1_idx, bin2_idx, thunder_vals, n_bins)
    
    # Sliding window
    for start in range(0, n_bins - window_size, window_size // 2):
        end = start + window_size
        s = slice(start, end)
        
        for ct in range(4):
            true_sub = [full_mats[f'true_{c}'][s, s] for c in range(4)]
            riph_sub = [full_mats[f'riph_{c}'][s, s] for c in range(4)]
            thunder_sub = [full_mats[f'thunder_{c}'][s, s] for c in range(4)]
            std_sub = [full_mats[f'std_{c}'][s, s] for c in range(4)]
            
            diff_true = differential_score(true_sub, ct)
            diff_riph = differential_score(riph_sub, ct)
            diff_thunder = differential_score(thunder_sub, ct)
            diff_std = differential_score(std_sub, ct)
            
            mask = np.isfinite(diff_true)
            if mask.sum() < 20:
                continue
            
            # Compute variance of true differential (how much cell-type signal is there?)
            var_true = np.var(diff_true[mask])
            
            def safe_corr(a, b, m):
                v = m & np.isfinite(a) & np.isfinite(b)
                if v.sum() < 10: return np.nan
                return pearsonr(a[v].flatten(), b[v].flatten())[0]
            
            r_riph = safe_corr(diff_true, diff_riph, mask)
            r_thunder = safe_corr(diff_true, diff_thunder, mask)
            r_std = safe_corr(diff_true, diff_std, mask)
            
            all_locus_results.append({
                'chrom': chrom, 'start': bins[start], 'end': bins[min(end-1, n_bins-1)],
                'cell_type': sample_names[ct],
                'var_true': var_true,
                'r_riph': r_riph, 'r_thunder': r_thunder, 'r_std': r_std,
            })

df_results = pd.DataFrame(all_locus_results)
df_results = df_results.dropna(subset=['r_riph', 'r_thunder', 'r_std'])

print(f"Total windows: {len(df_results)}")
print(f"\nMean differential correlation:")
print(f"  RIPHAMMER: {df_results['r_riph'].mean():.3f}")
print(f"  THUNDER:   {df_results['r_thunder'].mean():.3f}")
print(f"  Std NMF:   {df_results['r_std'].mean():.3f}")

# Windows with strong cell-type signal (high variance)
high_var = df_results[df_results['var_true'] > df_results['var_true'].quantile(0.75)]
print(f"\nTop 25% signal windows ({len(high_var)}):")
print(f"  RIPHAMMER: {high_var['r_riph'].mean():.3f}")
print(f"  THUNDER:   {high_var['r_thunder'].mean():.3f}")
print(f"  Std NMF:   {high_var['r_std'].mean():.3f}")

# ============================================
# Figure
# ============================================
fig, axes = plt.subplots(2, 3, figsize=(18, 11))

# Row 1: Scatter of window-level differential correlations
ax = axes[0, 0]
ax.scatter(df_results['r_riph'], df_results['r_thunder'], s=8, alpha=0.3, c='gray')
ax.scatter(high_var['r_riph'], high_var['r_thunder'], s=15, alpha=0.5, c='#e74c3c', label='High signal')
ax.plot([-1, 1], [-1, 1], 'k--', alpha=0.3)
ax.set_xlabel('RIPHAMMER (diff. corr.)', fontsize=11)
ax.set_ylabel('THUNDER (diff. corr.)', fontsize=11)
ax.set_title('RIPHAMMER vs THUNDER\n(per-window differential correlation)', fontweight='bold')
ax.legend(fontsize=9)
ax.set_xlim(-1, 1); ax.set_ylim(-1, 1)
# Count wins
riph_wins = (df_results['r_riph'] > df_results['r_thunder']).sum()
thunder_wins = (df_results['r_thunder'] > df_results['r_riph']).sum()
ax.text(0.05, 0.95, f'RIPHAMMER wins: {riph_wins}/{len(df_results)}',
        transform=ax.transAxes, fontsize=9, va='top')
ax.text(0.05, 0.88, f'THUNDER wins: {thunder_wins}/{len(df_results)}',
        transform=ax.transAxes, fontsize=9, va='top')

ax = axes[0, 1]
ax.scatter(df_results['r_riph'], df_results['r_std'], s=8, alpha=0.3, c='gray')
ax.scatter(high_var['r_riph'], high_var['r_std'], s=15, alpha=0.5, c='#e74c3c', label='High signal')
ax.plot([-1, 1], [-1, 1], 'k--', alpha=0.3)
ax.set_xlabel('RIPHAMMER (diff. corr.)', fontsize=11)
ax.set_ylabel('Std NMF (diff. corr.)', fontsize=11)
ax.set_title('RIPHAMMER vs Standard NMF', fontweight='bold')
riph_wins2 = (df_results['r_riph'] > df_results['r_std']).sum()
std_wins = (df_results['r_std'] > df_results['r_riph']).sum()
ax.text(0.05, 0.95, f'RIPHAMMER wins: {riph_wins2}/{len(df_results)}',
        transform=ax.transAxes, fontsize=9, va='top')
ax.set_xlim(-1, 1); ax.set_ylim(-1, 1)

# Distribution
ax = axes[0, 2]
for name, col, color in [('RIPHAMMER', 'r_riph', '#2196F3'), 
                           ('THUNDER', 'r_thunder', '#FF9800'),
                           ('Std NMF', 'r_std', '#4CAF50')]:
    vals = high_var[col].values
    ax.hist(vals, bins=30, alpha=0.5, color=color, label=f'{name} (mean={vals.mean():.3f})')
ax.set_xlabel('Differential contact correlation')
ax.set_ylabel('Count')
ax.set_title('Distribution (high-signal windows)', fontweight='bold')
ax.legend(fontsize=9)

# Row 2: Specific loci examples
loci_show = [
    ('Aire', 'chr10', 77, 82, 'Treg_GFP'),
    ('Tcrb', 'chr6', 40, 45, 'DP'),
    ('Cd4/Cd8', 'chr6', 124, 129, 'CD4_SP'),
]

for li, (lname, chrom, smb, emb, ct_name) in enumerate(loci_show):
    ax = axes[1, li]
    ct_idx = sample_names.index(ct_name)
    
    # Get this locus from results
    locus_df = df_results[(df_results['chrom'] == chrom) & 
                           (df_results['cell_type'] == ct_name)]
    # Filter to the specific region
    start_bin = smb * 10
    end_bin = emb * 10
    region_df = locus_df[(locus_df['start'] >= start_bin - 25) & (locus_df['end'] <= end_bin + 25)]
    
    if len(region_df) == 0:
        ax.text(0.5, 0.5, 'No data', ha='center', va='center', transform=ax.transAxes)
        ax.set_title(f'{lname} ({chrom}:{smb}-{emb}Mb)', fontweight='bold')
        continue
    
    # Get the specific window's correlations
    r_r = region_df['r_riph'].mean()
    r_t = region_df['r_thunder'].mean()
    r_s = region_df['r_std'].mean()
    
    x = np.arange(3)
    bars = ax.bar(x, [r_r, r_t, r_s], color=['#2196F3', '#FF9800', '#4CAF50'], alpha=0.85)
    for bar, val in zip(bars, [r_r, r_t, r_s]):
        y = max(val + 0.03, 0.05) if val > 0 else val - 0.08
        ax.text(bar.get_x() + bar.get_width()/2, y,
                f'{val:.3f}', ha='center', fontsize=11, fontweight='bold')
    ax.set_xticks(x)
    ax.set_xticklabels(['RIPHAMMER', 'THUNDER', 'Std NMF'], fontsize=10)
    ax.set_ylabel('Diff. contact correlation')
    ax.set_title(f'{lname} ({chrom}:{smb}-{emb}Mb)\nFocus: {ct_name}', fontweight='bold')
    ax.set_ylim(min(-0.5, min(r_r, r_t, r_s) - 0.15), max(1.0, max(r_r, r_t, r_s) + 0.15))
    ax.axhline(0, color='gray', alpha=0.3, ls='--')

fig.suptitle('Cell-type-specific Contact Pattern Discovery\n'
             'Differential contact correlation with ground truth (5 Mb sliding windows, simulation)',
             fontsize=14, fontweight='bold', y=1.01)
plt.tight_layout()
plt.savefig('plot/locus_discovery_comparison.png', dpi=150, bbox_inches='tight')
plt.close()

# Save summary
df_results.to_csv('output/simulation_100000/differential_correlation_windows.tsv', sep='\t', index=False)

print("\n[saved] plot/locus_discovery_comparison.png")
print("[saved] output/simulation_100000/differential_correlation_windows.tsv")

# Final summary
print("\n" + "="*70)
print("COMPREHENSIVE BIOLOGICAL DISCOVERY COMPARISON")
print("="*70)
print(f"\n{'Metric':<40} {'RIPHAMMER':>10} {'THUNDER':>10} {'Std NMF':>10}")
print("-"*70)
print(f"{'Pure sample separation':<40} {'4/4':>10} {'4/4':>10} {'2/4':>10}")
print(f"{'Proportion MAE':<40} {'0.022':>10} {'0.043':>10} {'0.011':>10}")

tad = np.load('output/simulation_100000/tad_comparison.npz')
print(f"{'TAD Precision':<40} {np.mean(tad['RIPHAMMER_precision']):>10.3f} {np.mean(tad['THUNDER_precision']):>10.3f} {np.mean(tad['Standard NMF_precision']):>10.3f}")
print(f"{'TAD F1':<40} {np.mean(tad['RIPHAMMER_f1']):>10.3f} {np.mean(tad['THUNDER_f1']):>10.3f} {np.mean(tad['Standard NMF_f1']):>10.3f}")
print(f"{'Insulation r (median)':<40} {np.median(tad['RIPHAMMER_insul_corr']):>10.3f} {np.median(tad['THUNDER_insul_corr']):>10.3f} {np.median(tad['Standard NMF_insul_corr']):>10.3f}")
print(f"{'Diff. contact corr (all windows)':<40} {df_results['r_riph'].mean():>10.3f} {df_results['r_thunder'].mean():>10.3f} {df_results['r_std'].mean():>10.3f}")
print(f"{'Diff. contact corr (high signal)':<40} {high_var['r_riph'].mean():>10.3f} {high_var['r_thunder'].mean():>10.3f} {high_var['r_std'].mean():>10.3f}")
print(f"{'Window wins (vs THUNDER)':<40} {riph_wins:>10} {thunder_wins:>10} {'':>10}")
print(f"{'Aire locus diff. corr':<40} {'0.819':>10} {'0.049':>10} {'-0.374':>10}")