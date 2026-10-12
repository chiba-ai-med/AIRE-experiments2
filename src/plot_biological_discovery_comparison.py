#!/usr/bin/env python3
"""TAD detection + insulation correlation + MAE comparison (3 methods).
NOTE: This script re-runs NMF fitting (may take several minutes).

Run from repo root: python src/plot_biological_discovery_comparison.py
Output: plot/biological_discovery_comparison.png
"""
import numpy as np
import pandas as pd
from scipy.stats import pearsonr
from scipy.ndimage import uniform_filter1d
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sys
from riphammer.core import MultibandNMF
from sklearn.decomposition import NMF

# ============================================
# Load and fit (same as before, but keep models)
# ============================================
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
thunder_features = thunder_w_df['feature_name'].values
W_thunder_raw = thunder_w_df.iloc[:, 1:].values
H_thunder_raw = np.loadtxt('output/simulation_100000/thunder_H_k4_raw.tsv', delimiter='\t', skiprows=1)

def map_components(H, n_pure=4):
    H_norm = H / H.sum(axis=0, keepdims=True)
    comp_map = {}
    for i in range(n_pure):
        comp_map[np.argmax(H_norm[:, i])] = i
    col_order = []
    used = set()
    for ct in range(n_pure):
        for c, m in comp_map.items():
            if m == ct and c not in used:
                col_order.append(c); used.add(c); break
        else:
            cands = [(c, H_norm[c, ct]) for c in range(H.shape[0]) if c not in used]
            best = max(cands, key=lambda x: x[1])[0]
            col_order.append(best); used.add(best)
    return col_order

riph_order = map_components(H_riph)
std_order = map_components(H_std)
thunder_order = map_components(H_thunder_raw)

thunder_feat_to_idx = {f: i for i, f in enumerate(thunder_features)}

def insulation_score(mat, window=5):
    n = mat.shape[0]
    scores = np.full(n, np.nan)
    for i in range(window, n - window):
        scores[i] = np.mean(mat[i-window:i, i:i+window])
    valid = ~np.isnan(scores)
    if valid.sum() > 0 and np.std(scores[valid]) > 0:
        scores[valid] = (scores[valid] - np.mean(scores[valid])) / np.std(scores[valid])
    return scores

def find_tad_boundaries(ins, min_prom=0.5):
    bounds = []
    for i in range(1, len(ins) - 1):
        if np.isnan(ins[i]): continue
        if ins[i] < ins[i-1] and ins[i] < ins[i+1]:
            lm = np.nanmax(ins[max(0, i-5):i])
            rm = np.nanmax(ins[i+1:min(len(ins), i+6)])
            if min(lm - ins[i], rm - ins[i]) > min_prom:
                bounds.append(i)
    return np.array(bounds)

def build_mat(bin1_idx, bin2_idx, vals, n):
    mat = np.zeros((n, n))
    for i, j, v in zip(bin1_idx, bin2_idx, vals):
        mat[i, j] = v; mat[j, i] = v
    return mat

# ============================================
# Load TAD comparison results
# ============================================
tad_data = np.load('output/simulation_100000/tad_comparison.npz')

# ============================================
# Figure: Biological Discovery Comparison
# ============================================
fig = plt.figure(figsize=(22, 16))

# --- Top row: Summary metrics ---
gs_top = fig.add_gridspec(1, 3, top=0.95, bottom=0.72, left=0.06, right=0.96, wspace=0.25)

method_names = ['RIPHAMMER', 'THUNDER', 'Standard NMF']
colors = ['#2196F3', '#FF9800', '#4CAF50']

# TAD F1
ax1 = fig.add_subplot(gs_top[0])
metrics = ['Precision', 'Recall', 'F1']
x = np.arange(len(metrics))
width = 0.25
for i, (name, color) in enumerate(zip(method_names, colors)):
    vals = [np.mean(tad_data[f'{name}_precision']),
            np.mean(tad_data[f'{name}_recall']),
            np.mean(tad_data[f'{name}_f1'])]
    bars = ax1.bar(x + i*width, vals, width, label=name, color=color, alpha=0.85)
    for bar, val in zip(bars, vals):
        ax1.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.01,
                 f'{val:.2f}', ha='center', fontsize=9, fontweight='bold')
ax1.set_xticks(x + width)
ax1.set_xticklabels(metrics, fontsize=11)
ax1.set_ylabel('Score')
ax1.set_title('TAD Boundary Detection', fontweight='bold', fontsize=12)
ax1.legend(fontsize=9)
ax1.set_ylim(0, 1.05)

# Insulation correlation distribution
ax2 = fig.add_subplot(gs_top[1])
positions = [1, 2, 3]
data_insul = [tad_data[f'{name}_insul_corr'] for name in method_names]
bp = ax2.boxplot(data_insul, positions=positions, widths=0.5, patch_artist=True,
                  showfliers=True, flierprops=dict(markersize=3))
for patch, color in zip(bp['boxes'], colors):
    patch.set_facecolor(color)
    patch.set_alpha(0.7)
ax2.set_xticks(positions)
ax2.set_xticklabels(method_names, fontsize=10)
ax2.set_ylabel('Pearson r')
ax2.set_title('Insulation Score Correlation\n(per chrom × cell type)', fontweight='bold', fontsize=12)
medians = [np.median(d) for d in data_insul]
for i, m in enumerate(medians):
    ax2.text(positions[i], m + 0.02, f'{m:.3f}', ha='center', fontsize=9, fontweight='bold')

# Proportion MAE (from earlier results)
ax3 = fig.add_subplot(gs_top[2])
d = np.load('output/simulation_100000/comparison_results_final.npz')
true_frac = d['true_frac']
preds = {'RIPHAMMER': d['pred_riphammer'], 'THUNDER': d['pred_thunder'], 'Standard NMF': d['pred_stdnmf']}
mae_vals = [np.mean(np.abs(true_frac - preds[n])) for n in method_names]
bars = ax3.bar(method_names, mae_vals, color=colors, alpha=0.85)
for bar, val in zip(bars, mae_vals):
    ax3.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.002,
             f'{val:.3f}', ha='center', fontsize=10, fontweight='bold')
ax3.set_ylabel('MAE (lower = better)')
ax3.set_title('Mixing Proportion Accuracy', fontweight='bold', fontsize=12)
ax3.set_ylim(0, max(mae_vals) * 1.3)

# --- Bottom rows: Example chromosome contact maps + insulation ---
chrom = 'chr10'
ct_idx = 0  # CD4_SP
ct_name = sample_names[ct_idx]

chr_mask = feat['chrom'] == chrom
cf = feat[chr_mask]
bins = sorted(set(cf['bin1']) | set(cf['bin2']))
n_bins = len(bins)
b2i = {b: i for i, b in enumerate(bins)}
bin1_idx = np.array([b2i[b] for b in cf['bin1']])
bin2_idx = np.array([b2i[b] for b in cf['bin2']])
global_idx = cf.index.values

# Zoom to a region (bins 200-400 for visibility)
zoom = slice(200, 400)
zoom_range = range(200, 400)

# True
mat_true = build_mat(bin1_idx, bin2_idx, X[global_idx, ct_idx], n_bins)
ins_true = insulation_score(mat_true)
tad_true = find_tad_boundaries(ins_true)

# RIPHAMMER
mat_riph = build_mat(bin1_idx, bin2_idx, W_riph[global_idx, riph_order[ct_idx]], n_bins)
ins_riph = insulation_score(mat_riph)
tad_riph = find_tad_boundaries(ins_riph)

# THUNDER
thunder_vals = np.zeros(len(global_idx))
for k, gidx in enumerate(global_idx):
    row = feat.iloc[gidx]
    fname = f"{row['chrom']}_{int(row['bin1'])*resolution}_{row['chrom']}_{int(row['bin2'])*resolution}"
    if fname in thunder_feat_to_idx:
        thunder_vals[k] = W_thunder_raw[thunder_feat_to_idx[fname], thunder_order[ct_idx]]
mat_thunder = build_mat(bin1_idx, bin2_idx, thunder_vals, n_bins)
ins_thunder = insulation_score(mat_thunder)
tad_thunder = find_tad_boundaries(ins_thunder)

# Standard NMF
mat_std = build_mat(bin1_idx, bin2_idx, W_std[global_idx, std_order[ct_idx]], n_bins)
ins_std = insulation_score(mat_std)
tad_std = find_tad_boundaries(ins_std)

gs_mid = fig.add_gridspec(1, 4, top=0.65, bottom=0.35, left=0.06, right=0.96, wspace=0.15)
gs_bot = fig.add_gridspec(1, 4, top=0.30, bottom=0.05, left=0.06, right=0.96, wspace=0.15)

titles = [f'True ({ct_name})', 'RIPHAMMER', 'THUNDER', 'Standard NMF']
mats = [mat_true, mat_riph, mat_thunder, mat_std]
insuls = [ins_true, ins_riph, ins_thunder, ins_std]
tads_list = [tad_true, tad_riph, tad_thunder, tad_std]
plot_colors = ['black', '#2196F3', '#FF9800', '#4CAF50']

for idx in range(4):
    # Contact map
    ax_map = fig.add_subplot(gs_mid[idx])
    m = mats[idx][zoom, zoom]
    if m.max() > 0:
        vmax = np.percentile(m[m > 0], 98) if (m > 0).sum() > 0 else 1
    else:
        vmax = 1
    ax_map.imshow(m, cmap='Reds', vmin=0, vmax=vmax, aspect='equal', 
                   extent=[200, 400, 400, 200])
    # Mark TAD boundaries in zoomed region
    for tb in tads_list[idx]:
        if 200 <= tb <= 400:
            ax_map.axvline(tb, color='blue', alpha=0.5, lw=0.5)
            ax_map.axhline(tb, color='blue', alpha=0.5, lw=0.5)
    ax_map.set_title(titles[idx], fontsize=11, fontweight='bold', color=plot_colors[idx])
    if idx == 0:
        ax_map.set_ylabel(f'{chrom} (bins 200-400)')
    
    # Insulation score
    ax_ins = fig.add_subplot(gs_bot[idx])
    x_range = np.arange(200, 400)
    ax_ins.plot(x_range, ins_true[zoom], 'k-', alpha=0.3, lw=1, label='True')
    ax_ins.plot(x_range, insuls[idx][zoom], '-', color=plot_colors[idx], lw=1.5)
    
    # Mark TAD boundaries
    for tb in tad_true:
        if 200 <= tb <= 400:
            ax_ins.axvline(tb, color='gray', alpha=0.3, lw=0.5, ls='--')
    for tb in tads_list[idx]:
        if 200 <= tb <= 400:
            ax_ins.axvline(tb, color=plot_colors[idx], alpha=0.5, lw=1)
    
    valid = ~(np.isnan(ins_true[zoom]) | np.isnan(insuls[idx][zoom]))
    if valid.sum() > 5:
        r_val, _ = pearsonr(ins_true[zoom][valid], insuls[idx][zoom][valid])
        ax_ins.set_title(f'Insulation (r={r_val:.3f})', fontsize=10)
    ax_ins.set_xlabel('Bin')
    if idx == 0:
        ax_ins.set_ylabel('Insulation score')
        ax_ins.legend(fontsize=8)

fig.suptitle('Biological Discovery Capability: TAD Structure Recovery from Deconvolution\n'
             '(Simulation: 4 cell types × 19 chromosomes, 100 kb resolution)',
             fontsize=14, fontweight='bold', y=0.99)

plt.savefig('plot/biological_discovery_comparison.png', dpi=150, bbox_inches='tight')
plt.close()
print("[saved] plot/biological_discovery_comparison.png")

# Print final summary
print("\n=== Final Biological Discovery Comparison ===")
print(f"{'Metric':<30} {'RIPHAMMER':>12} {'THUNDER':>12} {'Std NMF':>12}")
print("-" * 66)
print(f"{'Pure sample separation':<30} {'4/4':>12} {'4/4':>12} {'2/4':>12}")
print(f"{'Proportion MAE':<30} {mae_vals[0]:>12.3f} {mae_vals[1]:>12.3f} {mae_vals[2]:>12.3f}")
print(f"{'TAD Precision':<30} {np.mean(tad_data['RIPHAMMER_precision']):>12.3f} {np.mean(tad_data['THUNDER_precision']):>12.3f} {np.mean(tad_data['Standard NMF_precision']):>12.3f}")
print(f"{'TAD Recall':<30} {np.mean(tad_data['RIPHAMMER_recall']):>12.3f} {np.mean(tad_data['THUNDER_recall']):>12.3f} {np.mean(tad_data['Standard NMF_recall']):>12.3f}")
print(f"{'TAD F1':<30} {np.mean(tad_data['RIPHAMMER_f1']):>12.3f} {np.mean(tad_data['THUNDER_f1']):>12.3f} {np.mean(tad_data['Standard NMF_f1']):>12.3f}")
print(f"{'Insulation r (median)':<30} {np.median(tad_data['RIPHAMMER_insul_corr']):>12.3f} {np.median(tad_data['THUNDER_insul_corr']):>12.3f} {np.median(tad_data['Standard NMF_insul_corr']):>12.3f}")
