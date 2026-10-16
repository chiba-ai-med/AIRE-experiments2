#!/usr/bin/env python3
"""4-method comparison: RIPHAMMER vs THUNDER vs Standard NMF vs PCA.

Run from repo root: python src/plot_method_comparison_4way.py
Input: output/simulation_100000/{comparison_results.npz, thunder_H_fullgenome.tsv}
Output: plot/method_comparison_4way.png, output/simulation_100000/comparison_results_final.npz
"""
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.stats import pearsonr

# Load existing results
d = np.load('output/simulation_100000/comparison_results.npz', allow_pickle=True)
true_frac = d['true_frac']  # (6, 4)
pred_riphammer = d['pred_riphammer']
pred_stdnmf = d['pred_stdnmf']
pred_pca = d['pred_pca']

cell_types = ['CD4_SP', 'DN', 'DP', 'Treg_GFP']
mix_names = ['equal', 'CD4_dom', 'CD4_DN_dom', 'Treg_rare', 'gradient', 'extreme']

# Load full-genome THUNDER
thunder_H = pd.read_csv('output/simulation_100000/thunder_H_fullgenome.tsv', sep='\t')
thunder_mat = thunder_H.values  # 10 x 4

# Map components to cell types via pure samples
pure_block = thunder_mat[:4, :]
comp_to_ct = {}
for i in range(4):
    comp_to_ct[np.argmax(pure_block[i, :])] = i

col_order = [next(c for c, ct in comp_to_ct.items() if ct == i) for i in range(4)]
pred_thunder = thunder_mat[4:, col_order]

print("THUNDER full genome (reordered):")
print(pd.DataFrame(pred_thunder, columns=cell_types, index=mix_names).round(3))

# Metrics
def compute_metrics(true, pred):
    mae = np.mean(np.abs(true - pred))
    rmse = np.sqrt(np.mean((true - pred)**2))
    rs = []
    for ct in range(true.shape[1]):
        if np.std(true[:, ct]) > 1e-8 and np.std(pred[:, ct]) > 1e-8:
            r, _ = pearsonr(true[:, ct], pred[:, ct])
            rs.append(r)
    return mae, rmse, np.mean(rs) if rs else float('nan')

metrics = {
    'RIPHAMMER': compute_metrics(true_frac, pred_riphammer),
    'THUNDER': compute_metrics(true_frac, pred_thunder),
    'Standard NMF': compute_metrics(true_frac, pred_stdnmf),
    'PCA': compute_metrics(true_frac, pred_pca),
}

print("\n=== Final Comparison (full genome) ===")
print(f"{'Method':<15} {'Pure Sep':>8} {'MAE':>8} {'RMSE':>8} {'Mean r':>8}")
pure_seps = {'RIPHAMMER': '4/4', 'THUNDER': '4/4', 'Standard NMF': '2/4', 'PCA': '3/4'}
for name, (mae, rmse, mr) in metrics.items():
    print(f"{name:<15} {pure_seps[name]:>8} {mae:>8.4f} {rmse:>8.4f} {mr:>8.3f}")

# Save updated results
np.savez('output/simulation_100000/comparison_results_final.npz',
         pred_riphammer=pred_riphammer, pred_stdnmf=pred_stdnmf,
         pred_pca=pred_pca, pred_thunder=pred_thunder,
         true_frac=true_frac)

# ====== Final comparison figure ======
fig = plt.figure(figsize=(20, 12))
gs = fig.add_gridspec(2, 4, hspace=0.35, wspace=0.3)

methods_data = [
    ('RIPHAMMER', pred_riphammer, '4/4', '#2196F3'),
    ('THUNDER', pred_thunder, '4/4', '#FF9800'),
    ('Standard NMF\n(ablation)', pred_stdnmf, '2/4', '#4CAF50'),
    ('PCA\n(baseline)', pred_pca, '3/4', '#9C27B0'),
]
ct_colors = ['#e74c3c', '#3498db', '#2ecc71', '#f39c12']
ct_markers = ['o', 's', 'D', '^']

# Row 1: Scatter plots
for idx, (name, pred, ps, color) in enumerate(methods_data):
    ax = fig.add_subplot(gs[0, idx])
    mae, rmse, mr = metrics[name.split('\n')[0].strip('()')]

    for ct_idx, ct in enumerate(cell_types):
        ax.scatter(true_frac[:, ct_idx], pred[:, ct_idx],
                   s=100, alpha=0.8, label=ct, color=ct_colors[ct_idx],
                   marker=ct_markers[ct_idx], zorder=5, edgecolors='white', linewidth=0.5)
    ax.plot([0, 1], [0, 1], 'k--', alpha=0.3, lw=1)
    ax.set_xlim(-0.05, 1.05)
    ax.set_ylim(-0.05, 1.05)
    ax.set_xlabel('True fraction')
    if idx == 0:
        ax.set_ylabel('Predicted fraction')
    r_str = f'{mr:.3f}' if not np.isnan(mr) else 'N/A'
    ax.set_title(f'{name}\nMAE={mae:.3f}  r={r_str}\nPure separation: {ps}',
                 fontsize=11, fontweight='bold')
    ax.set_aspect('equal')
    ax.grid(True, alpha=0.2)
    if idx == 0:
        ax.legend(fontsize=8, loc='upper left')

# Row 2: Per-scenario MAE heatmap (left) + bar chart (right)
ax_heat = fig.add_subplot(gs[1, :2])
method_labels = ['RIPHAMMER', 'THUNDER', 'Std NMF', 'PCA']
preds = [pred_riphammer, pred_thunder, pred_stdnmf, pred_pca]
mae_matrix = np.array([
    [np.mean(np.abs(true_frac[i] - pred[i])) for pred in preds]
    for i in range(6)
])
im = ax_heat.imshow(mae_matrix, cmap='YlOrRd', aspect='auto', vmin=0, vmax=0.25)
ax_heat.set_xticks(range(4))
ax_heat.set_xticklabels(method_labels, fontsize=10, fontweight='bold')
ax_heat.set_yticks(range(6))
ax_heat.set_yticklabels(mix_names, fontsize=10)
for i in range(6):
    for j in range(4):
        ax_heat.text(j, i, f'{mae_matrix[i,j]:.3f}', ha='center', va='center', fontsize=10,
                     color='white' if mae_matrix[i,j] > 0.12 else 'black', fontweight='bold')
plt.colorbar(im, ax=ax_heat, label='MAE', shrink=0.8)
ax_heat.set_title('Per-scenario MAE (lower = better)', fontweight='bold', fontsize=12)

# Overall MAE bar
ax_bar = fig.add_subplot(gs[1, 2])
mae_vals = [metrics[m.split('\n')[0].strip('()')][0] for m, _, _, _ in methods_data]
bar_colors = [m[3] for m in methods_data]
bars = ax_bar.barh(range(3, -1, -1), mae_vals, color=bar_colors, alpha=0.85, edgecolor='white', height=0.6)
ax_bar.set_yticks(range(3, -1, -1))
ax_bar.set_yticklabels(method_labels, fontsize=10, fontweight='bold')
for i, (bar, val) in enumerate(zip(bars, mae_vals)):
    ax_bar.text(val + 0.003, bar.get_y() + bar.get_height()/2,
                f'{val:.3f}', va='center', fontsize=11, fontweight='bold')
ax_bar.set_xlabel('MAE')
ax_bar.set_title('Overall MAE', fontweight='bold', fontsize=12)
ax_bar.set_xlim(0, max(mae_vals) * 1.3)

# Correlation bar
ax_r = fig.add_subplot(gs[1, 3])
r_vals = [metrics[m.split('\n')[0].strip('()')][2] for m, _, _, _ in methods_data]
bars_r = ax_r.barh(range(3, -1, -1), r_vals, color=bar_colors, alpha=0.85, edgecolor='white', height=0.6)
ax_r.set_yticks(range(3, -1, -1))
ax_r.set_yticklabels(method_labels, fontsize=10, fontweight='bold')
for i, (bar, val) in enumerate(zip(bars_r, r_vals)):
    ax_r.text(val + 0.01 if val < 0.9 else val - 0.08, bar.get_y() + bar.get_height()/2,
              f'{val:.3f}', va='center', fontsize=11, fontweight='bold')
ax_r.set_xlabel('Mean Pearson r')
ax_r.set_title('Correlation with truth', fontweight='bold', fontsize=12)
ax_r.set_xlim(0, 1.1)

fig.suptitle('Simulation Benchmark: 4 cell types × 6 mixing scenarios (100 kb resolution)',
             fontsize=14, fontweight='bold', y=0.99)
plt.savefig('plot/method_comparison_4way.png', dpi=150, bbox_inches='tight')
plt.close()
print("\n[saved] plot/method_comparison_4way.png")
