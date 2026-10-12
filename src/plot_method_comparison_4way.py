#!/usr/bin/env python3
"""4-method comparison: RIPHAMMER vs THUNDER vs Standard NMF vs PCA.

Run from repo root: python src/plot_method_comparison_4way.py
Output: plot/method_comparison_4way.png
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
pred_riphammer = d['pred_riphammer']  # (6, 4)
pred_stdnmf = d['pred_stdnmf']  # (6, 4)
pred_pca = d['pred_pca']  # (6, 4)

# Load THUNDER
thunder_H = pd.read_csv('output/simulation_100000/thunder_H_final.tsv', sep='\t')
thunder_mat = thunder_H.values  # 10 x 4

# Pure sample separation (4/4 test)
pure_block = thunder_mat[:4, :]
comp_to_celltype = {}
for i in range(4):
    dom_comp = np.argmax(pure_block[i, :])
    comp_to_celltype[dom_comp] = i

thunder_col_order = []
for ct_idx in range(4):
    for comp_idx, mapped_ct in comp_to_celltype.items():
        if mapped_ct == ct_idx:
            thunder_col_order.append(comp_idx)
            break

pred_thunder = thunder_mat[4:, thunder_col_order]  # 6 x 4

# Check pure sample uniqueness for each method
# For RIPHAMMER, Std NMF, PCA - we need to check what happened with pure samples
# These pred arrays only have 6 mixed samples, pure sample check was done separately
# From context: RIPHAMMER 4/4, Std NMF 2/4, PCA 3/4, THUNDER 4/4

cell_types = ['CD4_SP', 'DN', 'DP', 'Treg_GFP']
mix_names = ['equal', 'CD4_dom', 'CD4_DN_dom', 'Treg_rare', 'gradient', 'extreme']

# Compute metrics
def compute_metrics(true, pred):
    mae = np.mean(np.abs(true - pred))
    rmse = np.sqrt(np.mean((true - pred)**2))
    # Per-cell-type correlation (across all mixing scenarios)
    rs = []
    for ct in range(true.shape[1]):
        if np.std(true[:, ct]) > 1e-8 and np.std(pred[:, ct]) > 1e-8:
            r, _ = pearsonr(true[:, ct], pred[:, ct])
            rs.append(r)
    mean_r = np.mean(rs) if rs else float('nan')
    return mae, rmse, mean_r

mae_r, rmse_r, r_r = compute_metrics(true_frac, pred_riphammer)
mae_t, rmse_t, r_t = compute_metrics(true_frac, pred_thunder)
mae_s, rmse_s, r_s = compute_metrics(true_frac, pred_stdnmf)
mae_p, rmse_p, r_p = compute_metrics(true_frac, pred_pca)

# Pure sample separation: unique dominant types in pure samples
pure_sep = {
    'RIPHAMMER': '4/4',
    'THUNDER': '4/4',
    'Standard NMF': '2/4',
    'PCA': '3/4',
}

print("=== Method Comparison Summary ===\n")
summary = pd.DataFrame({
    'Method': ['RIPHAMMER', 'THUNDER', 'Standard NMF', 'PCA'],
    'Pure Sep': ['4/4', '4/4', '2/4', '3/4'],
    'MAE': [mae_r, mae_t, mae_s, mae_p],
    'RMSE': [rmse_r, rmse_t, rmse_s, rmse_p],
    'Mean r': [r_r, r_t, r_s, r_p],
})
print(summary.to_string(index=False))

# Detailed per-scenario comparison
print("\n\n=== RIPHAMMER predicted vs true ===")
for i, mn in enumerate(mix_names):
    err = np.abs(true_frac[i] - pred_riphammer[i])
    print(f"  {mn}: pred={np.round(pred_riphammer[i],3)} true={np.round(true_frac[i],2)} MAE={err.mean():.4f}")

print("\n=== THUNDER predicted vs true ===")
for i, mn in enumerate(mix_names):
    err = np.abs(true_frac[i] - pred_thunder[i])
    print(f"  {mn}: pred={np.round(pred_thunder[i],3)} true={np.round(true_frac[i],2)} MAE={err.mean():.4f}")

# ====== Figure ======
fig = plt.figure(figsize=(22, 14))
gs = fig.add_gridspec(3, 4, hspace=0.4, wspace=0.3)

methods_data = [
    ('RIPHAMMER', pred_riphammer, mae_r, r_r, '4/4', '#2196F3'),
    ('THUNDER', pred_thunder, mae_t, r_t, '4/4', '#FF9800'),
    ('Standard NMF\n(ablation)', pred_stdnmf, mae_s, r_s, '2/4', '#4CAF50'),
    ('PCA\n(baseline)', pred_pca, mae_p, r_p, '3/4', '#9C27B0'),
]

ct_colors = ['#e74c3c', '#3498db', '#2ecc71', '#f39c12']
ct_markers = ['o', 's', 'D', '^']

# Row 1: Scatter (predicted vs true fraction)
for idx, (name, pred, mae, mr, ps, color) in enumerate(methods_data):
    ax = fig.add_subplot(gs[0, idx])
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

# Row 2: Stacked bars
for idx, (name, pred, mae, mr, ps, color) in enumerate(methods_data):
    ax = fig.add_subplot(gs[1, idx])
    x = np.arange(len(mix_names))
    width = 0.35
    
    # Predicted (left bar)
    bottom_p = np.zeros(len(mix_names))
    for ct_idx in range(4):
        ax.bar(x - width/2, pred[:, ct_idx], width, bottom=bottom_p,
               color=ct_colors[ct_idx], alpha=0.85, 
               label=cell_types[ct_idx] if idx==0 else None)
        bottom_p += pred[:, ct_idx]
    
    # True (right bar)  
    bottom_t = np.zeros(len(mix_names))
    for ct_idx in range(4):
        ax.bar(x + width/2, true_frac[:, ct_idx], width, bottom=bottom_t,
               color=ct_colors[ct_idx], alpha=0.3, edgecolor='black', linewidth=0.5)
        bottom_t += true_frac[:, ct_idx]
    
    ax.set_xticks(x)
    ax.set_xticklabels(mix_names, rotation=45, ha='right', fontsize=8)
    ax.set_ylim(0, 1.05)
    if idx == 0:
        ax.set_ylabel('Fraction')
    ax.set_title('Pred (solid) vs True (faded)', fontsize=9)

# Row 3: Summary bar chart
ax_bar = fig.add_subplot(gs[2, :2])
method_names = [m[0].replace('\n', ' ') for m in methods_data]
x = np.arange(len(method_names))
mae_vals = [m[2] for m in methods_data]
bar_colors = [m[5] for m in methods_data]
bars = ax_bar.bar(x, mae_vals, color=bar_colors, alpha=0.8, edgecolor='white', linewidth=2)
ax_bar.set_xticks(x)
ax_bar.set_xticklabels(method_names, fontsize=10)
ax_bar.set_ylabel('MAE (lower = better)')
ax_bar.set_title('Mean Absolute Error', fontweight='bold')
for bar, val in zip(bars, mae_vals):
    ax_bar.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.003,
                f'{val:.3f}', ha='center', fontsize=10, fontweight='bold')
ax_bar.set_ylim(0, max(mae_vals) * 1.2)

# Per-scenario MAE heatmap
ax_heat = fig.add_subplot(gs[2, 2:])
mae_matrix = np.array([
    [np.mean(np.abs(true_frac[i] - pred)) for pred in 
     [pred_riphammer[i], pred_thunder[i], pred_stdnmf[i], pred_pca[i]]]
    for i in range(6)
])
im = ax_heat.imshow(mae_matrix, cmap='YlOrRd', aspect='auto')
ax_heat.set_xticks(range(4))
ax_heat.set_xticklabels([m[0].replace('\n', ' ') for m in methods_data], fontsize=9)
ax_heat.set_yticks(range(6))
ax_heat.set_yticklabels(mix_names, fontsize=9)
for i in range(6):
    for j in range(4):
        ax_heat.text(j, i, f'{mae_matrix[i,j]:.3f}', ha='center', va='center', fontsize=9,
                     color='white' if mae_matrix[i,j] > 0.1 else 'black')
plt.colorbar(im, ax=ax_heat, label='MAE')
ax_heat.set_title('Per-scenario MAE', fontweight='bold')

fig.suptitle('Simulation Benchmark: 4 Methods × 6 Mixing Scenarios',
             fontsize=15, fontweight='bold', y=0.99)
plt.savefig('plot/method_comparison_4way.png', dpi=150, bbox_inches='tight')
plt.close()
print("\n[saved] plot/method_comparison_4way.png")
