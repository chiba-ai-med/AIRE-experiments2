#!/usr/bin/env python3
"""RIPHAMMER mixing fraction recovery scatter plot (4 cell types × 6 scenarios).

Run from repo root: python src/plot_simulation_fraction_recovery.py
Output: plot/simulation_fraction_recovery.png
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

H = np.loadtxt("output/nmf_simulation_100000/H_nmf_k4.tsv")
true_frac = np.loadtxt("output/simulation_100000/true_fractions.tsv", delimiter="\t", skiprows=1)

sample_names = ["CD4_SP", "DN", "DP", "Treg_GFP", 
                "equal", "CD4_dom", "CD4_DN_dom", "Treg_rare", "gradient", "extreme"]
mix_names = ["equal", "CD4_dom", "CD4_DN_dom", "Treg_rare", "gradient", "extreme"]
cell_names = ["CD4_SP", "DN", "DP", "Treg_GFP"]

# H columns for pure samples (first 4) and mixed (last 6)
H_pure = H[:, :4]   # (4 comps × 4 pure samples)
H_mix = H[:, 4:]    # (4 comps × 6 mixed samples)

# Predicted fraction: solve H_pure @ frac ≈ h_mix for each mixed sample
# Using least squares: frac = (H_pure^T H_pure)^{-1} H_pure^T h_mix
# With non-negativity constraint
from scipy.optimize import nnls

pred_fracs = []
for i in range(6):
    h = H_mix[:, i]
    frac, _ = nnls(H_pure, h)
    frac = frac / frac.sum()  # normalize to sum=1
    pred_fracs.append(frac)
pred_fracs = np.array(pred_fracs)

print("=== Mixing fraction recovery ===")
print(f"{'Scenario':<15} {'Cell type':<10} {'True':>6} {'Pred':>6} {'Error':>7}")
print("-" * 50)
for i, name in enumerate(mix_names):
    for j, cell in enumerate(cell_names):
        err = pred_fracs[i, j] - true_frac[i, j]
        print(f"{name:<15} {cell:<10} {true_frac[i,j]:6.3f} {pred_fracs[i,j]:6.3f} {err:+7.3f}")
    print()

# Overall metrics
errors = pred_fracs - true_frac
mae = np.mean(np.abs(errors))
rmse = np.sqrt(np.mean(errors**2))
# Per-scenario correlation
from scipy.stats import pearsonr
cors = []
for i in range(6):
    r, _ = pearsonr(true_frac[i], pred_fracs[i])
    cors.append(r)

print(f"Overall MAE:  {mae:.4f}")
print(f"Overall RMSE: {rmse:.4f}")
print(f"Per-scenario Pearson r: {[f'{r:.3f}' for r in cors]}")
print(f"Mean r: {np.mean(cors):.4f}")

# Plot
fig, axes = plt.subplots(1, 2, figsize=(14, 6))

# Scatter plot: true vs predicted
ax = axes[0]
colors = ["red", "green", "blue", "orange"]
for j, (cell, color) in enumerate(zip(cell_names, colors)):
    ax.scatter(true_frac[:, j], pred_fracs[:, j], c=color, s=80, label=cell, edgecolors="black", linewidth=0.5)
ax.plot([0, 1], [0, 1], "k--", alpha=0.5)
ax.set_xlabel("True fraction", fontsize=12)
ax.set_ylabel("Predicted fraction", fontsize=12)
ax.set_title(f"Mixing fraction recovery\nMAE={mae:.3f}, RMSE={rmse:.3f}", fontsize=13)
ax.legend(fontsize=10)
ax.set_xlim(-0.05, 1.05)
ax.set_ylim(-0.05, 1.05)

# Heatmap comparison
ax = axes[1]
combined = np.zeros((6, 8))
combined[:, :4] = true_frac
combined[:, 4:] = pred_fracs
im = ax.imshow(combined, cmap="YlOrRd", vmin=0, vmax=1, aspect="auto")
ax.set_xticks(range(8))
ax.set_xticklabels([f"T:{c}" for c in cell_names] + [f"P:{c}" for c in cell_names], 
                    fontsize=8, rotation=45, ha="right")
ax.set_yticks(range(6))
ax.set_yticklabels(mix_names, fontsize=9)
ax.set_title("True (T) vs Predicted (P)", fontsize=13)
for i in range(6):
    for j in range(8):
        ax.text(j, i, f"{combined[i,j]:.2f}", ha="center", va="center", fontsize=7)
ax.axvline(3.5, color="white", linewidth=2)
plt.colorbar(im, ax=ax, shrink=0.8)

plt.tight_layout()
plt.savefig("plot/simulation_fraction_recovery.png", dpi=200, bbox_inches="tight")
print("\n[saved] plot/simulation_fraction_recovery.png")
