# Reproducibility Guide for RIPHAMMER Plots

All 16 plots in `plot/` were generated during interactive analysis sessions.
This document maps each plot to its source script and input data.

## Script Inventory

All scripts are in `src/`. The `plot_*.py` scripts were extracted from the
conversation transcript and saved as standalone scripts. They should be run
from the repository root directory.

### Existing pipeline scripts (Snakemake-managed)
| Script | Purpose |
|---|---|
| `vectorize.py` | Convert mcool files to feature matrices |
| `plot_contact_maps.py` | Basic NMF component contact map visualization |

### THUNDER R scripts
| Script | Purpose | Runtime |
|---|---|---|
| `run_thunder_fullgenome.R` | Full THUNDER pipeline (3-step: NMF → feature selection → NMF) | ~28 min |
| `run_thunder_extract_W.R` | Extract W matrix via single-run NMF for comparison plots | ~5 min |

### Simulation benchmark plots
| Script | Output | Computation | Input Data |
|---|---|---|---|
| `plot_simulation_fraction_recovery.py` | simulation_fraction_recovery.png | Light (NNLS only) | nmf_simulation_100000/H_nmf_k4.tsv, true_fractions.tsv |
| `plot_method_comparison_simulation.py` | method_comparison_simulation.png | **Runs NMF + PCA** | X_bulk.npy, features.tsv, true_fractions.tsv |
| `plot_method_comparison_4way.py` | method_comparison_4way.png | Light (reads pre-computed) | comparison_results.npz, thunder_H_fullgenome.tsv |
| `plot_biological_discovery_comparison.py` | biological_discovery_comparison.png | **Runs RIPHAMMER + StdNMF** | X_bulk.npy, features.tsv, tad_comparison.npz, thunder_W_k4.tsv |
| `plot_locus_discovery_comparison.py` | locus_discovery_comparison.png | **Runs RIPHAMMER + StdNMF** | X_bulk.npy, features.tsv, thunder_W_k4.tsv |
| `plot_robustness_analysis.py` | robustness_analysis.png | **Heavy: multiple NMF runs** | X_bulk.npy, features.tsv, true_fractions.tsv |
| `plot_scalability_benchmark.py` | scalability_benchmark.png | **Timing benchmarks** | X_bulk.npy, features.tsv |
| `plot_noise_robustness.py` | noise_robustness.png | **Heavy: multiple NMF runs with noise** | X_bulk.npy, features.tsv, true_fractions.tsv |

### mTEC Aire locus plots (real data, 10kb)
| Script | Output | Input Data |
|---|---|---|
| `plot_aire_collaborator_summary.py` | aire_collaborator_summary.png | nmf_10000/W,H,cells,features |
| `plot_aire_locus_10kb_components.py` | aire_locus_10kb_components.png | nmf_10000/*, mtec_hic.mcool |
| `plot_aire_locus_10kb_deep.py` | aire_locus_10kb_deep.png | nmf_10000/*, tad_10000_k8/*, mtec_hic.mcool |
| `plot_aire_locus_10kb_selected.py` | aire_locus_10kb_selected.png | nmf_10000/*, mtec_hic.mcool |
| `plot_aire_locus_10kb_differential.py` | aire_locus_10kb_differential.png | nmf_10000/* |
| `plot_aire_locus_10kb_insulation.py` | aire_locus_10kb_insulation.png | nmf_10000/*, tad_10000_k8/* |

### Brain plots (real data, 100kb)
| Script | Output | Input Data |
|---|---|---|
| `plot_brain_dev_k2_100kb.py` | brain_dev_k2_100kb.png | nmf_brain_100000_k2/*, brain_hic_*.mcool |
| `plot_brain_dev_loci_100kb.py` | brain_dev_loci_100kb.png | nmf_brain_100000/*, brain_hic_*.mcool |

## Component-to-Cell-Type Mapping

### Simulation (4 cell types: CD4_SP, DN, DP, Treg_GFP)

All methods use pure-sample dominance mapping:

```python
# H matrix: (K components × N samples)
# Columns 0-3 are pure samples, columns 4-9 are mixtures
H_norm = H / H.sum(axis=0, keepdims=True)
for i in range(4):  # for each pure sample
    dominant_component = np.argmax(H_norm[:, i])
    comp_to_celltype[dominant_component] = i
```

"Pure sample separation" metric counts how many of the 4 pure samples
map to unique components. RIPHAMMER achieves 4/4 (perfect separation);
Standard NMF achieves only 2/4 because distance decay dominates.

When separation < 4/4 (e.g., under noise), MAE cannot be computed and
is reported as NaN. This is a limitation, not a bug.

### THUNDER mapping (same approach)

THUNDER's H matrix rows are NMF components. The pure sample columns
(first 4) each have ~1.0 in exactly one component. Column reordering:

```python
# thunder_mat: (10 samples × 4 components)
pure_block = thunder_mat[:4, :]
for i in range(4):
    comp_to_celltype[np.argmax(pure_block[i, :])] = i
```

### Real data (mTEC, 10kb, K=8)

No ground truth — mapping is based on H loading values and sample metadata:
- mTEC component identified as comp with `argmax(H_norm[:, mTEC_WT_col])` 
- For K=8 with 12 samples (1 mTEC + 11 thymocyte), comp1 = mTEC (H loading = 1.000)

## Differential Contact Correlation

The key discriminating metric between RIPHAMMER and THUNDER.

**Definition**: For each cell type t and feature (bin pair) f:
```
diff_true(f, t) = log2((X_pure[f, t] + eps) / (mean(X_pure[f, others]) + eps))
diff_pred(f, t) = log2((W_reordered[f, t] + eps) / (mean(W_reordered[f, others]) + eps))
```

Pearson correlation between diff_true and diff_pred computed per sliding window
(50 bins = 5 Mb, step 25 bins).

**Why THUNDER scores ~0**: THUNDER operates on raw (non-band-normalized) contact
values. The W matrix is dominated by distance decay, so differential patterns
between cell types are buried. RIPHAMMER's per-band L1 normalization removes
distance decay before NMF, enabling cell-type-specific patterns to emerge.

## Key Numerical Results

| Metric | RIPHAMMER | THUNDER | Std NMF | PCA |
|---|---|---|---|---|
| Pure separation | 4/4 | 4/4 | 2/4 | 3/4 |
| MAE | 0.022 | 0.043 | 0.011 | 0.193 |
| Mean r | 0.999 | 0.995 | 1.000 | 0.621 |
| TAD Precision | 0.867 | 0.802 | 0.679 | — |
| TAD F1 | 0.773 | 0.796 | 0.664 | — |
| Insulation r (median) | 0.954 | 0.902 | 0.862 | — |
| Diff contact corr (all windows) | 0.636 | -0.009 | 0.366 | — |
| Diff contact corr (Aire locus) | 0.819 | 0.049 | -0.374 | — |

Note: Std NMF has lowest MAE (0.011) but only 2/4 separation — it assigns
similar proportions to all samples because distance decay dominates.

## Required Python Packages

- riphammer (pip, from /home/koki/dev/RIPHAMMER or PyPI)
- numpy, pandas, matplotlib, scipy, sklearn
- cooler (for real data plots reading .mcool)
- cooltools (for insulation score in aire_locus_10kb_deep, brain_dev_loci)

## Required R Packages (THUNDER only)

- NMF (CRAN)
- thundeR (GitHub: brycerowland/thundeR)
- dplyr, readr
