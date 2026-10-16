# Reproducibility Guide

This document describes how to reproduce all RIPHAMMER validation analyses and figures.

## Prerequisites

### Software

| Tool | Version | Purpose |
|---|---|---|
| Snakemake | ≥ 8.10.0 | Pipeline orchestration |
| Python | ≥ 3.11 | NMF, vectorization, plotting |
| R | ≥ 4.2 | THUNDER comparison |
| RIPHAMMER | latest (`feat/single-matrix-and-h-norm` branch) | Multi-band NMF |

### Conda environments

Two conda environments are defined under `workflow/envs/`:

- **riphammer.yaml** — Python stack: numpy, pandas, cooler, cooltools, matplotlib, RIPHAMMER
- **thunder.yaml** — R stack: NMF, dplyr, readr (for THUNDER baseline)

### Input data

All input Hi-C data (mcool files) are listed in `workflow/config.yaml`:

- **mTEC**: `data/mtec/processed/mtec_hic.mcool`
- **Thymocyte subtypes** (11 samples): `data/4dn_thymocyte/*.mcool`
- **Brain** (ES, NPC, CN): `data/brain/processed/brain_hic_*.mcool`
- **Simulation**: `output/simulation_100000/` (pre-generated ground truth)

## Running the full pipeline

```bash
# From repository root
bash workflow/run.sh 4          # 4 cores (adjust as needed)

# Or directly with snakemake
snakemake --snakefile workflow/Snakefile --cores 4 --printshellcmds --keep-going --rerun-incomplete
```

### Estimated runtime (4 cores, 64 GB RAM)

| Stage | Time |
|---|---|
| Vectorize (100kb + 25kb + 10kb) | ~10 min |
| Vectorize brain | ~2 min |
| NMF fits (mTEC: 100kb, 25kb, 10kb × K=8) | ~3–8 h (10kb is slowest) |
| NMF fits (brain: K=2, K=3) | ~5 min |
| NMF fits (simulation: K=4) | ~1 min |
| TAD detection (all resolutions) | ~5–30 min |
| THUNDER (full genome + W extraction) | ~30 min |
| Comparison data generation | ~10 min |
| All plots | ~30 min (10kb plots are slowest due to 121M-row W matrix) |

Total: approximately 4–10 hours depending on hardware.

## Pipeline DAG

The Snakemake DAG contains 30 rules organized in 6 phases:

### Phase 0: Vectorization
- `vectorize` — mcool → feature matrices (mTEC + thymocyte, per resolution)
- `vectorize_brain` — mcool → feature matrices (brain, 100kb only)

### Phase 1: NMF fitting
- `fit` — Multi-band NMF on mTEC+thymocyte ({100kb, 25kb, 10kb} × K=8)
- `fit_brain` — Multi-band NMF on brain (100kb, K=2 and K=3)
- `fit_brain_k2` — Brain K=2 in separate output directory
- `fit_simulation` — Multi-band NMF on simulation (100kb, K=4)

### Phase 2: TAD detection
- `tad` — Insulation-based TAD calling per NMF component (mTEC+thymocyte)
- `tad_brain` — TAD calling for brain K=3

### Phase 3: THUNDER comparison
- `thunder_fullgenome` — Full-genome THUNDER pipeline (R)
- `thunder_extract_W` — Extract W/H matrices from THUNDER (R)

### Phase 4: Intermediate data
- `comparison_3way` — RIPHAMMER vs Standard NMF vs PCA metrics
- `tad_comparison` — TAD boundary comparison across methods

### Phase 5: Figures (19 plots)

| Plot file | Description | Figure |
|---|---|---|
| `simulation_fraction_recovery.png` | H matrix fraction recovery | Fig 2 |
| `method_comparison_simulation.png` | 3-method comparison | Fig 2 |
| `method_comparison_4way.png` | 4-method with THUNDER | Fig 3 |
| `biological_discovery_comparison.png` | TAD + insulation + MAE | Fig 3 |
| `locus_discovery_comparison.png` | Differential contact analysis | Fig 3 |
| `robustness_analysis.png` | K-sweep + seed robustness | Fig S |
| `scalability_benchmark.png` | Timing benchmarks | Fig S |
| `noise_robustness.png` | SNR robustness | Fig S |
| `contact_maps_{res}_k8.png` | Contact map overview (×3 res) | Fig 4 |
| `aire_collaborator_summary.png` | Aire locus summary | Fig 4 |
| `aire_locus_10kb_components.png` | Per-component Aire contacts | Fig 5 |
| `aire_locus_10kb_deep.png` | Aire deep analysis + TAD | Fig 5 |
| `aire_locus_10kb_differential.png` | Differential contact at Aire | Fig 5 |
| `aire_locus_10kb_insulation.png` | Insulation score comparison | Fig 5 |
| `aire_locus_10kb_selected.png` | Selected component maps | Fig 5 |
| `brain_dev_k2_100kb.png` | Brain ES vs differentiation | Fig 6 |
| `brain_dev_loci_100kb.png` | Brain developmental loci | Fig 6 |

## Key numerical results

### Simulation benchmark (100kb, K=4)

| Metric | RIPHAMMER | THUNDER | Std NMF | PCA |
|---|---|---|---|---|
| Pure separation | 4/4 | 4/4 | 2/4 | 3/4 |
| MAE | 0.022 | 0.043 | 0.011 | 0.193 |
| Mean correlation | 0.999 | 0.995 | 1.000 | 0.621 |
| TAD Precision | 0.892 | 0.838 | 0.750 | — |
| TAD F1 | 0.793 | 0.832 | 0.734 | — |
| Insulation correlation (median) | 0.953 | 0.902 | 0.862 | — |
| Diff contact correlation (all) | 0.631 | −0.009 | 0.366 | — |
| Diff contact correlation (Aire) | 0.819 | 0.049 | −0.374 | — |

### RIPHAMMER settings

All NMF fits use RIPHAMMER defaults:
- `mode='single'` (single-matrix NMF per distance band)
- `normalize_h=True` (H columns sum to 1)
- `normalize_bands=True` (distance-band normalization)
- `seed=42` for reproducibility

## File inventory

```
workflow/
  Snakefile             # 30-rule DAG (601 lines)
  config.yaml           # Sample paths, resolutions, parameters
  run.sh                # Wrapper script
  envs/
    riphammer.yaml      # Python conda env
    thunder.yaml        # R conda env

src/
  vectorize.py          # mcool → mTEC+thymocyte feature matrices
  vectorize_brain.py    # mcool → brain feature matrices
  compute_comparison_3way.py   # RIPHAMMER vs StdNMF vs PCA
  compute_tad_comparison.py    # TAD boundary comparison
  run_thunder_fullgenome.R     # THUNDER full-genome pipeline
  run_thunder_extract_W.R      # THUNDER W/H extraction
  plot_*.py (16 scripts)       # One script per figure

output/
  features_{res}/       # Vectorized feature matrices
  nmf_{res}/            # NMF W/H matrices (mTEC+thymocyte)
  nmf_brain_100000/     # Brain NMF results
  nmf_simulation_100000/  # Simulation NMF results
  tad_{res}_k8/         # TAD boundaries
  simulation_100000/    # Simulation data + comparison results

plot/                   # All generated figures (PNG, 150 dpi)
```
