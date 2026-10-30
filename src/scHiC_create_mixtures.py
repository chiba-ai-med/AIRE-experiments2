#!/usr/bin/env python3
"""Create synthetic bulk Hi-C mixtures from pseudo-bulk cell-type mcool files.

For each mixture, sample fractions from Dirichlet distribution and create
a weighted sum of pseudo-bulk contact matrices. Output is a multi-sample
mcool (one resolution) + ground truth fractions TSV.
"""
import argparse
import os

import cooler
import numpy as np
import pandas as pd


def load_pseudobulk_matrices(summary_path, pseudobulk_dir, resolution):
    """Load pseudo-bulk contact matrices as sparse pixel DataFrames."""
    summary = pd.read_csv(summary_path, sep="\t")
    celltypes = sorted(summary["celltype"].tolist())
    matrices = {}
    bins = None

    for _, row in summary.iterrows():
        ct = row["celltype"]
        mcool_path = row["mcool"]
        uri = f"{mcool_path}::resolutions/{resolution}"
        try:
            c = cooler.Cooler(uri)
            mat = c.pixels()[:]
            if bins is None:
                bins = c.bins()[:]
            matrices[ct] = mat
            print(f"  Loaded {ct}: {len(mat):,} pixels", flush=True)
        except Exception as e:
            print(f"  Warning: could not load {ct}: {e}", flush=True)

    return matrices, bins, sorted(matrices.keys())


def create_mixtures(matrices, bins, celltypes, n_mixtures, seed):
    """Create synthetic bulk mixtures with Dirichlet-sampled fractions.

    Each cell type's pseudo-bulk is normalized to equal total contacts
    before mixing, so the fraction directly controls the composition.
    """
    rng = np.random.RandomState(seed)
    K = len(celltypes)

    # Normalize each cell type to equal total contacts
    ct_totals = {}
    for ct in celltypes:
        ct_totals[ct] = matrices[ct]["count"].sum()
        print(f"  {ct} total contacts: {ct_totals[ct]:,.0f}", flush=True)

    target_total = min(ct_totals.values())
    norm_factors = {ct: target_total / ct_totals[ct] for ct in celltypes}
    print(f"  Normalizing to {target_total:,.0f} contacts per cell type", flush=True)

    fractions = rng.dirichlet(np.ones(K), size=n_mixtures)

    all_pixels = []
    sample_names = []

    for m in range(n_mixtures):
        frac = fractions[m]
        print(f"  Mixture {m}: " +
              ", ".join(f"{ct}={f:.3f}" for ct, f in zip(celltypes, frac)),
              flush=True)

        parts = []
        for i, ct in enumerate(celltypes):
            if frac[i] < 1e-6:
                continue
            mat = matrices[ct].copy()
            mat["count"] = mat["count"] * norm_factors[ct] * frac[i]
            parts.append(mat[["bin1_id", "bin2_id", "count"]])

        if parts:
            combined = pd.concat(parts, ignore_index=True)
            pixels_df = (combined.groupby(["bin1_id", "bin2_id"], sort=True)
                         ["count"].sum().reset_index())
        else:
            pixels_df = pd.DataFrame(columns=["bin1_id", "bin2_id", "count"])
        all_pixels.append(pixels_df)
        sample_names.append(f"mix_{m:03d}")

    frac_df = pd.DataFrame(fractions, columns=celltypes,
                           index=sample_names)
    return all_pixels, sample_names, frac_df, bins


def save_as_mcool(all_pixels, bins, sample_names, resolution, output_mcool):
    """Save mixtures as individual mcool files in a directory."""
    outdir = os.path.dirname(output_mcool)
    os.makedirs(outdir, exist_ok=True)

    mcool_paths = []
    for pixels, name in zip(all_pixels, sample_names):
        cool_path = os.path.join(outdir, f"{name}.cool")
        mcool_path = os.path.join(outdir, f"{name}.mcool")

        cooler.create_cooler(
            cool_path,
            bins=bins,
            pixels=pixels,
            ordered=True,
            columns=["count"],
            dtypes={"count": float},
        )
        os.system(
            f"cooler zoomify -r {resolution}N "
            f"-o {mcool_path} {cool_path} 2>/dev/null"
        )
        if os.path.exists(cool_path):
            os.remove(cool_path)
        mcool_paths.append(mcool_path)

    with open(output_mcool, "w") as f:
        for p in mcool_paths:
            f.write(p + "\n")

    return mcool_paths


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--pseudobulk-dir", required=True)
    parser.add_argument("--summary", required=True)
    parser.add_argument("--n-mixtures", type=int, default=20)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--resolution", type=int, default=100000)
    parser.add_argument("--output-mcool", required=True,
                        help="Output file listing mcool paths")
    parser.add_argument("--output-fractions", required=True)
    parser.add_argument("--output-samples", required=True)
    args = parser.parse_args()

    print("Loading pseudo-bulk matrices...", flush=True)
    matrices, bins, celltypes = load_pseudobulk_matrices(
        args.summary, args.pseudobulk_dir, args.resolution)
    print(f"Loaded {len(celltypes)} cell types: {celltypes}", flush=True)

    print(f"\nCreating {args.n_mixtures} synthetic mixtures...", flush=True)
    all_pixels, sample_names, frac_df, bins = create_mixtures(
        matrices, bins, celltypes, args.n_mixtures, args.seed)

    print("\nSaving mixtures as mcool...", flush=True)
    mcool_paths = save_as_mcool(
        all_pixels, bins, sample_names, args.resolution, args.output_mcool)

    frac_df.to_csv(args.output_fractions, sep="\t")
    print(f"Fractions saved to {args.output_fractions}", flush=True)

    with open(args.output_samples, "w") as f:
        for name, path in zip(sample_names, mcool_paths):
            f.write(f"{name}\tscHiC_mixture\n")
    print(f"Sample names saved to {args.output_samples}", flush=True)

    print(f"\nTrue fractions:\n{frac_df.to_string()}", flush=True)


if __name__ == "__main__":
    main()
