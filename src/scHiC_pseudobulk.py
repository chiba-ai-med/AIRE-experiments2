#!/usr/bin/env python3
"""Aggregate single-cell Hi-C contacts by cell type into pseudo-bulk mcool files.

Memory-efficient version: bins contacts on-the-fly into Counter dicts
instead of storing raw (chr, pos) tuples. At 100kb resolution this reduces
peak RAM from ~90GB to ~2-5GB for 500M pairs.
"""
import argparse
import collections
import gzip
import json
import os
import re

import cooler
import numpy as np
import pandas as pd


def load_annotations(path):
    df = pd.read_csv(path, sep="\t")
    return dict(zip(df["barcode"], df["celltype"]))


def load_pairs_spec(path):
    with open(path) as f:
        return json.load(f)


def extract_barcode(fields, spec):
    src = spec["barcode_source"]
    bc_col = spec.get("barcode_col")

    if src == "column" and bc_col is not None:
        return fields[bc_col] if len(fields) > bc_col else None
    elif src == "readID_prefix":
        m = re.match(r"^([ACGTN]+)", fields[0])
        return m.group(1) if m else None
    elif src == "readID_suffix":
        return fields[0].split(":")[-1] if ":" in fields[0] else None
    return None


def get_chromsizes(genome):
    try:
        import bioframe
        return bioframe.fetch_chromsizes(genome)
    except Exception:
        pass
    import urllib.request
    url = f"https://hgdownload.cse.ucsc.edu/goldenpath/{genome}/bigZips/{genome}.chrom.sizes"
    path = f"/tmp/{genome}.chrom.sizes"
    if not os.path.exists(path):
        urllib.request.urlretrieve(url, path)
    sizes = pd.read_csv(path, sep="\t", header=None, names=["name", "length"])
    sizes = sizes[~sizes["name"].str.contains("_")]
    return sizes.set_index("name")["length"]


def build_pseudobulk_mcools(pairs_files, spec, bc_to_ct, outdir,
                            resolution, chromsizes_name, min_cells,
                            min_contacts, threads):
    chromsizes = get_chromsizes(chromsizes_name)
    autosomes = [c for c in chromsizes.index if re.match(r"^chr\d+$", c)]
    chromsizes = chromsizes[autosomes]
    bins = cooler.binnify(chromsizes, resolution)

    chrom_offsets = {}
    offset = 0
    for chrom in chromsizes.index:
        chrom_offsets[chrom] = offset
        offset += int(np.ceil(chromsizes[chrom] / resolution))

    ct_counts = collections.Counter(bc_to_ct.values())
    valid_cts = {ct for ct, n in ct_counts.items() if n >= min_cells}
    print(f"Cell types with >= {min_cells} cells: {sorted(valid_cts)}", flush=True)

    # Accumulate directly into per-celltype pixel counters (memory-efficient)
    ct_pixels = {ct: collections.Counter() for ct in valid_cts}
    ct_n_contacts = collections.Counter()

    chr1_col = spec.get("chr1_col", 0)
    pos1_col = spec.get("pos1_col", 1)
    chr2_col = spec.get("chr2_col", 2)
    pos2_col = spec.get("pos2_col", 3)
    min_cols = max(chr1_col, pos1_col, chr2_col, pos2_col,
                   spec.get("barcode_col") or 0) + 1

    total_count = 0
    total_matched = 0

    for pairs_path in pairs_files:
        print(f"  Streaming {os.path.basename(pairs_path)}...", flush=True)
        with gzip.open(pairs_path, "rt") as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                fields = line.rstrip("\n").split("\t")
                if len(fields) < min_cols:
                    continue

                total_count += 1

                bc = extract_barcode(fields, spec)
                if bc is None:
                    continue
                ct = bc_to_ct.get(bc)
                if ct is None or ct not in valid_cts:
                    continue

                chr1 = fields[chr1_col].replace("mm10_", "")
                chr2 = fields[chr2_col].replace("mm10_", "")
                if chr1 not in chrom_offsets or chr2 not in chrom_offsets:
                    continue

                pos1 = int(fields[pos1_col])
                pos2 = int(fields[pos2_col])
                idx1 = chrom_offsets[chr1] + pos1 // resolution
                idx2 = chrom_offsets[chr2] + pos2 // resolution
                if idx1 > idx2:
                    idx1, idx2 = idx2, idx1

                ct_pixels[ct][(idx1, idx2)] += 1
                ct_n_contacts[ct] += 1
                total_matched += 1

                if total_count % 10_000_000 == 0:
                    print(f"    processed {total_count:,} pairs, "
                          f"matched {total_matched:,}", flush=True)

        print(f"  File done. Running total: {total_count:,} processed, "
              f"{total_matched:,} matched", flush=True)

    print(f"\nTotal: {total_count:,} pairs, {total_matched:,} matched",
          flush=True)

    summary_rows = []
    for ct in sorted(ct_pixels.keys()):
        pixels = ct_pixels[ct]
        n_contacts = ct_n_contacts[ct]
        n_cells = ct_counts[ct]

        if n_contacts < min_contacts:
            print(f"  Skipping {ct}: {n_contacts:,} contacts < {min_contacts:,}",
                  flush=True)
            continue

        print(f"  Writing cool for {ct}: {n_cells} cells, "
              f"{n_contacts:,} contacts, {len(pixels):,} unique pixels",
              flush=True)

        keys = sorted(pixels.keys())
        pixels_df = pd.DataFrame({
            "bin1_id": [k[0] for k in keys],
            "bin2_id": [k[1] for k in keys],
            "count": [pixels[k] for k in keys],
        })

        cool_path = os.path.join(outdir, f"{ct}.cool")
        mcool_path = os.path.join(outdir, f"{ct}.mcool")

        cooler.create_cooler(
            cool_path, bins=bins, pixels=pixels_df,
            ordered=True, columns=["count"],
        )
        print(f"  Building mcool for {ct}...", flush=True)
        os.system(
            f"cooler zoomify -r {resolution}N "
            f"-o {mcool_path} {cool_path} 2>/dev/null"
        )
        if os.path.exists(cool_path):
            os.remove(cool_path)

        summary_rows.append({
            "celltype": ct,
            "n_cells": n_cells,
            "n_contacts": n_contacts,
            "mcool": mcool_path,
        })

        # Free memory after writing
        del ct_pixels[ct]
        print(f"  Done: {ct}", flush=True)

    summary = pd.DataFrame(summary_rows)
    return summary


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--pairs-files", nargs="+", required=True)
    parser.add_argument("--annotations", required=True)
    parser.add_argument("--pairs-spec", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--resolution", type=int, required=True)
    parser.add_argument("--chromsizes", default="mm10")
    parser.add_argument("--min-cells", type=int, default=50)
    parser.add_argument("--min-contacts", type=int, default=100000)
    parser.add_argument("--threads", type=int, default=4)
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    bc_to_ct = load_annotations(args.annotations)
    spec = load_pairs_spec(args.pairs_spec)

    print(f"Loaded {len(bc_to_ct)} barcode->celltype mappings", flush=True)

    summary = build_pseudobulk_mcools(
        args.pairs_files, spec, bc_to_ct, args.outdir,
        args.resolution, args.chromsizes,
        args.min_cells, args.min_contacts, args.threads,
    )

    summary_path = os.path.join(args.outdir, "summary.tsv")
    summary.to_csv(summary_path, sep="\t", index=False)
    print(f"\nSummary saved to {summary_path}", flush=True)
    print(summary.to_string(index=False), flush=True)


if __name__ == "__main__":
    main()
