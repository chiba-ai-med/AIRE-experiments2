#!/usr/bin/env python3
"""Merge pseudo-bulk mcool files by cell-type grouping.

Reads existing per-celltype pseudo-bulk mcool files and merges (sums)
their contact matrices according to a provided mapping, producing
new merged mcool files and a new summary.tsv.
"""
import argparse
import json
import os
import tempfile

import cooler
import pandas as pd


def merge_coolers(mcool_paths, resolution, output_mcool):
    """Sum contact matrices from multiple mcool files into one."""
    all_pixels = []
    bins = None

    for path in mcool_paths:
        uri = f"{path}::resolutions/{resolution}"
        c = cooler.Cooler(uri)
        px = c.pixels()[:]
        if bins is None:
            bins = c.bins()[:]
        all_pixels.append(px)

    combined = pd.concat(all_pixels, ignore_index=True)
    merged = (combined.groupby(["bin1_id", "bin2_id"], sort=True)
              ["count"].sum().reset_index())

    with tempfile.NamedTemporaryFile(suffix=".cool", delete=False) as tmp:
        tmp_cool = tmp.name

    cooler.create_cooler(
        tmp_cool,
        bins=bins,
        pixels=merged,
        ordered=True,
        columns=["count"],
        dtypes={"count": float},
    )
    os.system(
        f"cooler zoomify -r {resolution}N "
        f"-o {output_mcool} {tmp_cool} 2>/dev/null"
    )
    os.remove(tmp_cool)
    return len(merged)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--summary", required=True,
                        help="Original pseudobulk summary.tsv")
    parser.add_argument("--merge-map", required=True,
                        help="JSON file: {merged_name: [orig_ct, ...]}")
    parser.add_argument("--resolution", type=int, default=100000)
    parser.add_argument("--output-dir", required=True)
    args = parser.parse_args()

    with open(args.merge_map) as f:
        merge_map = json.load(f)

    summary = pd.read_csv(args.summary, sep="\t")
    ct_to_mcool = dict(zip(summary["celltype"], summary["mcool"]))
    ct_to_ncells = dict(zip(summary["celltype"], summary["n_cells"]))
    ct_to_ncontacts = dict(zip(summary["celltype"], summary["n_contacts"]))

    os.makedirs(args.output_dir, exist_ok=True)
    rows = []

    for merged_name, orig_cts in sorted(merge_map.items()):
        print(f"\nMerging {merged_name}: {orig_cts}", flush=True)
        mcool_paths = [ct_to_mcool[ct] for ct in orig_cts]
        n_cells = sum(ct_to_ncells[ct] for ct in orig_cts)
        n_contacts = sum(ct_to_ncontacts[ct] for ct in orig_cts)

        output_mcool = os.path.join(args.output_dir, f"{merged_name}.mcool")
        n_pixels = merge_coolers(mcool_paths, args.resolution, output_mcool)

        print(f"  {len(orig_cts)} types, {n_cells} cells, "
              f"{n_contacts:,} contacts, {n_pixels:,} pixels", flush=True)
        rows.append({
            "celltype": merged_name,
            "n_cells": n_cells,
            "n_contacts": n_contacts,
            "mcool": output_mcool,
        })

    out_summary = os.path.join(args.output_dir, "summary.tsv")
    pd.DataFrame(rows).to_csv(out_summary, sep="\t", index=False)
    print(f"\nSaved {out_summary}", flush=True)


if __name__ == "__main__":
    main()
