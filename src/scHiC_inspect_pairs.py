#!/usr/bin/env python3
"""Inspect a .pairs.gz file header to determine column layout and barcode encoding.

Supports two formats:
- 4DN standard: #columns header, readID chr1 pos1 chr2 pos2 strand1 strand2 [barcode]
- GAGE-seq: no header, chr1 pos1 chr2 pos2 barcode (5 columns)
- Droplet Hi-C: 4DN standard with barcode in readID or extra column
"""
import argparse
import gzip
import json
import re
import sys


def inspect(pairs_path, n_data_lines=100):
    columns = None
    header_lines = []
    data_lines = []

    with gzip.open(pairs_path, "rt") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith("#"):
                header_lines.append(line)
                if line.startswith("#columns:"):
                    columns = line.split(":")[1].strip().split()
            else:
                data_lines.append(line)
                if len(data_lines) >= n_data_lines:
                    break

    if not data_lines:
        return {"error": "No data lines found"}

    first = data_lines[0]
    n_cols = len(first.split("\t"))

    barcode_source = "unknown"
    barcode_col = None

    if columns is not None:
        # 4DN standard pairs format with header
        pass
    elif n_cols == 5:
        # GAGE-seq format: chr1 pos1 chr2 pos2 barcode
        columns = ["chr1", "pos1", "chr2", "pos2", "barcode"]
        barcode_source = "column"
        barcode_col = 4
    elif n_cols >= 7:
        columns = (["readID", "chr1", "pos1", "chr2", "pos2",
                     "strand1", "strand2"] +
                   [f"col{i}" for i in range(7, n_cols)])
    else:
        columns = [f"col{i}" for i in range(n_cols)]

    # Auto-detect barcode in various positions
    if barcode_source == "unknown" and data_lines:
        fields = data_lines[0].split("\t")

        # Check last column for barcode-like pattern
        last = fields[-1] if fields else ""
        if re.match(r"^[A-Z]\d+,[A-Z]\d+$", last):
            barcode_source = "column"
            barcode_col = len(fields) - 1
        elif re.match(r"^[ACGTN]{8,}$", last):
            barcode_source = "column"
            barcode_col = len(fields) - 1

        # Check readID (first column) for barcode prefix
        if barcode_source == "unknown":
            read_id = fields[0]
            if re.match(r"^[ACGTN]{8,}", read_id):
                barcode_source = "readID_prefix"
            elif ":" in read_id:
                parts = read_id.split(":")
                if re.match(r"^[ACGTN]{8,}$", parts[-1]):
                    barcode_source = "readID_suffix"

    def find_col(names, default):
        for n in names:
            if n in columns:
                return columns.index(n)
        return default

    spec = {
        "columns": columns,
        "n_columns": n_cols,
        "has_header": len(header_lines) > 0,
        "barcode_source": barcode_source,
        "barcode_col": barcode_col,
        "chr1_col": find_col(["chr1", "chrom1"], 0),
        "pos1_col": find_col(["pos1"], 1),
        "chr2_col": find_col(["chr2", "chrom2"], 2),
        "pos2_col": find_col(["pos2"], 3),
        "header_lines": header_lines[:10],
        "sample_data": data_lines[:5],
    }

    return spec


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--pairs", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    spec = inspect(args.pairs)
    with open(args.output, "w") as f:
        json.dump(spec, f, indent=2)

    print(f"Columns ({spec['n_columns']}): {spec['columns']}", flush=True)
    print(f"Has header: {spec['has_header']}", flush=True)
    print(f"Barcode source: {spec['barcode_source']}", flush=True)
    if spec["barcode_col"] is not None:
        print(f"Barcode column index: {spec['barcode_col']}", flush=True)
    print(f"Coordinate columns: chr1={spec['chr1_col']}, pos1={spec['pos1_col']}, "
          f"chr2={spec['chr2_col']}, pos2={spec['pos2_col']}", flush=True)
    print(f"Sample data:", flush=True)
    for line in spec["sample_data"][:3]:
        print(f"  {line}", flush=True)


if __name__ == "__main__":
    main()
