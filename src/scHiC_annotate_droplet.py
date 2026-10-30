#!/usr/bin/env python3
"""Parse Droplet Hi-C metadata to extract cell-type annotations.

The metadata file is a TSV (despite .xls extension) with columns:
  Unnamed:0 (lib:barcode), rnatype (celltype), score, bc (barcode),
  total, mapped, cis, cis_1kb-, cis_1kb+, cis_10kb+, trans, sample (library)
"""
import argparse
import pandas as pd


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--metadata", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    df = pd.read_csv(args.metadata, sep="\t")

    print(f"Metadata shape: {df.shape}", flush=True)
    print(f"Columns: {list(df.columns)}", flush=True)

    out = pd.DataFrame({
        "barcode": df["bc"],
        "celltype": df["rnatype"],
        "library": df["sample"],
    })
    out = out.dropna(subset=["barcode", "celltype"])

    dup_bc = out["barcode"][out["barcode"].duplicated(keep=False)].unique()
    if len(dup_bc) > 0:
        print(f"Dropping {len(dup_bc)} ambiguous barcodes "
              f"({out['barcode'].isin(dup_bc).sum()} cells) shared across libraries",
              flush=True)
        out = out[~out["barcode"].isin(dup_bc)]

    print(f"\nCell type counts:", flush=True)
    for ct, cnt in out["celltype"].value_counts().items():
        print(f"  {ct}: {cnt}", flush=True)

    print(f"\nLibrary counts:", flush=True)
    for lib, cnt in out["library"].value_counts().items():
        print(f"  {lib}: {cnt}", flush=True)

    out[["barcode", "celltype"]].to_csv(args.output, sep="\t", index=False)
    print(f"\nSaved {len(out)} annotations to {args.output}", flush=True)


if __name__ == "__main__":
    main()
