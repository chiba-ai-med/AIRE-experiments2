#!/usr/bin/env python3
"""Derive cell-type annotations from GAGE-seq RNA component.

GAGE-seq RNA files are per-read TSVs (no header):
  chr  barcode_pair  gene_ensembl_id

We build a barcode x gene count matrix, then cluster with PCA + KMeans
and annotate using known cortex marker genes.
"""
import argparse
import collections
import gzip
import os
import sys

import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix
from sklearn.decomposition import TruncatedSVD
from sklearn.cluster import KMeans
from sklearn.preprocessing import normalize


def load_rna_counts(rna_paths):
    """Build barcode x gene count matrix from per-read RNA files."""
    counts = collections.defaultdict(lambda: collections.defaultdict(int))
    total_reads = 0

    for path in rna_paths:
        print(f"Loading {os.path.basename(path)}...", flush=True)
        n = 0
        opener = gzip.open if path.endswith(".gz") else open
        with opener(path, "rt") as fh:
            for line in fh:
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 3:
                    continue
                barcode = fields[1]
                gene = fields[2]
                counts[barcode][gene] += 1
                n += 1
        print(f"  {n:,} reads", flush=True)
        total_reads += n

    barcodes = sorted(counts.keys())
    all_genes = set()
    for bc_counts in counts.values():
        all_genes.update(bc_counts.keys())
    genes = sorted(all_genes)

    print(f"Total: {total_reads:,} reads, {len(barcodes)} barcodes, "
          f"{len(genes)} genes", flush=True)

    rows, cols, vals = [], [], []
    gene_idx = {g: i for i, g in enumerate(genes)}
    for i, bc in enumerate(barcodes):
        for gene, count in counts[bc].items():
            rows.append(i)
            cols.append(gene_idx[gene])
            vals.append(count)

    mat = csr_matrix((vals, (rows, cols)),
                     shape=(len(barcodes), len(genes)), dtype=np.float32)
    return mat, barcodes, genes


def cluster_and_annotate(mat, barcodes, genes, n_clusters=8, n_pcs=30):
    """PCA + KMeans clustering, annotate via marker genes."""
    gene_counts = np.diff(mat.indptr)
    cell_mask = gene_counts >= 100
    mat_filt = mat[cell_mask]
    barcodes_filt = [b for b, m in zip(barcodes, cell_mask) if m]

    gene_cell_counts = (mat_filt > 0).sum(axis=0).A1
    gene_mask = gene_cell_counts >= 5
    mat_filt = mat_filt[:, gene_mask]
    genes_filt = [g for g, m in zip(genes, gene_mask) if m]

    print(f"After filtering: {mat_filt.shape[0]} cells, {mat_filt.shape[1]} genes",
          flush=True)

    row_sums = mat_filt.sum(axis=1).A1
    row_sums[row_sums == 0] = 1
    mat_norm = mat_filt.multiply(1e4 / row_sums[:, None]).tocsr()
    mat_norm.data = np.log1p(mat_norm.data)

    n_comps = min(n_pcs, mat_norm.shape[0] - 1, mat_norm.shape[1] - 1)
    svd = TruncatedSVD(n_components=n_comps, random_state=42)
    X_pca = svd.fit_transform(mat_norm)
    X_pca = normalize(X_pca)

    print(f"PCA: {n_comps} components, explained var={svd.explained_variance_ratio_.sum():.3f}",
          flush=True)

    km = KMeans(n_clusters=n_clusters, random_state=42, n_init=10)
    labels = km.fit_predict(X_pca)

    markers_ensembl = {
        "Excitatory": ["ENSMUSG00000070570", "ENSMUSG00000030067"],
        "Inhibitory": ["ENSMUSG00000070880", "ENSMUSG00000020263"],
        "Astrocyte": ["ENSMUSG00000024411", "ENSMUSG00000020932"],
        "Oligodendrocyte": ["ENSMUSG00000041607", "ENSMUSG00000032517"],
        "OPC": ["ENSMUSG00000029231"],
        "Microglia": ["ENSMUSG00000052336", "ENSMUSG00000024397"],
    }

    gene_to_idx = {g: i for i, g in enumerate(genes_filt)}
    mat_dense_markers = {}
    for ct, ensembl_ids in markers_ensembl.items():
        cols = [gene_to_idx[g] for g in ensembl_ids if g in gene_to_idx]
        if cols:
            mat_dense_markers[ct] = mat_norm[:, cols].toarray()

    cluster_labels = {}
    for cl in range(n_clusters):
        mask = labels == cl
        n_cells = mask.sum()
        best_type = f"Cluster_{cl}"
        best_score = -1

        for ct, expr in mat_dense_markers.items():
            mean_expr = np.mean(expr[mask])
            if mean_expr > best_score:
                best_score = mean_expr
                best_type = ct

        cluster_labels[cl] = best_type
        print(f"  Cluster {cl} ({n_cells} cells) -> {best_type} "
              f"(score={best_score:.3f})", flush=True)

    celltypes = [cluster_labels[l] for l in labels]
    result = pd.DataFrame({
        "barcode": barcodes_filt,
        "celltype": celltypes,
    })
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--rna-files", nargs="+", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    mat, barcodes, genes = load_rna_counts(args.rna_files)
    annotations = cluster_and_annotate(mat, barcodes, genes)

    print(f"\nCell type counts:", flush=True)
    for ct, cnt in annotations["celltype"].value_counts().items():
        print(f"  {ct}: {cnt}", flush=True)

    annotations.to_csv(args.output, sep="\t", index=False)
    print(f"\nSaved {len(annotations)} annotations to {args.output}", flush=True)


if __name__ == "__main__":
    main()
