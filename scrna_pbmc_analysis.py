#!/usr/bin/env python3
"""
Single-cell RNA-seq analysis of 3k PBMCs (10x Genomics).

Downloads the real, public 10x Genomics PBMC3k dataset and runs a standard
Scanpy workflow end to end:

    QC  ->  normalisation  ->  HVG selection  ->  PCA  ->  neighbours
        ->  UMAP  ->  Leiden clustering  ->  marker genes  ->  cell-type labels

All figures and result tables are written to figures/ and results/ so the
analysis is fully reproducible: delete those folders, re-run, get them back.

Run locally:
    pip install -r requirements.txt
    python scrna_pbmc_analysis.py
"""
import os
import warnings
from pathlib import Path

warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

FIG_DIR = Path("figures")
RES_DIR = Path("results")
FIG_DIR.mkdir(exist_ok=True)
RES_DIR.mkdir(exist_ok=True)

sc.settings.verbosity = 1
sc.settings.figdir = FIG_DIR
sc.settings.set_figure_params(dpi=150, facecolor="white")

# Canonical marker genes for the major PBMC populations. Used to give each
# Leiden cluster a human-readable label instead of just a number.
PBMC_MARKERS = {
    "CD4 T": ["IL7R", "CD3D", "CD3E"],
    "CD8 T": ["CD8A", "CD8B", "GZMK"],
    "B": ["MS4A1", "CD79A", "CD79B"],
    "NK": ["GNLY", "NKG7", "KLRD1"],
    "CD14 Monocyte": ["CD14", "LYZ", "S100A8"],
    "FCGR3A Monocyte": ["FCGR3A", "MS4A7"],
    "Dendritic": ["FCER1A", "CST3"],
    "Platelet": ["PPBP"],
}


def label_clusters(adata: sc.AnnData) -> dict:
    """Assign a cell-type label to each Leiden cluster by mean marker expression."""
    labels = {}
    for cluster in adata.obs["leiden"].cat.categories:
        mask = adata.obs["leiden"] == cluster
        best_type, best_score = "Unknown", -np.inf
        for cell_type, genes in PBMC_MARKERS.items():
            present = [g for g in genes if g in adata.raw.var_names]
            if not present:
                continue
            expr = adata.raw[mask, present].X
            score = float(np.asarray(expr.mean()))
            if score > best_score:
                best_type, best_score = cell_type, score
        labels[cluster] = best_type
    return labels


def main() -> None:
    print("Downloading real 10x PBMC3k dataset...")
    adata = sc.datasets.pbmc3k()
    print(f"Raw data: {adata.n_obs} cells x {adata.n_vars} genes")

    # ---- Quality control -------------------------------------------------
    adata.var_names_make_unique()
    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=3)
    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], percent_top=None, inplace=True)

    sc.pl.violin(
        adata, ["n_genes_by_counts", "total_counts", "pct_counts_mt"],
        jitter=0.4, multi_panel=True, show=False, save="_qc.png",
    )

    n_before = adata.n_obs
    adata = adata[adata.obs.n_genes_by_counts < 2500].copy()
    adata = adata[adata.obs.pct_counts_mt < 5].copy()
    print(f"QC: kept {adata.n_obs}/{n_before} cells after filtering")

    # ---- Normalisation + feature selection -------------------------------
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
    adata.raw = adata
    adata = adata[:, adata.var.highly_variable].copy()
    sc.pp.regress_out(adata, ["total_counts", "pct_counts_mt"])
    sc.pp.scale(adata, max_value=10)

    # ---- Dimensionality reduction + clustering ---------------------------
    sc.tl.pca(adata, svd_solver="arpack", n_comps=40)
    sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
    sc.tl.umap(adata)
    sc.tl.leiden(adata, resolution=0.5)
    print(f"Found {adata.obs.leiden.nunique()} clusters")

    # ---- Marker genes + cell-type labels ---------------------------------
    sc.tl.rank_genes_groups(adata, "leiden", method="wilcoxon")
    labels = label_clusters(adata)
    adata.obs["cell_type"] = adata.obs["leiden"].map(labels).astype("category")

    # ---- Figures ---------------------------------------------------------
    sc.pl.umap(adata, color="leiden", legend_loc="on data",
               title="Leiden clusters", show=False, save="_clusters.png")
    sc.pl.umap(adata, color="cell_type",
               title="Annotated cell types (PBMC)", show=False, save="_celltypes.png")

    marker_flat = [g for genes in PBMC_MARKERS.values() for g in genes
                   if g in adata.raw.var_names]
    sc.pl.dotplot(adata, marker_flat, groupby="cell_type", show=False,
                  save="_markers.png")

    # ---- Result tables ---------------------------------------------------
    top = sc.get.rank_genes_groups_df(adata, group=None)
    top.to_csv(RES_DIR / "cluster_markers.csv", index=False)

    counts = adata.obs["cell_type"].value_counts()
    counts.to_csv(RES_DIR / "cell_type_counts.csv", header=["n_cells"])

    with open(RES_DIR / "summary.txt", "w") as fh:
        fh.write("Single-cell RNA-seq analysis of 3k PBMCs (10x Genomics)\n")
        fh.write("=" * 55 + "\n\n")
        fh.write(f"Cells after QC : {adata.n_obs}\n")
        fh.write(f"Genes (HVG)    : {adata.n_vars}\n")
        fh.write(f"Clusters       : {adata.obs.leiden.nunique()}\n\n")
        fh.write("Cell types identified:\n")
        for ct, n in counts.items():
            fh.write(f"  {ct:<18} {n} cells\n")

    print("Done. Figures in figures/, tables in results/.")


if __name__ == "__main__":
    main()
