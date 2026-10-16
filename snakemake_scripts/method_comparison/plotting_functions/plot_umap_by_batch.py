#!/usr/bin/env python
"""
Plot UMAP colored by Leiden cluster, one panel per batch, where only
cells belonging to that batch are colored (all other cells shown in gray).

Uses scanpy's built-in `mask_obs` parameter (scanpy >= 1.10) if available;
otherwise falls back to a manual matplotlib implementation that does the
same thing.

Usage:
    python plot_umap_by_batch.py <path_to_h5ad> [--outdir OUTDIR]
"""

import argparse
import os
from packaging.version import Version

import numpy as np
import scanpy as sc
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def plot_with_mask_obs(adata, outdir):
    """Use scanpy's native mask_obs parameter (scanpy >= 1.10)."""
    sc.settings.figdir = outdir
    for batch in adata.obs["batch"].cat.categories:
        safe_name = str(batch).replace("/", "_")
        mask = (adata.obs["batch"] == batch).values
        sc.pl.umap(
            adata,
            color="leiden",
            mask_obs=mask,
            save=f"_{safe_name}.png",
            title=f"{batch} (n={mask.sum()})",
            frameon=False,
            show=False,
        )
        # scanpy prefixes with "umap_batch_leiden", so also emit a matching pdf
        sc.pl.umap(
            adata,
            color="leiden",
            mask_obs=mask,
            save=f"_{safe_name}.pdf",
            title=f"{batch} (n={mask.sum()})",
            frameon=False,
            show=False,
        )
        print(f"{batch}: saved via mask_obs")


def plot_manual(adata, outdir):
    """Fallback for scanpy < 1.10: manual matplotlib implementation."""
    os.makedirs(outdir, exist_ok=True)

    umap_coords = adata.obsm["X_umap"]
    leiden_cats = adata.obs["leiden"].cat.categories
    leiden_colors = adata.uns.get("leiden_colors", None)
    if leiden_colors is None:
        # fall back to default matplotlib tab20 cycling if not present
        cmap = plt.get_cmap("tab20")
        leiden_colors = [cmap(i % 20) for i in range(len(leiden_cats))]

    for batch in adata.obs["batch"].cat.categories:
        safe_name = str(batch).replace("/", "_")
        mask = (adata.obs["batch"] == batch).values

        fig, ax = plt.subplots(figsize=(5, 5))

        # background: all cells, light gray
        ax.scatter(
            umap_coords[~mask, 0], umap_coords[~mask, 1],
            c="lightgray", s=2, linewidths=0, alpha=0.3,
        )

        # foreground: this batch's cells, colored by leiden
        for j, cat in enumerate(leiden_cats):
            cat_mask = mask & (adata.obs["leiden"] == cat).values
            if cat_mask.sum() == 0:
                continue
            ax.scatter(
                umap_coords[cat_mask, 0], umap_coords[cat_mask, 1],
                c=[leiden_colors[j]], s=3, linewidths=0, label=cat,
            )

        ax.set_title(f"{batch} (n={mask.sum()})", fontsize=10)
        ax.set_xticks([])
        ax.set_yticks([])
        fig.tight_layout()

        fig.savefig(os.path.join(outdir, f"umap_{safe_name}.png"), dpi=200)
        fig.savefig(os.path.join(outdir, f"umap_{safe_name}.pdf"))
        plt.close(fig)

        size = os.path.getsize(os.path.join(outdir, f"umap_{safe_name}.png"))
        print(f"{batch}: saved manually, size={size} bytes")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("h5ad_path", help="Path to input .h5ad file")
    parser.add_argument(
        "--outdir", default="umap_by_batch",
        help="Output directory for plots (default: umap_by_batch)",
    )
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    print(f"scanpy version: {sc.__version__}")
    adata = sc.read_h5ad(args.h5ad_path)

    if "batch" not in adata.obs.columns:
        raise ValueError("adata.obs must contain a 'batch' column")
    if "leiden" not in adata.obs.columns:
        raise ValueError("adata.obs must contain a 'leiden' column")
    if not str(adata.obs["batch"].dtype) == "category":
        adata.obs["batch"] = adata.obs["batch"].astype("category")

    if Version(sc.__version__) >= Version("1.10"):
        print("Using native mask_obs support")
        plot_with_mask_obs(adata, args.outdir)
    else:
        print("scanpy < 1.10 detected, using manual fallback")
        plot_manual(adata, args.outdir)

    print(f"Done. Plots saved to: {args.outdir}")


if __name__ == "__main__":
    main()
