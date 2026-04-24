#!/usr/bin/env python3
"""
plot_umap_genes.py

Plot UMAP colored by gene expression for a list of genes from a config file.
Each gene gets its own 2x2 inch PDF with gene name and FlyBase ID in the title.

Usage:
    python plot_umap_genes.py --h5ad <path.h5ad> --config <genes.tsv> --outdir <output_dir>
    
    python scripts/method_comparison/plot_custom_umap.py \
        --h5ad results/integrated/integrated.h5ad \
        --config results/pseudotime_genes/wolbachia_infection/tradeseq_inputs/custom_genes.csv \
        --outdir results/custom_umaps/ 

Config file format (TSV or CSV, with header):
    gene_name   flybase_id
    sxl         FBgn0264270
    msl-2       FBgn0005640
    ...
"""

import argparse
import os
import sys
import warnings

import pandas as pd
import scanpy as sc
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

warnings.filterwarnings("ignore")


def parse_args():
    parser = argparse.ArgumentParser(description="Plot UMAPs colored by gene expression.")
    parser.add_argument("--h5ad", required=True, help="Path to input .h5ad file")
    parser.add_argument("--config", required=True,
                        help="TSV/CSV config with columns: Gene, FlyBaseId")
    parser.add_argument("--outdir", default="umap_gene_plots",
                        help="Output directory for PDFs (default: umap_gene_plots)")
    parser.add_argument("--umap-key", default="X_umap",
                        help="obsm key for UMAP coordinates (default: X_umap)")
    parser.add_argument("--layer", default=None,
                        help="Layer to use for expression (default: None = adata.X)")
    parser.add_argument("--cmap", default="magma",
                        help="Colormap for expression (default: magma)")
    parser.add_argument("--vmax", default="p99",
                        help="vmax for color scale: 'p99', 'p95', or a float (default: p99)")
    return parser.parse_args()


def load_config(config_path):
    """Load gene config file. Accepts TSV or CSV with gene and flybase_id columns.

    Handles:
    - Space-padded CSV headers/values (e.g. 'Gene, FlyBaseId')
    - Gene names containing spaces, dashes, or special characters
    - Flexible column name aliases
    """
    sep = "\t" if config_path.endswith(".tsv") else ","
    df = pd.read_csv(config_path, sep=sep, skipinitialspace=True)

    # Strip whitespace from all column names, normalize to lowercase_underscore
    df.columns = df.columns.str.strip().str.lower().str.replace(r"\s+", "_", regex=True)

    # Strip whitespace from all string values
    for col in df.select_dtypes(include="object").columns:
        df[col] = df[col].str.strip()

    # Accept flexible column names
    name_aliases = ["gene_name", "gene", "name", "symbol"]
    fb_aliases = ["flybase_id", "flybaseid", "fbgn", "flybase"]

    name_col = next((c for c in name_aliases if c in df.columns), None)
    fb_col = next((c for c in fb_aliases if c in df.columns), None)

    if name_col is None or fb_col is None:
        sys.exit(
            f"Config must have gene name and FlyBase ID columns.\n"
            f"  Recognized gene name columns: {name_aliases}\n"
            f"  Recognized FlyBase ID columns: {fb_aliases}\n"
            f"  Found columns: {list(df.columns)}"
        )

    return df[[name_col, fb_col]].rename(
        columns={name_col: "gene_name", fb_col: "flybase_id"}
    ).dropna()


def resolve_vmax(expr_values, vmax_arg):
    """Return a numeric vmax from string or float argument."""
    if isinstance(vmax_arg, str):
        if vmax_arg == "p99":
            return float(expr_values.quantile(0.99))
        elif vmax_arg == "p95":
            return float(expr_values.quantile(0.95))
        else:
            return float(vmax_arg)
    return float(vmax_arg)


def safe_filename(name):
    """Replace characters unsafe for filenames with underscores."""
    import re
    return re.sub(r"[^\w\-.]", "_", name)


def plot_gene(adata, gene_name, flybase_id, outdir, umap_key, layer, cmap, vmax_arg):
    """Generate and save a 2x2 inch PDF UMAP plot for one gene."""
    import scipy.sparse as sp

    # Check gene exists
    if gene_name not in adata.var_names:
        print(f"  [WARN] '{gene_name}' not found in adata.var_names — skipping.")
        return False

    # Extract expression
    if layer is not None:
        if layer not in adata.layers:
            print(f"  [WARN] Layer '{layer}' not found — using adata.X.")
            expr = adata[:, gene_name].X
        else:
            expr = adata[:, gene_name].layers[layer]
    else:
        expr = adata[:, gene_name].X

    # Flatten to 1D
    if sp.issparse(expr):
        expr = expr.toarray().flatten()
    else:
        expr = expr.flatten()

    expr_series = pd.Series(expr)
    vmax = resolve_vmax(expr_series, vmax_arg)

    # UMAP coordinates
    if umap_key not in adata.obsm:
        sys.exit(
            f"UMAP key '{umap_key}' not found in adata.obsm. "
            f"Available: {list(adata.obsm.keys())}"
        )
    umap = adata.obsm[umap_key]
    x, y = umap[:, 0], umap[:, 1]

    # Plot
    fig, ax = plt.subplots(figsize=(2, 2))

    sc_plot = ax.scatter(
        x, y,
        c=expr_series,
        cmap=cmap,
        vmin=0,
        vmax=vmax,
        s=0.5,
        linewidths=0,
        rasterized=True
    )

    title = f"{gene_name}\n{flybase_id}"
    ax.set_title(title, fontsize=5, pad=2)
    ax.set_xlabel("UMAP 1", fontsize=4)
    ax.set_ylabel("UMAP 2", fontsize=4)
    ax.tick_params(labelsize=3, length=2, pad=1)
    ax.set_aspect("equal", adjustable="datalim")

    cbar = fig.colorbar(sc_plot, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label("Expression", fontsize=4)
    cbar.ax.tick_params(labelsize=3)

    plt.tight_layout(pad=0.3)

    # Save
    safe_name = safe_filename(gene_name)
    out_path = os.path.join(outdir, f"{safe_name}_{flybase_id}.pdf")
    fig.savefig(out_path, format="pdf", bbox_inches="tight", dpi=300)
    plt.close(fig)

    print(f"  Saved: {out_path}")
    return True


def main():
    args = parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    print(f"Loading h5ad: {args.h5ad}")
    adata = sc.read_h5ad(args.h5ad)
    print(f"  {adata.n_obs} cells x {adata.n_vars} genes")

    print(f"Loading config: {args.config}")
    genes_df = load_config(args.config)
    print(f"  {len(genes_df)} genes to plot")

    n_plotted = 0
    n_skipped = 0
    for _, row in genes_df.iterrows():
        gene = str(row["gene_name"]).strip()
        fbid = str(row["flybase_id"]).strip()
        print(f"  Plotting {gene} ({fbid})")
        success = plot_gene(
            adata, gene, fbid,
            outdir=args.outdir,
            umap_key=args.umap_key,
            layer=args.layer,
            cmap=args.cmap,
            vmax_arg=args.vmax
        )
        if success:
            n_plotted += 1
        else:
            n_skipped += 1

    print(f"\nDone. {n_plotted} plots saved to '{args.outdir}/', {n_skipped} skipped.")


if __name__ == "__main__":
    main()
