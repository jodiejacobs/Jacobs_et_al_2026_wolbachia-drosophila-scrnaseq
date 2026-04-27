#!/usr/bin/env python3
"""
plot_custom_umap.py

Plot UMAP colored by gene expression for a list of genes from a config file.
Each gene gets its own 2x2 inch PDF with gene name and FlyBase ID in the title.

Supports adata objects where var_names are either:
  - Gene symbols (e.g. "Rab7")
  - FlyBase IDs (e.g. "FBgn0015795")
  - Any other index, with FlyBase IDs stored in a adata.var column

Use --use-raw to pull expression from adata.raw (all genes before HVG filtering).
Use --filter-obs to subset cells before plotting (e.g. --filter-obs method==pipseq).
UMAP coordinates and cell metadata always come from the main adata object.

Usage:
    python plot_custom_umap.py --h5ad <path.h5ad> --config <genes.csv> --outdir <output_dir>
    python plot_custom_umap.py --h5ad <path.h5ad> --config <genes.csv> --outdir <output_dir> --use-raw

    python scripts/method_comparison/plot_custom_umap.py \
        --h5ad results/integrated/integrated.h5ad \
        --config results/pseudotime_genes/wolbachia_infection/tradeseq_inputs/custom_genes.csv \
        --outdir results/custom_umaps/pipseq/ \
        --use-raw \
        --cmap viridis \
        --filter-obs method==pipseq

Config file format (TSV or CSV, with header):
    Gene, FlyBaseId
    Rab7, FBgn0015795
    sxl, FBgn0264270
    ...
"""

import argparse
import os
import re
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
    parser.add_argument("--use-raw", action="store_true",
                        help="Pull expression from adata.raw (pre-HVG all-gene matrix)")
    parser.add_argument("--layer", default=None,
                        help="Layer to use for expression (ignored if --use-raw; default: adata.X)")
    parser.add_argument("--cmap", default="magma",
                        help="Colormap for expression (default: magma)")
    parser.add_argument("--vmax", default="p99",
                        help="vmax for color scale: 'p99', 'p95', or a float (default: p99)")
    parser.add_argument("--filter-obs", default=None, metavar="COL==VALUE",
                        help="Subset cells by an obs column before plotting. "
                             "Format: column==value (e.g. method==pipseq). "
                             "Case-insensitive value matching.")
    return parser.parse_args()


def parse_filter(filter_str):
    """Parse 'col==value' into (col, value). Exits on bad format."""
    if filter_str is None:
        return None, None
    if "==" not in filter_str:
        sys.exit(f"ERROR: --filter-obs must be in the format col==value, got: '{filter_str}'")
    col, val = filter_str.split("==", 1)
    return col.strip(), val.strip()


def apply_filter(adata, col, val):
    """Subset adata to cells where obs[col] == val (case-insensitive)."""
    if col not in adata.obs.columns:
        sys.exit(
            f"ERROR: obs column '{col}' not found.\n"
            f"  Available columns: {list(adata.obs.columns)}"
        )
    obs_vals = adata.obs[col].astype(str).str.strip().str.lower()
    mask = obs_vals == val.lower()
    n_match = mask.sum()
    if n_match == 0:
        unique_vals = adata.obs[col].astype(str).unique().tolist()
        sys.exit(
            f"ERROR: No cells found where {col}=='{val}'.\n"
            f"  Unique values in '{col}': {unique_vals}"
        )
    subset = adata[mask].copy()
    print(f"  Filtered to {n_match} cells where {col}=='{val}' (dropped {adata.n_obs - n_match})")
    return subset


def load_config(config_path):
    """Load gene config. Accepts TSV or CSV with flexible column names.
    Handles space-padded headers like 'Gene, FlyBaseId'.
    """
    sep = "\t" if config_path.endswith(".tsv") else ","
    df = pd.read_csv(config_path, sep=sep, skipinitialspace=True)

    # Normalize column names
    df.columns = df.columns.str.strip().str.lower().str.replace(r"\s+", "_", regex=True)

    # Strip whitespace from all string values
    for col in df.select_dtypes(include="object").columns:
        df[col] = df[col].str.strip()

    name_aliases = ["gene_name", "gene", "name", "symbol"]
    fb_aliases   = ["flybase_id", "flybaseid", "fbgn", "flybase"]

    name_col = next((c for c in name_aliases if c in df.columns), None)
    fb_col   = next((c for c in fb_aliases   if c in df.columns), None)

    if name_col is None or fb_col is None:
        sys.exit(
            f"Config must have gene name and FlyBase ID columns.\n"
            f"  Recognized gene name columns : {name_aliases}\n"
            f"  Recognized FlyBase ID columns: {fb_aliases}\n"
            f"  Found columns: {list(df.columns)}"
        )

    return (
        df[[name_col, fb_col]]
        .rename(columns={name_col: "gene_name", fb_col: "flybase_id"})
        .dropna()
    )


def build_fbgn_lookup(var_names, var_df):
    """Return a dict mapping FlyBase ID -> var_name, or None if not applicable.

    Works on any var_names + var DataFrame (adata or adata.raw).
    Strategy (in order):
    1. var_names are FBgn IDs directly
    2. A column in var_df contains FBgn IDs
    3. No FBgn lookup available (fall back to gene symbol matching)
    """
    fbgn_pattern = re.compile(r"^FBgn\d+$")

    # Strategy 1: var_names are FBgn IDs
    if sum(bool(fbgn_pattern.match(v)) for v in var_names[:50]) > 25:
        print("  Detected: var_names are FlyBase IDs -- will look up genes by FBgn ID")
        return {v: v for v in var_names}

    # Strategy 2: a var column holds FBgn IDs
    for col in var_df.columns:
        sample = var_df[col].dropna().astype(str).head(50)
        if sum(bool(fbgn_pattern.match(v)) for v in sample) > 25:
            print(f"  Detected: FlyBase IDs in var['{col}'] -- building lookup")
            return dict(zip(var_df[col].astype(str), var_names))

    print("  No FBgn column detected -- will look up genes directly by symbol in var_names")
    return None


def resolve_vmax(expr_values, vmax_arg):
    """Return a numeric vmax."""
    if isinstance(vmax_arg, str):
        if vmax_arg == "p99":
            return float(expr_values.quantile(0.99))
        elif vmax_arg == "p95":
            return float(expr_values.quantile(0.95))
        else:
            return float(vmax_arg)
    return float(vmax_arg)


def safe_filename(name):
    """Replace filename-unsafe characters with underscores."""
    return re.sub(r"[^\w\-.]", "_", name)


def plot_gene(adata, expr_source, gene_name, flybase_id, var_key, outdir,
              umap_key, layer, use_raw, cmap, vmax_arg):
    """Generate and save a 2x2 inch PDF UMAP for one gene.

    adata:       main AnnData (for UMAP coords); may already be subset
    expr_source: adata or adata.raw (for expression)
    var_key:     key in expr_source.var_names to slice
    """
    import scipy.sparse as sp

    var_names = list(expr_source.var_names)
    if var_key not in var_names:
        print(f"  [WARN] '{var_key}' not found in var_names -- skipping.")
        return False

    # Extract expression
    if use_raw:
        idx = var_names.index(var_key)
        expr = expr_source.X[:, idx]
    elif layer is not None and layer in adata.layers:
        expr = adata[:, var_key].layers[layer]
    else:
        if layer is not None:
            print(f"  [WARN] Layer '{layer}' not found -- using adata.X.")
        expr = adata[:, var_key].X

    # Flatten to 1D
    if sp.issparse(expr):
        expr = expr.toarray().flatten()
    else:
        expr = expr.flatten()

    expr_series = pd.Series(expr)
    vmax = resolve_vmax(expr_series, vmax_arg)

    # UMAP coordinates always from main adata (already subset if filtered)
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
        vmax=max(vmax, 1e-6),  # avoid vmax=0 for all-zero genes
        s=0.5,
        linewidths=0,
        rasterized=True
    )

    ax.set_title(f"{gene_name}\n{flybase_id}", fontsize=5, pad=2)
    ax.set_xlabel("UMAP 1", fontsize=4)
    ax.set_ylabel("UMAP 2", fontsize=4)
    ax.tick_params(labelsize=3, length=2, pad=1)
    ax.set_aspect("equal", adjustable="datalim")

    cbar = fig.colorbar(sc_plot, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label("Expression", fontsize=4)
    cbar.ax.tick_params(labelsize=3)

    plt.tight_layout(pad=0.3)

    out_path = os.path.join(outdir, f"{safe_filename(gene_name)}_{flybase_id}.pdf")
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

    # Apply obs filter if requested (subsets cells, preserves UMAP coords)
    filter_col, filter_val = parse_filter(args.filter_obs)
    if filter_col is not None:
        adata = apply_filter(adata, filter_col, filter_val)

    # Decide expression source
    if args.use_raw:
        if adata.raw is None:
            sys.exit("ERROR: --use-raw specified but adata.raw is None.")
        expr_source = adata.raw
        print(f"  Using adata.raw: {adata.raw.X.shape[0]} cells x {adata.raw.X.shape[1]} genes")
        var_names = list(adata.raw.var_names)
        var_df    = adata.raw.var
    else:
        expr_source = adata
        var_names = list(adata.var_names)
        var_df    = adata.var

    print(f"  var_names sample: {var_names[:5]}")
    fbgn_lookup = build_fbgn_lookup(var_names, var_df)

    print(f"Loading config: {args.config}")
    genes_df = load_config(args.config)
    print(f"  {len(genes_df)} genes to plot")

    n_plotted = 0
    n_skipped = 0

    for _, row in genes_df.iterrows():
        gene = str(row["gene_name"]).strip()
        fbid = str(row["flybase_id"]).strip()

        # Resolve the actual var_names key to use for slicing
        if fbgn_lookup is not None:
            var_key = fbgn_lookup.get(fbid)
            if var_key is None:
                print(f"  [WARN] FlyBase ID '{fbid}' ({gene}) not found in expression source -- skipping.")
                n_skipped += 1
                continue
        else:
            var_key = gene

        print(f"  Plotting {gene} ({fbid})  [var_key={var_key}]")
        success = plot_gene(
            adata, expr_source, gene, fbid, var_key,
            outdir=args.outdir,
            umap_key=args.umap_key,
            layer=args.layer,
            use_raw=args.use_raw,
            cmap=args.cmap,
            vmax_arg=args.vmax,
        )
        if success:
            n_plotted += 1
        else:
            n_skipped += 1

    print(f"\nDone. {n_plotted} plots saved to '{args.outdir}/', {n_skipped} skipped.")


if __name__ == "__main__":
    main()