"""
plot_gene_by_cellcycle.py
=========================
Violin + ridge plots of gene expression by Cyclum cell-cycle phase.

Usage
-----
python plot_gene_by_cellcycle.py --gene FBgn0000404
python plot_gene_by_cellcycle.py --gene FBgn0000404 --out CycA
"""

import argparse
import sys
import numpy as np
import matplotlib.pyplot as plt
import scipy.sparse
import scanpy as sc


# ── colours ───────────────────────────────────────────────────────────────────

CC_ORDER  = ["g0/g1", "s", "g2/m"]
CC_COLORS = {"g0/g1": "#ea546f", "s": "#7bcdca", "g2/m": "#2d9bb4"}

DEFAULT_H5AD = (
    "/private/groups/russelllab/jodie/scRNAseq/"
    "Jacobs_et_al_2026_wolbachia-drosophila-scrnaseq/"
    "cell_culture_system/results/integrated/"
    "integrated_uninfected_with_cellcycle.h5ad"
)


# ── helpers ───────────────────────────────────────────────────────────────────

def find_gene(adata, gene_id):
    if gene_id in adata.var_names:
        return gene_id
    if adata.raw is not None and gene_id in adata.raw.var_names:
        return gene_id
    # case-insensitive fallback
    for vnames in ([adata.var_names] +
                   ([adata.raw.var_names] if adata.raw is not None else [])):
        lmap = {v.lower(): v for v in vnames}
        if gene_id.lower() in lmap:
            return lmap[gene_id.lower()]
    return None


def get_expression(adata, gene_name):
    """Pull log-normalised expression from adata.raw.X (already log-normed)."""
    if adata.raw is not None and gene_name in adata.raw.var_names:
        idx = list(adata.raw.var_names).index(gene_name)
        x   = adata.raw.X[:, idx]
        if scipy.sparse.issparse(x):
            x = np.asarray(x.todense()).flatten()
        return x.astype(np.float32)
    # fallback: adata.X
    idx = list(adata.var_names).index(gene_name)
    x   = adata.X[:, idx]
    if scipy.sparse.issparse(x):
        x = np.asarray(x.todense()).flatten()
    return x.astype(np.float32)


# ── plots ─────────────────────────────────────────────────────────────────────

def make_violin(expr, stages, gene, out_path):
    observed = [s for s in CC_ORDER if s in stages.unique()]
    observed += [s for s in stages.unique() if s not in CC_ORDER]
    colors   = [CC_COLORS.get(s, "#aaaaaa") for s in observed]
    data     = [expr[stages == s] for s in observed]

    # stats
    print(f"\n{gene} expression by CC phase:")
    for s, d in zip(observed, data):
        print(f"  {s:8s}  n={len(d)}  median={np.median(d):.3f}  "
              f"pct_expressing={100*(d>0).mean():.1f}%")

    fig = plt.figure(figsize=(1.15, 0.9))
    ax  = fig.add_axes([0, 0, 1, 1])

    # strip behind violins
    rng = np.random.default_rng(42)
    for i, (d, c) in enumerate(zip(data, colors), 1):
        jitter = rng.uniform(-0.15, 0.15, size=len(d))
        ax.scatter(np.full(len(d), i) + jitter, d,
                   color=c, s=0.3, alpha=0.3, linewidths=0,
                   rasterized=True, zorder=2)

    parts = ax.violinplot(data, positions=range(1, len(observed)+1),
                          widths=0.6, showmedians=True, showextrema=False)
    for body, c in zip(parts["bodies"], colors):
        body.set_facecolor(c); body.set_edgecolor(c)
        body.set_alpha(0.7); body.set_zorder(3)
    parts["cmedians"].set_color("white")
    parts["cmedians"].set_linewidth(1.2)
    parts["cmedians"].set_zorder(4)

    ax.set_xlim(0.3, len(observed) + 0.7)
    ax.set_xticks(range(1, len(observed)+1))
    ax.set_xticklabels([s.upper() for s in observed], fontsize=6)
    ax.set_ylabel("log-norm expression", fontsize=6)
    ax.set_title(gene, fontsize=7, fontweight="bold")
    ax.tick_params(labelsize=6, length=2)
    ax.spines[["top", "right"]].set_visible(False)

    plt.savefig(out_path, bbox_inches="tight", dpi=300)
    plt.close()
    print(f"Saved -> {out_path}")


def make_ridgeplot(expr, pseudotime, stages, gene, out_path, n_bins=50):
    """Mean expression per pseudotime bin, one ridge per CC phase."""
    from scipy.ndimage import gaussian_filter1d

    observed = [s for s in CC_ORDER if s in stages.unique()]
    observed += [s for s in stages.unique() if s not in CC_ORDER]

    pt_min, pt_max = float(pseudotime.min()), float(pseudotime.max())
    bins        = np.linspace(pt_min, pt_max, n_bins + 1)
    bin_centers = (bins[:-1] + bins[1:]) / 2

    fig, axes = plt.subplots(len(observed), 1,
                             figsize=(2.0, 0.55 * len(observed)),
                             sharex=True)
    if len(observed) == 1:
        axes = [axes]

    for ax, stage in zip(axes, observed):
        mask     = (stages == stage).values
        bin_idx  = np.clip(np.digitize(pseudotime[mask], bins) - 1, 0, n_bins-1)
        expr_s   = expr[mask]
        bin_mean = np.array([
            expr_s[bin_idx == b].mean() if (bin_idx == b).any() else 0.0
            for b in range(n_bins)
        ])
        smoothed = gaussian_filter1d(bin_mean, sigma=2)

        color = CC_COLORS.get(stage, "#aaaaaa")
        ax.fill_between(bin_centers, smoothed, alpha=0.75, color=color)
        ax.plot(bin_centers, smoothed, color=color, lw=0.8)
        ax.set_xlim(pt_min, pt_max)
        ax.set_ylim(bottom=0)
        ax.set_yticks([])
        ax.spines[["top", "right", "left"]].set_visible(False)
        ax.set_ylabel(stage.upper(), fontsize=6, rotation=0,
                      labelpad=2, va="center")
        ax.yaxis.set_label_position("right")

    tick_vals = np.linspace(pt_min, pt_max, 5)
    axes[-1].set_xticks(tick_vals)
    axes[-1].set_xticklabels([f"{v:.1f}" for v in tick_vals], fontsize=6)
    axes[-1].set_xlabel("Cyclum pseudotime", fontsize=6)
    axes[0].set_title(gene, fontsize=7, fontweight="bold")

    plt.tight_layout(h_pad=0.1)
    plt.savefig(out_path, bbox_inches="tight", dpi=300)
    plt.close()
    print(f"Saved -> {out_path}")


# ── main ──────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--gene",           required=True)
    parser.add_argument("--h5ad",           default=DEFAULT_H5AD)
    parser.add_argument("--stage_col",      default="cyclum_stage")
    parser.add_argument("--pseudotime_col", default="cyclum_pseudotime")
    parser.add_argument("--out",            default=None,
                        help="Output stem; _violin.pdf and _ridge.pdf appended")
    args = parser.parse_args()

    print(f"Loading {args.h5ad} ...")
    adata = sc.read_h5ad(args.h5ad)

    gene_name = find_gene(adata, args.gene)
    if gene_name is None:
        print(f"ERROR: '{args.gene}' not found. First 10 var_names: {list(adata.var_names[:10])}")
        sys.exit(1)
    print(f"  Gene: {gene_name}")

    stages = adata.obs[args.stage_col].astype(str).str.strip().str.lower()
    expr   = get_expression(adata, gene_name)
    stem   = (args.out or args.gene).rstrip(".pdf")

    make_violin(expr, stages, gene_name, f"{stem}_violin.pdf")

    if args.pseudotime_col in adata.obs.columns:
        pseudotime = adata.obs[args.pseudotime_col].values.astype(np.float32)
        make_ridgeplot(expr, pseudotime, stages, gene_name, f"{stem}_ridge.pdf")
    else:
        print(f"WARNING: '{args.pseudotime_col}' not in obs, skipping ridge plot.")


if __name__ == "__main__":
    main()