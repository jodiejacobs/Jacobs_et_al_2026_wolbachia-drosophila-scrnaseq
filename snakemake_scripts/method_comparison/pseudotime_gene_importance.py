#!/usr/bin/env python3
"""
pseudotime_gene_importance.py

Identifies genes important across SCEPTIC pseudotime using:
  1. Spearman correlation (fast, exploratory)
  2. tradeSeq GAM via embedded R subprocess (rigorous, publication-ready)
  3. cNMF program–pseudotime correlation (module-level)
  4. Binned differential expression (early vs late)
  5. Pseudotime-ordered heatmap of top dynamic genes

Usage:
    python pseudotime_gene_importance.py \\
        --h5ad results/sceptic/adata_with_pseudotime.h5ad \\
        --outdir results/pseudotime_genes/ \\
        [--skip-tradeseq]  # skip if R/tradeSeq not installed

Requirements (conda/mamba):
    mamba install -c conda-forge -c bioconda rpy2 bioconductor-tradeseq
    pip install statsmodels adjusttext

Author: Jodie Jacobs et al. 2026
"""

import argparse
import os
import sys
import subprocess
import tempfile
import warnings
warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
import scanpy as sc
from scipy import stats
from scipy.sparse import issparse
from statsmodels.stats.multitest import multipletests

# ─────────────────────────────────────────────────────────────────────────────
# CONFIG — edit if your column names differ
# ─────────────────────────────────────────────────────────────────────────────

PSEUDOTIME_KEY = "sceptic_pseudotime"
WEIGHTS_KEY    = "sceptic_weights"   # set None if not present
CLUSTER_KEY    = "leiden"
NMF_PREFIX     = "Program_"         # prefix for cNMF usage columns in .obs

# Empirical cluster order by median pseudotime (early → late)
# Derived from: adata.obs.groupby("leiden")["sceptic_pseudotime"].median().sort_values()
# Cluster  2  median=0.010   ← earliest
# Cluster  8  median=0.053
# Cluster  6  median=0.067
# Cluster  0  median=0.144
# Cluster  1  median=0.685
# Cluster 12  median=1.001
# Cluster 10  median=1.057
# Cluster  9  median=6.671
# Cluster 11  median=12.557
# Cluster 13  median=28.006
# Cluster  4  median=55.996
# Cluster 14  median=62.684
# Cluster  3  median=96.820
# Cluster  5  median=99.980
# Cluster  7  median=99.998   ← latest
ORDERED_CLUSTERS = ['2','8','6','0','1','12','10','9','11','13','4','14','3','5','7']

PADJ_THRESH    = 0.05
RHO_THRESH     = 0.25               # min |Spearman rho| to call a gene dynamic
MIN_EXPR_CELLS = 50                 # min cells with >0 counts to test a gene
N_TOP_GENES    = 200                # genes to carry through to downstream steps
N_KNOTS        = 6                  # tradeSeq GAM knots

# ─────────────────────────────────────────────────────────────────────────────
# HELPERS
# ─────────────────────────────────────────────────────────────────────────────

def get_counts(adata):
    """
    Return (count_matrix, gene_names) — cells x genes float32 + matching gene list.
    Priority: layers['counts'] > layers['raw'] > adata.raw.X > adata.X
    Accepts float matrices that are integer-valued (common after Scanpy QC).
    Rejects log-normalised matrices (negatives, non-trivial decimals).
    adata.raw may have more genes than adata.X — gene_names reflects this.
    """
    def _looks_like_counts(mat):
        s = mat[:min(500, mat.shape[0]), :min(500, mat.shape[1])]
        if issparse(s): s = s.toarray()
        if s.min() < 0:
            return False
        return np.allclose(s, np.round(s), atol=1e-3)

    raw_var = adata.raw.var_names if adata.raw is not None else None
    raw_X   = adata.raw.X         if adata.raw is not None else None

    candidates = [
        ("layers['counts']", adata.layers.get("counts"), adata.var_names),
        ("layers['raw']",    adata.layers.get("raw"),    adata.var_names),
        ("adata.raw.X",        raw_X,                      raw_var),
        ("adata.X",            adata.X,                    adata.var_names),
    ]

    for src_name, mat, var_names in candidates:
        if mat is None:
            continue
        if issparse(mat): mat = mat.toarray()
        if _looks_like_counts(mat):
            print(f"  Raw counts source: {src_name}  shape={mat.shape}  "
                  f"max={mat.max():.0f}  genes={len(var_names)}")
            return mat.astype(np.float32), list(var_names)

    # ── fallback: try reversing log1p on adata.raw.X ────────────────────────
    # Happens when sc.pp.log1p() was called before adata.raw was set.
    # expm1 + round recovers integer counts if normalisation was log1p(X/total*10000).
    # Note: this recovers relative counts (CPM-scaled), not original UMI counts,
    # which is acceptable for tradeSeq/Spearman but not UMI-level analyses.
    if raw_X is not None:
        if issparse(raw_X): raw_X = raw_X.toarray()
        recovered = np.expm1(raw_X)
        recovered_rounded = np.round(recovered)
        residual = np.abs(recovered - recovered_rounded).mean()
        print(f"  ⚠️  No integer matrix found. Attempting expm1(adata.raw.X) ...")
        print(f"     mean residual after rounding: {residual:.4f}  "
              f"max: {np.abs(recovered - recovered_rounded).max():.4f}")
        print(f"     recovered range: {recovered.min():.1f} – {recovered.max():.1f}")
        print(f"     Using expm1-recovered counts (CPM-scaled, not raw UMI).")
        return recovered_rounded.astype(np.float32), list(raw_var)

    raise ValueError(
        "Could not find or recover a count matrix.\n"
        f"  adata.X range: {adata.X.min():.2f} to {adata.X.max():.2f}\n"
        "  adata.raw.X is also log-normalised and expm1 recovery failed.\n"
        "  Fix: re-run pipeline storing raw UMI counts in adata.layers['counts'] "
        "before any normalisation."
    )


def filter_genes(counts, gene_names, min_cells=MIN_EXPR_CELLS):
    """Keep only genes expressed in >= min_cells cells."""
    gene_names = np.array(gene_names)
    expressed = (counts > 0).sum(axis=0) >= min_cells
    print(f"  Genes after min_cells filter ({min_cells}): {expressed.sum()} / {counts.shape[1]}")
    return counts[:, expressed], gene_names[expressed].tolist()


def save_csv(df, path, msg=""):
    df.to_csv(path)
    print(f"  Saved {msg}: {path}  ({len(df)} rows)")

# ─────────────────────────────────────────────────────────────────────────────
# STEP 1 — SPEARMAN CORRELATION
# ─────────────────────────────────────────────────────────────────────────────

def run_spearman(counts, gene_names, pseudotime, outdir):
    print("\n[1/5] Spearman correlation with pseudotime ...")
    pt = pseudotime.values
    n_cells, n_genes = counts.shape
    print(f"  Vectorized Spearman: {n_cells} cells × {n_genes} genes ...")

    # Vectorized Spearman: rank both pt and all genes at once, then correlate.
    # This is O(n_genes * n_cells) but in numpy — ~100x faster than a Python loop.
    from scipy.stats import rankdata

    pt_ranks = rankdata(pt).astype(np.float32)                  # (n_cells,)
    pt_ranks -= pt_ranks.mean()
    pt_std = pt_ranks.std()

    # Rank each gene column; work in chunks to keep memory under control
    CHUNK = 2000
    rhos = np.empty(n_genes, dtype=np.float32)

    for start in range(0, n_genes, CHUNK):
        end = min(start + CHUNK, n_genes)
        chunk = counts[:, start:end].astype(np.float32)          # (n_cells, chunk)
        ranked = np.apply_along_axis(rankdata, 0, chunk)         # rank each gene
        ranked -= ranked.mean(axis=0)
        gene_stds = ranked.std(axis=0)
        # Pearson on ranks = Spearman
        with np.errstate(invalid='ignore', divide='ignore'):
            rhos[start:end] = np.where(
                gene_stds > 0,
                (pt_ranks @ ranked) / (n_cells * pt_std * gene_stds),
                0.0
            )
        if (start // CHUNK) % 5 == 0:
            print(f"    {end}/{n_genes} genes done ...")

    # p-values: use t-distribution approximation (same as scipy.stats.spearmanr)
    # t = rho * sqrt((n-2) / (1 - rho^2))
    with np.errstate(invalid='ignore', divide='ignore'):
        t_stat = rhos * np.sqrt((n_cells - 2) / np.maximum(1 - rhos**2, 1e-14))
    from scipy.stats import t as t_dist
    pvals = 2 * t_dist.sf(np.abs(t_stat), df=n_cells - 2)
    rhos_list  = rhos.tolist()
    pvals_list = pvals.tolist()

    padj = multipletests(pvals_list, method="fdr_bh")[1]
    df = pd.DataFrame({
        "gene":    gene_names,
        "rho":     rhos_list,
        "pval":    pvals_list,
        "padj":    padj,
        "abs_rho": np.abs(rhos),
    }).set_index("gene").sort_values("abs_rho", ascending=False)

    sig = df[(df.padj < PADJ_THRESH) & (df.abs_rho >= RHO_THRESH)]
    print(f"  Significant (padj<{PADJ_THRESH}, |rho|>={RHO_THRESH}): {len(sig)} genes")

    save_csv(df,  os.path.join(outdir, "spearman_all.csv"),  "all Spearman results")
    save_csv(sig, os.path.join(outdir, "spearman_sig.csv"),  "significant Spearman genes")

    # Volcano-style rho vs -log10(padj)
    fig, ax = plt.subplots(figsize=(7, 5))
    colors = ["#c0392b" if (r > RHO_THRESH and p < PADJ_THRESH)
              else "#2980b9" if (r < -RHO_THRESH and p < PADJ_THRESH)
              else "#bdc3c7"
              for r, p in zip(df.rho, df.padj)]
    ax.scatter(df.rho, -np.log10(df.padj + 1e-300), c=colors, s=4, alpha=0.6, rasterized=True)
    ax.axvline( RHO_THRESH, ls="--", lw=0.8, c="gray")
    ax.axvline(-RHO_THRESH, ls="--", lw=0.8, c="gray")
    ax.axhline(-np.log10(PADJ_THRESH), ls="--", lw=0.8, c="gray")
    ax.set_xlabel("Spearman ρ with pseudotime")
    ax.set_ylabel("−log₁₀(padj)")
    ax.set_title("Gene–pseudotime correlation")
    plt.tight_layout()
    fig.savefig(os.path.join(outdir, "spearman_volcano.pdf"), dpi=150)
    plt.close()

    return sig.head(N_TOP_GENES).index.tolist(), df


# ─────────────────────────────────────────────────────────────────────────────
# STEP 2 — tradeSeq via embedded R script
# ─────────────────────────────────────────────────────────────────────────────

TRADESEQ_R = r"""
# ─── tradeSeq single-trajectory analysis ───────────────────────────────────
suppressPackageStartupMessages({
    library(tradeSeq)
    library(BiocParallel)
    library(Matrix)
})

args      <- commandArgs(trailingOnly = TRUE)
rds_file  <- args[1]   # RDS: list(counts=dgCMatrix genes×cells, pt=numeric, genes=char)
out_dir   <- args[2]
n_knots   <- as.integer(args[3])
n_workers <- as.integer(args[4])

# ── load data ────────────────────────────────────────────────────────────────
cat("  Loading counts from RDS ...
")
obj    <- readRDS(rds_file)
counts <- obj$counts   # sparse genes x cells (dgCMatrix)
pt_vec <- obj$pt       # named numeric pseudotime vector (cells)

# align cell order
shared  <- intersect(colnames(counts), names(pt_vec))
counts  <- counts[, shared]
pt_vec  <- pt_vec[shared]

pt_mat <- matrix(pt_vec, ncol = 1, dimnames = list(shared, "pseudotime"))
wt_mat <- matrix(1,      nrow = length(shared), ncol = 1,
                 dimnames = list(shared, "w1"))

cat(sprintf("  Genes: %d   Cells: %d   Workers: %d
",
            nrow(counts), ncol(counts), n_workers))

# ── fit GAMs ─────────────────────────────────────────────────────────────────
cat(sprintf("  Fitting GAMs (nknots=%d) ...
", n_knots))
BPPARAM <- MulticoreParam(workers = n_workers, progressbar = TRUE)

sce <- fitGAM(counts      = counts,
              pseudotime  = pt_mat,
              cellWeights = wt_mat,
              nknots      = n_knots,
              parallel    = TRUE,
              BPPARAM     = BPPARAM,
              verbose     = FALSE)

# save SCE so tests can be re-run without re-fitting
saveRDS(sce, file.path(out_dir, "tradeseq_sce.rds"))
cat("  Saved SCE to tradeseq_sce.rds
")

# ── association test ──────────────────────────────────────────────────────────
cat("  Running associationTest ...
")
assoc      <- as.data.frame(associationTest(sce, lineages = FALSE))
assoc$gene <- rownames(assoc)
assoc$padj <- p.adjust(assoc$pvalue, method = "BH")
assoc      <- assoc[order(assoc$waldStat, decreasing = TRUE), ]
write.csv(assoc, file.path(out_dir, "tradeseq_association.csv"))
cat(sprintf("  Association test: %d significant genes (padj<0.05)
",
            sum(assoc$padj < 0.05, na.rm = TRUE)))

# ── start vs end test ─────────────────────────────────────────────────────────
cat("  Running startVsEndTest ...
")
sve      <- as.data.frame(startVsEndTest(sce))
sve$gene <- rownames(sve)
sve$padj <- p.adjust(sve$pvalue, method = "BH")
sve      <- sve[order(sve$waldStat, decreasing = TRUE), ]
write.csv(sve, file.path(out_dir, "tradeseq_startvsend.csv"))

# ── smooth predictions for top dynamic genes ──────────────────────────────────
cat("  Predicting smooth curves for top genes ...
")
top_genes <- head(assoc$gene[assoc$padj < 0.05], 300)
if (length(top_genes) > 0) {
    yhat <- predictSmooth(sce, gene = top_genes, nPoints = 100, tidy = FALSE)
    write.csv(yhat, file.path(out_dir, "tradeseq_smooth_predictions.csv"))
    cat(sprintf("  Saved smooth predictions for %d genes
", length(top_genes)))
}

cat("  tradeSeq complete.
")
"""


def prepare_tradeseq_inputs(counts, gene_names, cell_names, pseudotime, outdir,
                             spearman_df=None):
    """
    Write tradeSeq input CSVs (counts + pseudotime) for later use with
    run_tradeseq.R. Does NOT run R — call this script first, then run:

        Rscript scripts/method_comparison/run_tradeseq.R \
            --counts  {outdir}/tradeseq_inputs/counts_genesXcells.csv \
            --pt      {outdir}/tradeseq_inputs/pseudotime.csv \
            --outdir  {outdir} \
            --nknots  6 \
            --nworkers 16
    """
    print("\n[2/5] Preparing tradeSeq inputs ...")

    ts_dir = os.path.join(outdir, "tradeseq_inputs")
    os.makedirs(ts_dir, exist_ok=True)

    # ── Pre-filter genes using Spearman results ──────────────────────────────
    # Fitting GAMs on all 17k genes is prohibitively slow.
    # Use Spearman sig genes + top candidates by |rho|, capped at 5000.
    gene_names_arr = np.array(gene_names)
    if spearman_df is not None:
        sig_genes  = set(spearman_df[spearman_df.padj < PADJ_THRESH].index)
        top_genes  = set(spearman_df.nlargest(3000, "abs_rho").index)
        keep_set   = (sig_genes | top_genes) & set(gene_names)
        keep_idx   = np.array([i for i, g in enumerate(gene_names) if g in keep_set])
        if len(keep_idx) == 0:
            keep_idx = np.arange(min(5000, len(gene_names)))
        keep_idx      = keep_idx[:5000]
        counts_ts     = counts[:, keep_idx]
        gene_names_ts = gene_names_arr[keep_idx].tolist()
        print(f"  Pre-filtered to {len(gene_names_ts)} genes "
              f"({len(sig_genes)} Spearman sig + top 3000 by |rho|, cap=5000)")
    else:
        counts_ts     = counts
        gene_names_ts = gene_names
        print(f"  No Spearman filter — using all {len(gene_names_ts)} genes")

    # ── Write counts CSV (genes × cells) ────────────────────────────────────
    cnt_path = os.path.join(ts_dir, "counts_genesXcells.csv")
    pt_path  = os.path.join(ts_dir, "pseudotime.csv")

    print(f"  Writing counts: {len(gene_names_ts)} genes × {counts_ts.shape[0]} cells ...")
    pd.DataFrame(
        counts_ts.T.astype(int),
        index=gene_names_ts,
        columns=cell_names,
    ).to_csv(cnt_path)

    pd.DataFrame(
        {"pseudotime": pseudotime.values},
        index=pseudotime.index,
    ).to_csv(pt_path)

    print(f"  counts  -> {cnt_path}  ({os.path.getsize(cnt_path) / 1e6:.0f} MB)")
    print(f"  pt      -> {pt_path}")
    print(f"\n  To run tradeSeq:")
    print(f"    Rscript scripts/method_comparison/run_tradeseq.R \\")
    print(f"        --counts  {cnt_path} \\")
    print(f"        --pt      {pt_path} \\")
    print(f"        --outdir  {outdir} \\")
    print(f"        --nknots  6 \\")
    print(f"        --nworkers 16")

    # ── Load tradeSeq results if they already exist ──────────────────────────
    assoc_path = os.path.join(outdir, "tradeseq_association.csv")
    sve_path   = os.path.join(outdir, "tradeseq_startvsend.csv")

    if os.path.exists(assoc_path) and os.path.exists(sve_path):
        print(f"\n  Found existing tradeSeq results — loading ...")
        assoc_df = pd.read_csv(assoc_path, index_col=0)
        sve_df   = pd.read_csv(sve_path,   index_col=0)
        sig_assoc = assoc_df[assoc_df.padj < PADJ_THRESH]
        print(f"  Association test sig genes: {len(sig_assoc)}")
        return assoc_df, sve_df
    else:
        print(f"\n  No tradeSeq results found yet — run the Rscript above first.")
        return None, None


# ─────────────────────────────────────────────────────────────────────────────
# STEP 3 — cNMF program × pseudotime correlation
# ─────────────────────────────────────────────────────────────────────────────

def run_cnmf_correlation(adata, pseudotime, program_dir, outdir):
    """
    Correlates each cNMF program's per-cell usage with pseudotime.
    Then pulls top genes for dynamic programs from Program_N_genes.csv files.

    program_dir: directory containing Program_0_genes.csv ... Program_N_genes.csv
                 (output of nmf_programs rule)
    """
    print("\n[3/5] cNMF program–pseudotime correlation ...")

    # Find program usage columns in adata.obs
    prog_cols = [c for c in adata.obs.columns if c.startswith(NMF_PREFIX)]
    if not prog_cols:
        print(f"  ⚠️  No columns with prefix '{NMF_PREFIX}' found in adata.obs — skipping.")
        return [], None

    usage = adata.obs[prog_cols].copy()
    pt    = pseudotime.reindex(usage.index)

    rows = []
    for prog in prog_cols:
        r, p = stats.spearmanr(pt, usage[prog])
        rows.append({"program": prog, "rho": r, "pval": p})

    prog_df = pd.DataFrame(rows).set_index("program")
    prog_df["padj"]    = multipletests(prog_df.pval, method="fdr_bh")[1]
    prog_df["abs_rho"] = prog_df.rho.abs()
    prog_df = prog_df.sort_values("abs_rho", ascending=False)

    save_csv(prog_df, os.path.join(outdir, "cnmf_program_pseudotime_corr.csv"),
             "cNMF program correlations")

    # Bar chart of program correlations
    fig, ax = plt.subplots(figsize=(8, 4))
    colors = ["#c0392b" if r > 0 else "#2980b9" for r in prog_df.rho]
    ax.bar(prog_df.index, prog_df.rho, color=colors)
    ax.axhline(0, color="black", lw=0.8)
    ax.set_xticklabels(prog_df.index, rotation=45, ha="right", fontsize=8)
    ax.set_ylabel("Spearman ρ with pseudotime")
    ax.set_title("cNMF program correlation with pseudotime")
    plt.tight_layout()
    fig.savefig(os.path.join(outdir, "cnmf_program_pseudotime.pdf"), dpi=150)
    plt.close()

    # Identify dynamic programs
    dynamic_progs = prog_df[prog_df.padj < PADJ_THRESH].index.tolist()
    if not dynamic_progs:
        dynamic_progs = prog_df.head(3).index.tolist()
        print(f"  No programs pass padj threshold — using top 3 by |rho|")
    print(f"  Dynamic programs: {dynamic_progs}")

    # Load top genes from Program_N_genes.csv files
    # Program_N_genes.csv expected to have a column of gene names (first col or 'gene')
    top_cnmf_genes = []
    if program_dir and os.path.isdir(program_dir):
        for prog in dynamic_progs:
            # prog is e.g. "Program_3" — map to "Program_3_genes.csv"
            csv_path = os.path.join(program_dir, f"{prog}_genes.csv")
            if os.path.exists(csv_path):
                df = pd.read_csv(csv_path)
                # gene names are in first column, whatever it's called
                genes = df.iloc[:, 0].dropna().tolist()
                genes = [g for g in genes if not g.startswith("GQX67")]
                top_cnmf_genes.extend(genes)
                print(f"    {prog}: {len(genes)} genes from {csv_path}")
            else:
                print(f"    ⚠️  {csv_path} not found — skipping {prog}")
        top_cnmf_genes = list(dict.fromkeys(top_cnmf_genes))  # dedup, preserve order
        print(f"  Total unique cNMF genes from dynamic programs: {len(top_cnmf_genes)}")
    else:
        print("  No --program-dir provided — skipping gene list extraction.")

    return top_cnmf_genes, prog_df


# ─────────────────────────────────────────────────────────────────────────────
# STEP 4 — BINNED DIFFERENTIAL EXPRESSION (early vs late)
# ─────────────────────────────────────────────────────────────────────────────

def run_binned_de(adata, pseudotime, outdir):
    """
    Bin pseudotime into quartiles, run Wilcoxon early (Q1) vs late (Q4).
    """
    print("\n[4/5] Binned DE (early Q1 vs late Q4) ...")

    adata.obs["pt_quartile"] = pd.qcut(
        pseudotime.reindex(adata.obs_names),
        q=4,
        labels=["Q1_early", "Q2", "Q3", "Q4_late"]
    )

    # subset to Q1 and Q4 only
    sub = adata[adata.obs.pt_quartile.isin(["Q1_early", "Q4_late"])].copy()
    sc.tl.rank_genes_groups(
        sub,
        groupby="pt_quartile",
        groups=["Q4_late"],
        reference="Q1_early",
        method="wilcoxon",
        key_added="pt_de",
        tie_correct=True,
    )

    de = sc.get.rank_genes_groups_df(sub, group="Q4_late", key="pt_de")
    de = de.rename(columns={"names": "gene"}).set_index("gene")
    de["padj"] = multipletests(de.pvals, method="fdr_bh")[1]
    de_sig = de[(de.padj < PADJ_THRESH) & (de.scores.abs() > 1)]
    de_sig = de_sig.sort_values("scores", ascending=False)

    save_csv(de,     os.path.join(outdir, "binned_de_all.csv"),  "all binned DE")
    save_csv(de_sig, os.path.join(outdir, "binned_de_sig.csv"),  "sig binned DE genes")
    print(f"  Sig upregulated in late: {(de_sig.scores > 0).sum()}")
    print(f"  Sig downregulated in late: {(de_sig.scores < 0).sum()}")

    # Dot plot of top 10 up/down
    top_up   = de_sig[de_sig.scores > 0].head(10).index.tolist()
    top_down = de_sig[de_sig.scores < 0].tail(10).index.tolist()
    show_genes = top_up + top_down
    if show_genes:
        sc.pl.dotplot(sub, var_names=show_genes, groupby="pt_quartile",
                      show=False, save=False)
        plt.savefig(os.path.join(outdir, "binned_de_dotplot.pdf"),
                    bbox_inches="tight", dpi=150)
        plt.close()

    return de_sig.index.tolist(), de


# ─────────────────────────────────────────────────────────────────────────────
# STEP 5 — PSEUDOTIME HEATMAP
# ─────────────────────────────────────────────────────────────────────────────

def make_pseudotime_heatmap(adata, pseudotime, gene_lists, smooth_path, outdir):
    """
    Heatmap of top dynamic genes with cells ordered by pseudotime.
    If tradeSeq smooth predictions exist, use those (cleaner signal).
    Otherwise, use log-normalised expression from adata.
    """
    print("\n[5/5] Pseudotime heatmap ...")

    # Collect and prioritise genes: tradeSeq > Spearman > DE > cNMF
    all_genes = []
    for gl in gene_lists:
        all_genes.extend(gl)
    # dedup, preserve priority order
    seen, unique_genes = set(), []
    for g in all_genes:
        if g not in seen:
            seen.add(g)
            unique_genes.append(g)

    # Keep only genes present in adata
    avail = set(adata.var_names)
    plot_genes = [g for g in unique_genes if g in avail][:100]

    if len(plot_genes) == 0:
        print("  ⚠️  No genes to plot — skipping heatmap.")
        return

    print(f"  Plotting {len(plot_genes)} genes")

    if smooth_path and os.path.exists(smooth_path):
        # tradeSeq smooth: rows = genes, cols = pseudotime points
        smooth = pd.read_csv(smooth_path, index_col=0)
        smooth_genes = [g for g in plot_genes if g in smooth.index]
        if smooth_genes:
            mat = smooth.loc[smooth_genes]
            mat_scaled = pd.DataFrame(
                (mat.values - mat.values.mean(axis=1, keepdims=True))
                / (mat.values.std(axis=1, keepdims=True) + 1e-9),
                index=mat.index,
                columns=mat.columns,
            )
            fig, ax = plt.subplots(figsize=(12, max(4, len(smooth_genes) * 0.12)))
            sns.heatmap(mat_scaled, ax=ax,
                        cmap="RdBu_r", center=0,
                        xticklabels=False, yticklabels=(len(smooth_genes) < 80),
                        rasterized=True, cbar_kws={"label": "Scaled expression"})
            ax.set_xlabel("Pseudotime →")
            ax.set_title("Dynamic genes (tradeSeq smooth)")
            plt.tight_layout()
            fig.savefig(os.path.join(outdir, "heatmap_tradeseq_smooth.pdf"),
                        dpi=150, bbox_inches="tight")
            plt.close()
            print(f"  Saved tradeSeq smooth heatmap")

    # Cell-level heatmap ordered by pseudotime
    pt_order = pseudotime.sort_values().index
    sub = adata[pt_order, :]
    sc.pp.normalize_total(sub, target_sum=1e4, inplace=True)
    sc.pp.log1p(sub)

    avail_plot = [g for g in plot_genes if g in sub.var_names]
    X = sub[:, avail_plot].X
    if issparse(X): X = X.toarray()

    # Smooth with rolling window to reduce noise
    window = max(1, X.shape[0] // 100)
    X_smooth = pd.DataFrame(X).rolling(window, center=True, min_periods=1).mean().values
    X_scaled = (X_smooth - X_smooth.mean(axis=0)) / (X_smooth.std(axis=0) + 1e-9)

    fig, ax = plt.subplots(figsize=(14, max(5, len(avail_plot) * 0.14)))
    sns.heatmap(X_scaled.T, ax=ax,
                cmap="RdBu_r", center=0,
                xticklabels=False,
                yticklabels=(len(avail_plot) < 80),
                rasterized=True,
                cbar_kws={"label": "Scaled expression (rolling avg)"})
    ax.set_xlabel("Cells ordered by pseudotime →")
    ax.set_title(f"Top {len(avail_plot)} dynamic genes")
    plt.tight_layout()
    fig.savefig(os.path.join(outdir, "heatmap_cells_by_pseudotime.pdf"),
                dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  Saved cell-level heatmap")


# ─────────────────────────────────────────────────────────────────────────────
# SUMMARY TABLE
# ─────────────────────────────────────────────────────────────────────────────

def make_summary_table(spearman_df, tradeseq_assoc, de_df, cnmf_genes, outdir):
    """
    Merge all results into a single ranked summary table.
    Gene gets a score = number of methods that call it significant.
    """
    print("\n  Building summary table ...")

    all_genes = set()
    method_sets = {}

    if spearman_df is not None:
        sp_sig = set(spearman_df[(spearman_df.padj < PADJ_THRESH) &
                                  (spearman_df.abs_rho >= RHO_THRESH)].index)
        method_sets["spearman"] = sp_sig
        all_genes |= sp_sig

    if tradeseq_assoc is not None:
        ts_sig = set(tradeseq_assoc[tradeseq_assoc.padj < PADJ_THRESH].index)
        method_sets["tradeseq"] = ts_sig
        all_genes |= ts_sig

    if de_df is not None:
        de_sig = set(de_df[(de_df.padj < PADJ_THRESH) &
                            (de_df.scores.abs() > 1)].index)
        method_sets["binned_de"] = de_sig
        all_genes |= de_sig

    if cnmf_genes:
        method_sets["cnmf"] = set(cnmf_genes)
        all_genes |= set(cnmf_genes)

    rows = []
    for g in all_genes:
        row = {"gene": g}
        for m, s in method_sets.items():
            row[m] = g in s
        row["n_methods"] = sum(row[m] for m in method_sets)

        if spearman_df is not None and g in spearman_df.index:
            row["spearman_rho"]  = spearman_df.loc[g, "rho"]
            row["spearman_padj"] = spearman_df.loc[g, "padj"]

        if tradeseq_assoc is not None and g in tradeseq_assoc.index:
            row["tradeseq_waldStat"] = tradeseq_assoc.loc[g, "waldStat"]
            row["tradeseq_padj"]     = tradeseq_assoc.loc[g, "padj"]

        rows.append(row)

    summary = pd.DataFrame(rows).set_index("gene")
    summary = summary.sort_values("n_methods", ascending=False)
    save_csv(summary, os.path.join(outdir, "summary_dynamic_genes.csv"),
             "summary table")

    # Upset-style bar chart of method overlaps
    from itertools import combinations
    method_names = list(method_sets.keys())
    print(f"\n  Method overlap:")
    for i, m in enumerate(method_names):
        n = len(method_sets[m])
        print(f"    {m}: {n} genes")
    for a, b in combinations(method_names, 2):
        overlap = method_sets[a] & method_sets[b]
        print(f"    {a} ∩ {b}: {len(overlap)}")

    return summary


# ─────────────────────────────────────────────────────────────────────────────
# MAIN
# ─────────────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(description="Pseudotime gene importance pipeline")
    parser.add_argument("--h5ad",           required=True,  help="Input h5ad file")
    parser.add_argument("--outdir",         default="results/pseudotime_genes/",
                        help="Output directory")
    parser.add_argument("--program-dir",    default=None,
                        help="Directory containing Program_N_genes.csv files from NMF step")
    parser.add_argument("--skip-tradeseq",  action="store_true",
                        help="Skip tradeSeq (if R/tradeSeq not installed)")
    parser.add_argument("--pseudotime-key", default=PSEUDOTIME_KEY,
                        help=f"obs column for pseudotime (default: {PSEUDOTIME_KEY})")
    parser.add_argument("--n-knots",        type=int, default=N_KNOTS,
                        help=f"tradeSeq GAM knots (default: {N_KNOTS})")
    parser.add_argument("--n-threads",      type=int, default=8,
                        help="Threads for tradeSeq BiocParallel (default: 8, match SLURM --cpus)")
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    # ── load ────────────────────────────────────────────────────────────────
    print(f"\nLoading: {args.h5ad}")
    adata = sc.read_h5ad(args.h5ad)
    print(f"  Cells: {adata.n_obs}  Genes: {adata.n_vars}")

    # ── validate pseudotime ──────────────────────────────────────────────────
    if args.pseudotime_key not in adata.obs.columns:
        print(f"\n⚠️  Column '{args.pseudotime_key}' not found in adata.obs")
        print(f"  Available obs columns containing 'pseudo' or 'sceptic':")
        cands = [c for c in adata.obs.columns
                 if any(x in c.lower() for x in ["pseudo", "sceptic", "palantir", "dpt"])]
        print(f"  {cands}")
        sys.exit(1)

    pseudotime = adata.obs[args.pseudotime_key].dropna()
    print(f"  Pseudotime range: {pseudotime.min():.3f} – {pseudotime.max():.3f}")
    print(f"  Cells with pseudotime: {len(pseudotime)} / {adata.n_obs}")

    # Work only on cells with pseudotime
    adata_pt = adata[pseudotime.index].copy()

    # ── raw counts ───────────────────────────────────────────────────────────
    counts, gene_names = get_counts(adata_pt)
    counts_filt, gene_names = filter_genes(counts, gene_names)

    # Remove bacterial (Wolbachia) transcripts — prefixed with GQX67
    n_before = len(gene_names)
    gene_names_arr = np.array(gene_names)
    host_mask = np.array([not g.startswith("GQX67") for g in gene_names])
    counts_filt = counts_filt[:, host_mask]
    gene_names  = gene_names_arr[host_mask].tolist()
    print(f"  Removed {n_before - len(gene_names)} bacterial transcripts (GQX67*), "
          f"{len(gene_names)} host genes remaining")

    # ── run analyses ─────────────────────────────────────────────────────────
    spearman_genes, spearman_df = run_spearman(
        counts_filt, gene_names, pseudotime.reindex(adata_pt.obs_names), args.outdir
    )

    tradeseq_assoc, tradeseq_sve = None, None
    tradeseq_genes = []
    if not args.skip_tradeseq:
        tradeseq_assoc, tradeseq_sve = prepare_tradeseq_inputs(
            counts_filt, gene_names,
            adata_pt.obs_names.tolist(),
            pseudotime.reindex(adata_pt.obs_names),
            args.outdir,
            spearman_df = spearman_df,
        )
        if tradeseq_assoc is not None:
            tradeseq_genes = (tradeseq_assoc[tradeseq_assoc.padj < PADJ_THRESH]
                              .head(N_TOP_GENES).index.tolist())

    cnmf_genes, prog_df = run_cnmf_correlation(
        adata_pt,
        pseudotime.reindex(adata_pt.obs_names),
        args.program_dir,
        args.outdir
    )

    de_genes, de_df = run_binned_de(adata_pt, pseudotime.reindex(adata_pt.obs_names),
                                    args.outdir)

    smooth_path = os.path.join(args.outdir, "tradeseq_smooth_predictions.csv")
    make_pseudotime_heatmap(
        adata_pt,
        pseudotime.reindex(adata_pt.obs_names),
        gene_lists=[tradeseq_genes, spearman_genes, de_genes, cnmf_genes],
        smooth_path=smooth_path,
        outdir=args.outdir
    )

    summary = make_summary_table(spearman_df, tradeseq_assoc, de_df, cnmf_genes,
                                 args.outdir)

    print(f"\n✓ Done. Top dynamic genes (by method count):")
    print(summary.head(20).to_string())
    print(f"\nAll results in: {args.outdir}")


if __name__ == "__main__":
    main()
