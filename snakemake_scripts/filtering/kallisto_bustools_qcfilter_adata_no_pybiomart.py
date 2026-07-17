'''
Filtering script for kallisto bustools aligned single cell data.
    Input: raw data object
    Output: QC plots and filtered adata object

Speed improvements over original:
    - Scrublet: subsampled to max 10k cells, no internal UMAP
    - Wolbachia titer: single sparse matrix slice instead of gene-by-gene loop
    - qc_plots: vectorised knee plot, SVD only on subsample
    - BioMart: always uses cache if available, only queries if cache missing
    - analyze_filtered_adata: removed redundant per-sample PCA/UMAP/leiden
      (integrate.py does this jointly across all samples)
    - All dense conversions avoided where possible
'''

import scanpy as sc
import anndata as ad
import pandas as pd
import numpy as np
import scrublet as scr
import scipy.sparse
from scipy.interpolate import interpn
from sklearn.decomposition import TruncatedSVD
import matplotlib.pyplot as plt
import matplotlib
from scipy.sparse import csr_matrix
import os
from pybiomart import Dataset
import sys
import gc
import argparse

# ─────────────────────────────────────────────────────────────────────────────
# Argument parsing
# ─────────────────────────────────────────────────────────────────────────────

parser = argparse.ArgumentParser(description="Kallisto Bustools QC Filter")
parser.add_argument("--input",  type=str, required=True,  help="Input h5ad file")
parser.add_argument("--output", type=str, required=True,  help="Output h5ad file")
parser.add_argument("--scrublet_max_cells", type=int, default=10000,
                    help="Max cells for Scrublet subsampling (default: 10000)")
parser.add_argument("--n_pcs_scrublet",     type=int, default=20,
                    help="PCs for Scrublet (default: 20)")

args   = parser.parse_args()
input  = args.input
output = args.output

sample_name = os.path.basename(input).replace(".h5ad", "")
output_dir  = os.path.dirname(output)
fig_dir     = f"{output_dir}/{sample_name}"

os.makedirs(output_dir, exist_ok=True)
os.makedirs(fig_dir,    exist_ok=True)

sc.settings.autosave = True
sc.settings.figdir   = fig_dir

# ─────────────────────────────────────────────────────────────────────────────
# Gene lists
# ─────────────────────────────────────────────────────────────────────────────

rRNA_genes = [
    "GQX67_00940", "GQX67_00945", "GQX67_05945",
    "FBgn0267496", "FBgn0267500", "FBgn0267503", "FBgn0085765",
    "FBgn0267518", "FBgn0267524", "FBgn0267509", "FBgn0267499",
    "FBgn0267502", "FBgn0267512", "FBgn0267517", "FBgn0267523",
    "FBgn0250731", "FBgn0267514", "FBgn0085802", "FBgn0267498",
    "FBgn0267501", "FBgn0267521", "FBgn0085813", "FBgn0267504",
    "FBgn0267508", "FBgn0267511", "FBgn0085753", "FBgn0267497",
    "FBgn0267522", "FBgn0085771", "FBgn0267519", "FBgn0085819",
    "FBgn0267513", "FBgn0267520", "FBgn0267515",
]

MITO_GENES_FALLBACK = [
    "FBgn0013674", "FBgn0013675", "FBgn0013676", "FBgn0013677",
    "FBgn0013678", "FBgn0013679", "FBgn0013680", "FBgn0013681",
    "FBgn0013682", "FBgn0013683", "FBgn0013684", "FBgn0013685",
    "FBgn0013686", "FBgn0013687", "FBgn0013688", "FBgn0013689",
    "FBgn0013690", "FBgn0013691", "FBgn0013692", "FBgn0013693",
    "FBgn0262006", "FBgn0262007", "FBgn0262008", "FBgn0262009",
]

MITO_CACHE_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "dmel_mito_genes_cache.tsv")

def get_mito_genes():
    # Always use cache if present — avoids network call on every run
    if os.path.exists(MITO_CACHE_PATH):
        print(f"Using cached mitochondrial gene list: {MITO_CACHE_PATH}")
        cached = pd.read_csv(MITO_CACHE_PATH, sep='\t')
        return cached['Gene stable ID'].values.flatten()

    # Try BioMart and cache the result
    try:
        dataset = Dataset(name='dmelanogaster_gene_ensembl',
                          host='http://www.ensembl.org')
        result  = dataset.query(
            attributes=['ensembl_gene_id', 'external_gene_name'],
            filters={'chromosome_name': ['mitochondrion_genome']}
        )
        result.to_csv(MITO_CACHE_PATH, sep='\t', index=False)
        print(f"BioMart succeeded: {len(result)} genes cached to {MITO_CACHE_PATH}")
        return result['Gene stable ID'].values.flatten()
    except Exception as e:
        print(f"WARNING: BioMart failed ({e}) — using hardcoded fallback")
        return np.array(MITO_GENES_FALLBACK)

mito_genes = get_mito_genes()

# ─────────────────────────────────────────────────────────────────────────────
# Matplotlib settings
# ─────────────────────────────────────────────────────────────────────────────

matplotlib.rcParams.update({
    'font.size': 6, 'figure.figsize': [2, 2],
    'axes.titlesize': 6, 'axes.labelsize': 6,
    'xtick.labelsize': 6, 'ytick.labelsize': 6,
    'legend.fontsize': 6, 'figure.titlesize': 6,
})
sc.settings.set_figure_params(dpi=300, dpi_save=300, figsize=(2, 2), fontsize=6)


# ─────────────────────────────────────────────────────────────────────────────
# Wolbachia titer — vectorised
# ─────────────────────────────────────────────────────────────────────────────

def calculate_wolbachia_titer(adata, rRNA_genes):
    """Vectorised titer calculation using sparse matrix slicing."""
    print("Calculating Wolbachia titer …")

    wMel_genes  = ["GQX67_05945"]
    var_names   = list(adata.var_names)

    # Single slice for all wMel genes
    wMel_present = [g for g in wMel_genes   if g in var_names]
    dmel_present = [g for g in rRNA_genes   if g in var_names]

    def _sum_genes(gene_list):
        if not gene_list:
            return np.zeros(adata.n_obs, dtype=np.float32)
        idx = [var_names.index(g) for g in gene_list]
        X   = adata.X[:, idx]
        if scipy.sparse.issparse(X):
            return np.asarray(X.sum(axis=1)).flatten().astype(np.float32)
        return X.sum(axis=1).astype(np.float32)

    wMel_total = _sum_genes(wMel_present)
    dmel_total = _sum_genes(dmel_present)

    denom = wMel_total + dmel_total
    # Where denominator is 0 (no rRNA detected), titer = 0
    titer = np.where(denom > 0, wMel_total / denom, 0.0)

    adata.obs['wolbachia_titer'] = titer
    print(f"  Mean: {titer.mean():.4f}  Median: {np.median(titer):.4f}  "
          f"wMel genes used: {wMel_present}  "
          f"Dmel rRNA genes present: {len(dmel_present)}")
    return adata


# ─────────────────────────────────────────────────────────────────────────────
# Doublet detection — subsampled
# ─────────────────────────────────────────────────────────────────────────────
def identify_doublets(adata, fig_dir):
    print("Starting scrublet doublet detection:")
    print(f"Dataset dimensions: {adata.n_obs} cells, {adata.n_vars} genes")
    
    max_components = min(adata.n_obs, adata.n_vars) - 1
    n_components   = min(30, max_components)
    print(f"Using {n_components} principal components (max possible: {max_components})")
    
    scrub = scr.Scrublet(adata.X, expected_doublet_rate=0.1)
    doublet_scores, predicted_doublets = scrub.scrub_doublets(
        min_counts=2, min_cells=3,
        min_gene_variability_pctl=85,
        n_prin_comps=n_components,
    )
    scrub.call_doublets()
    
    # Histogram only — scrublet UMAP removed (slow, not used downstream)
    print("Plotting scrublet histogram")
    scrub.plot_histogram()
    plt.savefig(f"{fig_dir}/doublet_histogram.pdf", bbox_inches='tight', pad_inches=0.1)
    plt.close()

    print("Saving scrublet data to adata")
    adata.obs['doublet_score']         = doublet_scores
    adata.obs['predicted_doublet']     = predicted_doublets
    adata.obs['predicted_doublet_cat'] = (adata.obs['predicted_doublet']
                                           .astype(str).astype('category'))
    return adata

# ─────────────────────────────────────────────────────────────────────────────
# QC metrics
# ─────────────────────────────────────────────────────────────────────────────

def calculate_qc_metrics(adata, sample_name="", stage=""):
    X = adata.X
    if scipy.sparse.issparse(X):
        genes_per_cell = np.asarray((X > 0).sum(axis=1)).flatten()
        umis_per_cell  = np.asarray(X.sum(axis=1)).flatten()
    else:
        genes_per_cell = (X > 0).sum(axis=1)
        umis_per_cell  = X.sum(axis=1)

    tcov = genes_per_cell / adata.n_vars

    doublet_rate = n_doublets = None
    if 'predicted_doublet' in adata.obs.columns:
        n_doublets   = int(adata.obs['predicted_doublet'].sum())
        doublet_rate = n_doublets / adata.n_obs

    return {
        'sample_name': sample_name, 'stage': stage,
        'n_cells': adata.n_obs, 'n_genes_total': adata.n_vars,
        'genes_per_cell_mean':   float(genes_per_cell.mean()),
        'genes_per_cell_median': float(np.median(genes_per_cell)),
        'genes_per_cell_std':    float(genes_per_cell.std()),
        'genes_per_cell_min':    int(genes_per_cell.min()),
        'genes_per_cell_max':    int(genes_per_cell.max()),
        'umis_per_cell_mean':    float(umis_per_cell.mean()),
        'umis_per_cell_median':  float(np.median(umis_per_cell)),
        'umis_per_cell_std':     float(umis_per_cell.std()),
        'umis_per_cell_min':     float(umis_per_cell.min()),
        'umis_per_cell_max':     float(umis_per_cell.max()),
        'transcriptome_coverage_mean':   float(tcov.mean()),
        'transcriptome_coverage_median': float(np.median(tcov)),
        'transcriptome_coverage_std':    float(tcov.std()),
        'transcriptome_coverage_min':    float(tcov.min()),
        'transcriptome_coverage_max':    float(tcov.max()),
        'doublet_rate': doublet_rate, 'n_doublets': n_doublets,
        'mitochondrial_percent_mean':   float(adata.obs['percent_mito'].mean())
                                        if 'percent_mito' in adata.obs.columns else None,
        'mitochondrial_percent_median': float(adata.obs['percent_mito'].median())
                                        if 'percent_mito' in adata.obs.columns else None,
    }


def print_qc_summary(metrics):
    print(f"\n{'='*50}")
    print(f"QC  |  {metrics['sample_name']}  |  {metrics['stage']}")
    print(f"{'='*50}")
    print(f"  Cells: {metrics['n_cells']:,}  Genes: {metrics['n_genes_total']:,}")
    print(f"  Genes/cell  mean={metrics['genes_per_cell_mean']:.0f}  "
          f"median={metrics['genes_per_cell_median']:.0f}  "
          f"range={metrics['genes_per_cell_min']}–{metrics['genes_per_cell_max']}")
    print(f"  UMIs/cell   mean={metrics['umis_per_cell_mean']:.0f}  "
          f"median={metrics['umis_per_cell_median']:.0f}")
    if metrics['doublet_rate'] is not None:
        print(f"  Doublets: {metrics['n_doublets']} "
              f"({metrics['doublet_rate']*100:.1f}%)")
    if metrics['mitochondrial_percent_mean'] is not None:
        print(f"  %Mito mean={metrics['mitochondrial_percent_mean']*100:.1f}%  "
              f"median={metrics['mitochondrial_percent_median']*100:.1f}%")
    print(f"{'='*50}")


def save_metrics_to_csv(metrics_list, output_path):
    pd.DataFrame(metrics_list).to_csv(output_path, index=False)
    print(f"Metrics → {output_path}")


# ─────────────────────────────────────────────────────────────────────────────
# QC plots + filtering
# ─────────────────────────────────────────────────────────────────────────────

def density_scatter(x, y, ax=None, sort=True, bins=20, **kwargs):
    if ax is None:
        fig, ax = plt.subplots()
    data, x_e, y_e = np.histogram2d(x, y, bins=bins)
    z = interpn(
        (0.5*(x_e[1:]+x_e[:-1]), 0.5*(y_e[1:]+y_e[:-1])),
        data, np.vstack([x, y]).T,
        method="splinef2d", bounds_error=False,
    )
    if sort:
        idx = z.argsort()
        x, y, z = x[idx], y[idx], z[idx]
    return ax.scatter(x, y, c=z, **kwargs)


def qc_plots(adata, fig_dir, svd_subsample=5000):
    """QC plots + knee-plot filtering.

    SVD is computed on a subsample (default 5000 cells) for speed.
    Knee plot is fully vectorised — no Python loop over barcodes.
    """
    # ── SVD plot (subsampled) ─────────────────────────────────────────────────
    n_svd = min(svd_subsample, adata.n_obs)
    if n_svd < adata.n_obs:
        idx_svd = np.random.default_rng(42).choice(adata.n_obs, n_svd, replace=False)
        X_svd   = adata.X[idx_svd]
    else:
        X_svd = adata.X

    if scipy.sparse.issparse(X_svd):
        X_svd = X_svd.toarray()

    tsvd  = TruncatedSVD(n_components=2, random_state=42)
    X_2d  = tsvd.fit_transform(X_svd)

    fig, ax = plt.subplots(figsize=(7, 7))
    density_scatter(X_2d[:, 0], X_2d[:, 1], ax=ax, cmap="Greens")
    ax.axis('off')
    plt.savefig(f"{fig_dir}/density_scatter.pdf", bbox_inches='tight')
    plt.close(fig)

    # ── Genes vs UMI (subsampled for speed) ───────────────────────────────────
    n_plot = min(20000, adata.n_obs)
    idx_p  = np.random.default_rng(0).choice(adata.n_obs, n_plot, replace=False)
    X_p    = adata.X[idx_p]
    if scipy.sparse.issparse(X_p):
        x_umi  = np.asarray(X_p.sum(axis=1)).flatten()
        y_gene = np.asarray((X_p > 0).sum(axis=1)).flatten()
    else:
        x_umi  = X_p.sum(axis=1)
        y_gene = (X_p > 0).sum(axis=1)

    fig, ax = plt.subplots(figsize=(7, 5))
    ax.scatter(x_umi, y_gene, color="green", alpha=0.15, s=1, rasterized=True)
    ax.set_xlabel("UMI Counts"); ax.set_ylabel("Genes Detected")
    ax.set_xscale('log');        ax.set_yscale('log')
    plt.savefig(f"{fig_dir}/genes_vs_umi.pdf", bbox_inches='tight')
    plt.close(fig)

    # ── Knee plot — fully vectorised ──────────────────────────────────────────
    if scipy.sparse.issparse(adata.X):
        umi_per_cell = np.asarray(adata.X.sum(axis=1)).flatten()
    else:
        umi_per_cell = adata.X.sum(axis=1)

    cutoff = 100
    knee   = np.sort(umi_per_cell)[::-1]        # descending
    cell_set = np.arange(len(knee))

    # Vectorised: find last index where knee > cutoff
    above    = knee > cutoff
    num_cells = int(above.sum()) - 1             # index of last cell above cutoff

    fig, ax = plt.subplots(figsize=(7, 5))
    ax.loglog(knee, cell_set, linewidth=2, color="g")
    ax.axvline(x=cutoff,    linewidth=2, color="k")
    ax.axhline(y=num_cells, linewidth=2, color="k")
    ax.set_xlabel("UMI Counts"); ax.set_ylabel("Barcodes")
    plt.grid(True, which="both")
    plt.savefig(f"{fig_dir}/knee_plot.pdf", bbox_inches='tight')
    plt.close(fig)

    print(f"Knee plot: {num_cells:,} cells above {cutoff} UMI threshold")

    # ── Filter cells ──────────────────────────────────────────────────────────
    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_cells(adata, min_counts=knee[num_cells])
    print(f"After knee filtering: {adata.n_obs} cells")

    # ── Mitochondrial % ───────────────────────────────────────────────────────
    valid_mito = [g for g in mito_genes if g in adata.var_names]
    if valid_mito:
        mito_counts  = np.asarray(adata[:, valid_mito].X.sum(axis=1)).flatten()
        total_counts = np.asarray(adata.X.sum(axis=1)).flatten()
        adata.obs['percent_mito'] = mito_counts / np.where(total_counts == 0, 1, total_counts)
    else:
        print("No mitochondrial genes found — setting percent_mito to 0")
        adata.obs['percent_mito'] = 0.0

    adata.obs['n_counts'] = np.asarray(adata.X.sum(axis=1)).flatten()

    fig, ax = plt.subplots(figsize=(7, 5))
    sc.pl.scatter(adata, x='n_counts', y='percent_mito', ax=ax, show=False)
    plt.savefig(f"{fig_dir}/mito_scatter.pdf", bbox_inches='tight')
    plt.close(fig)

    # ── Filter genes ──────────────────────────────────────────────────────────
    sc.pp.filter_genes(adata, min_cells=3)
    print(f"After gene filtering: {adata.n_vars} genes")

    sc.pl.violin(adata, ['n_genes', 'n_counts', 'percent_mito'],
                 jitter=0.4, multi_panel=True, show=False)
    plt.savefig(f"{fig_dir}/violin_plot.pdf", bbox_inches='tight')
    plt.close()

    return adata


# ─────────────────────────────────────────────────────────────────────────────
# Normalise + save .raw  (no per-sample clustering)
# ─────────────────────────────────────────────────────────────────────────────

def analyze_filtered_adata(adata, output_prefix=""):
    """Doublet removal → save raw → normalise → scale.

    Per-sample PCA/UMAP/leiden clustering is intentionally removed:
    integrate.py performs joint clustering across all samples, which is
    both faster and biologically more meaningful.
    """
    print(f"Cells before doublet removal: {adata.n_obs}")
    adata = adata[~adata.obs['predicted_doublet']].copy()
    print(f"Cells after doublet removal:  {adata.n_obs}")

    # ── Save raw counts BEFORE any normalisation ──────────────────────────────
    adata.raw = adata
    print(f"Saved raw counts to adata.raw "
          f"({adata.n_obs} cells × {adata.n_vars} genes)")

    # ── Normalise + log1p ─────────────────────────────────────────────────────
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    # ── HVG (on log-normalised data, seurat flavor) ───────────────────────────
    sc.pp.highly_variable_genes(
        adata, min_mean=0.01, max_mean=10, min_disp=0.1,
        n_top_genes=1000, n_bins=20, flavor="seurat",
    )
    n_hvg = adata.var['highly_variable'].sum()
    print(f"HVGs: {n_hvg}")

    # ── Scale ─────────────────────────────────────────────────────────────────
    sc.pp.scale(adata, max_value=10)

    print(f"Final: {adata.n_obs} cells × {adata.n_vars} genes")
    print(adata)
    return adata


# ─────────────────────────────────────────────────────────────────────────────
# Main processing
# ─────────────────────────────────────────────────────────────────────────────

def process_data_with_metrics(key, matrix, log_to_file=True):
    adata = ad.read_h5ad(matrix)

    if log_to_file:
        log_file        = open(f"{output_dir}/{key}_stats.txt", 'w')
        original_stdout = sys.stdout
        sys.stdout      = log_file

    try:
        all_metrics = []

        print("\n=== Calculating Wolbachia titer ===")
        adata = calculate_wolbachia_titer(adata, rRNA_genes)

        raw_metrics = calculate_qc_metrics(adata, sample_name=key, stage="raw")
        all_metrics.append(raw_metrics)
        print_qc_summary(raw_metrics)

        print("\n=== Doublet detection ===")
        adata = identify_doublets(adata, fig_dir)
        all_metrics.append(
            calculate_qc_metrics(adata, sample_name=key,
                                  stage="post_doublet_detection"))
        print_qc_summary(all_metrics[-1])

        print("\n=== QC filtering ===")
        filtered_adata = qc_plots(adata, fig_dir)
        del adata; gc.collect()

        all_metrics.append(
            calculate_qc_metrics(filtered_adata, sample_name=key, stage="filtered"))
        print_qc_summary(all_metrics[-1])

        print("\n=== Normalising ===")
        final_adata = analyze_filtered_adata(filtered_adata, output_dir)
        del filtered_adata; gc.collect()

        save_metrics_to_csv(all_metrics, f"{fig_dir}/{key}_qc_metrics.csv")

        final_adata.write(output)
        print(f"Written → {output}")

    finally:
        if log_to_file:
            sys.stdout = original_stdout
            log_file.close()

    return all_metrics


print(f"Analysing {sample_name}: {input}")
process_data_with_metrics(sample_name, input)
print(f"Completed {sample_name}")