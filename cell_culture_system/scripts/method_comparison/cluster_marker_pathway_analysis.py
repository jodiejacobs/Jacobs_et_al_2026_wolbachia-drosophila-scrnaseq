'''
Comprehensive cluster analysis: transcriptional activity, marker genes, and pathway enrichment
Uses FlyEnrichr API for automated pathway analysis with FlyBase annotations

Key features:
  - DE run on .raw (log-normalised counts), NOT scaled .X
  - Per-cluster adaptive marker thresholds (targets 30-150 markers per cluster)
  - Background gene set passed to FlyEnrichr (all detected genes in dataset)
  - Significance filter (adj_p < 0.1) applied before combined_score sorting
  - mt/ribo/cell_cycle/bacterial/TE genes excluded from enrichment input
  - Dot plot + network plot visualisation for enrichment results
  - Sparse-cluster rescue for clusters with 0 adaptive markers (broad/ubiquitous clusters):
      A. Pairwise DE vs nearest-neighbour clusters in Harmony PCA space (primary)
      D. Tau-like specificity score across all clusters (fallback if A yields < 5 genes)
    Rescue markers are flagged with rescue_strategy in the output CSVs.
  - Outputs: *_markers_cluster_fbgn_pval.csv and *_background_genes.csv
'''

import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import argparse
import os
import requests
import time
import gzip
from scipy.stats import kruskal
from matplotlib.patches import Patch
import scipy.sparse
from scipy.sparse import issparse

try:
    import networkx as nx
    NETWORKX_AVAILABLE = True
except ImportError:
    NETWORKX_AVAILABLE = False


# ─────────────────────────────────────────────────────────────────────────────
# Gene symbol mapping
# ─────────────────────────────────────────────────────────────────────────────

def load_fbgn_to_symbol_mapping(mapping_file):
    """Load FBgn -> gene symbol from transcripts_to_genes.txt"""
    print("\n" + "="*60)
    print("LOADING GENE SYMBOL MAPPING")
    print("="*60)
    print(f"  File: {mapping_file}")

    fbgn_to_symbol = {}
    try:
        opener, mode = (gzip.open, 'rt') if mapping_file.endswith('.gz') else (open, 'r')
        with opener(mapping_file, mode) as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) < 3:
                    continue
                fbgn_id, gene_symbol = parts[1], parts[2]
                if fbgn_id.startswith('FBgn'):
                    fbgn_to_symbol[fbgn_id] = gene_symbol
        print(f"  Loaded {len(fbgn_to_symbol):,} mappings")
        for fbgn, sym in list(fbgn_to_symbol.items())[:5]:
            print(f"    {fbgn} -> {sym}")
        return fbgn_to_symbol
    except FileNotFoundError:
        print(f"  ERROR: File not found: {mapping_file}")
        return None
    except Exception as e:
        print(f"  ERROR: {e}")
        import traceback; traceback.print_exc()
        return None


def symbols_from_fbgn(fbgn_list, fbgn_to_symbol):
    """Convert a list of FBgn IDs to symbols, return (symbols, n_unmapped)."""
    symbols, unmapped = [], []
    for g in fbgn_list:
        sym = fbgn_to_symbol.get(g)
        if sym:
            symbols.append(sym)
        else:
            unmapped.append(g)
    return symbols, len(unmapped)


def _is_te(var_names_series):
    """
    Return boolean mask for transposable element IDs.
    Covers FBti* IDs and *_transposable_element names.
    """
    return (
        var_names_series.str.startswith('FBti') |
        var_names_series.str.contains('transposable_element', regex=False)
    )


# ─────────────────────────────────────────────────────────────────────────────
# Transcriptional activity
# ─────────────────────────────────────────────────────────────────────────────

def plot_transcriptional_activity(adata, output_dir, sample_name):
    print("\n" + "="*60)
    print("TRANSCRIPTIONAL ACTIVITY BY CLUSTER")
    print("="*60)

    required = ['n_counts', 'n_genes', 'leiden']
    missing = [c for c in required if c not in adata.obs.columns]
    if missing:
        print(f"  ERROR: Missing columns: {missing}")
        return None

    clusters = sorted(adata.obs['leiden'].unique(),
                      key=lambda x: int(x) if str(x).isdigit() else x)

    summary = adata.obs.groupby('leiden')[['n_counts', 'n_genes']].agg(
        ['mean', 'median', 'std', 'min', 'max'])
    print(summary)
    summary.to_csv(os.path.join(output_dir,
                                f'{sample_name}_transcriptional_activity_summary.csv'))

    groups_counts = [adata.obs[adata.obs['leiden'] == c]['n_counts'].values for c in clusters]
    groups_genes  = [adata.obs[adata.obs['leiden'] == c]['n_genes'].values  for c in clusters]
    h_counts, p_counts = kruskal(*groups_counts)
    h_genes,  p_genes  = kruskal(*groups_genes)

    n, k = len(adata.obs), len(clusters)
    eta_counts = (h_counts - k + 1) / (n - k)
    eta_genes  = (h_genes  - k + 1) / (n - k)

    def eta_label(e):
        if e < 0.01: return "negligible"
        if e < 0.06: return "small"
        if e < 0.14: return "medium"
        return "large"

    print(f"\nn_counts: H={h_counts:.2f}  p={p_counts:.2e}  "
          f"η²={eta_counts:.4f} ({eta_label(eta_counts)})")
    print(f"n_genes:  H={h_genes:.2f}  p={p_genes:.2e}  "
          f"η²={eta_genes:.4f} ({eta_label(eta_genes)})")

    if 'leiden_colors' in adata.uns:
        palette = dict(zip(clusters, adata.uns['leiden_colors']))
    else:
        cmap = plt.colormaps.get_cmap('tab20')
        palette = {c: cmap(i % 20) for i, c in enumerate(clusters)}

    fig, axes = plt.subplots(1, 2, figsize=(16, 6))
    for ax, metric, h, eta, label in [
        (axes[0], 'n_counts', h_counts, eta_counts, 'Total UMI Counts'),
        (axes[1], 'n_genes',  h_genes,  eta_genes,  'Genes Detected'),
    ]:
        box_data = [adata.obs[adata.obs['leiden'] == c][metric].values for c in clusters]
        bp = ax.boxplot(box_data, labels=clusters, patch_artist=True,
                        showfliers=False, widths=0.6)
        for patch, c in zip(bp['boxes'], clusters):
            patch.set_facecolor(palette[c]); patch.set_alpha(0.7)
        ax.set_xlabel('Leiden Cluster', fontsize=12)
        ax.set_ylabel(label, fontsize=12)
        ax.set_title(f'{label} by Cluster\nη²={eta:.4f}  H={h:.2f}', fontsize=13)
        ax.grid(axis='y', alpha=0.3)

    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, f'{sample_name}_transcriptional_activity.pdf'),
                dpi=150, bbox_inches='tight')
    plt.close()

    if 'X_umap' in adata.obsm:
        sc.pl.umap(adata, color=['leiden', 'n_counts', 'n_genes'],
                   save=f'_{sample_name}_transcriptional_activity.pdf',
                   cmap='viridis', ncols=3)

    return dict(h_counts=h_counts, p_counts=p_counts,
                h_genes=h_genes,   p_genes=p_genes,
                eta_counts=eta_counts, eta_genes=eta_genes)


# ─────────────────────────────────────────────────────────────────────────────
# Sparse-cluster rescue: Strategy A (pairwise) + Strategy D (specificity)
# ─────────────────────────────────────────────────────────────────────────────

def _nearest_neighbour_clusters(adata_de, target_cl, n_neighbours=3):
    """
    Return the n_neighbours clusters closest to target_cl in mean Harmony PCA space.
    Falls back to raw PCA, then X_umap, then random selection if none available.
    """
    # Choose embedding: prefer Harmony PCA, then PCA, then UMAP
    if 'X_pca_harmony' in adata_de.obsm:
        coords_key = 'X_pca_harmony'
    elif 'X_pca' in adata_de.obsm:
        coords_key = 'X_pca'
    elif 'X_umap' in adata_de.obsm:
        coords_key = 'X_umap'
    else:
        coords_key = None

    all_clusters = [c for c in adata_de.obs['leiden'].unique() if c != target_cl]

    if coords_key is None or len(all_clusters) == 0:
        return all_clusters[:n_neighbours]

    coords = adata_de.obsm[coords_key]
    target_mask = adata_de.obs['leiden'] == target_cl
    target_centroid = coords[target_mask].mean(axis=0)

    dists = {}
    for cl in all_clusters:
        mask = adata_de.obs['leiden'] == cl
        centroid = coords[mask].mean(axis=0)
        dists[cl] = np.linalg.norm(target_centroid - centroid)

    neighbours = sorted(dists, key=dists.get)[:n_neighbours]
    print(f"    Nearest neighbours (by {coords_key}): {neighbours}")
    return neighbours


def _pairwise_de(adata_de, target_cl, compare_cls, method='wilcoxon',
                 log2fc_min=0.25, pct_in_min=0.05, adj_p_max=0.1):
    """
    Strategy A: run DE of target_cl vs each neighbour cluster individually,
    then take the union of genes that pass filters in >= 1 comparison.
    Returns a DataFrame with the same columns as the main marker_df, plus
    'rescue_strategy' = 'pairwise_DE' and 'reference_clusters'.
    """
    passing = {}  # gene -> best row across comparisons

    for ref_cl in compare_cls:
        subset_mask = adata_de.obs['leiden'].isin([target_cl, ref_cl])
        adata_sub = adata_de[subset_mask].copy()

        # Need >= 2 cells in each group
        n_target = (adata_sub.obs['leiden'] == target_cl).sum()
        n_ref    = (adata_sub.obs['leiden'] == ref_cl).sum()
        if n_target < 2 or n_ref < 2:
            print(f"    Skipping vs cluster {ref_cl}: too few cells "
                  f"(target={n_target}, ref={n_ref})")
            continue

        try:
            sc.tl.rank_genes_groups(
                adata_sub, 'leiden', groups=[target_cl],
                reference=ref_cl, method=method,
                key_added='rgg_pairwise',
                tie_correct=True, rankby_abs=False, pts=True,
            )
        except Exception as e:
            print(f"    WARNING: pairwise DE vs {ref_cl} failed: {e}")
            continue

        res = adata_sub.uns['rgg_pairwise']
        for i in range(len(res['names'][target_cl])):
            gene = res['names'][target_cl][i]
            lfc  = res['logfoldchanges'][target_cl][i]
            padj = res['pvals_adj'][target_cl][i]
            pval = res['pvals'][target_cl][i]
            scr  = res['scores'][target_cl][i]
            pct_in   = res['pts'][target_cl][i]   if 'pts'      in res else np.nan
            pct_rest = res['pts_rest'][target_cl][i] if 'pts_rest' in res else np.nan

            if pd.isna(lfc) or np.isinf(lfc):
                continue
            if lfc < log2fc_min or pct_in < pct_in_min or padj >= adj_p_max:
                continue

            # Keep best (lowest padj) row per gene across comparisons
            if gene not in passing or padj < passing[gene]['pval_adj']:
                passing[gene] = dict(
                    cluster          = target_cl,
                    gene             = gene,
                    log2fc           = lfc,
                    pval             = pval,
                    pval_adj         = padj,
                    score            = scr,
                    pct_in           = pct_in,
                    pct_rest         = pct_rest,
                    rescue_strategy  = 'pairwise_DE',
                    reference_clusters = str(compare_cls),
                )

    if not passing:
        return pd.DataFrame()

    df = pd.DataFrame(list(passing.values()))
    df['pct_ratio'] = df['pct_in'] / (df['pct_rest'] + 1e-9)
    return df.sort_values('score', ascending=False)


def _specificity_score(adata_de, target_cl, top_n=50, pct_in_min=0.1):
    """
    Strategy D: tau-like specificity score.

    For each gene, compute mean expression per cluster (log-normalised space).
    Tau = (1 - x_i / x_max) summed across clusters, normalised.
    High tau = expressed specifically in few clusters.

    We then filter to genes where target_cl has the *highest* mean and
    tau >= 0.5, returning the top_n by tau * mean_in_target.

    Returns a DataFrame with the same core columns as marker_df, plus
    'rescue_strategy' = 'specificity_tau' and 'tau' column.
    """
    # Dense mean expression per cluster per gene
    X = adata_de.X
    if scipy.sparse.issparse(X):
        X = X.toarray()

    genes   = adata_de.var_names.tolist()
    clusters = adata_de.obs['leiden'].unique().tolist()

    # Build (n_clusters x n_genes) mean matrix
    means = np.zeros((len(clusters), len(genes)), dtype=np.float32)
    for i, cl in enumerate(clusters):
        mask = (adata_de.obs['leiden'] == cl).values
        means[i] = X[mask].mean(axis=0)

    # Tau: ranges 0 (ubiquitous) -> 1 (perfectly specific)
    x_max = means.max(axis=0)
    x_max[x_max == 0] = 1  # avoid div by zero for unexpressed genes
    n = len(clusters)
    tau = ((1 - means / x_max).sum(axis=0)) / (n - 1) if n > 1 else np.zeros(len(genes))

    target_idx = clusters.index(target_cl)
    mean_target = means[target_idx]

    # pct expressed in target cluster
    target_mask = (adata_de.obs['leiden'] == target_cl).values
    X_target = X[target_mask]
    pct_in = (X_target > 0).mean(axis=0)

    # Filter: target must have highest mean AND tau >= 0.5 AND pct_in >= threshold
    argmax_cluster = means.argmax(axis=0)
    keep = (
        (argmax_cluster == target_idx) &
        (tau >= 0.5) &
        (pct_in >= pct_in_min)
    )

    if keep.sum() == 0:
        # Relax: just highest mean in target + tau >= 0.3
        keep = (argmax_cluster == target_idx) & (tau >= 0.3) & (pct_in >= pct_in_min)

    if keep.sum() == 0:
        return pd.DataFrame()

    kept_genes   = [genes[i]       for i in range(len(genes)) if keep[i]]
    kept_tau     = [tau[i]         for i in range(len(genes)) if keep[i]]
    kept_mean    = [mean_target[i] for i in range(len(genes)) if keep[i]]
    kept_pct     = [pct_in[i]      for i in range(len(genes)) if keep[i]]

    # Score = tau * mean (ranks genes that are both specific and well-expressed)
    spec_score = np.array(kept_tau) * np.array(kept_mean)
    order = np.argsort(spec_score)[::-1][:top_n]

    rows = []
    for idx in order:
        rows.append(dict(
            cluster           = target_cl,
            gene              = kept_genes[idx],
            log2fc            = np.log2(kept_mean[idx] + 1),  # pseudo log2FC vs 0
            pval              = np.nan,
            pval_adj          = np.nan,
            score             = spec_score[idx],
            pct_in            = kept_pct[idx],
            pct_rest          = np.nan,
            pct_ratio         = np.nan,
            tau               = kept_tau[idx],
            rescue_strategy   = 'specificity_tau',
            reference_clusters= 'all',
        ))

    return pd.DataFrame(rows)


def rescue_sparse_clusters(adata_de, marker_df_filtered, all_marker_df,
                           output_dir, sample_name, method='wilcoxon',
                           n_neighbours=3):
    """
    For any cluster with 0 markers after adaptive filtering, attempt rescue:

    Strategy A — pairwise DE vs n_neighbours nearest clusters in PCA space.
                 Filters: log2FC >= 0.5, pct_in >= 0.05, adj_p < 0.1.
                 Takes union of genes significant in >= 1 pairwise comparison.

    Strategy D — Tau-like specificity score across all clusters.
                 Selects genes where this cluster has the highest mean expression
                 and tau >= 0.5 (relaxed to 0.3 if needed), pct_in >= 0.1.

    A is tried first. D is used as fallback only if A yields < 5 genes.
    Both can contribute: if A gives >= 5 genes it is used alone; otherwise
    D supplements A (or replaces it entirely).

    Rescued markers are tagged with 'rescue_strategy' and 'reference_clusters'
    columns. Non-rescued markers get rescue_strategy = 'adaptive'.

    Returns the merged marker DataFrame and a rescue summary dict.
    """
    all_clusters = sorted(
        adata_de.obs['leiden'].unique(),
        key=lambda x: int(x) if str(x).isdigit() else x,
    )
    sparse_clusters = [
        cl for cl in all_clusters
        if cl not in marker_df_filtered['cluster'].values
        or len(marker_df_filtered[marker_df_filtered['cluster'] == cl]) == 0
    ]

    if not sparse_clusters:
        print("  No sparse clusters — rescue not needed.")
        marker_df_filtered['rescue_strategy']   = 'adaptive'
        marker_df_filtered['reference_clusters'] = 'all_others'
        marker_df_filtered['tau']                = np.nan
        return marker_df_filtered, {}

    print(f"\n  Sparse clusters (0 adaptive markers): {sparse_clusters}")
    rescue_summary = {}
    rescue_parts   = []

    for cl in sparse_clusters:
        print(f"\n  {'─'*50}")
        print(f"  Rescuing cluster {cl} ...")

        # ── Strategy A: pairwise DE vs nearest neighbours ─────────────────────
        neighbours = _nearest_neighbour_clusters(adata_de, cl, n_neighbours=n_neighbours)
        print(f"  [A] Pairwise DE vs {neighbours} ...")
        a_df = _pairwise_de(
            adata_de, cl, neighbours, method=method,
            log2fc_min=0.5, pct_in_min=0.05, adj_p_max=0.1,
        )
        print(f"  [A] → {len(a_df)} genes")

        if len(a_df) >= 5:
            rescue_parts.append(a_df)
            rescue_summary[cl] = {'strategy': 'A_pairwise', 'n_genes': len(a_df),
                                  'neighbours': neighbours}
            print(f"  ✓ Strategy A succeeded ({len(a_df)} genes)")
            continue

        # ── Strategy D: tau-like specificity score ────────────────────────────
        print(f"  [D] Tau specificity score (A gave only {len(a_df)} genes) ...")
        d_df = _specificity_score(adata_de, cl, top_n=50, pct_in_min=0.1)
        print(f"  [D] → {len(d_df)} genes")

        if len(a_df) == 0 and len(d_df) == 0:
            print(f"  ✗ Both strategies failed for cluster {cl} — "
                  f"this cluster may be genuinely transcriptionally indistinct.")
            rescue_summary[cl] = {'strategy': 'none', 'n_genes': 0, 'neighbours': neighbours}
            continue

        # Combine A and D if A gave > 0 genes, otherwise use D alone
        if len(a_df) > 0 and len(d_df) > 0:
            combined = pd.concat([a_df, d_df], ignore_index=True)
            combined = combined.drop_duplicates(subset='gene', keep='first')
            strategy_label = 'A_pairwise+D_tau'
        elif len(d_df) > 0:
            combined = d_df
            strategy_label = 'D_tau'
        else:
            combined = a_df
            strategy_label = 'A_pairwise'

        rescue_parts.append(combined)
        rescue_summary[cl] = {'strategy': strategy_label,
                               'n_genes': len(combined),
                               'neighbours': neighbours}
        print(f"  ✓ Rescue: {strategy_label} ({len(combined)} genes)")

    # Tag original adaptive markers
    marker_df_filtered = marker_df_filtered.copy()
    marker_df_filtered['rescue_strategy']    = 'adaptive'
    marker_df_filtered['reference_clusters'] = 'all_others'
    if 'tau' not in marker_df_filtered.columns:
        marker_df_filtered['tau'] = np.nan

    if rescue_parts:
        rescue_df = pd.concat(rescue_parts, ignore_index=True)
        # Ensure all columns align
        for col in ['tau', 'reference_clusters']:
            if col not in rescue_df.columns:
                rescue_df[col] = np.nan

        merged = pd.concat([marker_df_filtered, rescue_df], ignore_index=True)
    else:
        merged = marker_df_filtered

    # Print rescue summary
    print(f"\n  {'─'*50}")
    print(f"  RESCUE SUMMARY")
    for cl, info in rescue_summary.items():
        print(f"    Cluster {cl}: strategy={info['strategy']}  "
              f"n_genes={info['n_genes']}  neighbours={info.get('neighbours','')}")

    # Save rescue-specific outputs
    if rescue_parts:
        rescue_df_all = pd.concat(rescue_parts, ignore_index=True)
        rescue_df_all.to_csv(
            os.path.join(output_dir, f'{sample_name}_markers_rescued.csv'), index=False)
        print(f"\n  Saved: {sample_name}_markers_rescued.csv")

    return merged, rescue_summary


# ─────────────────────────────────────────────────────────────────────────────
# Marker genes — per-cluster adaptive thresholds
# ─────────────────────────────────────────────────────────────────────────────

def find_marker_genes(adata, output_dir, sample_name, method='wilcoxon',
                      log2fc_min=0.25, pct_min=0.05, pct_ratio_min=1.0):
    """
    Find cluster marker genes using per-cluster adaptive thresholds.

    DE is run on adata.raw (log-normalised counts). Scaled .X is NOT used.

    Adaptive thresholds:
      Progressively stricter filters are applied per cluster, targeting
      30-150 markers. This prevents large/distinct clusters from flooding
      the enrichment analysis with thousands of generic markers, while
      preserving sensitivity for small/similar clusters.

    The initial log2fc_min/pct_min/pct_ratio_min args set the floor
    (loosest threshold). Adaptive tightening starts from there.

    Bacterial genes (GQX67_) and transposable elements (FBti*, *_transposable_element)
    are excluded from DE so they don't dilute the enrichment foreground.

    Outputs
    -------
    *_markers_all.csv              : all DE results pre-filter
    *_markers_filtered.csv         : adaptive-filtered markers (full columns)
    *_markers_top50.csv            : top 50 per cluster
    *_marker_thresholds.csv        : per-cluster threshold log
    *_markers_cluster_fbgn_pval.csv: cluster / flybase_id / p_value / adj_p_value
    """
    print("\n" + "="*60)
    print("DIFFERENTIAL GENE EXPRESSION")
    print("="*60)
    print(f"  Method: {method}")
    print(f"  Adaptive thresholds: targeting 30-150 markers per cluster")

    # ── Use .raw for DE ───────────────────────────────────────────────────────
    if adata.raw is not None:
        print("  Using adata.raw for DE (log-normalised, pre-scale counts)")
        adata_de = adata.raw.to_adata()
        adata_de.X = scipy.sparse.csr_matrix(adata_de.X)
        adata_de.obs['leiden'] = adata.obs['leiden'].values
    else:
        print("  WARNING: adata.raw is None — using adata.X")
        adata_de = adata.copy()
        adata_de.X = scipy.sparse.csr_matrix(adata_de.X)

    # Filter out Wolbachia transcripts
    n_before_bact = adata_de.n_vars
    adata_de = adata_de[:, ~adata_de.var_names.str.startswith('GQX67_')].copy()
    print(f"  Removed {n_before_bact - adata_de.n_vars:,} bacterial transcripts "
          f"({adata_de.n_vars:,} remaining)")

    # Filter out transposable elements (FBti* and *_transposable_element)
    # These cannot be mapped to FlyEnrichr GO terms and dilute marker gene lists
    n_before_te = adata_de.n_vars
    te_mask = _is_te(pd.Series(adata_de.var_names))
    adata_de = adata_de[:, ~te_mask.values].copy()
    print(f"  Removed {n_before_te - adata_de.n_vars:,} transposable elements "
          f"(FBti* / *_transposable_element) ({adata_de.n_vars:,} remaining)")

    # Filter to genes expressed in >= 3 cells
    n_cells_per_gene = np.array((adata_de.X > 0).sum(axis=0)).flatten()
    n_before = adata_de.n_vars
    adata_de = adata_de[:, n_cells_per_gene >= 3].copy()
    print(f"  Genes after expression filter: {adata_de.n_vars:,} / {n_before:,}")

    sc.settings.n_jobs = -1
    sc.tl.rank_genes_groups(
        adata_de, 'leiden', method=method,
        key_added='rank_genes_groups',
        tie_correct=True,
        rankby_abs=False,
        pts=True,
    )

    adata.uns['rank_genes_groups'] = adata_de.uns['rank_genes_groups']
    result = adata.uns['rank_genes_groups']
    groups = result['names'].dtype.names

    # Collect all DE results
    rows = []
    for grp in groups:
        for i in range(len(result['names'][grp])):
            lfc  = result['logfoldchanges'][grp][i]
            padj = result['pvals_adj'][grp][i]
            if pd.isna(lfc) or np.isinf(lfc):
                continue
            pct_in   = result['pts'][grp][i]      if 'pts'      in result else np.nan
            pct_rest = result['pts_rest'][grp][i]  if 'pts_rest' in result else np.nan
            rows.append(dict(
                cluster  = grp,
                gene     = result['names'][grp][i],
                log2fc   = lfc,
                pval     = result['pvals'][grp][i],
                pval_adj = padj,
                score    = result['scores'][grp][i],
                pct_in   = pct_in,
                pct_rest = pct_rest,
            ))

    marker_df = pd.DataFrame(rows)
    marker_df['pct_ratio'] = marker_df['pct_in'] / (marker_df['pct_rest'] + 1e-9)
    print(f"  Raw DE results: {len(marker_df):,} gene×cluster entries")

    # Save all results before filtering
    marker_df.to_csv(os.path.join(output_dir, f'{sample_name}_markers_all.csv'), index=False)

    # ── Per-cluster adaptive thresholds ──────────────────────────────────────
    # Progressively tighten filters per cluster, targeting 30-150 markers.
    # If a cluster cannot reach 30 markers at any threshold, keep the loosest.
    TARGET_MIN = 30
    TARGET_MAX = 150
    THRESHOLDS = [
        (0.25, 0.03, 1.0),   # loosest floor
        (0.50, 0.05, 1.2),
        (0.75, 0.08, 1.3),
        (1.00, 0.10, 1.5),
        (1.25, 0.12, 1.75),
        (1.50, 0.15, 2.0),
        (2.00, 0.20, 2.5),
    ]

    adaptive_parts = []
    threshold_log  = []

    for cl in sorted(marker_df['cluster'].unique(),
                     key=lambda x: int(x) if str(x).isdigit() else x):
        cm = marker_df[marker_df['cluster'] == cl].copy()
        best_filt   = None
        best_thresh = THRESHOLDS[0]

        for lfc, pct_min_t, pct_rat in THRESHOLDS:
            filt = cm[
                (cm['log2fc']    >= lfc) &
                (cm['pval_adj']   < 0.1) &
                (cm['pct_in']    >= pct_min_t) &
                (cm['pct_ratio'] >= pct_rat)
            ]
            # Keep track of loosest result
            if best_filt is None:
                best_filt   = filt
                best_thresh = (lfc, pct_min_t, pct_rat)
            # Accept first threshold that hits the target range
            if TARGET_MIN <= len(filt) <= TARGET_MAX:
                best_filt   = filt
                best_thresh = (lfc, pct_min_t, pct_rat)
                break
            # Update best if we're above the minimum
            if len(filt) >= TARGET_MIN:
                best_filt   = filt
                best_thresh = (lfc, pct_min_t, pct_rat)

        # Hard cap at TARGET_MAX by score to prevent over-representation
        if len(best_filt) > TARGET_MAX:
            best_filt = best_filt.nlargest(TARGET_MAX, 'score')

        adaptive_parts.append(best_filt)
        threshold_log.append({
            'cluster':       cl,
            'log2fc_thresh': best_thresh[0],
            'pct_thresh':    best_thresh[1],
            'ratio_thresh':  best_thresh[2],
            'n_markers':     len(best_filt),
        })

    marker_df_filtered = pd.concat(adaptive_parts, ignore_index=True)
    thresh_df = pd.DataFrame(threshold_log)

    print(f"\n  Adaptive thresholds per cluster:")
    print(f"  {'Cluster':>10} {'log2fc':>8} {'pct':>6} {'ratio':>7} {'n_markers':>10}")
    for _, row in thresh_df.iterrows():
        flag = ' ✓' if TARGET_MIN <= row.n_markers <= TARGET_MAX else ' ⚠'
        print(f"  {int(row.cluster):>10} {row.log2fc_thresh:>8.2f} "
              f"{row.pct_thresh:>6.2f} {row.ratio_thresh:>7.2f} "
              f"{int(row.n_markers):>10}{flag}")

    print(f"\n  Total adaptive markers: {len(marker_df_filtered):,}")

    # ── Sparse-cluster rescue (A: pairwise DE, D: tau specificity) ───────────
    print("\n" + "="*60)
    print("SPARSE CLUSTER RESCUE")
    print("="*60)
    marker_df_filtered, rescue_summary = rescue_sparse_clusters(
        adata_de, marker_df_filtered, marker_df,
        output_dir, sample_name, method=method,
    )

    # Print top 5 per cluster (flag rescued clusters)
    print("\n  Top 5 markers per cluster (by score):")
    all_clusters_sorted = sorted(
        marker_df_filtered['cluster'].unique(),
        key=lambda x: int(x) if str(x).isdigit() else x,
    )
    for cl in all_clusters_sorted:
        sub = marker_df_filtered[marker_df_filtered['cluster'] == cl].nlargest(5, 'score')
        strat = sub['rescue_strategy'].iloc[0] if len(sub) > 0 else 'adaptive'
        tag = f' [rescued: {strat}]' if strat != 'adaptive' else ''
        print(f"\n  Cluster {cl}{tag}:")
        if len(sub) == 0:
            print("    (no markers)")
            continue
        for _, row in sub.iterrows():
            tau_str  = f"  tau={row['tau']:.2f}" if not pd.isna(row.get('tau', np.nan)) else ''
            padj_val = row['pval_adj']
            padj_str = 'n/a' if pd.isna(padj_val) else f"{padj_val:.2e}"
            print(f"    {row['gene']:<22} log2FC={row['log2fc']:>6.2f}  "
                  f"pct_in={row['pct_in']:.2f}  pval_adj={padj_str}{tau_str}")

    # ── Save filtered results ─────────────────────────────────────────────────
    marker_df_filtered.to_csv(
        os.path.join(output_dir, f'{sample_name}_markers_filtered.csv'), index=False)
    marker_df_filtered.groupby('cluster').head(50).to_csv(
        os.path.join(output_dir, f'{sample_name}_markers_top50.csv'), index=False)
    thresh_df.to_csv(
        os.path.join(output_dir, f'{sample_name}_marker_thresholds.csv'), index=False)

    # ── Slim export: cluster / FBgn ID / p-value / adj p-value / rescue info ─
    export_cols = ['cluster', 'gene', 'pval', 'pval_adj', 'rescue_strategy']
    if 'tau' in marker_df_filtered.columns:
        export_cols.append('tau')
    marker_export = (
        marker_df_filtered[export_cols]
        .rename(columns={
            'gene':     'flybase_id',
            'pval':     'p_value',
            'pval_adj': 'adj_p_value',
        })
        .sort_values(['cluster', 'adj_p_value'])
    )
    marker_export_path = os.path.join(
        output_dir, f'{sample_name}_markers_cluster_fbgn_pval.csv')
    marker_export.to_csv(marker_export_path, index=False)
    print(f"\n  Saved marker export ({len(marker_export):,} rows): "
          f"{sample_name}_markers_cluster_fbgn_pval.csv")

    # Scanpy plots
    sc.pl.rank_genes_groups(adata, n_genes=25,
                            save=f'_{sample_name}_ranked_genes.pdf',
                            key='rank_genes_groups')
    try:
        sc.pl.rank_genes_groups_heatmap(adata, n_genes=10,
                                        save=f'_{sample_name}_top10_heatmap.pdf',
                                        show_gene_labels=True, cmap='viridis',
                                        key='rank_genes_groups')
        sc.pl.rank_genes_groups_dotplot(adata, n_genes=5,
                                        save=f'_{sample_name}_top5_dotplot.pdf',
                                        key='rank_genes_groups')
    except Exception as e:
        print(f"  WARNING: Some marker plots failed: {e}")

    return marker_df_filtered


# ─────────────────────────────────────────────────────────────────────────────
# FlyEnrichr
# ─────────────────────────────────────────────────────────────────────────────

def _submit_to_flyenrichr(gene_symbols, description="gene_list"):
    url = 'https://maayanlab.cloud/FlyEnrichr/addList'
    genes_str = '\n'.join(str(g).strip() for g in gene_symbols if g)
    payload = {'list': (None, genes_str), 'description': (None, description)}
    resp = requests.post(url, files=payload, timeout=30)
    resp.raise_for_status()
    data = resp.json()
    if 'userListId' not in data:
        raise ValueError(f"No userListId: {data}")
    return data['userListId']


def flyenrichr_analysis(gene_symbols, background_symbols,
                        gene_set_library='GO_Biological_Process_2018',
                        description="gene_list"):
    BASE = 'https://maayanlab.cloud/FlyEnrichr'
    try:
        fg_id = _submit_to_flyenrichr(gene_symbols, description)
        bg_id = _submit_to_flyenrichr(background_symbols, f"{description}_background")
        url = (f"{BASE}/enrich?userListId={fg_id}"
               f"&backgroundType={gene_set_library}"
               f"&backgroundListId={bg_id}")
        resp = requests.get(url, timeout=30)
        resp.raise_for_status()
        data = resp.json()
        if gene_set_library not in data:
            return None
        results = []
        for entry in data[gene_set_library]:
            results.append({
                'term':           entry[1],
                'p_value':        entry[2],
                'z_score':        entry[3],
                'combined_score': entry[4],
                'genes':          entry[5],
                'adj_p_value':    entry[6],
            })
        return pd.DataFrame(results) if results else None
    except Exception as e:
        print(f"    ERROR ({gene_set_library}): {e}")
        return None


# ─────────────────────────────────────────────────────────────────────────────
# Background gene set
# ─────────────────────────────────────────────────────────────────────────────

def build_background(adata, fbgn_to_symbol, output_dir, sample_name, min_cells=3):
    """
    Return (gene_symbols, background_fbgn) for all genes detected in >= min_cells cells.

    Excluded from background (not in FlyEnrichr, not meaningful):
      - Wolbachia genes (GQX67_*)
      - Transposable elements (FBti*, *_transposable_element)

    mt / ribo / cell_cycle genes are intentionally KEPT in the background
    so that enrichment p-values are correctly calibrated against the full
    expressed transcriptome.

    Also writes *_background_genes.csv (flybase_id, gene_symbol) to output_dir.
    """
    var_names = adata.raw.var_names if adata.raw is not None else adata.var_names
    X = adata.raw.X if adata.raw is not None else adata.X
    if scipy.sparse.issparse(X):
        n_cells_per_gene = np.array((X > 0).sum(axis=0)).flatten()
    else:
        n_cells_per_gene = (X > 0).sum(axis=0)

    expressed_mask = n_cells_per_gene >= min_cells

    var_series = pd.Series(var_names)

    # Exclude bacterial genes and transposable elements
    bact_mask = var_series.str.startswith('GQX67_').values
    te_mask   = _is_te(var_series).values

    keep = expressed_mask & ~bact_mask & ~te_mask
    background_fbgn = var_names[keep].tolist()
    symbols, n_unmapped = symbols_from_fbgn(background_fbgn, fbgn_to_symbol)

    print(f"  Background: {len(background_fbgn):,} genes → {len(symbols):,} symbols "
          f"({n_unmapped:,} unmapped, {int(bact_mask.sum()):,} bacterial excluded, "
          f"{int(te_mask.sum()):,} TEs excluded)")

    # ── Write background gene file ────────────────────────────────────────────
    bg_df = pd.DataFrame({
        'flybase_id':  background_fbgn,
        'gene_symbol': [fbgn_to_symbol.get(g, '') for g in background_fbgn],
    })
    bg_path = os.path.join(output_dir, f'{sample_name}_background_genes.csv')
    bg_df.to_csv(bg_path, index=False)
    print(f"  Saved background genes ({len(bg_df):,} rows): "
          f"{sample_name}_background_genes.csv")

    return symbols, background_fbgn


# ─────────────────────────────────────────────────────────────────────────────
# Enrichment network plot
# ─────────────────────────────────────────────────────────────────────────────

def plot_enrichment_network(sig_df, output_dir, sample_name,
                            jaccard_thresh=0.3, top_n_per_cluster=15,
                            min_adj_p=0.1):
    """
    Enrichment map network plot for GO results.

    Nodes  : GO terms (top N significant per cluster)
    Edges  : Jaccard similarity of gene sets >= jaccard_thresh
    Color  : cluster with strongest significance for that term
    Size   : -log10(adj p-value)
    Layout : spring layout (Fruchterman-Reingold)
    """
    if not NETWORKX_AVAILABLE:
        print("  ⚠️  networkx not installed — skipping network plot")
        print("     mamba install -c conda-forge networkx")
        return

    print("\n  Building GO enrichment network ...")

    go_sig = sig_df[
        sig_df['library'].str.startswith('GO_') &
        (sig_df['adj_p_value'] < min_adj_p)
    ].copy()

    if len(go_sig) == 0:
        print("  No significant GO terms — skipping network")
        return

    go_sig['neg_log10p'] = -np.log10(go_sig['adj_p_value'] + 1e-300)
    go_sig['term_short'] = go_sig['term'].apply(lambda x: x.split('(')[0][:45].rstrip())

    def _parse_genes(g):
        if isinstance(g, list): return set(g)
        if isinstance(g, str):  return set(g.replace(';', ',').split(','))
        return set()

    go_sig['gene_set'] = go_sig['genes'].apply(_parse_genes)

    # Top N terms per cluster by significance
    top_terms = (go_sig.sort_values('adj_p_value')
                       .groupby('cluster')
                       .head(top_n_per_cluster)['term'].unique())
    go_plot = go_sig[go_sig['term'].isin(top_terms)].copy()

    # One row per term: use cluster where it is most significant
    term_best = (go_plot.sort_values('adj_p_value')
                        .drop_duplicates('term')
                        [['term', 'term_short', 'neg_log10p', 'cluster', 'gene_set']]
                        .set_index('term'))

    terms = term_best.index.tolist()
    print(f"  Terms in network: {len(terms)}")
    if len(terms) < 3:
        print("  Too few terms — skipping network")
        return

    # Build graph
    G = nx.Graph()
    for term in terms:
        row = term_best.loc[term]
        G.add_node(term,
                   label      = row['term_short'],
                   cluster    = str(row['cluster']),
                   neg_log10p = float(row['neg_log10p']),
                   gene_set   = row['gene_set'])

    n_edges = 0
    for i, t1 in enumerate(terms):
        for t2 in terms[i+1:]:
            g1 = term_best.loc[t1, 'gene_set']
            g2 = term_best.loc[t2, 'gene_set']
            union = g1 | g2
            if not union:
                continue
            j = len(g1 & g2) / len(union)
            if j >= jaccard_thresh:
                G.add_edge(t1, t2, weight=j)
                n_edges += 1

    print(f"  Edges (Jaccard ≥ {jaccard_thresh}): {n_edges}")

    # Remove isolates if enough connected nodes remain
    connected = [n for n in G.nodes if G.degree(n) > 0]
    if len(connected) >= 5:
        G.remove_nodes_from([n for n in list(G.nodes) if G.degree(n) == 0])
        print(f"  Nodes after removing isolates: {G.number_of_nodes()}")

    # Relax threshold if too sparse
    if G.number_of_nodes() < 3:
        print("  Too few connected nodes — lowering Jaccard threshold to 0.1")
        jaccard_thresh = 0.1
        for i, t1 in enumerate(terms):
            for t2 in terms[i+1:]:
                if G.has_edge(t1, t2):
                    continue
                g1 = term_best.loc[t1, 'gene_set'] if t1 in term_best.index else set()
                g2 = term_best.loc[t2, 'gene_set'] if t2 in term_best.index else set()
                union = g1 | g2
                if not union:
                    continue
                j = len(g1 & g2) / len(union)
                if j >= jaccard_thresh:
                    G.add_edge(t1, t2, weight=j)
        print(f"  Edges after relaxation: {G.number_of_edges()}")

    np.random.seed(42)
    pos = (nx.spring_layout(G, weight='weight', k=2.5, iterations=100, seed=42)
           if G.number_of_edges() > 0 else nx.circular_layout(G))

    cmap = plt.cm.get_cmap('tab20')
    cluster_list = sorted(set(nx.get_node_attributes(G, 'cluster').values()),
                          key=lambda x: int(x) if str(x).isdigit() else x)
    c_colors     = {c: cmap(i % 20) for i, c in enumerate(cluster_list)}
    node_colors  = [c_colors[G.nodes[n]['cluster']] for n in G.nodes]
    node_sizes   = [G.nodes[n]['neg_log10p'] * 80 + 100 for n in G.nodes]
    edge_weights = [G[u][v]['weight'] * 3 for u, v in G.edges]
    labels       = {n: G.nodes[n]['label'] for n in G.nodes}

    # ── Plot 1: full network ──────────────────────────────────────────────────
    fig, ax = plt.subplots(figsize=(16, 14))
    nx.draw_networkx_edges(G, pos, ax=ax, width=edge_weights,
                           alpha=0.35, edge_color='#888888')
    nx.draw_networkx_nodes(G, pos, ax=ax, node_color=node_colors,
                           node_size=node_sizes, alpha=0.9,
                           linewidths=0.8, edgecolors='white')
    nx.draw_networkx_labels(G, pos, labels=labels, ax=ax,
                            font_size=6.5, font_color='black', font_weight='bold',
                            bbox=dict(boxstyle='round,pad=0.2', facecolor='white',
                                      alpha=0.6, edgecolor='none'))
    legend_els = [plt.scatter([], [], c=[c_colors[c]], s=80,
                              label=f'Cluster {c}', alpha=0.9)
                  for c in cluster_list]
    for sig_val, lbl in [(5, 'p=0.01'), (10, 'p=1e-10'), (20, 'p=1e-20')]:
        legend_els.append(plt.scatter([], [], c='gray',
                                      s=sig_val * 80 + 100, alpha=0.6, label=lbl))
    ax.legend(handles=legend_els, loc='upper left', bbox_to_anchor=(1.01, 1),
              fontsize=9, title='Cluster / Significance', title_fontsize=9,
              framealpha=0.8)
    ax.set_title(f'GO Enrichment Network — {sample_name}\n'
                 f'Node color = cluster  |  Node size = −log₁₀(p)  |  '
                 f'Edge = Jaccard ≥ {jaccard_thresh}',
                 fontsize=12, fontweight='bold')
    ax.axis('off')
    plt.tight_layout()
    plt.savefig(f"{output_dir}/{sample_name}_GO_enrichment_network.pdf",
                dpi=150, bbox_inches='tight')
    plt.close()
    print(f"    Saved: GO_enrichment_network.pdf")

    # ── Plot 2: per-cluster highlight panels ──────────────────────────────────
    ncols = min(3, len(cluster_list))
    nrows = int(np.ceil(len(cluster_list) / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(ncols * 7, nrows * 6))
    axes = np.array(axes).flatten()

    for i, cl in enumerate(cluster_list):
        ax = axes[i]
        hi_nodes = [n for n in G.nodes if G.nodes[n]['cluster'] == cl]
        lo_nodes = [n for n in G.nodes if G.nodes[n]['cluster'] != cl]
        hi_sizes = [G.nodes[n]['neg_log10p'] * 100 + 120 for n in hi_nodes]

        nx.draw_networkx_edges(G, pos, ax=ax,
                               width=[G[u][v]['weight'] * 2 for u, v in G.edges],
                               alpha=0.2, edge_color='#aaaaaa')
        nx.draw_networkx_nodes(G, pos, nodelist=lo_nodes, ax=ax,
                               node_color='#dddddd', node_size=60, alpha=0.5)
        nx.draw_networkx_nodes(G, pos, nodelist=hi_nodes, ax=ax,
                               node_color=[c_colors[cl]] * len(hi_nodes),
                               node_size=hi_sizes, alpha=0.95,
                               linewidths=1.0, edgecolors='white')
        nx.draw_networkx_labels(G, pos,
                                labels={n: G.nodes[n]['label'] for n in hi_nodes},
                                ax=ax, font_size=6, font_weight='bold',
                                bbox=dict(boxstyle='round,pad=0.2', facecolor='white',
                                          alpha=0.7, edgecolor='none'))
        ax.set_title(f'Cluster {cl}  ({len(hi_nodes)} terms)',
                     fontsize=10, fontweight='bold', color=c_colors[cl])
        ax.axis('off')

    for j in range(i + 1, len(axes)):
        axes[j].axis('off')

    plt.suptitle(f'GO Network — per-cluster highlight — {sample_name}',
                 fontsize=12, fontweight='bold')
    plt.tight_layout()
    plt.savefig(f"{output_dir}/{sample_name}_GO_network_per_cluster.pdf",
                dpi=150, bbox_inches='tight')
    plt.close()
    print(f"    Saved: GO_network_per_cluster.pdf")


# ─────────────────────────────────────────────────────────────────────────────
# Enrichment visualisation — bar plots + heatmap
# ─────────────────────────────────────────────────────────────────────────────

def _plot_enrichment(sig_df, output_dir, sample_name, clusters, combine_go=True):
    print("\n  Generating enrichment plots ...")

    category_colors = {
        'Biological_Process': '#e74c3c',
        'Molecular_Function': '#3498db',
        'Cellular_Component': '#2ecc71',
    }

    if combine_go:
        go_sig = sig_df[sig_df['library'].str.startswith('GO_')].copy()
        go_sig['go_category'] = (go_sig['library']
                                 .str.replace('GO_', '')
                                 .str.replace('_2018', ''))

        if len(go_sig) > 0:
            n = len(clusters)
            fig, axes = plt.subplots(n, 1, figsize=(16, 5 * n))
            if n == 1: axes = [axes]

            for idx, cluster in enumerate(clusters):
                ax = axes[idx]
                cgo = (go_sig[go_sig['cluster'] == cluster]
                       .sort_values('adj_p_value').head(15))

                if len(cgo) == 0:
                    ax.text(0.5, 0.5, f'No significant GO terms\nCluster {cluster}',
                            ha='center', va='center', transform=ax.transAxes)
                    ax.axis('off')
                    continue

                cgo = cgo.copy()
                cgo['term_short'] = cgo['term'].apply(
                    lambda x: x.split('(')[0][:55].rstrip())
                cgo['-log10p'] = -np.log10(cgo['adj_p_value'] + 1e-300)
                y_pos = range(len(cgo))
                colors = [category_colors.get(cat, '#888888') for cat in cgo['go_category']]
                ax.barh(y_pos, cgo['-log10p'], color=colors, alpha=0.75)
                ax.set_yticks(y_pos)
                ax.set_yticklabels(cgo['term_short'], fontsize=9)
                ax.set_xlabel('-log10(adj p-value)', fontsize=11)
                ax.set_title(f'Cluster {cluster} — Top GO Terms',
                             fontsize=12, fontweight='bold')
                ax.invert_yaxis()
                ax.grid(axis='x', alpha=0.3)

                if idx == 0:
                    legend_els = [Patch(facecolor=v, label=k.replace('_', ' '), alpha=0.75)
                                  for k, v in category_colors.items()]
                    ax.legend(handles=legend_els, loc='lower right', fontsize=9)

            plt.tight_layout()
            plt.savefig(f"{output_dir}/{sample_name}_GO_combined_barplots.pdf",
                        dpi=150, bbox_inches='tight')
            plt.close()
            print(f"    Saved: GO_combined_barplots.pdf")

            # Cross-cluster heatmap
            top_terms = go_sig.nsmallest(60, 'adj_p_value')['term'].unique()[:30]
            hdata = []
            for term in top_terms:
                row = {'term': term.split('(')[0][:50].rstrip()}
                for cl in clusters:
                    sub = go_sig[(go_sig['cluster'] == cl) & (go_sig['term'] == term)]
                    row[str(cl)] = (min(-np.log10(sub.iloc[0]['adj_p_value'] + 1e-300), 50)
                                    if len(sub) > 0 else 0)
                hdata.append(row)

            hdf = pd.DataFrame(hdata).set_index('term')
            fig, ax = plt.subplots(figsize=(max(10, len(clusters)),
                                            max(8, len(top_terms) * 0.35)))
            sns.heatmap(hdf, cmap='YlOrRd', ax=ax,
                        cbar_kws={'label': '-log10(adj p-value)'},
                        linewidths=0.3, linecolor='lightgray')
            ax.set_xlabel('Leiden Cluster', fontsize=12)
            ax.set_ylabel('GO Term', fontsize=12)
            ax.set_title('Top GO Terms Across All Clusters', fontsize=13)
            plt.tight_layout()
            plt.savefig(f"{output_dir}/{sample_name}_GO_combined_heatmap.pdf",
                        dpi=150, bbox_inches='tight')
            plt.close()
            print(f"    Saved: GO_combined_heatmap.pdf")

    # Per-library bar plots
    for lib in sig_df['library'].unique():
        lib_sig = sig_df[sig_df['library'] == lib]
        if len(lib_sig) == 0:
            continue
        n = len(clusters)
        fig, axes = plt.subplots(n, 1, figsize=(14, 4 * n))
        if n == 1: axes = [axes]

        for idx, cluster in enumerate(clusters):
            ax = axes[idx]
            top = (lib_sig[lib_sig['cluster'] == cluster]
                   .sort_values('adj_p_value').head(10).copy())
            if len(top) == 0:
                ax.text(0.5, 0.5, f'No significant terms\nCluster {cluster}',
                        ha='center', va='center', transform=ax.transAxes)
                ax.axis('off')
                continue
            top['term_short'] = top['term'].apply(lambda x: x.split('(')[0][:50].rstrip())
            top['-log10p'] = -np.log10(top['adj_p_value'] + 1e-300)
            y_pos = range(len(top))
            colors = ['#d62728' if p < 0.01 else '#ff7f0e' for p in top['adj_p_value']]
            ax.barh(y_pos, top['-log10p'], color=colors, alpha=0.7)
            ax.set_yticks(y_pos)
            ax.set_yticklabels(top['term_short'], fontsize=9)
            ax.set_xlabel('-log10(adj p-value)', fontsize=11)
            ax.set_title(f'Cluster {cluster}', fontsize=12, fontweight='bold')
            ax.invert_yaxis()
            ax.grid(axis='x', alpha=0.3)

        lib_short = lib.replace('_2018', '').replace('_2019', '')
        plt.tight_layout()
        plt.savefig(f"{output_dir}/{sample_name}_enrichment_{lib_short}.pdf",
                    dpi=150, bbox_inches='tight')
        plt.close()
        print(f"    Saved: enrichment_{lib_short}.pdf")


# ─────────────────────────────────────────────────────────────────────────────
# Enrichment per cluster
# ─────────────────────────────────────────────────────────────────────────────

def enrichment_analysis_per_cluster(adata, marker_df, fbgn_to_symbol,
                                    output_dir, sample_name, combine_go=True):
    print("\n" + "="*60)
    print("PATHWAY ENRICHMENT ANALYSIS (FlyEnrichr)")
    print("="*60)

    if fbgn_to_symbol is None:
        print("  ERROR: No gene symbol mapping — cannot run enrichment")
        return None

    print("\n  Building background gene set ...")
    background_symbols, _ = build_background(
        adata, fbgn_to_symbol, output_dir, sample_name)
    if len(background_symbols) < 100:
        print("  WARNING: Very small background — check mapping file")

    libraries = [
        'GO_Biological_Process_2018',
        'GO_Molecular_Function_2018',
        'GO_Cellular_Component_2018',
        'KEGG_2019',
        'WikiPathways_2018',
    ]

    clusters = sorted(marker_df['cluster'].unique(),
                      key=lambda x: int(x) if str(x).isdigit() else x)
    all_results = []

    for cluster in clusters:
        print(f"\n{'='*50}\n  Cluster {cluster}\n{'='*50}")
        cmarkers = (marker_df[marker_df['cluster'] == cluster]
                    .sort_values('score', ascending=False))

        if len(cmarkers) < 5:
            print(f"  Skipping: only {len(cmarkers)} markers (need ≥ 5)")
            continue

        print(f"  Using {len(cmarkers)} markers for enrichment")
        genes_fbgn = cmarkers['gene'].tolist()
        genes_symbols, n_unmapped = symbols_from_fbgn(genes_fbgn, fbgn_to_symbol)
        print(f"  {len(genes_fbgn)} FBgn IDs → {len(genes_symbols)} symbols "
              f"({n_unmapped} unmapped)")

        if len(genes_symbols) < 5:
            print(f"  Skipping: too few mapped symbols")
            continue

        for lib in libraries:
            print(f"  [{lib}] ", end='', flush=True)
            result_df = flyenrichr_analysis(
                genes_symbols, background_symbols,
                gene_set_library=lib,
                description=f"cluster_{cluster}",
            )
            if result_df is not None and len(result_df) > 0:
                sig = result_df[result_df['adj_p_value'] < 0.1]
                result_df['cluster'] = cluster
                result_df['library'] = lib
                print(f"{len(result_df)} terms, {len(sig)} significant (adj_p<0.1)")
                all_results.append(result_df)
            else:
                print("no results")
            time.sleep(0.5)

    if not all_results:
        print("\n  No enrichment results obtained")
        return None

    combined_df = pd.concat(all_results, ignore_index=True)
    combined_df.to_csv(f"{output_dir}/{sample_name}_flyenrichr_all_results.csv", index=False)

    sig_df = combined_df[combined_df['adj_p_value'] < 0.1].copy()
    sig_df.to_csv(f"{output_dir}/{sample_name}_flyenrichr_significant.csv", index=False)
    print(f"\n  Total terms: {len(combined_df)}")
    print(f"  Significant (adj_p < 0.1): {len(sig_df)}")

    if combine_go:
        go_sig = sig_df[sig_df['library'].str.startswith('GO_')].copy()
        if len(go_sig) > 0:
            go_sig.sort_values(['cluster', 'adj_p_value']).to_csv(
                f"{output_dir}/{sample_name}_GO_combined_significant.csv", index=False)
            go_sig.groupby('cluster').head(20).to_csv(
                f"{output_dir}/{sample_name}_GO_combined_top20_per_cluster.csv", index=False)
            print(f"  Significant GO terms: {len(go_sig)}")

    # Summary to stdout
    print(f"\n{'='*60}\nENRICHMENT SUMMARY\n{'='*60}")
    for lib in libraries:
        lib_sig = sig_df[sig_df['library'] == lib]
        if len(lib_sig) == 0:
            continue
        print(f"\n{lib}")
        for cluster in clusters:
            top = lib_sig[lib_sig['cluster'] == cluster].sort_values('adj_p_value').head(3)
            if len(top) == 0:
                continue
            print(f"  Cluster {cluster}:")
            for _, row in top.iterrows():
                term = row['term'].split('(')[0][:60]
                print(f"    {term:<60} p_adj={row['adj_p_value']:.2e}")

    _plot_enrichment(sig_df, output_dir, sample_name, clusters, combine_go=combine_go)
    plot_enrichment_network(sig_df, output_dir, sample_name)

    return combined_df


# ─────────────────────────────────────────────────────────────────────────────
# Main
# ─────────────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description='Cluster analysis with marker genes and FlyEnrichr GO enrichment',
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument('--input',   '-i', required=True,  help='Path to integrated .h5ad')
    parser.add_argument('--output',  '-o', default='cluster_analysis', help='Output directory')
    parser.add_argument('--sample',  '-s', default='sample', help='Sample name prefix')
    parser.add_argument('--mapping', '-map', required=True,
                        help='transcripts_to_genes.txt (FBgn -> symbol)')
    parser.add_argument('--method',  '-m', default='wilcoxon',
                        choices=['wilcoxon', 't-test', 'logreg'])
    parser.add_argument('--log2fc-min', type=float, default=0.25,
                        help='Floor log2FC for adaptive threshold sweep (default: 0.25)')
    parser.add_argument('--pct-min',   type=float, default=0.05,
                        help='Floor pct_in for adaptive threshold sweep (default: 0.05)')
    parser.add_argument('--pct-ratio', type=float, default=1.0,
                        help='Floor pct_ratio for adaptive threshold sweep (default: 1.0)')
    parser.add_argument('--skip-enrichment', action='store_true')
    parser.add_argument('--no-combine-go',   action='store_true')

    args = parser.parse_args()
    os.makedirs(args.output, exist_ok=True)
    sc.settings.figdir = args.output

    print("="*60)
    print("LOADING DATA")
    print("="*60)
    adata = sc.read_h5ad(args.input)
    print(f"  Cells: {adata.n_obs:,}  Genes: {adata.n_vars:,}  "
          f"Clusters: {adata.obs['leiden'].nunique()}")
    if adata.raw is None:
        print("  WARNING: adata.raw is None — DE will use adata.X")
    else:
        print(f"  adata.raw: {adata.raw.n_vars:,} genes (will be used for DE)")

    fbgn_to_symbol = load_fbgn_to_symbol_mapping(args.mapping)
    if fbgn_to_symbol is None and not args.skip_enrichment:
        print("\nERROR: Could not load gene mapping. Use --skip-enrichment to skip GO.")
        return

    plot_transcriptional_activity(adata, args.output, args.sample)

    marker_df = find_marker_genes(
        adata, args.output, args.sample,
        method=args.method,
        log2fc_min=args.log2fc_min,
        pct_min=args.pct_min,
        pct_ratio_min=args.pct_ratio,
    )

    if not args.skip_enrichment:
        enrichment_analysis_per_cluster(
            adata, marker_df, fbgn_to_symbol,
            args.output, args.sample,
            combine_go=not args.no_combine_go,
        )

    print("\n" + "="*60)
    print("DONE")
    print("="*60)
    print(f"  Results: {args.output}/")


if __name__ == '__main__':
    main()