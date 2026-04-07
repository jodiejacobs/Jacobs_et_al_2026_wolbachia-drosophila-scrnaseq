'''
Comprehensive cluster analysis: transcriptional activity, marker genes, and pathway enrichment
Uses FlyEnrichr API for automated pathway analysis with FlyBase annotations

Key fixes vs previous version:
  - DE run on .raw (log-normalised counts), NOT scaled .X
  - Marker filtering uses log2fc > 1.0 + pct_in >= 0.10 + pct_ratio >= 1.5
  - Background gene set passed to FlyEnrichr (all detected genes in dataset)
  - Significance filter (adj_p < 0.05) applied before combined_score sorting
  - mt/ribo/cell_cycle genes excluded from enrichment input
'''

import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import argparse
import os
import requests
import json
import time
import gzip
from scipy.stats import kruskal
from matplotlib.patches import Patch
import scipy.sparse
from scipy.sparse import issparse


# ─────────────────────────────────────────────────────────────────────────────
# Gene symbol mapping
# ─────────────────────────────────────────────────────────────────────────────

def load_fbgn_to_symbol_mapping(mapping_file):
    """
    Load FBgn -> gene symbol from transcripts_to_genes.txt
    Format: transcript_id  FBgn_id  gene_symbol  ...
    """
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
        sample = list(fbgn_to_symbol.items())[:5]
        for fbgn, sym in sample:
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
    summary.to_csv(os.path.join(output_dir, f'{sample_name}_transcriptional_activity_summary.csv'))

    # Kruskal-Wallis
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

    print(f"\nn_counts: H={h_counts:.2f}  p={p_counts:.2e}  η²={eta_counts:.4f} ({eta_label(eta_counts)})")
    print(f"n_genes:  H={h_genes:.2f}  p={p_genes:.2e}  η²={eta_genes:.4f} ({eta_label(eta_genes)})")

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
# Marker genes  (FIX: use .raw, strict filters)
# ─────────────────────────────────────────────────────────────────────────────

def find_marker_genes(adata, output_dir, sample_name, method='wilcoxon',
                      log2fc_min=1.0, pct_min=0.10, pct_ratio_min=1.5):
    """
    Find cluster marker genes.

    Critical: DE is run on adata.raw (log-normalised counts).
    Scaled .X is NOT appropriate for DE — it distorts fold changes.

    Filters applied:
      log2fc     >= log2fc_min   (default 1.0 = 2× fold change)
      pct_in     >= pct_min      (default 10% of cells in cluster express gene)
      pct_ratio  >= pct_ratio_min (expressed proportionally more in cluster vs rest)
      pval_adj   < 0.05
    """
    print("\n" + "="*60)
    print("DIFFERENTIAL GENE EXPRESSION")
    print("="*60)
    print(f"  Method: {method}")
    print(f"  Filters: log2fc>={log2fc_min}, pct_in>={pct_min}, "
          f"pct_ratio>={pct_ratio_min}, pval_adj<0.05")

    # ── Use .raw for DE ───────────────────────────────────────────────────────
    if adata.raw is not None:
        print("  Using adata.raw for DE (log-normalised, pre-scale counts)")
        adata_de = adata.raw.to_adata()
        adata_de.X = scipy.sparse.csr_matrix(adata_de.X)  # ensure sparse
        adata_de.obs['leiden'] = adata.obs['leiden'].values
    else:
        print("  WARNING: adata.raw is None — using adata.X (check this is "
              "log-normalised, not scaled!)")
        adata_de = adata.copy()
        adata_de.X = scipy.sparse.csr_matrix(adata_de.X)  # ensure sparse

    # Filter out Wolbachia transcripts (GQX67_ prefix)
    n_before_bact = adata_de.n_vars
    adata_de = adata_de[:, ~adata_de.var_names.str.startswith('GQX67_')].copy()
    print(f"  Genes after Wolbachia filter: {adata_de.n_vars:,} / {n_before_bact:,} "
          f"({n_before_bact - adata_de.n_vars:,} bacterial transcripts removed)")

    # Filter to genes expressed in >= 3 cells (matches build_background threshold)
    n_cells_per_gene = np.array((adata_de.X > 0).sum(axis=0)).flatten()
    n_before = adata_de.n_vars
    adata_de = adata_de[:, n_cells_per_gene >= 3].copy()
    print(f"  Genes after expression filter: {adata_de.n_vars:,} / {n_before:,}")

    sc.settings.n_jobs = -1  # use all available cores
    sc.tl.rank_genes_groups(
        adata_de, 'leiden', method=method,
        key_added='rank_genes_groups',
        tie_correct=True,
        rankby_abs=False,
        pts=True,
    )

    # Copy result back so plotting functions work on original adata
    adata.uns['rank_genes_groups'] = adata_de.uns['rank_genes_groups']

    result = adata.uns['rank_genes_groups']
    groups = result['names'].dtype.names

    rows = []
    for grp in groups:
        n = len(result['names'][grp])
        for i in range(n):
            lfc  = result['logfoldchanges'][grp][i]
            padj = result['pvals_adj'][grp][i]
            if pd.isna(lfc) or np.isinf(lfc):
                continue
            pct_in   = result['pts'][grp][i]      if 'pts'      in result else np.nan
            pct_rest = result['pts_rest'][grp][i]  if 'pts_rest' in result else np.nan
            rows.append(dict(
                cluster=grp,
                gene=result['names'][grp][i],
                log2fc=lfc,
                pval=result['pvals'][grp][i],
                pval_adj=padj,
                score=result['scores'][grp][i],
                pct_in=pct_in,
                pct_rest=pct_rest,
            ))

    marker_df = pd.DataFrame(rows)
    print(f"  Raw DE results: {len(marker_df)} gene×cluster entries")

    # ── Strict filtering ──────────────────────────────────────────────────────
    mask = (
        (marker_df['log2fc']  >= log2fc_min) &
        (marker_df['pval_adj'] < 0.05)
    )
    if not marker_df['pct_in'].isna().all():
        mask &= (marker_df['pct_in'] >= pct_min)
        pct_ratio = marker_df['pct_in'] / (marker_df['pct_rest'] + 1e-9)
        mask &= (pct_ratio >= pct_ratio_min)

    marker_df_filtered = marker_df[mask].copy()
    print(f"  After filtering:  {len(marker_df_filtered)} entries")

    # Print top 10 per cluster
    print("\n  Top 10 markers per cluster:")
    for grp in groups:
        sub = marker_df_filtered[marker_df_filtered['cluster'] == grp].head(10)
        print(f"\n  Cluster {grp}:")
        if len(sub) == 0:
            print("    (no markers passed filter)")
            continue
        for _, row in sub.iterrows():
            print(f"    {row['gene']:<22} log2FC={row['log2fc']:>6.2f}  "
                  f"pct_in={row['pct_in']:.2f}  pval_adj={row['pval_adj']:.2e}")

    # Save
    marker_df.to_csv(os.path.join(output_dir, f'{sample_name}_markers_all.csv'), index=False)
    marker_df_filtered.to_csv(
        os.path.join(output_dir, f'{sample_name}_markers_filtered.csv'), index=False)
    marker_df_filtered.groupby('cluster').head(50).to_csv(
        os.path.join(output_dir, f'{sample_name}_markers_top50.csv'), index=False)

    # Plots (top 10 from filtered set)
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
# FlyEnrichr  (FIX: pass background gene list)
# ─────────────────────────────────────────────────────────────────────────────

def _submit_to_flyenrichr(gene_symbols, description="gene_list"):
    """Submit a gene list to FlyEnrichr and return userListId."""
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
    """
    Run FlyEnrichr enrichment with a proper background gene set.

    Parameters
    ----------
    gene_symbols      : list of str  — foreground (marker genes)
    background_symbols: list of str  — background (all genes detected in dataset)
    gene_set_library  : str
    description       : str

    Notes
    -----
    FlyEnrichr supports a background via a second list submission.
    We submit both foreground and background, then query enrich with
    the background userListId as the reference.
    """
    BASE = 'https://maayanlab.cloud/FlyEnrichr'

    try:
        # Submit foreground
        fg_id = _submit_to_flyenrichr(gene_symbols, description)

        # Submit background
        bg_id = _submit_to_flyenrichr(background_symbols, f"{description}_background")

        # Query enrichment with background
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
                'term':          entry[1],
                'p_value':       entry[2],
                'z_score':       entry[3],
                'combined_score': entry[4],
                'genes':         entry[5],
                'adj_p_value':   entry[6],
            })

        return pd.DataFrame(results) if results else None

    except Exception as e:
        print(f"    ERROR ({gene_set_library}): {e}")
        return None


# ─────────────────────────────────────────────────────────────────────────────
# Build background gene set
# ─────────────────────────────────────────────────────────────────────────────

def build_background(adata, fbgn_to_symbol, min_cells=3):
    """
    Return gene symbols for all genes detected in >= min_cells cells,
    excluding mt / ribo / cell_cycle / Wolbachia genes.
    """
    var_names = adata.raw.var_names if adata.raw is not None else adata.var_names

    if adata.raw is not None:
        X = adata.raw.X
    else:
        X = adata.X
    import scipy.sparse
    if scipy.sparse.issparse(X):
        n_cells_per_gene = np.array((X > 0).sum(axis=0)).flatten()
    else:
        n_cells_per_gene = (X > 0).sum(axis=0)

    expressed_mask = n_cells_per_gene >= min_cells

    # Build exclude mask — initialise first, then OR in each category
    exclude = np.zeros(len(var_names), dtype=bool)

    # Exclude QC gene classes flagged in adata.var
    for flag in ('mt', 'ribo', 'cell_cycle'):
        if flag in adata.var.columns:
            in_raw = pd.Series(
                adata.var[flag].values,
                index=adata.var_names,
            ).reindex(var_names).fillna(False).values.astype(bool)
            exclude |= in_raw

    # Exclude Wolbachia transcripts
    bact_mask = pd.Series(var_names).str.startswith('GQX67_').values
    exclude |= bact_mask

    keep = expressed_mask & ~exclude
    background_fbgn = var_names[keep].tolist()

    symbols, n_unmapped = symbols_from_fbgn(background_fbgn, fbgn_to_symbol)
    print(f"  Background: {len(background_fbgn)} genes → {len(symbols)} symbols "
          f"({n_unmapped} unmapped, {bact_mask.sum()} bacterial excluded)")
    return symbols

# ─────────────────────────────────────────────────────────────────────────────
# Enrichment per cluster
# ─────────────────────────────────────────────────────────────────────────────

def enrichment_analysis_per_cluster(adata, marker_df, fbgn_to_symbol,
                                    output_dir, sample_name,
                                    top_n=200, combine_go=True):
    print("\n" + "="*60)
    print("PATHWAY ENRICHMENT ANALYSIS (FlyEnrichr)")
    print("="*60)

    if fbgn_to_symbol is None:
        print("  ERROR: No gene symbol mapping — cannot run enrichment")
        return None

    # Build background once
    print("\n  Building background gene set …")
    background_symbols = build_background(adata, fbgn_to_symbol)
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

        # ── Select marker genes ───────────────────────────────────────────────
        # Already filtered to log2fc>=1, pct_in>=0.1, pval_adj<0.05 by find_marker_genes
        # Take top_n by score
        cmarkers = (marker_df[marker_df['cluster'] == cluster]
                    .sort_values('score', ascending=False)
                    .head(top_n))

        if len(cmarkers) < 5:
            print(f"  Skipping: only {len(cmarkers)} markers (need ≥ 5)")
            continue

        genes_fbgn = cmarkers['gene'].tolist()
        genes_symbols, n_unmapped = symbols_from_fbgn(genes_fbgn, fbgn_to_symbol)
        print(f"  Markers: {len(genes_fbgn)} FBgn IDs → {len(genes_symbols)} symbols "
              f"({n_unmapped} unmapped)")

        if len(genes_symbols) < 5:
            print(f"  Skipping: too few mapped symbols")
            continue

        # Remove any background exclusions from foreground too
        genes_symbols = [g for g in genes_symbols if g in set(background_symbols)]
        print(f"  Foreground after background intersection: {len(genes_symbols)}")

        for lib in libraries:
            print(f"  [{lib}] ", end='', flush=True)
            result_df = flyenrichr_analysis(
                genes_symbols, background_symbols,
                gene_set_library=lib,
                description=f"cluster_{cluster}",
            )
            if result_df is not None and len(result_df) > 0:
                # Filter to significant terms only before saving
                sig = result_df[result_df['adj_p_value'] < 0.05]
                result_df['cluster'] = cluster
                result_df['library'] = lib
                sig_count = len(sig)
                print(f"{len(result_df)} terms, {sig_count} significant")
                all_results.append(result_df)
            else:
                print("no results")
            time.sleep(0.5)   # be polite to the API

    if not all_results:
        print("\n  No enrichment results obtained")
        return None

    combined_df = pd.concat(all_results, ignore_index=True)
    combined_df.to_csv(f"{output_dir}/{sample_name}_flyenrichr_all_results.csv", index=False)

    # Significant subset
    sig_df = combined_df[combined_df['adj_p_value'] < 0.05].copy()
    sig_df.to_csv(f"{output_dir}/{sample_name}_flyenrichr_significant.csv", index=False)
    print(f"\n  Total terms: {len(combined_df)}")
    print(f"  Significant (adj_p < 0.05): {len(sig_df)}")

    # Combined GO
    if combine_go:
        go_sig = sig_df[sig_df['library'].str.startswith('GO_')].copy()
        if len(go_sig) > 0:
            go_sig = go_sig.sort_values(['cluster', 'adj_p_value'])
            go_sig.to_csv(f"{output_dir}/{sample_name}_GO_combined_significant.csv", index=False)
            go_sig.groupby('cluster').head(20).to_csv(
                f"{output_dir}/{sample_name}_GO_combined_top20_per_cluster.csv", index=False)
            print(f"  Significant GO terms: {len(go_sig)}")

    # Summary printout
    print(f"\n{'='*60}\nENRICHMENT SUMMARY\n{'='*60}")
    for lib in libraries:
        lib_sig = sig_df[sig_df['library'] == lib]
        if len(lib_sig) == 0:
            continue
        print(f"\n{lib}")
        for cluster in clusters:
            top = (lib_sig[lib_sig['cluster'] == cluster]
                   .sort_values('adj_p_value').head(3))
            if len(top) == 0:
                continue
            print(f"  Cluster {cluster}:")
            for _, row in top.iterrows():
                term = row['term'].split('(')[0][:60]
                print(f"    {term:<60} p_adj={row['adj_p_value']:.2e}")

    # Visualisations
    _plot_enrichment(sig_df, output_dir, sample_name, clusters, combine_go=combine_go)

    return combined_df


# ─────────────────────────────────────────────────────────────────────────────
# Visualisation
# ─────────────────────────────────────────────────────────────────────────────

def _plot_enrichment(sig_df, output_dir, sample_name, clusters, combine_go=True):
    print("\n  Generating enrichment plots …")

    category_colors = {
        'Biological_Process': '#e74c3c',
        'Molecular_Function': '#3498db',
        'Cellular_Component': '#2ecc71',
    }

    if combine_go:
        go_sig = sig_df[sig_df['library'].str.startswith('GO_')].copy()
        go_sig['go_category'] = go_sig['library'].str.replace('GO_', '').str.replace('_2018', '')

        if len(go_sig) > 0:
            # Per-cluster combined GO barplot
            n = len(clusters)
            fig, axes = plt.subplots(n, 1, figsize=(16, 5 * n))
            if n == 1: axes = [axes]

            for idx, cluster in enumerate(clusters):
                ax = axes[idx]
                # Sort by adj_p_value (ascending), not combined_score
                cgo = (go_sig[go_sig['cluster'] == cluster]
                       .sort_values('adj_p_value')
                       .head(15))

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
                colors = [category_colors[cat] for cat in cgo['go_category']]
                ax.barh(y_pos, cgo['-log10p'], color=colors, alpha=0.75)
                ax.set_yticks(y_pos)
                ax.set_yticklabels(cgo['term_short'], fontsize=9)
                ax.set_xlabel('-log10(adj p-value)', fontsize=11)
                ax.set_title(f'Cluster {cluster} — Top GO Terms',
                             fontsize=12, fontweight='bold')
                ax.invert_yaxis()
                ax.grid(axis='x', alpha=0.3)

                if idx == 0:
                    legend_els = [
                        Patch(facecolor=v, label=k.replace('_', ' '), alpha=0.75)
                        for k, v in category_colors.items()
                    ]
                    ax.legend(handles=legend_els, loc='lower right', fontsize=9)

            plt.tight_layout()
            plt.savefig(f"{output_dir}/{sample_name}_GO_combined_barplots.pdf",
                        dpi=150, bbox_inches='tight')
            plt.close()
            print(f"    Saved: GO_combined_barplots.pdf")

            # Cross-cluster GO heatmap
            top_terms = (go_sig.nsmallest(60, 'adj_p_value')['term'].unique()[:30])
            hdata = []
            for term in top_terms:
                row = {'term': term.split('(')[0][:50].rstrip()}
                for cl in clusters:
                    sub = go_sig[(go_sig['cluster'] == cl) & (go_sig['term'] == term)]
                    row[str(cl)] = (min(-np.log10(sub.iloc[0]['adj_p_value'] + 1e-300), 50)
                                    if len(sub) > 0 else 0)
                hdata.append(row)

            hdf = pd.DataFrame(hdata).set_index('term')

            fig, ax = plt.subplots(figsize=(max(10, len(clusters)), max(8, len(top_terms) * 0.35)))
            sns.heatmap(hdf, cmap='YlOrRd', ax=ax,
                        cbar_kws={'label': '-log10(adj p-value)'},
                        linewidths=0.3, linecolor='lightgray')
            ax.set_xlabel('Leiden Cluster', fontsize=12)
            ax.set_ylabel('GO Term', fontsize=12)
            ax.set_title('Top GO Terms Across All Clusters\n(sorted by significance)',
                         fontsize=13)
            plt.tight_layout()
            plt.savefig(f"{output_dir}/{sample_name}_GO_combined_heatmap.pdf",
                        dpi=150, bbox_inches='tight')
            plt.close()
            print(f"    Saved: GO_combined_heatmap.pdf")

    # Individual library barplots
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
                   .sort_values('adj_p_value')
                   .head(10)
                   .copy())

            if len(top) == 0:
                ax.text(0.5, 0.5, f'No significant terms\nCluster {cluster}',
                        ha='center', va='center', transform=ax.transAxes)
                ax.axis('off')
                continue

            top['term_short'] = top['term'].apply(lambda x: x.split('(')[0][:50].rstrip())
            top['-log10p'] = -np.log10(top['adj_p_value'] + 1e-300)
            y_pos = range(len(top))
            colors = ['#d62728' if p < 0.01 else '#ff7f0e' for p in top['adj_p_value']]
            bars = ax.barh(y_pos, top['-log10p'], color=colors, alpha=0.7)
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
# Main
# ─────────────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description='Cluster analysis with marker genes and FlyEnrichr GO enrichment',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='''
Examples:
  python cluster_marker_pathway_analysis.py \\
      --input integrated.h5ad \\
      --output results/cluster_analysis \\
      --sample wolbachia_infection \\
      --mapping /path/to/transcripts_to_genes.txt

  # Adjust stringency
  python cluster_marker_pathway_analysis.py \\
      --input integrated.h5ad --output results \\
      --sample test --mapping t2g.txt \\
      --log2fc-min 0.5 --pct-min 0.05

  # Skip enrichment (faster, markers only)
  python cluster_marker_pathway_analysis.py \\
      --input integrated.h5ad --output results \\
      --sample test --mapping t2g.txt \\
      --skip-enrichment
        '''
    )

    parser.add_argument('--input',  '-i', required=True,  help='Path to integrated .h5ad')
    parser.add_argument('--output', '-o', default='cluster_analysis', help='Output directory')
    parser.add_argument('--sample', '-s', default='sample', help='Sample name prefix')
    parser.add_argument('--mapping', '-map', required=True,
                        help='transcripts_to_genes.txt (FBgn -> symbol)')
    parser.add_argument('--method', '-m', default='wilcoxon',
                        choices=['wilcoxon', 't-test', 'logreg'])
    parser.add_argument('--log2fc-min',  type=float, default=0.5,
                        help='Minimum log2FC for markers (default: 0.5 = 1.41x)')
    parser.add_argument('--pct-min',    type=float, default=0.10,
                        help='Min fraction of cluster cells expressing marker (default: 0.10)')
    parser.add_argument('--pct-ratio',  type=float, default=1.2,
                        help='Min pct_in/pct_rest ratio (default: 1.2)')
    parser.add_argument('--top-n',      type=int,   default=200,
                        help='Max markers per cluster for enrichment (default: 200)')
    parser.add_argument('--skip-enrichment', action='store_true')
    parser.add_argument('--no-combine-go',   action='store_true',
                        help='Skip combined GO visualisation')

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
        print("  WARNING: adata.raw is None — DE will use adata.X (verify it is "
              "log-normalised, not scaled)")
    else:
        print(f"  adata.raw: {adata.raw.n_vars:,} genes (will be used for DE)")

    fbgn_to_symbol = load_fbgn_to_symbol_mapping(args.mapping)
    if fbgn_to_symbol is None and not args.skip_enrichment:
        print("\nERROR: Could not load gene mapping. Use --skip-enrichment to skip GO.")
        return

    # Step 1 — transcriptional activity
    plot_transcriptional_activity(adata, args.output, args.sample)

    # Step 2 — marker genes
    marker_df = find_marker_genes(
        adata, args.output, args.sample,
        method=args.method,
        log2fc_min=args.log2fc_min,
        pct_min=args.pct_min,
        pct_ratio_min=args.pct_ratio,
    )

    # Step 3 — enrichment
    if not args.skip_enrichment:
        enrichment_analysis_per_cluster(
            adata, marker_df, fbgn_to_symbol,
            args.output, args.sample,
            top_n=args.top_n,
            combine_go=not args.no_combine_go,
        )

    print("\n" + "="*60)
    print("DONE")
    print("="*60)
    print(f"  Results: {args.output}/")


if __name__ == '__main__':
    main()
