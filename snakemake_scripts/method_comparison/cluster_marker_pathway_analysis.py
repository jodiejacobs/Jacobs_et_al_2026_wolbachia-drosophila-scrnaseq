'''
Comprehensive cluster analysis: transcriptional activity, marker genes, and pathway enrichment
Uses FlyEnrichr API for automated pathway analysis with FlyBase annotations

Key features:
  - DE run on adata.X (log-normalised counts), NOT scaled or raw
  - Background gene set passed to FlyEnrichr (all detected genes in dataset)
  - Top 50 upregulated markers (log2fc > 0) per cluster passed to FlyEnrichr
  - Bacterial/TE genes excluded from DE and enrichment input
  - Dot plot + network plot visualisation for enrichment results
  - Scanpy marker gene plots use FBgn IDs directly (symbol conversion for FlyEnrichr only)
  - Workaround for scanpy 1.10.x rankby_abs bug: re-ranks by signed score after extraction
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
                   show=False,
                   cmap='viridis', ncols=3)

    return dict(h_counts=h_counts, p_counts=p_counts,
                h_genes=h_genes,   p_genes=p_genes,
                eta_counts=eta_counts, eta_genes=eta_genes)


# ─────────────────────────────────────────────────────────────────────────────
# Marker genes
# ─────────────────────────────────────────────────────────────────────────────

def find_marker_genes(adata, output_dir, sample_name, method='wilcoxon'):
    """
    Find cluster marker genes.

    DE is run on adata.X (log-normalised, NOT scaled). use_raw=False.

    Gene filters applied before DE:
      - Bacterial genes (var_names starting with 'G')
      - FBtr transcripts (var_names starting with 'FBtr')
      - Transposable elements (FBti*, *_transposable_element)

    Workaround for scanpy 1.10.x rankby_abs bug:
      After extraction, genes are re-ranked by signed score so that
      upregulated genes sort to the top per cluster.
    """
    print("\n" + "="*60)
    print("DIFFERENTIAL GENE EXPRESSION")
    print("="*60)
    print(f"  Method: {method}")
    print(f"  Using adata.X (log-normalised, use_raw=False)")
    print(f"  Scanpy version: {sc.__version__} — applying signed score re-rank workaround")

    sc.settings.figdir = output_dir
    os.makedirs(output_dir, exist_ok=True)

    # ── Gene filtering ────────────────────────────────────────────────────────
    var_series = pd.Series(adata.var_names)

    bact_mask = var_series.str.startswith('G').values
    fbtr_mask = var_series.str.startswith('FBtr').values
    te_mask   = _is_te(var_series).values
    exclude   = bact_mask | fbtr_mask | te_mask

    adata_de = adata[:, ~exclude].copy()
    print(f"  Removed {bact_mask.sum():,} bacterial genes (G*), "
          f"{fbtr_mask.sum():,} FBtr transcripts, "
          f"{te_mask.sum():,} transposable elements")
    print(f"  Genes remaining for DE: {adata_de.n_vars:,}")

    # Filter to genes expressed in >= 3 cells
    X = adata_de.X
    if issparse(X):
        n_cells_per_gene = np.array((X > 0).sum(axis=0)).flatten()
    else:
        n_cells_per_gene = (X > 0).sum(axis=0)

    n_before = adata_de.n_vars
    adata_de = adata_de[:, n_cells_per_gene >= 3].copy()
    print(f"  Genes after expression filter (>=3 cells): "
          f"{adata_de.n_vars:,} / {n_before:,}")

    # ── Run DE ────────────────────────────────────────────────────────────────
    sc.settings.n_jobs = -1
    sc.tl.rank_genes_groups(
        adata_de,
        groupby='leiden',
        method=method,
        use_raw=False,
        key_added='rank_genes_groups',
        tie_correct=True,
        pts=True,
    )

    # Store results back on original adata
    adata.uns['rank_genes_groups'] = adata_de.uns['rank_genes_groups']

    # ── Scanpy marker plots ───────────────────────────────────────────────────
    print("\n  Saving scanpy marker gene plots ...")
    print(f"  Plots will be saved to: {output_dir}/")

    sc.pl.rank_genes_groups(
        adata_de,
        n_genes=25,
        save=f'_{sample_name}_ranked_genes.pdf',
        key='rank_genes_groups',
        show=False,
    )
    print(f"  Saved: ranked_genes_{sample_name}.pdf")

    sc.pl.rank_genes_groups_heatmap(
        adata_de,
        n_genes=10,
        save=f'_{sample_name}_top10_heatmap.pdf',
        show_gene_labels=True,
        cmap='viridis',
        key='rank_genes_groups',
        show=False,
    )
    print(f"  Saved: top10_heatmap_{sample_name}.pdf")

    sc.pl.rank_genes_groups_dotplot(
        adata_de,
        n_genes=5,
        save=f'_{sample_name}_top5_dotplot.pdf',
        key='rank_genes_groups',
        groupby='leiden',
        color_map='viridis',
        show=False,
    )
    print(f"  Saved: top5_dotplot_{sample_name}.pdf")

    sc.pl.rank_genes_groups_matrixplot(
        adata_de,
        n_genes=5,
        save=f'_{sample_name}_top5_matrixplot.pdf',
        key='rank_genes_groups',
        groupby='leiden',
        cmap='viridis',
        show=False,
    )
    print(f"  Saved: top5_matrixplot_{sample_name}.pdf")

    # ── Extract full DE results ───────────────────────────────────────────────
    result = adata_de.uns['rank_genes_groups']
    groups = result['names'].dtype.names

    rows = []
    for grp in groups:
        for i in range(len(result['names'][grp])):
            lfc  = result['logfoldchanges'][grp][i]
            padj = result['pvals_adj'][grp][i]
            if pd.isna(lfc) or np.isinf(lfc):
                continue
            pct_in   = result['pts'][grp][i]     if 'pts'      in result else np.nan
            pct_rest = result['pts_rest'][grp][i] if 'pts_rest' in result else np.nan
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

    # ── Workaround for scanpy 1.10.x rankby_abs bug ───────────────────────────
    # Re-rank by signed score so upregulated genes sort to top per cluster
    marker_df['signed_score'] = np.where(
        marker_df['log2fc'] > 0,
         marker_df['score'],
        -marker_df['score'],
    )
    marker_df = marker_df.sort_values(
        ['cluster', 'signed_score'], ascending=[True, False]
    ).drop(columns='signed_score')

    print(f"\n  Raw DE results: {len(marker_df):,} gene×cluster entries")

    # Print top 5 upregulated per cluster
    print("\n  Top 5 upregulated markers per cluster (by signed score):")
    for cl in sorted(marker_df['cluster'].unique(),
                     key=lambda x: int(x) if str(x).isdigit() else x):
        sub = marker_df[
            (marker_df['cluster'] == cl) &
            (marker_df['log2fc']   > 0)
        ].head(5)
        print(f"\n  Cluster {cl}:")
        if len(sub) == 0:
            print("    (no upregulated markers found)")
            continue
        for _, row in sub.iterrows():
            print(f"    {row['gene']:<22} log2FC={row['log2fc']:>6.2f}  "
                  f"pct_in={row['pct_in']:.2f}  pval_adj={row['pval_adj']:.2e}")

    # Save results
    marker_df.to_csv(
        os.path.join(output_dir, f'{sample_name}_markers_all.csv'), index=False)
    marker_df[marker_df['log2fc'] > 0].groupby('cluster').head(50).to_csv(
        os.path.join(output_dir, f'{sample_name}_markers_top50.csv'), index=False)
    print(f"\n  Saved: {sample_name}_markers_all.csv")
    print(f"  Saved: {sample_name}_markers_top50.csv")

    return marker_df


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

def build_background(adata, fbgn_to_symbol, min_cells=3):
    """
    Return gene symbols for all genes detected in >= min_cells cells.
    Excludes bacterial (G*), FBtr transcripts, and transposable elements.
    mt/ribo/cell_cycle genes are kept so enrichment p-values are correctly
    calibrated against the full expressed transcriptome.
    """
    var_series = pd.Series(adata.var_names)
    X = adata.X
    if issparse(X):
        n_cells_per_gene = np.array((X > 0).sum(axis=0)).flatten()
    else:
        n_cells_per_gene = (X > 0).sum(axis=0)

    expressed_mask = n_cells_per_gene >= min_cells
    bact_mask = var_series.str.startswith('G').values
    fbtr_mask = var_series.str.startswith('FBtr').values
    te_mask   = _is_te(var_series).values

    keep = expressed_mask & ~bact_mask & ~fbtr_mask & ~te_mask
    background_fbgn = adata.var_names[keep].tolist()
    symbols, n_unmapped = symbols_from_fbgn(background_fbgn, fbgn_to_symbol)
    print(f"  Background: {len(background_fbgn):,} FBgn IDs → {len(symbols):,} symbols "
          f"({n_unmapped:,} unmapped, {int(bact_mask.sum()):,} bacterial excluded, "
          f"{int(fbtr_mask.sum()):,} FBtr excluded, {int(te_mask.sum()):,} TEs excluded)")
    return symbols


# ─────────────────────────────────────────────────────────────────────────────
# Enrichment network plot
# ─────────────────────────────────────────────────────────────────────────────

def plot_enrichment_network(sig_df, output_dir, sample_name,
                            jaccard_thresh=0.3, top_n_per_cluster=15,
                            min_adj_p=0.1):
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

    top_terms = (go_sig.sort_values('adj_p_value')
                       .groupby('cluster')
                       .head(top_n_per_cluster)['term'].unique())
    go_plot = go_sig[go_sig['term'].isin(top_terms)].copy()

    term_best = (go_plot.sort_values('adj_p_value')
                        .drop_duplicates('term')
                        [['term', 'term_short', 'neg_log10p', 'cluster', 'gene_set']]
                        .set_index('term'))

    terms = term_best.index.tolist()
    print(f"  Terms in network: {len(terms)}")
    if len(terms) < 3:
        print("  Too few terms — skipping network")
        return

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

    connected = [n for n in G.nodes if G.degree(n) > 0]
    if len(connected) >= 5:
        G.remove_nodes_from([n for n in list(G.nodes) if G.degree(n) == 0])
        print(f"  Nodes after removing isolates: {G.number_of_nodes()}")

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

    nx.write_graphml(G, f"{output_dir}/{sample_name}_GO_enrichment_network.graphml")
    print(f"    Saved: GO_enrichment_network.graphml")


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
    background_symbols = build_background(adata, fbgn_to_symbol)
    if len(background_symbols) < 100:
        print("  WARNING: Very small background — check mapping file")

    libraries = [
        'GO_Biological_Process_2018',
        'WikiPathways_2018',
    ]

    clusters = sorted(marker_df['cluster'].unique(),
                      key=lambda x: int(x) if str(x).isdigit() else x)
    all_results = []

    for cluster in clusters:
        print(f"\n{'='*50}\n  Cluster {cluster}\n{'='*50}")

        # Top 50 upregulated markers per cluster
        cmarkers = (
            marker_df[
                (marker_df['cluster'] == cluster) &
                (marker_df['log2fc']   > 0)
            ]
            .head(50)
        )

        print(f"  Top 50 upregulated markers (log2fc>0): {len(cmarkers)}")

        if len(cmarkers) < 5:
            print(f"  Skipping: too few markers (need ≥ 5)")
            continue

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
    parser.add_argument('--skip-enrichment', action='store_true')
    parser.add_argument('--no-combine-go',   action='store_true')

    args = parser.parse_args()
    os.makedirs(args.output, exist_ok=True)

    sc.settings.figdir = args.output
    print(f"  sc.settings.figdir = {sc.settings.figdir}")

    print("="*60)
    print("LOADING DATA")
    print("="*60)
    adata = sc.read_h5ad(args.input)
    print(f"  Cells: {adata.n_obs:,}  Genes: {adata.n_vars:,}  "
          f"Clusters: {adata.obs['leiden'].nunique()}")
    print(f"  adata.X will be used for DE (use_raw=False)")

    fbgn_to_symbol = load_fbgn_to_symbol_mapping(args.mapping)
    if fbgn_to_symbol is None and not args.skip_enrichment:
        print("\nERROR: Could not load gene mapping. Use --skip-enrichment to skip GO.")
        return

    plot_transcriptional_activity(adata, args.output, args.sample)

    marker_df = find_marker_genes(
        adata, args.output, args.sample,
        method=args.method,
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