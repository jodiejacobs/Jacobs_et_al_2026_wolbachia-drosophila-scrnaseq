'''
Cell cycle annotation using Scanpy with Drosophila-specific markers
and cluster-cell cycle association analysis

# Reference (leiden column, use --ctrl-only if you want)
python scanpy_cellcycle_analysis.py --input integrated_reference.h5ad \
    --output results/cellcycle/ref --sample ref --run-scoring --save-output

# Query (leiden_ref column, do NOT use --ctrl-only)
python scanpy_cellcycle_analysis.py --input integrated_query.h5ad \
    --output results/cellcycle/query --sample query --run-scoring

# Combined, ctrl only
python scanpy_cellcycle_analysis.py --input integrated_combined.h5ad \
    --output results/cellcycle/combined_ctrl --sample combined_ctrl \
    --ctrl-only --run-scoring

# Infection cell cycle analysis (JW18DOX vs rest)
python scanpy_cellcycle_analysis.py --input integrated_combined.h5ad \
    --output results/cellcycle/infection --sample infection \
    --run-scoring --infection-analysis

# Custom infected label
python scanpy_cellcycle_analysis.py --input integrated_combined.h5ad \
    --output results/cellcycle/infection --sample infection \
    --run-scoring --infection-analysis --infected-label MyInfectedGroup

'''
import scanpy as sc 
import argparse
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import chi2_contingency, kruskal

# Drosophila cell cycle genes (FlyBase IDs)
# Mapping from gene symbols to FlyBase IDs
FLYBASE_CELL_CYCLE_GENES = {
    # S phase genes
    'Pcna': 'FBgn0005655',
    'RPA1': 'FBgn0015806', 
    'RPA2': 'FBgn0034898',
    'pol-alpha1': 'FBgn0011230',
    'DNApol-alpha60': 'FBgn0015278',
    'DNApol-delta': 'FBgn0019624',
    'RnrL': 'FBgn0020369',
    'RnrS': 'FBgn0261933',
    'Mcm2': 'FBgn0020651',
    'Mcm3': 'FBgn0020652',
    'Mcm5': 'FBgn0015929',
    'Mcm6': 'FBgn0032435',
    'Mcm7': 'FBgn0015308',
    'E2f1': 'FBgn0011766',
    'E2f2': 'FBgn0262656',
    'CycE': 'FBgn0010382',
    'Cdk2': 'FBgn0010314',
    'Dp': 'FBgn0000499',
    'Rbf': 'FBgn0015799',
    'Rbf2': 'FBgn0028396',
    'Orc1': 'FBgn0015270',
    'Orc2': 'FBgn0015714',
    'Orc6': 'FBgn0025926',
    'Rrp1': 'FBgn0003257',
    # G2/M genes
    'CycA': 'FBgn0010114',
    'CycB': 'FBgn0010113',
    'CycB3': 'FBgn0011577',
    'Cdk1': 'FBgn0004107',
    'stg': 'FBgn0003525',
    'polo': 'FBgn0003124',
    'aurA': 'FBgn0025564',
    'aurB': 'FBgn0025948',
    'Nek2': 'FBgn0027548',
    'Pbl': 'FBgn0005619',
    'Wee1': 'FBgn0011739',
    'myt': 'FBgn0002863',
    'BubR1': 'FBgn0024822',
    'Mad2': 'FBgn0002610',
    'Cdc20': 'FBgn0010309',
    'APC2': 'FBgn0261823',
    'APC10': 'FBgn0036449',
}

# Define S and G2M gene sets by FlyBase ID
S_GENES_FBGN = [
    'FBgn0005655',  # Pcna
    'FBgn0015806',  # RPA1
    'FBgn0034898',  # RPA2
    'FBgn0011230',  # pol-alpha1
    'FBgn0015278',  # DNApol-alpha60
    'FBgn0019624',  # DNApol-delta
    'FBgn0020369',  # RnrL
    'FBgn0261933',  # RnrS
    'FBgn0020651',  # Mcm2
    'FBgn0020652',  # Mcm3
    'FBgn0015929',  # Mcm5
    'FBgn0032435',  # Mcm6
    'FBgn0015308',  # Mcm7
    'FBgn0011766',  # E2f1
    'FBgn0262656',  # E2f2
    'FBgn0010382',  # CycE
    'FBgn0010314',  # Cdk2
    'FBgn0000499',  # Dp
    'FBgn0015799',  # Rbf
    'FBgn0028396',  # Rbf2
    'FBgn0015270',  # Orc1
    'FBgn0015714',  # Orc2
    'FBgn0025926',  # Orc6
    'FBgn0003257',  # Rrp1
]

G2M_GENES_FBGN = [
    'FBgn0010114',  # CycA
    'FBgn0010113',  # CycB
    'FBgn0011577',  # CycB3
    'FBgn0004107',  # Cdk1
    'FBgn0003525',  # stg
    'FBgn0003124',  # polo
    'FBgn0025564',  # aurA
    'FBgn0025948',  # aurB
    'FBgn0027548',  # Nek2
    'FBgn0005619',  # Pbl
    'FBgn0011739',  # Wee1
    'FBgn0002863',  # myt
    'FBgn0024822',  # BubR1
    'FBgn0002610',  # Mad2
    'FBgn0010309',  # Cdc20
    'FBgn0261823',  # APC2
    'FBgn0036449',  # APC10
]

# Create reverse mapping for reporting
FBGN_TO_SYMBOL = {v: k for k, v in FLYBASE_CELL_CYCLE_GENES.items()}


def _get_leiden_col(adata):
    """
    Return the leiden column name present in adata.obs.
    Prefers 'leiden_ref' (query/combined from integrate.py) over 'leiden' (reference).
    Raises a clear error if neither is found.
    """
    for col in ('leiden_ref', 'leiden'):
        if col in adata.obs.columns:
            return col
    raise KeyError(
        "No leiden column found in adata.obs. "
        "Expected 'leiden_ref' (query/combined) or 'leiden' (reference). "
        f"Available columns: {list(adata.obs.columns)}"
    )


def _ensure_lognorm(adata):
    """
    sc.tl.score_genes_cell_cycle expects log-normalised counts.
    If the matrix looks like raw counts (max > 20 is a rough heuristic),
    normalise and log-transform in-place.
    """
    import scipy.sparse
    X = adata.X
    mat_max = X.max() if not scipy.sparse.issparse(X) else X.max()
    if hasattr(mat_max, 'A1'):      # sparse scalar
        mat_max = float(mat_max.A1[0])
    else:
        mat_max = float(mat_max)

    if mat_max > 20:
        print("WARNING: Data appears to be raw counts (max value = "
              f"{mat_max:.1f}). Running normalize_total + log1p before scoring.")
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
    else:
        print(f"Data appears log-normalised (max value = {mat_max:.3f}). Skipping normalisation.")


def check_gene_names(adata):
    """Check what format the gene names are in."""
    print("\nChecking gene naming format...")
    sample_genes = list(adata.var_names[:10])
    print(f"Sample gene names: {sample_genes}")

    fbgn_count   = sum(1 for g in adata.var_names if str(g).startswith('FBgn'))
    symbol_count = sum(1 for g in adata.var_names if not str(g).startswith('FBgn'))

    print(f"\nGenes starting with 'FBgn': {fbgn_count}")
    print(f"Genes not starting with 'FBgn': {symbol_count}")

    if fbgn_count > symbol_count:
        print("-> Detected FlyBase ID format")
        return 'flybase'
    else:
        print("-> Detected gene symbol format")
        return 'symbol'


def score_cell_cycle_scanpy(adata, output_dir, sample_name, s_genes=None, g2m_genes=None):
    """Score cell cycle using Scanpy with Drosophila genes."""
    print("\n" + "=" * 60)
    print("SCANPY CELL CYCLE SCORING (DROSOPHILA)")
    print("=" * 60)

    check_gene_names(adata)

    # Fix 1: ensure data is log-normalised before scoring
    _ensure_lognorm(adata)

    if s_genes is None:
        s_genes = S_GENES_FBGN
    if g2m_genes is None:
        g2m_genes = G2M_GENES_FBGN

    s_genes_present   = [g for g in s_genes   if g in adata.var_names]
    g2m_genes_present = [g for g in g2m_genes if g in adata.var_names]

    print(f"\nS phase genes: {len(s_genes_present)}/{len(s_genes)} found")
    if s_genes_present:
        print("  Present (showing symbols): " +
              ', '.join([FBGN_TO_SYMBOL.get(g, g) for g in s_genes_present[:10]]) +
              (f"... (+{len(s_genes_present)-10} more)" if len(s_genes_present) > 10 else ""))
    else:
        print("  None found")

    print(f"\nG2/M phase genes: {len(g2m_genes_present)}/{len(g2m_genes)} found")
    if g2m_genes_present:
        print("  Present (showing symbols): " +
              ', '.join([FBGN_TO_SYMBOL.get(g, g) for g in g2m_genes_present[:10]]) +
              (f"... (+{len(g2m_genes_present)-10} more)" if len(g2m_genes_present) > 10 else ""))
    else:
        print("  None found")

    if len(s_genes_present) < 3 or len(g2m_genes_present) < 3:
        print("\nWARNING: Very few marker genes found. Results may not be reliable.")
        missing_s   = [FBGN_TO_SYMBOL.get(g, g) for g in s_genes   if g not in adata.var_names]
        missing_g2m = [FBGN_TO_SYMBOL.get(g, g) for g in g2m_genes if g not in adata.var_names]
        print(f"Missing S genes:   {', '.join(missing_s[:20])}")
        print(f"Missing G2/M genes: {', '.join(missing_g2m[:20])}")
        if len(s_genes_present) == 0 and len(g2m_genes_present) == 0:
            print("\nERROR: No cell cycle genes found! Check gene annotations.")
            return None

    print("\nScoring cell cycle phases...")
    sc.tl.score_genes_cell_cycle(adata, s_genes=s_genes_present, g2m_genes=g2m_genes_present)

    print("\nCell cycle phase distribution:")
    phase_counts = adata.obs['phase'].value_counts()
    for phase in ['G1', 'S', 'G2M']:
        count = phase_counts.get(phase, 0)
        print(f"  {phase}: {count} cells ({count / adata.n_obs * 100:.1f}%)")

    create_scanpy_plots(adata, output_dir, sample_name, s_genes_present, g2m_genes_present)
    return adata


def create_scanpy_plots(adata, output_dir, sample_name, s_genes, g2m_genes):
    """Create diagnostic plots for Scanpy cell cycle scoring."""
    os.makedirs(output_dir, exist_ok=True)

    colors = {'G1': '#FF6B6B', 'S': '#4ECDC4', 'G2M': '#45B7D1'}

    # 1. Phase distribution overview
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))

    phase_counts = adata.obs['phase'].value_counts()
    phase_counts.plot(kind='bar', ax=axes[0],
                      color=[colors.get(p, 'gray') for p in phase_counts.index])
    axes[0].set_xlabel('Cell Cycle Phase')
    axes[0].set_ylabel('Number of Cells')
    axes[0].set_title('Cell Cycle Phase Distribution')
    axes[0].set_xticklabels(axes[0].get_xticklabels(), rotation=45)

    for phase, color in colors.items():
        mask = adata.obs['phase'] == phase
        axes[1].scatter(adata.obs.loc[mask, 'S_score'],
                        adata.obs.loc[mask, 'G2M_score'],
                        c=color, label=phase, alpha=0.5, s=10)
    axes[1].set_xlabel('S Score')
    axes[1].set_ylabel('G2M Score')
    axes[1].set_title('Cell Cycle Scores')
    axes[1].legend()
    axes[1].axhline(y=0, color='k', linestyle='--', alpha=0.3)
    axes[1].axvline(x=0, color='k', linestyle='--', alpha=0.3)

    axes[2].hist(adata.obs['S_score'],   bins=50, alpha=0.5, label='S score',   color=colors['S'])
    axes[2].hist(adata.obs['G2M_score'], bins=50, alpha=0.5, label='G2M score', color=colors['G2M'])
    axes[2].set_xlabel('Score')
    axes[2].set_ylabel('Number of Cells')
    axes[2].set_title('Score Distributions')
    axes[2].legend()
    axes[2].axvline(x=0, color='k', linestyle='--', alpha=0.3)

    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, f'{sample_name}_cellcycle_overview.pdf'),
                dpi=300, bbox_inches='tight')
    plt.close()

    # 2. Marker gene expression heatmap
    import scipy.sparse
    all_markers = list(set(s_genes + g2m_genes))
    if all_markers:
        adata_sorted = adata[adata.obs.sort_values('phase').index, :]
        X = adata_sorted[:, all_markers].X
        if scipy.sparse.issparse(X):
            X = X.toarray()
        marker_expr = pd.DataFrame(X, index=adata_sorted.obs_names, columns=all_markers).T
        marker_expr_norm = marker_expr.apply(lambda x: (x - x.mean()) / (x.std() + 1e-10), axis=1)
        gene_labels = [FBGN_TO_SYMBOL.get(g, g) for g in all_markers]

        fig, ax = plt.subplots(figsize=(12, max(8, len(all_markers) * 0.3)))
        sns.heatmap(marker_expr_norm, cmap='RdBu_r', center=0,
                    cbar_kws={'label': 'Z-score'},
                    xticklabels=False, yticklabels=gene_labels,
                    vmin=-2, vmax=2, ax=ax)
        ax.set_xlabel('Cells (sorted by phase)')
        ax.set_ylabel('Cell Cycle Marker Genes')
        ax.set_title(f'Cell cycle marker expression ({len(all_markers)} genes)\n{sample_name}')
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, f'{sample_name}_marker_heatmap.pdf'),
                    dpi=300, bbox_inches='tight')
        plt.close()

    print(f"\nCell cycle plots saved to {output_dir}")


def analyze_cluster_cellcycle_association(adata, fig_dir, sample):
    """Test and visualize association between Leiden clusters and cell cycle."""
    print("\n" + "=" * 60)
    print("CLUSTER - CELL CYCLE ASSOCIATION ANALYSIS")
    print("=" * 60)

    # Fix 1: handle both 'leiden' (reference) and 'leiden_ref' (query/combined)
    try:
        leiden_col = _get_leiden_col(adata)
    except KeyError as e:
        print(f"ERROR: {e}")
        return None

    print(f"Using leiden column: '{leiden_col}'")

    if 'phase' not in adata.obs.columns:
        print("ERROR: No 'phase' found in adata.obs. Run with --run-scoring first.")
        return None

    os.makedirs(fig_dir, exist_ok=True)
    # Set figdir before any sc.pl calls
    sc.settings.figdir = fig_dir

    clusters = sorted(adata.obs[leiden_col].unique())
    cmap = plt.cm.get_cmap('tab20')
    leiden_colors = [cmap(i % 20) for i in range(len(clusters))]

    # ── 1. Chi-square test ────────────────────────────────────────────────────
    print("\n" + "=" * 60)
    print("1. CHI-SQUARE TEST: Cluster vs Cell Cycle Stage")
    print("=" * 60)

    contingency = pd.crosstab(adata.obs[leiden_col], adata.obs['phase'])
    chi2, p_value, dof, expected = chi2_contingency(contingency)

    print(f"chi2 = {chi2:.2f}")
    print(f"degrees of freedom = {dof}")
    print(f"p-value = {p_value:.2e}")
    print(f"\nConclusion: Clusters are "
          f"{'SIGNIFICANTLY' if p_value < 0.05 else 'NOT significantly'} "
          f"associated with cell cycle stage")

    n = contingency.sum().sum()
    cramers_v = np.sqrt(chi2 / (n * (min(contingency.shape) - 1)))
    print(f"Cramer's V = {cramers_v:.3f}")
    if cramers_v < 0.1:       effect = "negligible"
    elif cramers_v < 0.3:     effect = "weak"
    elif cramers_v < 0.5:     effect = "moderate"
    else:                     effect = "strong"
    print(f"Effect size: {effect}")

    # ── 2. Heatmap ────────────────────────────────────────────────────────────
    contingency_norm = contingency.div(contingency.sum(axis=1), axis=0) * 100
    print("\nCell cycle stage distribution by cluster (%):")
    print(contingency_norm.round(1))

    fig, ax = plt.subplots(figsize=(12, 8))
    sns.heatmap(contingency_norm, annot=True, fmt='.1f', cmap='YlOrRd', ax=ax,
                cbar_kws={'label': '% of cells in cluster'})
    ax.set_xlabel('Cell Cycle Phase')
    ax.set_ylabel(f'Leiden Cluster ({leiden_col})')
    ax.set_title(f"Cell cycle phase distribution by cluster\n"
                 f"chi2={chi2:.2f}, p={p_value:.2e}, Cramer's V={cramers_v:.3f}")
    plt.tight_layout()
    plt.savefig(os.path.join(fig_dir, f'heatmap_cluster_cellcycle_{sample}.pdf'))
    plt.close()

    # ── 3. Stacked bar ────────────────────────────────────────────────────────
    fig, ax = plt.subplots(figsize=(14, 7))
    contingency_norm.plot(kind='bar', stacked=True, ax=ax, width=0.8)
    ax.set_xlabel(f'Leiden Cluster ({leiden_col})', fontsize=12)
    ax.set_ylabel('Percentage of cells', fontsize=12)
    ax.set_title('Cell cycle phase composition by cluster', fontsize=14)
    ax.legend(title='Cell Cycle Phase', bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.xticks(rotation=0)
    plt.tight_layout()
    plt.savefig(os.path.join(fig_dir, f'barplot_cluster_cellcycle_{sample}.pdf'))
    plt.close()

    # ── 4. Enrichment summary ─────────────────────────────────────────────────
    print("\n" + "=" * 60)
    print("2. CELL CYCLE ENRICHED CLUSTERS")
    print("=" * 60)

    dominant_phase  = contingency_norm.idxmax(axis=1)
    max_percentage  = contingency_norm.max(axis=1)

    print(f"\n{'Cluster':<10} {'Dominant Phase':<15} {'Percentage':<12} Status")
    print("-" * 60)
    for cluster in clusters:
        phase = dominant_phase[cluster]
        pct   = max_percentage[cluster]
        enrichment = ("STRONGLY ENRICHED" if pct > 50
                      else "ENRICHED" if pct > 40
                      else "Mixed")
        print(f"{cluster:<10} {phase:<15} {pct:>6.1f}%      {enrichment}")

    # ── 5. Cell cycle scores by cluster ──────────────────────────────────────
    print("\n" + "=" * 60)
    print("3. CELL CYCLE SCORES BY CLUSTER")
    print("=" * 60)

    score_by_cluster = adata.obs.groupby(leiden_col)[['S_score', 'G2M_score']].agg(['mean', 'std'])
    print(score_by_cluster.round(3))

    groups_s   = [adata.obs[adata.obs[leiden_col] == c]['S_score'].dropna().values   for c in clusters]
    groups_g2m = [adata.obs[adata.obs[leiden_col] == c]['G2M_score'].dropna().values for c in clusters]
    h_stat_s,   p_value_s   = kruskal(*groups_s)
    h_stat_g2m, p_value_g2m = kruskal(*groups_g2m)

    print(f"\nKruskal-Wallis S_score:   H={h_stat_s:.2f},   p={p_value_s:.2e}")
    print(f"Kruskal-Wallis G2M_score: H={h_stat_g2m:.2f}, p={p_value_g2m:.2e}")

    # Violin plots
    fig, axes = plt.subplots(1, 2, figsize=(16, 6))
    sc.pl.violin(adata, 'S_score',   groupby=leiden_col, ax=axes[0], show=False, rotation=0)
    axes[0].set_title(f'S Score by Cluster\nKruskal-Wallis H={h_stat_s:.2f}, p={p_value_s:.2e}')
    sc.pl.violin(adata, 'G2M_score', groupby=leiden_col, ax=axes[1], show=False, rotation=0)
    axes[1].set_title(f'G2M Score by Cluster\nKruskal-Wallis H={h_stat_g2m:.2f}, p={p_value_g2m:.2e}')
    plt.tight_layout()
    plt.savefig(os.path.join(fig_dir, f'violin_scores_by_cluster_{sample}.pdf'))
    plt.close()

    # Bar plots for mean scores
    fig, axes = plt.subplots(1, 2, figsize=(16, 6))
    for ax, score_key, title in zip(
        axes,
        ['S_score', 'G2M_score'],
        ['Mean S Score by Cluster', 'Mean G2M Score by Cluster'],
    ):
        means = score_by_cluster[score_key]['mean']
        stds  = score_by_cluster[score_key]['std']
        ax.bar(range(len(means)), means, yerr=stds,
               color=leiden_colors, alpha=0.7, capsize=5)
        ax.set_xlabel(f'Leiden Cluster ({leiden_col})')
        ax.set_ylabel(f'Mean {score_key}')
        ax.set_title(title)
        ax.set_xticks(range(len(means)))
        ax.set_xticklabels(means.index)
        ax.axhline(y=0, color='k', linestyle='--', alpha=0.3)
    plt.tight_layout()
    plt.savefig(os.path.join(fig_dir, f'barplot_scores_by_cluster_{sample}.pdf'))
    plt.close()

    # ── 6. UMAPs ──────────────────────────────────────────────────────────────
    if 'X_umap' in adata.obsm:
        fig, axes = plt.subplots(1, 2, figsize=(16, 6))
        sc.pl.umap(adata, color=leiden_col, ax=axes[0], show=False,
                   title=f'Leiden Clusters ({leiden_col})', frameon=False, legend_loc='on data')
        sc.pl.umap(adata, color='phase', ax=axes[1], show=False,
                   title='Cell Cycle Phase', frameon=False)
        plt.tight_layout()
        plt.savefig(os.path.join(fig_dir, f'umap_cluster_vs_cellcycle_{sample}.pdf'))
        plt.close()

        sc.pl.umap(adata, color=leiden_col,
                   save=f'_{sample}_leiden.pdf',
                   title=f'Leiden Clusters ({leiden_col})', legend_loc='on data')
        sc.pl.umap(adata, color='phase',
                   save=f'_{sample}_cellcycle_phase.pdf',
                   title='Cell Cycle Phase')
        sc.pl.umap(adata, color=['S_score', 'G2M_score'],
                   save=f'_{sample}_cellcycle_scores.pdf',
                   cmap='viridis')

    # ── 7. Summary table ──────────────────────────────────────────────────────
    print("\n" + "=" * 60)
    print("4. COMPREHENSIVE SUMMARY TABLE")
    print("=" * 60)

    summary_df = pd.DataFrame({
        'Cluster':          clusters,
        'N_Cells':          [contingency.loc[c].sum()                             for c in clusters],
        'Dominant_Phase':   [dominant_phase[c]                                    for c in clusters],
        'Phase_Percentage': [max_percentage[c]                                    for c in clusters],
        'Mean_S_Score':     [score_by_cluster.loc[c, ('S_score',   'mean')]       for c in clusters],
        'Mean_G2M_Score':   [score_by_cluster.loc[c, ('G2M_score', 'mean')]       for c in clusters],
    })
    print(summary_df.to_string(index=False))

    summary_df.to_csv(os.path.join(fig_dir, f'cluster_cellcycle_summary_{sample}.csv'),   index=False)
    contingency.to_csv(os.path.join(fig_dir,      f'contingency_table_counts_{sample}.csv'))
    contingency_norm.to_csv(os.path.join(fig_dir, f'contingency_table_percentages_{sample}.csv'))

    stats_df = pd.DataFrame({
        'Test':       ['Chi-square', "Cramer's V",
                       'Kruskal-Wallis (S_score)', 'Kruskal-Wallis (G2M_score)'],
        'Statistic':  [f'{chi2:.6e}',     f'{cramers_v:.6e}',
                       f'{h_stat_s:.6e}', f'{h_stat_g2m:.6e}'],
        'P-value':    [f'{p_value:.6e}',  'NA',
                       f'{p_value_s:.6e}', f'{p_value_g2m:.6e}'],
        'Interpretation': [
            'Significant association'     if p_value   < 0.05 else 'No significant association',
            f'{effect.capitalize()} effect size',
            'S scores differ'             if p_value_s < 0.05 else 'No difference in S scores',
            'G2M scores differ'           if p_value_g2m < 0.05 else 'No difference in G2M scores',
        ],
    })
    stats_df.to_csv(os.path.join(fig_dir, f'statistical_tests_{sample}.csv'), index=False)

    print("\n" + "=" * 60)
    print("ANALYSIS COMPLETE")
    print("=" * 60)
    print(f"Output directory: {fig_dir}")

    return {
        'chi2':           chi2,
        'chi2_pvalue':    p_value,
        'cramers_v':      cramers_v,
        'kw_s_pvalue':    p_value_s,
        'kw_g2m_pvalue':  p_value_g2m,
        'summary':        summary_df,
    }


def analyze_infection_cellcycle(adata, fig_dir, sample,
                                 infection_col='cell_line',
                                 infected_label='JW18DOX'):
    """
    Test whether Wolbachia strain/infection (infected_label) shifts cell cycle
    phase proportions relative to all other groups in infection_col.

    Key change from v2: infection_col is now a parameter (default: 'cell_line')
    so comparisons like JW18DOX vs JW18wMel work out of the box, as do
    comparisons using 'treatment', 'bio_condition', or any other obs column.

    Statistical approach
    --------------------
    1. Overall chi-square across all groups x phases.
    2. Pairwise chi-square: infected_label vs each other group, Bonferroni-corrected.
    3. Mann-Whitney U on S_score and G2M_score: infected_label vs each other group,
       Bonferroni-corrected.

    Outputs
    -------
    Figures
      infection_cellcycle_bars_{sample}.pdf    – stacked bar + Δ% panel
      infection_cellcycle_heatmap_{sample}.pdf – phase-% heatmap (all groups)
      infection_scores_violin_{sample}.pdf     – S/G2M score violins by group
    CSVs
      infection_contingency_counts_{sample}.csv
      infection_contingency_pct_{sample}.csv
      infection_pairwise_chisq_{sample}.csv
      infection_pairwise_scores_{sample}.csv
      infection_overall_stats_{sample}.csv
    """
    from scipy.stats import chi2_contingency, mannwhitneyu

    print("\n" + "=" * 60)
    print("INFECTION vs CELL CYCLE ANALYSIS")
    print(f"  Column        : '{infection_col}'")
    print(f"  Infected label: '{infected_label}'")
    print("=" * 60)

    # ── Validation ────────────────────────────────────────────────────────────
    if 'phase' not in adata.obs.columns:
        print("ERROR: No 'phase' column found. Run --run-scoring first.")
        return None

    if infection_col not in adata.obs.columns:
        print(f"ERROR: '{infection_col}' not found in adata.obs.")
        print(f"Available columns: {list(adata.obs.columns)}")
        return None

    os.makedirs(fig_dir, exist_ok=True)

    groups = sorted(adata.obs[infection_col].astype(str).unique())
    print(f"\nGroups in '{infection_col}': {groups}")

    if infected_label not in groups:
        print(f"\nERROR: '{infected_label}' not found in '{infection_col}'.")
        print(f"Available values: {groups}")
        print(f"Hint: use --infection-col and --infected-label to match your metadata.")
        return None

    other_groups  = [g for g in groups if g != infected_label]
    n_comparisons = len(other_groups)
    bonf_alpha    = 0.05 / n_comparisons if n_comparisons > 0 else 0.05

    # ── Contingency table: group x phase ──────────────────────────────────────
    contingency = pd.crosstab(
        adata.obs[infection_col].astype(str), adata.obs['phase']
    )
    for ph in ['G1', 'S', 'G2M']:
        if ph not in contingency.columns:
            contingency[ph] = 0
    contingency     = contingency[['G1', 'S', 'G2M']]
    contingency_pct = contingency.div(contingency.sum(axis=1), axis=0) * 100

    print(f"\nCell cycle distribution per '{infection_col}' group (%):")
    print(contingency_pct.round(2))

    # ── 1. Overall chi-square ─────────────────────────────────────────────────
    chi2_all, p_all, dof_all, _ = chi2_contingency(contingency)
    n_all     = contingency.sum().sum()
    cramers_v = np.sqrt(chi2_all / (n_all * (min(contingency.shape) - 1)))
    if   cramers_v < 0.1: effect = "negligible"
    elif cramers_v < 0.3: effect = "weak"
    elif cramers_v < 0.5: effect = "moderate"
    else:                 effect = "strong"

    print(f"\nOverall chi-square (all groups × phase):")
    print(f"  chi2={chi2_all:.4f}, dof={dof_all}, p={p_all:.2e}")
    print(f"  Cramer's V={cramers_v:.3f} ({effect} effect)")

    # ── 2. Pairwise chi-square: infected vs each other ────────────────────────
    print(f"\nPairwise chi-square: '{infected_label}' vs each other group")
    print(f"Bonferroni threshold: p < {bonf_alpha:.4f}  ({n_comparisons} comparisons)")
    print(f"\n{'Comparison':<40} {'chi2':>8} {'p-value':>12} {'Bonf.sig':>9} {'Cramer V':>10}")
    print("-" * 85)

    pairwise_chisq = []
    for other in other_groups:
        sub = contingency.loc[[infected_label, other], :]
        sub = sub.loc[:, sub.sum(axis=0) > 0]
        chi2_pw, p_pw, _, _ = chi2_contingency(sub)
        n_pw  = sub.sum().sum()
        cv_pw = np.sqrt(chi2_pw / (n_pw * (min(sub.shape) - 1))) if min(sub.shape) > 1 else 0.0
        sig   = "YES" if p_pw < bonf_alpha else "no"
        label = f"{infected_label} vs {other}"
        print(f"{label:<40} {chi2_pw:>8.3f} {p_pw:>12.2e} {sig:>9} {cv_pw:>10.3f}")
        pairwise_chisq.append({
            'comparison':       label,
            'infected_label':   infected_label,
            'other_group':      other,
            'chi2':             chi2_pw,
            'p_value':          p_pw,
            'bonf_significant': sig == "YES",
            'cramers_v':        cv_pw,
        })
    pairwise_chisq_df = pd.DataFrame(pairwise_chisq)

    # ── 3. Mann-Whitney U on continuous scores ────────────────────────────────
    print(f"\nMann-Whitney U: '{infected_label}' vs each other group")
    print(f"(S_score and G2M_score, Bonferroni threshold p < {bonf_alpha:.4f})")
    print(f"\n{'Comparison':<40} {'S U-stat':>10} {'S p':>12} {'S sig':>7}  "
          f"{'G2M U-stat':>10} {'G2M p':>12} {'G2M sig':>7}")
    print("-" * 105)

    mask_infected = adata.obs[infection_col].astype(str) == infected_label
    infected_S    = adata.obs.loc[mask_infected, 'S_score'].dropna().values
    infected_G2M  = adata.obs.loc[mask_infected, 'G2M_score'].dropna().values

    pairwise_scores = []
    for other in other_groups:
        mask_other = adata.obs[infection_col].astype(str) == other
        other_S    = adata.obs.loc[mask_other, 'S_score'].dropna().values
        other_G2M  = adata.obs.loc[mask_other, 'G2M_score'].dropna().values

        u_s,   p_s   = mannwhitneyu(infected_S,   other_S,   alternative='two-sided')
        u_g2m, p_g2m = mannwhitneyu(infected_G2M, other_G2M, alternative='two-sided')
        sig_s   = "YES" if p_s   < bonf_alpha else "no"
        sig_g2m = "YES" if p_g2m < bonf_alpha else "no"
        label   = f"{infected_label} vs {other}"
        print(f"{label:<40} {u_s:>10.1f} {p_s:>12.2e} {sig_s:>7}  "
              f"{u_g2m:>10.1f} {p_g2m:>12.2e} {sig_g2m:>7}")
        pairwise_scores.append({
            'comparison':         label,
            'infected_label':     infected_label,
            'other_group':        other,
            'S_score_U':          u_s,
            'S_score_p':          p_s,
            'S_score_bonf_sig':   sig_s == "YES",
            'G2M_score_U':        u_g2m,
            'G2M_score_p':        p_g2m,
            'G2M_score_bonf_sig': sig_g2m == "YES",
        })
    pairwise_scores_df = pd.DataFrame(pairwise_scores)

    # ── 4. Stacked bar + Δ% panel ─────────────────────────────────────────────
    phase_colors = {'G1': '#FF6B6B', 'S': '#4ECDC4', 'G2M': '#45B7D1'}
    ordered      = [infected_label] + other_groups
    pct_ordered  = contingency_pct.loc[ordered]

    fig, axes = plt.subplots(1, 2, figsize=(16, 6))

    pct_ordered.plot(
        kind='bar', stacked=True, ax=axes[0], width=0.7,
        color=[phase_colors.get(p, 'gray') for p in pct_ordered.columns],
    )
    axes[0].set_xlabel(infection_col, fontsize=12)
    axes[0].set_ylabel('% of cells', fontsize=12)
    axes[0].set_title(f'Cell cycle phase composition\nby {infection_col}', fontsize=13)
    axes[0].legend(title='Phase', bbox_to_anchor=(1.01, 1), loc='upper left')
    axes[0].set_xticklabels(axes[0].get_xticklabels(), rotation=30, ha='right')
    for lbl in axes[0].get_xticklabels():
        if lbl.get_text() == infected_label:
            lbl.set_fontweight('bold')
            lbl.set_color('darkred')

    infected_pct = contingency_pct.loc[infected_label]
    others_mean  = contingency_pct.loc[other_groups].mean()
    delta        = infected_pct - others_mean
    bar_colors   = ['#2ecc71' if d >= 0 else '#e74c3c' for d in delta]
    axes[1].bar(delta.index, delta.values, color=bar_colors, edgecolor='black', linewidth=0.8)
    axes[1].axhline(0, color='black', linewidth=0.8, linestyle='--')
    axes[1].set_xlabel('Cell Cycle Phase', fontsize=12)
    axes[1].set_ylabel(f'Δ% ({infected_label} − mean of others)', fontsize=11)
    axes[1].set_title(f'Phase enrichment: {infected_label}\nvs mean of other groups', fontsize=13)
    for i, (ph, val) in enumerate(delta.items()):
        axes[1].text(
            i, val + (0.4 if val >= 0 else -0.9),
            f'{val:+.1f}%', ha='center', va='bottom', fontsize=11, fontweight='bold',
        )

    plt.suptitle(
        f"Wolbachia strain vs cell cycle  ({infection_col})  |  "
        f"Overall chi2={chi2_all:.2f}, p={p_all:.2e}, Cramer's V={cramers_v:.3f}",
        fontsize=11, y=1.01,
    )
    plt.tight_layout()
    plt.savefig(os.path.join(fig_dir, f'infection_cellcycle_bars_{sample}.pdf'),
                bbox_inches='tight')
    plt.close()

    # ── 5. Phase-% heatmap ────────────────────────────────────────────────────
    fig, ax = plt.subplots(figsize=(7, max(4, len(groups) * 0.55 + 2)))
    sns.heatmap(
        pct_ordered, annot=True, fmt='.1f', cmap='YlOrRd', ax=ax,
        linewidths=0.5, cbar_kws={'label': '% of cells'},
    )
    ax.set_xlabel('Cell Cycle Phase', fontsize=12)
    ax.set_ylabel(infection_col, fontsize=12)
    ax.set_title(
        f'Cell cycle phase % by {infection_col}\n'
        f"Overall chi2={chi2_all:.2f}, p={p_all:.2e}, Cramer's V={cramers_v:.3f}",
        fontsize=11,
    )
    for lbl in ax.get_yticklabels():
        if lbl.get_text() == infected_label:
            lbl.set_fontweight('bold')
            lbl.set_color('darkred')
    plt.tight_layout()
    plt.savefig(os.path.join(fig_dir, f'infection_cellcycle_heatmap_{sample}.pdf'),
                bbox_inches='tight')
    plt.close()

    # ── 6. Score violins by group ─────────────────────────────────────────────
    sc.settings.figdir = fig_dir
    if 'S_score' in adata.obs.columns and 'G2M_score' in adata.obs.columns:
        # Temporarily set a clean categorical for plotting without mutating original
        adata.obs['_infection_plot_col'] = pd.Categorical(
            adata.obs[infection_col].astype(str), categories=ordered, ordered=True
        )
        fig, axes = plt.subplots(1, 2, figsize=(max(10, len(groups) * 2), 6))
        sc.pl.violin(adata, 'S_score',   groupby='_infection_plot_col',
                     ax=axes[0], show=False, rotation=30, order=ordered)
        axes[0].set_title(f'S Score by {infection_col}\n({infected_label} highlighted)')
        axes[0].set_xlabel(infection_col)
        axes[0].axhline(0, color='k', linestyle='--', alpha=0.4)
        for lbl in axes[0].get_xticklabels():
            if lbl.get_text() == infected_label:
                lbl.set_color('darkred'); lbl.set_fontweight('bold')

        sc.pl.violin(adata, 'G2M_score', groupby='_infection_plot_col',
                     ax=axes[1], show=False, rotation=30, order=ordered)
        axes[1].set_title(f'G2M Score by {infection_col}\n({infected_label} highlighted)')
        axes[1].set_xlabel(infection_col)
        axes[1].axhline(0, color='k', linestyle='--', alpha=0.4)
        for lbl in axes[1].get_xticklabels():
            if lbl.get_text() == infected_label:
                lbl.set_color('darkred'); lbl.set_fontweight('bold')

        plt.tight_layout()
        plt.savefig(os.path.join(fig_dir, f'infection_scores_violin_{sample}.pdf'),
                    bbox_inches='tight')
        plt.close()
        # Clean up temp column
        del adata.obs['_infection_plot_col']

    # ── 7. Save CSVs ──────────────────────────────────────────────────────────
    contingency.to_csv(
        os.path.join(fig_dir, f'infection_contingency_counts_{sample}.csv'))
    contingency_pct.to_csv(
        os.path.join(fig_dir, f'infection_contingency_pct_{sample}.csv'))
    pairwise_chisq_df.to_csv(
        os.path.join(fig_dir, f'infection_pairwise_chisq_{sample}.csv'), index=False)
    pairwise_scores_df.to_csv(
        os.path.join(fig_dir, f'infection_pairwise_scores_{sample}.csv'), index=False)

    overall_stats = pd.DataFrame([{
        'test':            'chi-square (all groups)',
        'infection_col':   infection_col,
        'infected_label':  infected_label,
        'chi2':            chi2_all,
        'dof':             dof_all,
        'p_value':         p_all,
        'cramers_v':       cramers_v,
        'effect_size':     effect,
        'n_comparisons':   n_comparisons,
        'bonf_threshold':  bonf_alpha,
    }])
    overall_stats.to_csv(
        os.path.join(fig_dir, f'infection_overall_stats_{sample}.csv'), index=False)

    print(f"\nInfection cell cycle outputs saved to: {fig_dir}")

    return {
        'chi2':            chi2_all,
        'p_value':         p_all,
        'cramers_v':       cramers_v,
        'pairwise_chisq':  pairwise_chisq_df,
        'pairwise_scores': pairwise_scores_df,
        'contingency_pct': contingency_pct,
        'delta_vs_mean':   delta,
    }

def analyze_titer(adata, fig_dir, sample, titer_col='titer'):
    """
    Three-part titer analysis:
      A. Titer vs Leiden cluster  – does titer differ across clusters?
      B. Titer vs Cell cycle phase – does titer differ across phases?
      C. Titer vs continuous cell cycle scores (Spearman correlation)
      D. Three-way summary: mean titer per cluster x phase heatmap

    Requires:
      adata.obs[titer_col]  – numeric titer per cell (log-scale recommended)
      adata.obs['phase']    – from score_cell_cycle_scanpy
      adata.obs[leiden_col] – 'leiden_ref' or 'leiden'
    """
    from scipy.stats import kruskal, spearmanr, mannwhitneyu
    from itertools import combinations
    import warnings

    print("\n" + "=" * 60)
    print("TITER ANALYSIS")
    print(f"Titer column: '{titer_col}'")
    print("=" * 60)

    # ── Validation ────────────────────────────────────────────────────────────
    if titer_col not in adata.obs.columns:
        print(f"ERROR: '{titer_col}' not found in adata.obs.")
        print(f"Available numeric columns: "
              f"{[c for c in adata.obs.columns if pd.api.types.is_numeric_dtype(adata.obs[c])]}")
        return None

    if 'phase' not in adata.obs.columns:
        print("ERROR: No 'phase' column. Run --run-scoring first.")
        return None

    try:
        leiden_col = _get_leiden_col(adata)
    except KeyError as e:
        print(f"ERROR: {e}")
        return None

    os.makedirs(fig_dir, exist_ok=True)
    sc.settings.figdir = fig_dir

    # ── Sanity check: show titer distribution by cell_line if present ────────
    for ref_col in ('cell_line', 'treatment'):
        if ref_col in adata.obs.columns:
            print(f"\nTiter distribution by '{ref_col}':")
            desc = adata.obs.groupby(ref_col, observed=True)[titer_col].describe()
            print(desc.round(4))
            frac = adata.obs.groupby(ref_col, observed=True)[titer_col].apply(
                lambda x: (x > 0).mean()
            )
            print(f"\nFraction titer > 0 per '{ref_col}':")
            print(frac.round(4))
            print()
            break   # only print once for the first matching column

    # Work on non-NaN titer cells only
    obs = adata.obs[[titer_col, 'phase', leiden_col, 'S_score', 'G2M_score']].dropna(
        subset=[titer_col]
    ).copy()
    obs[titer_col] = pd.to_numeric(obs[titer_col], errors='coerce')
    obs = obs.dropna(subset=[titer_col])

    n_titer = len(obs)
    n_total = adata.n_obs
    print(f"\nCells with valid titer: {n_titer} / {n_total} "
          f"({n_titer / n_total * 100:.1f}%)")
    print(f"Titer range: {obs[titer_col].min():.3f} – {obs[titer_col].max():.3f}  "
          f"(median {obs[titer_col].median():.3f})")

    clusters = sorted(obs[leiden_col].unique())
    phases   = [p for p in ['G1', 'S', 'G2M'] if p in obs['phase'].unique()]

    # Colour palettes
    cmap_cluster = plt.cm.get_cmap('tab20')
    cluster_colors = {c: cmap_cluster(i % 20) for i, c in enumerate(clusters)}
    phase_colors   = {'G1': '#FF6B6B', 'S': '#4ECDC4', 'G2M': '#45B7D1'}

    results = {}

    # ══════════════════════════════════════════════════════════════════════════
    # A. TITER vs CLUSTER
    # ══════════════════════════════════════════════════════════════════════════
    print("\n" + "─" * 60)
    print("A. TITER vs LEIDEN CLUSTER")
    print("─" * 60)

    groups_by_cluster = [obs.loc[obs[leiden_col] == c, titer_col].values for c in clusters]
    # Drop clusters with < 3 cells (can't do stats)
    valid_clusters  = [c for c, g in zip(clusters, groups_by_cluster) if len(g) >= 3]
    valid_groups_cl = [g for g in groups_by_cluster if len(g) >= 3]

    h_cl, p_cl = kruskal(*valid_groups_cl)
    print(f"Kruskal-Wallis: H={h_cl:.3f}, p={p_cl:.2e}")
    print(f"({'SIGNIFICANT' if p_cl < 0.05 else 'not significant'} at α=0.05)")

    # Mean titer per cluster
    cluster_stats = obs.groupby(leiden_col)[titer_col].agg(['mean', 'median', 'std', 'count'])
    print(f"\nMean titer per cluster:\n{cluster_stats.round(3)}")

    # Pairwise Mann-Whitney (Bonferroni)
    pairs_cl     = list(combinations(valid_clusters, 2))
    bonf_cl      = 0.05 / len(pairs_cl) if pairs_cl else 0.05
    pairwise_cl  = []
    for a, b in pairs_cl:
        ga = obs.loc[obs[leiden_col] == a, titer_col].values
        gb = obs.loc[obs[leiden_col] == b, titer_col].values
        u, p = mannwhitneyu(ga, gb, alternative='two-sided')
        pairwise_cl.append({
            'cluster_A': a, 'cluster_B': b,
            'U': u, 'p_value': p,
            'bonf_significant': p < bonf_cl,
        })
    pairwise_cl_df = pd.DataFrame(pairwise_cl)
    n_sig_cl = pairwise_cl_df['bonf_significant'].sum() if not pairwise_cl_df.empty else 0
    print(f"\nPairwise MW-U (Bonferroni α={bonf_cl:.4f}): "
          f"{n_sig_cl}/{len(pairs_cl)} pairs significant")

    results['kruskal_cluster_H']  = h_cl
    results['kruskal_cluster_p']  = p_cl
    results['cluster_stats']      = cluster_stats
    results['pairwise_cluster']   = pairwise_cl_df

    # ── Figure A: violin + swarm of titer by cluster ──────────────────────────
    fig, axes = plt.subplots(1, 2, figsize=(max(14, len(clusters) * 0.9 + 4), 6))

    # Violin
    data_for_violin = [obs.loc[obs[leiden_col] == c, titer_col].values for c in clusters]
    vp = axes[0].violinplot(data_for_violin, positions=range(len(clusters)),
                            showmedians=True, showextrema=False)
    for i, body in enumerate(vp['bodies']):
        body.set_facecolor(cluster_colors[clusters[i]])
        body.set_alpha(0.7)
    vp['cmedians'].set_color('black')
    axes[0].set_xticks(range(len(clusters)))
    axes[0].set_xticklabels(clusters)
    axes[0].set_xlabel(f'Leiden Cluster ({leiden_col})', fontsize=11)
    axes[0].set_ylabel(f'{titer_col}', fontsize=11)
    axes[0].set_title(
        f'Titer by cluster\nKruskal-Wallis H={h_cl:.2f}, p={p_cl:.2e}', fontsize=12)
    axes[0].axhline(obs[titer_col].median(), color='grey', linestyle='--',
                    alpha=0.5, label='overall median')
    axes[0].legend(fontsize=9)

    # Bar: mean ± SD
    means = cluster_stats['mean']
    stds  = cluster_stats['std']
    axes[1].bar(range(len(clusters)), means,
                yerr=stds, capsize=4,
                color=[cluster_colors[c] for c in clusters], alpha=0.8, edgecolor='black')
    axes[1].set_xticks(range(len(clusters)))
    axes[1].set_xticklabels(clusters)
    axes[1].set_xlabel(f'Leiden Cluster ({leiden_col})', fontsize=11)
    axes[1].set_ylabel(f'Mean {titer_col} ± SD', fontsize=11)
    axes[1].set_title('Mean titer per cluster', fontsize=12)
    axes[1].axhline(obs[titer_col].mean(), color='grey', linestyle='--',
                    alpha=0.5, label='overall mean')
    axes[1].legend(fontsize=9)

    plt.tight_layout()
    plt.savefig(os.path.join(fig_dir, f'titer_vs_cluster_{sample}.pdf'), bbox_inches='tight')
    plt.close()

    # ── UMAP coloured by titer (if available) ─────────────────────────────────
    if 'X_umap' in adata.obsm and n_titer > 0:
        sc.pl.umap(adata, color=titer_col, cmap='viridis',
                   save=f'_{sample}_titer.pdf',
                   title=f'{titer_col} (UMAP)')

    # ══════════════════════════════════════════════════════════════════════════
    # B. TITER vs CELL CYCLE PHASE
    # ══════════════════════════════════════════════════════════════════════════
    print("\n" + "─" * 60)
    print("B. TITER vs CELL CYCLE PHASE")
    print("─" * 60)

    groups_by_phase = [obs.loc[obs['phase'] == p, titer_col].values for p in phases]
    valid_phases    = [p for p, g in zip(phases, groups_by_phase) if len(g) >= 3]
    valid_groups_ph = [g for p, g in zip(phases, groups_by_phase) if len(g) >= 3]

    h_ph, p_ph = kruskal(*valid_groups_ph)
    print(f"Kruskal-Wallis: H={h_ph:.3f}, p={p_ph:.2e}")
    print(f"({'SIGNIFICANT' if p_ph < 0.05 else 'not significant'} at α=0.05)")

    phase_stats = obs.groupby('phase')[titer_col].agg(['mean', 'median', 'std', 'count'])
    print(f"\nMean titer per phase:\n{phase_stats.round(3)}")

    # Pairwise
    pairs_ph    = list(combinations(valid_phases, 2))
    bonf_ph     = 0.05 / len(pairs_ph) if pairs_ph else 0.05
    pairwise_ph = []
    for a, b in pairs_ph:
        ga = obs.loc[obs['phase'] == a, titer_col].values
        gb = obs.loc[obs['phase'] == b, titer_col].values
        u, p = mannwhitneyu(ga, gb, alternative='two-sided')
        pairwise_ph.append({
            'phase_A': a, 'phase_B': b,
            'U': u, 'p_value': p,
            'bonf_significant': p < bonf_ph,
        })
    pairwise_ph_df = pd.DataFrame(pairwise_ph)
    n_sig_ph = pairwise_ph_df['bonf_significant'].sum() if not pairwise_ph_df.empty else 0
    print(f"\nPairwise MW-U (Bonferroni α={bonf_ph:.4f}): "
          f"{n_sig_ph}/{len(pairs_ph)} pairs significant")

    results['kruskal_phase_H'] = h_ph
    results['kruskal_phase_p'] = p_ph
    results['phase_stats']     = phase_stats
    results['pairwise_phase']  = pairwise_ph_df

    # ── Figure B: violin + bar for titer by phase ─────────────────────────────
    fig, axes = plt.subplots(1, 2, figsize=(10, 6))

    data_by_phase = [obs.loc[obs['phase'] == p, titer_col].values for p in phases]
    vp2 = axes[0].violinplot(data_by_phase, positions=range(len(phases)),
                             showmedians=True, showextrema=False)
    for i, body in enumerate(vp2['bodies']):
        body.set_facecolor(phase_colors.get(phases[i], 'grey'))
        body.set_alpha(0.7)
    vp2['cmedians'].set_color('black')
    axes[0].set_xticks(range(len(phases)))
    axes[0].set_xticklabels(phases)
    axes[0].set_xlabel('Cell Cycle Phase', fontsize=11)
    axes[0].set_ylabel(f'{titer_col}', fontsize=11)
    axes[0].set_title(
        f'Titer by phase\nKruskal-Wallis H={h_ph:.2f}, p={p_ph:.2e}', fontsize=12)

    # Annotate pairwise significance above violins
    y_max = obs[titer_col].max()
    y_range = obs[titer_col].max() - obs[titer_col].min()
    step = y_range * 0.07
    for k, row in pairwise_ph_df.iterrows():
        if row['bonf_significant']:
            i = phases.index(row['phase_A'])
            j = phases.index(row['phase_B'])
            y = y_max + step * (k + 1)
            axes[0].plot([i, j], [y, y], color='black', linewidth=1)
            axes[0].text((i + j) / 2, y + step * 0.2, '*', ha='center',
                         fontsize=12, fontweight='bold')

    axes[1].bar(range(len(phases)),
                [phase_stats.loc[p, 'mean'] for p in phases],
                yerr=[phase_stats.loc[p, 'std'] for p in phases],
                color=[phase_colors.get(p, 'grey') for p in phases],
                capsize=5, alpha=0.8, edgecolor='black')
    axes[1].set_xticks(range(len(phases)))
    axes[1].set_xticklabels(phases)
    axes[1].set_xlabel('Cell Cycle Phase', fontsize=11)
    axes[1].set_ylabel(f'Mean {titer_col} ± SD', fontsize=11)
    axes[1].set_title('Mean titer per phase', fontsize=12)

    plt.tight_layout()
    plt.savefig(os.path.join(fig_dir, f'titer_vs_phase_{sample}.pdf'), bbox_inches='tight')
    plt.close()

    # ══════════════════════════════════════════════════════════════════════════
    # C. TITER vs CONTINUOUS CELL CYCLE SCORES (Spearman)
    # ══════════════════════════════════════════════════════════════════════════
    print("\n" + "─" * 60)
    print("C. TITER vs CONTINUOUS CELL CYCLE SCORES (Spearman)")
    print("─" * 60)

    score_corrs = {}
    for score_col in ['S_score', 'G2M_score']:
        if score_col not in obs.columns:
            continue
        valid = obs[[titer_col, score_col]].dropna()
        r, p = spearmanr(valid[titer_col], valid[score_col])
        print(f"  {score_col}: rho={r:.4f}, p={p:.2e} "
              f"({'sig' if p < 0.05 else 'ns'})")
        score_corrs[score_col] = {'rho': r, 'p': p}

    results['score_correlations'] = score_corrs

    # ── Figure C: scatter titer vs S_score and G2M_score ─────────────────────
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    for ax, score_col in zip(axes, ['S_score', 'G2M_score']):
        if score_col not in obs.columns:
            ax.set_visible(False)
            continue
        valid = obs[[titer_col, score_col, 'phase']].dropna()
        for ph in phases:
            mask = valid['phase'] == ph
            ax.scatter(valid.loc[mask, titer_col],
                       valid.loc[mask, score_col],
                       c=phase_colors.get(ph, 'grey'), label=ph,
                       alpha=0.4, s=8, rasterized=True)
        # Regression line
        from numpy.polynomial.polynomial import polyfit as npfit
        x_ = valid[titer_col].values
        y_ = valid[score_col].values
        if len(x_) > 2:
            coefs = np.polyfit(x_, y_, 1)
            x_line = np.linspace(x_.min(), x_.max(), 200)
            ax.plot(x_line, np.polyval(coefs, x_line), 'k--', linewidth=1.5)
        corr = score_corrs.get(score_col, {})
        ax.set_xlabel(f'{titer_col}', fontsize=11)
        ax.set_ylabel(score_col, fontsize=11)
        ax.set_title(
            f'{titer_col} vs {score_col}\n'
            f"Spearman rho={corr.get('rho', 0):.3f}, p={corr.get('p', 1):.2e}",
            fontsize=12,
        )
        ax.legend(title='Phase', fontsize=9)
        ax.axhline(0, color='grey', linestyle='--', alpha=0.4)

    plt.tight_layout()
    plt.savefig(os.path.join(fig_dir, f'titer_vs_scores_{sample}.pdf'), bbox_inches='tight')
    plt.close()

    # ══════════════════════════════════════════════════════════════════════════
    # D. THREE-WAY SUMMARY: mean titer per cluster x phase heatmap
    # ══════════════════════════════════════════════════════════════════════════
    print("\n" + "─" * 60)
    print("D. THREE-WAY SUMMARY: mean titer per cluster × phase")
    print("─" * 60)

    three_way = obs.groupby([leiden_col, 'phase'])[titer_col].mean().unstack('phase')
    # Ensure phase column order
    for ph in ['G1', 'S', 'G2M']:
        if ph not in three_way.columns:
            three_way[ph] = np.nan
    three_way = three_way[['G1', 'S', 'G2M']]
    print(three_way.round(3))

    fig, axes = plt.subplots(1, 2, figsize=(16, max(6, len(clusters) * 0.5 + 3)))

    # Heatmap
    sns.heatmap(three_way, annot=True, fmt='.2f', cmap='viridis', ax=axes[0],
                cbar_kws={'label': f'Mean {titer_col}'},
                linewidths=0.5)
    axes[0].set_xlabel('Cell Cycle Phase', fontsize=11)
    axes[0].set_ylabel(f'Leiden Cluster ({leiden_col})', fontsize=11)
    axes[0].set_title(f'Mean {titer_col} per cluster × phase', fontsize=12)

    # Grouped bar: clusters on x, hue = phase
    x   = np.arange(len(clusters))
    w   = 0.25
    for j, ph in enumerate(['G1', 'S', 'G2M']):
        if ph not in three_way.columns:
            continue
        vals = [three_way.loc[c, ph] if c in three_way.index else np.nan
                for c in clusters]
        axes[1].bar(x + j * w, vals, width=w,
                    color=phase_colors.get(ph, 'grey'), label=ph,
                    alpha=0.8, edgecolor='black')
    axes[1].set_xticks(x + w)
    axes[1].set_xticklabels(clusters)
    axes[1].set_xlabel(f'Leiden Cluster ({leiden_col})', fontsize=11)
    axes[1].set_ylabel(f'Mean {titer_col}', fontsize=11)
    axes[1].set_title('Mean titer per cluster, split by phase', fontsize=12)
    axes[1].legend(title='Phase', fontsize=9)

    plt.tight_layout()
    plt.savefig(os.path.join(fig_dir, f'titer_cluster_phase_threeway_{sample}.pdf'),
                bbox_inches='tight')
    plt.close()

    # ── Save CSVs ──────────────────────────────────────────────────────────────
    cluster_stats.to_csv(os.path.join(fig_dir, f'titer_cluster_stats_{sample}.csv'))
    phase_stats.to_csv(  os.path.join(fig_dir, f'titer_phase_stats_{sample}.csv'))
    pairwise_cl_df.to_csv(os.path.join(fig_dir, f'titer_pairwise_cluster_{sample}.csv'), index=False)
    pairwise_ph_df.to_csv(os.path.join(fig_dir, f'titer_pairwise_phase_{sample}.csv'),   index=False)
    three_way.to_csv(     os.path.join(fig_dir, f'titer_cluster_phase_threeway_{sample}.csv'))

    corr_df = pd.DataFrame([
        {'score': k, 'spearman_rho': v['rho'], 'p_value': v['p']}
        for k, v in score_corrs.items()
    ])
    corr_df.to_csv(os.path.join(fig_dir, f'titer_score_correlations_{sample}.csv'), index=False)

    print(f"\nTiter analysis outputs saved to: {fig_dir}")
    print(f"  titer_vs_cluster_{sample}.pdf")
    print(f"  titer_vs_phase_{sample}.pdf")
    print(f"  titer_vs_scores_{sample}.pdf")
    print(f"  titer_cluster_phase_threeway_{sample}.pdf")

    return results



def analyze_titer_regression(adata, fig_dir, sample,
                              titer_col='wolbachia_titer',
                              leiden_col=None):
    """
    OLS regression: wolbachia_titer ~ C(phase) + C(leiden)

    Tests whether cell cycle phase predicts titer after controlling for
    cluster identity (and vice versa). This separates the independent
    contributions of cluster and phase to titer variation.

    Outputs
    -------
    Figures
      titer_regression_coefs_{sample}.pdf  – coefficient plot with 95% CIs
    CSVs
      titer_regression_summary_{sample}.csv – full coefficient table
      titer_regression_stats_{sample}.csv   – model-level stats (R², F, p)
    """
    try:
        import statsmodels.formula.api as smf
    except ImportError:
        print("ERROR: statsmodels not installed. Run: pip install statsmodels")
        return None

    print("\n" + "=" * 60)
    print("TITER REGRESSION: titer ~ phase + cluster")
    print("=" * 60)

    if 'phase' not in adata.obs.columns:
        print("ERROR: No 'phase' column. Run --run-scoring first.")
        return None
    if titer_col not in adata.obs.columns:
        print(f"ERROR: '{titer_col}' not found in adata.obs.")
        return None

    # Resolve leiden column
    if leiden_col is None:
        try:
            leiden_col = _get_leiden_col(adata)
        except KeyError as e:
            print(f"ERROR: {e}")
            return None

    os.makedirs(fig_dir, exist_ok=True)

    # Build regression dataframe
    reg_df = adata.obs[[titer_col, 'phase', leiden_col]].copy()
    reg_df[titer_col]  = pd.to_numeric(reg_df[titer_col], errors='coerce')
    reg_df[leiden_col] = reg_df[leiden_col].astype('category')
    reg_df = reg_df.dropna()

    n_cells = len(reg_df)
    print(f"\nCells in regression: {n_cells}")
    print(f"Titer median: {reg_df[titer_col].median():.4f}")
    print(f"Phase counts:\n{reg_df['phase'].value_counts().to_string()}")
    print(f"Cluster counts:\n{reg_df[leiden_col].value_counts().sort_index().to_string()}")

    # ── Fit model ─────────────────────────────────────────────────────────────
    formula = f"{titer_col} ~ C(phase) + C({leiden_col})"
    print(f"\nFormula: {formula}")
    print(f"Reference levels: phase=G1, {leiden_col}={sorted(reg_df[leiden_col].unique())[0]}")

    model  = smf.ols(formula, data=reg_df).fit()

    print("\n" + model.summary().as_text())

    # ── Clean coefficient table ───────────────────────────────────────────────
    coef_df = pd.DataFrame({
        'term':    model.params.index,
        'coef':    model.params.values,
        'se':      model.bse.values,
        'ci_low':  model.conf_int()[0].values,
        'ci_high': model.conf_int()[1].values,
        't':       model.tvalues.values,
        'p':       model.pvalues.values,
    })
    coef_df['sig'] = pd.cut(
        coef_df['p'],
        bins=[-1, 0.001, 0.01, 0.05, 1],
        labels=['***', '**', '*', 'ns'],
    )

    # Separate phase and cluster terms for printing
    phase_rows   = coef_df[coef_df['term'].str.contains('phase',   case=False)]
    cluster_rows = coef_df[coef_df['term'].str.contains(leiden_col, case=False)]

    print("\n" + "─" * 60)
    print("PHASE COEFFICIENTS (vs G1 reference)")
    print("─" * 60)
    print(phase_rows[['term', 'coef', 'se', 'ci_low', 'ci_high', 'p', 'sig']].to_string(index=False))

    print("\n" + "─" * 60)
    print(f"CLUSTER COEFFICIENTS (vs cluster {sorted(reg_df[leiden_col].unique())[0]} reference)")
    print("─" * 60)
    print(cluster_rows[['term', 'coef', 'se', 'ci_low', 'ci_high', 'p', 'sig']].to_string(index=False))

    print("\n" + "─" * 60)
    print("MODEL FIT")
    print("─" * 60)
    print(f"  R²        = {model.rsquared:.4f}")
    print(f"  Adj. R²   = {model.rsquared_adj:.4f}")
    print(f"  F-stat    = {model.fvalue:.2f}")
    print(f"  F p-value = {model.f_pvalue:.2e}")
    print(f"  AIC       = {model.aic:.1f}")
    print(f"  n cells   = {n_cells}")

    # ── Interpretation note ────────────────────────────────────────────────────
    phase_sig = (phase_rows['p'] < 0.05).any()
    max_phase_effect = phase_rows['coef'].abs().max() if not phase_rows.empty else 0
    print("\n" + "─" * 60)
    print("INTERPRETATION")
    print("─" * 60)
    if phase_sig:
        print(f"  Phase IS a significant predictor of titer after controlling")
        print(f"  for cluster (max |coef| = {max_phase_effect:.4f}).")
    else:
        print(f"  Phase is NOT a significant predictor of titer after controlling")
        print(f"  for cluster. The titer~phase association seen in the univariate")
        print(f"  analysis is likely driven by cluster composition differences.")
    print(f"  Model explains {model.rsquared * 100:.2f}% of titer variance.")

    # ── Coefficient plot ───────────────────────────────────────────────────────
    plot_rows = coef_df[coef_df['term'] != 'Intercept'].copy()
    # Shorten term labels for readability
    plot_rows['label'] = (
        plot_rows['term']
        .str.replace("C(phase)[T.", "phase: ", regex=False)
        .str.replace("C(" + leiden_col + ")[T.", "cluster: ", regex=False)
        .str.replace("]", "", regex=False)
    )
    # Color by term type
    plot_rows['color'] = plot_rows['term'].apply(
        lambda x: '#4ECDC4' if 'phase' in x.lower() else '#FF6B6B'
    )

    fig, ax = plt.subplots(figsize=(8, max(5, len(plot_rows) * 0.45 + 1)))
    y_pos = range(len(plot_rows))
    ax.barh(
        list(y_pos), plot_rows['coef'].values,
        xerr=[(plot_rows['coef'] - plot_rows['ci_low']).values,
              (plot_rows['ci_high'] - plot_rows['coef']).values],
        color=plot_rows['color'].values, alpha=0.8,
        edgecolor='black', linewidth=0.6, capsize=3, height=0.6,
    )
    ax.axvline(0, color='black', linewidth=1, linestyle='--', alpha=0.6)
    ax.set_yticks(list(y_pos))
    ax.set_yticklabels(plot_rows['label'].values, fontsize=9)
    ax.set_xlabel(f'Coefficient (effect on {titer_col})', fontsize=11)
    ax.set_title(
        f'OLS: {titer_col} ~ phase + cluster\n'
        f'R²={model.rsquared:.3f}, F p={model.f_pvalue:.2e}, n={n_cells}',
        fontsize=12,
    )
    # Add significance stars
    x_max = plot_rows['ci_high'].abs().max() * 1.15
    for i, (_, row) in enumerate(plot_rows.iterrows()):
        if row['sig'] != 'ns':
            ax.text(x_max * 0.98, i, row['sig'],
                    ha='right', va='center', fontsize=9, color='black')
    # Legend
    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor='#4ECDC4', label='Phase term', alpha=0.8),
        Patch(facecolor='#FF6B6B', label='Cluster term', alpha=0.8),
    ]
    ax.legend(handles=legend_elements, loc='lower right', fontsize=9)
    plt.tight_layout()
    plt.savefig(os.path.join(fig_dir, f'titer_regression_coefs_{sample}.pdf'),
                bbox_inches='tight')
    plt.close()

    # ── Save outputs ──────────────────────────────────────────────────────────
    coef_df.to_csv(
        os.path.join(fig_dir, f'titer_regression_summary_{sample}.csv'), index=False)

    stats_df = pd.DataFrame([{
        'formula':       formula,
        'n_cells':       n_cells,
        'r_squared':     model.rsquared,
        'adj_r_squared': model.rsquared_adj,
        'f_stat':        model.fvalue,
        'f_pvalue':      model.f_pvalue,
        'aic':           model.aic,
        'phase_sig':     phase_sig,
        'max_phase_coef': max_phase_effect,
    }])
    stats_df.to_csv(
        os.path.join(fig_dir, f'titer_regression_stats_{sample}.csv'), index=False)

    print(f"\nRegression outputs saved to: {fig_dir}")
    print(f"  titer_regression_coefs_{sample}.pdf")
    print(f"  titer_regression_summary_{sample}.csv")
    print(f"  titer_regression_stats_{sample}.csv")

    return {
        'model':       model,
        'coef_df':     coef_df,
        'r_squared':   model.rsquared,
        'f_pvalue':    model.f_pvalue,
        'phase_sig':   phase_sig,
    }


def main():
    parser = argparse.ArgumentParser(
        description='Cell cycle annotation using Scanpy with Drosophila genes',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='''
Examples:
  # Reference object
  python scanpy_cellcycle_analysis.py \\
      --input integrated_reference.h5ad \\
      --output results/cellcycle/ref --sample ref --run-scoring --save-output

  # JW18DOX vs JW18wMel (cell_line column)
  python scanpy_cellcycle_analysis.py \\
      --input integrated.h5ad \\
      --output results/cellcycle/strain \\
      --sample strain --run-scoring \\
      --infection-analysis \\
      --infection-col cell_line \\
      --infected-label JW18DOX

  # JW18DOX-SV-D1 vs others using bio_condition
  python scanpy_cellcycle_analysis.py \\
      --input integrated.h5ad \\
      --output results/cellcycle/biocond \\
      --sample biocond --run-scoring \\
      --infection-analysis \\
      --infection-col bio_condition \\
      --infected-label JW18DOX-SV-D1

  # Titer analysis restricted to infected cells only (JW18wMel)
  # Use --subset-col/--subset-val to restrict adata before all analyses
  python scanpy_cellcycle_analysis.py \\
      --input integrated.h5ad \\
      --output results/cellcycle/titer_wMel \\
      --sample wMel \\
      --subset-col cell_line --subset-val JW18wMel \\
      --titer-analysis --titer-col wolbachia_titer

  # All analyses at once (infection comparison uses full data,
  # titer analysis uses wMel-only subset — run separately)
  python scanpy_cellcycle_analysis.py \\
      --input integrated.h5ad \\
      --output results/cellcycle/full \\
      --sample full --run-scoring \\
      --infection-analysis --infection-col cell_line --infected-label JW18wMel \\
      --titer-analysis --titer-col wolbachia_titer
        '''
    )

    parser.add_argument('--input',  '-i', required=True,
                        help='Path to h5ad file')
    parser.add_argument('--output', '-o', default='cellcycle_analysis',
                        help='Output directory (default: cellcycle_analysis)')
    parser.add_argument('--sample', '-s', default='sample',
                        help='Sample name for output files (default: sample)')
    parser.add_argument('--run-scoring', action='store_true',
                        help='Run Scanpy cell cycle scoring')
    parser.add_argument('--save-output', action='store_true',
                        help='Save updated h5ad with cell cycle annotations')
    parser.add_argument('--output-h5ad', default=None,
                        help='Path for saved h5ad (default: input with _with_cellcycle suffix)')    
    parser.add_argument('--ctrl-only', action='store_true',
                        help='Restrict to treatment==Ctrl cells only')
    parser.add_argument('--subset-col', default=None,
                        help="obs column to subset on before all analyses "
                             "(e.g. 'cell_line'). Must be paired with --subset-val.")
    parser.add_argument('--subset-val', default=None,
                        help="Value to keep in --subset-col "
                             "(e.g. 'JW18wMel'). Applied after --ctrl-only.")

    # Infection analysis args
    parser.add_argument('--infection-analysis', action='store_true',
                        help='Run infection vs cell cycle analysis')
    parser.add_argument('--infection-col', default='cell_line',
                        help="obs column to group by for infection analysis "
                             "(default: 'cell_line'). "
                             "e.g. 'cell_line' for JW18DOX vs JW18wMel, "
                             "'treatment' for Ctrl vs SV, "
                             "'bio_condition' for fine-grained comparisons.")
    parser.add_argument('--infected-label', default='JW18DOX',
                        help="Label in --infection-col identifying the infected group "
                             "(default: 'JW18DOX').")

    # Titer analysis args
    parser.add_argument('--titer-analysis', action='store_true',
                        help='Run titer vs cluster/phase/score analysis')
    parser.add_argument('--titer-col', default='wolbachia_titer',
                        help="adata.obs column for per-cell titer "
                             "(default: 'wolbachia_titer').")

    args = parser.parse_args()

    os.makedirs(args.output, exist_ok=True)

    print(f"Loading data from {args.input}...")
    adata = sc.read_h5ad(args.input)

    if args.ctrl_only:
        if 'treatment' not in adata.obs.columns:
            print("WARNING: --ctrl-only requested but 'treatment' column not found. Skipping.")
        else:
            before = adata.n_obs
            adata  = adata[adata.obs['treatment'] == 'Ctrl'].copy()
            print(f"Filtered to Ctrl cells: {before} -> {adata.n_obs} cells")
            if adata.n_obs == 0:
                raise ValueError("No cells remaining after --ctrl-only filter.")

    # ── Generic subset (e.g. restrict to JW18wMel for titer analysis) ──────────
    if args.subset_col is not None and args.subset_val is not None:
        if args.subset_col not in adata.obs.columns:
            raise ValueError(
                f"--subset-col '{args.subset_col}' not found in adata.obs. "
                f"Available: {list(adata.obs.columns)}"
            )
        before = adata.n_obs
        adata = adata[adata.obs[args.subset_col].astype(str) == args.subset_val].copy()
        print(f"\nSubset to {args.subset_col}=='{args.subset_val}': "
              f"{before} -> {adata.n_obs} cells")
        if adata.n_obs == 0:
            raise ValueError(
                f"No cells remaining after subsetting {args.subset_col}=='{args.subset_val}'. "
                f"Available values: {sorted(adata.obs[args.subset_col].unique())}"
            )
    elif (args.subset_col is None) != (args.subset_val is None):
        raise ValueError("--subset-col and --subset-val must be provided together.")

    print(f"\nLoaded AnnData: {adata.n_obs} cells, {adata.n_vars} genes")
    print(f"obs columns: {list(adata.obs.columns)}")

    # Print a quick value_counts for the infection column so user can verify
    if args.infection_analysis:
        if args.infection_col in adata.obs.columns:
            print(f"\nValues in '{args.infection_col}':")
            print(adata.obs[args.infection_col].value_counts().to_string())
        else:
            print(f"\nWARNING: --infection-col '{args.infection_col}' not found in obs.")

    has_scoring = all(c in adata.obs.columns for c in ['phase', 'S_score', 'G2M_score'])

    if args.run_scoring or not has_scoring:
        result = score_cell_cycle_scanpy(adata, args.output, args.sample)
        if result is None:
            print("\nERROR: Cell cycle scoring failed. Exiting.")
            return
    if args.save_output:
        out_path = args.output_h5ad or args.input.replace('.h5ad', '_with_cellcycle.h5ad')
        os.makedirs(os.path.dirname(out_path) or '.', exist_ok=True)
        print(f"\nSaving updated h5ad to {out_path}")
        adata.write_h5ad(out_path)
        print(f"Saved: {adata.n_obs} cells with columns: "
            f"{[c for c in ['phase','S_score','G2M_score'] if c in adata.obs.columns]}")
    else:
        print("\nCell cycle scoring already present, skipping scoring step.")
        print(f"Phase distribution: {adata.obs['phase'].value_counts().to_dict()}")

    # ── Standard cluster-phase analysis ──────────────────────────────────────
    results = analyze_cluster_cellcycle_association(adata, args.output, args.sample)

    if results:
        print(f"\n{'=' * 60}")
        print("KEY FINDINGS — CLUSTER vs CELL CYCLE")
        print(f"{'=' * 60}")
        print(f"Chi-square: chi2={results['chi2']:.2f}, p={results['chi2_pvalue']:.2e}")
        print(f"Cramer's V: {results['cramers_v']:.3f}")
        print(f"Kruskal-Wallis S score:   p={results['kw_s_pvalue']:.2e}")
        print(f"Kruskal-Wallis G2M score: p={results['kw_g2m_pvalue']:.2e}")

    # ── Infection analysis ────────────────────────────────────────────────────
    if args.infection_analysis:
        infection_results = analyze_infection_cellcycle(
            adata,
            fig_dir        = args.output,
            sample         = args.sample,
            infection_col  = args.infection_col,
            infected_label = args.infected_label,
        )
        if infection_results:
            print(f"\n{'=' * 60}")
            print("KEY FINDINGS — INFECTION vs CELL CYCLE")
            print(f"  Column: '{args.infection_col}', Label: '{args.infected_label}'")
            print(f"{'=' * 60}")
            print(f"Overall chi2={infection_results['chi2']:.2f}, "
                  f"p={infection_results['p_value']:.2e}, "
                  f"Cramer's V={infection_results['cramers_v']:.3f}")
            print(f"\nPhase Δ% ({args.infected_label} − mean of others):")
            for ph, val in infection_results['delta_vs_mean'].items():
                print(f"  {ph}: {val:+.2f}%")
            sig_chisq = infection_results['pairwise_chisq']['bonf_significant'].sum()
            sig_score = (infection_results['pairwise_scores']['S_score_bonf_sig'] |
                         infection_results['pairwise_scores']['G2M_score_bonf_sig']).sum()
            print(f"\nPairwise chi-square: {sig_chisq}/"
                  f"{len(infection_results['pairwise_chisq'])} Bonf. significant")
            print(f"Pairwise score MW-U: {sig_score}/"
                  f"{len(infection_results['pairwise_scores'])} Bonf. significant (S or G2M)")

    # ── Titer analysis ────────────────────────────────────────────────────────
    if args.titer_analysis:
        titer_results = analyze_titer(
            adata,
            fig_dir   = args.output,
            sample    = args.sample,
            titer_col = args.titer_col,
        )
        if titer_results:
            print(f"\n{'=' * 60}")
            print("KEY FINDINGS — TITER")
            print(f"{'=' * 60}")
            print(f"Titer vs cluster: Kruskal-Wallis H={titer_results['kruskal_cluster_H']:.2f}, "
                  f"p={titer_results['kruskal_cluster_p']:.2e}")
            print(f"Titer vs phase:   Kruskal-Wallis H={titer_results['kruskal_phase_H']:.2f}, "
                  f"p={titer_results['kruskal_phase_p']:.2e}")
            for score, corr in titer_results['score_correlations'].items():
                print(f"Titer vs {score}: Spearman rho={corr['rho']:.3f}, p={corr['p']:.2e}")

        # Regression: titer ~ phase + cluster (controls for confounding)
        reg_results = analyze_titer_regression(
            adata,
            fig_dir   = args.output,
            sample    = args.sample,
            titer_col = args.titer_col,
        )
        if reg_results:
            print(f"\n{'=' * 60}")
            print("KEY FINDINGS — TITER REGRESSION")
            print(f"{'=' * 60}")
            print(f"R² = {reg_results['r_squared']:.4f}")
            print(f"F-test p = {reg_results['f_pvalue']:.2e}")
            print(f"Phase significant after controlling for cluster: "
                  f"{'YES' if reg_results['phase_sig'] else 'NO'}")


if __name__ == "__main__":
    main()