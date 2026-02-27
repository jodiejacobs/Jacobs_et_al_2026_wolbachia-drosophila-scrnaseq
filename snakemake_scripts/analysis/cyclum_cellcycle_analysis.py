'''
Cyclum cell cycle analysis with:
  1. Cyclum pseudotime-based phase assignment (g0/g1, s, g2/m)
  2. Validation of cyclum phases using Drosophila FlyBase marker gene expression
     - S-phase markers should be enriched in cyclum 's' cells
     - G2/M markers should be enriched in cyclum 'g2/m' cells
     - Continuous S_score / G2M_score computed ONLY for validation, not reclassification
  3. Leiden cluster ~ cyclum cell cycle phase association

Usage:
  python cyclum_cellcycle_analysis.py \
      --input filtered.h5ad \
      --output results/cyclum_cellcycle/sample1 \
      --sample sample1 \
      --save-h5ad
'''

import cyclum
import cyclum.models
import cyclum.tuning
import cyclum.illustration
import scanpy as sc
import argparse
import os
from sklearn.neighbors import NearestNeighbors
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import kruskal, mannwhitneyu
from itertools import combinations
import scipy.sparse

# ── FlyBase cell cycle gene sets ──────────────────────────────────────────────
FLYBASE_CELL_CYCLE_GENES = {
    # S phase
    'Pcna':           'FBgn0005655',
    'RPA1':           'FBgn0015806',
    'RPA2':           'FBgn0034898',
    'pol-alpha1':     'FBgn0011230',
    'DNApol-alpha60': 'FBgn0015278',
    'DNApol-delta':   'FBgn0019624',
    'RnrL':           'FBgn0020369',
    'RnrS':           'FBgn0261933',
    'Mcm2':           'FBgn0020651',
    'Mcm3':           'FBgn0020652',
    'Mcm5':           'FBgn0015929',
    'Mcm6':           'FBgn0032435',
    'Mcm7':           'FBgn0015308',
    'E2f1':           'FBgn0011766',
    'E2f2':           'FBgn0262656',
    'CycE':           'FBgn0010382',
    'Cdk2':           'FBgn0010314',
    'Dp':             'FBgn0000499',
    'Rbf':            'FBgn0015799',
    'Rbf2':           'FBgn0028396',
    'Orc1':           'FBgn0015270',
    'Orc2':           'FBgn0015714',
    'Orc6':           'FBgn0025926',
    'Rrp1':           'FBgn0003257',
    # G2/M
    'CycA':           'FBgn0010114',
    'CycB':           'FBgn0010113',
    'CycB3':          'FBgn0011577',
    'Cdk1':           'FBgn0004107',
    'stg':            'FBgn0003525',
    'polo':           'FBgn0003124',
    'aurA':           'FBgn0025564',
    'aurB':           'FBgn0025948',
    'Nek2':           'FBgn0027548',
    'Pbl':            'FBgn0005619',
    'Wee1':           'FBgn0011739',
    'myt':            'FBgn0002863',
    'BubR1':          'FBgn0024822',
    'Mad2':           'FBgn0002610',
    'Cdc20':          'FBgn0010309',
    'APC2':           'FBgn0261823',
    'APC10':          'FBgn0036449',
}

S_GENES_FBGN = [
    'FBgn0005655', 'FBgn0015806', 'FBgn0034898', 'FBgn0011230',
    'FBgn0015278', 'FBgn0019624', 'FBgn0020369', 'FBgn0261933',
    'FBgn0020651', 'FBgn0020652', 'FBgn0015929', 'FBgn0032435',
    'FBgn0015308', 'FBgn0011766', 'FBgn0262656', 'FBgn0010382',
    'FBgn0010314', 'FBgn0000499', 'FBgn0015799', 'FBgn0028396',
    'FBgn0015270', 'FBgn0015714', 'FBgn0025926', 'FBgn0003257',
]

G2M_GENES_FBGN = [
    'FBgn0010114', 'FBgn0010113', 'FBgn0011577', 'FBgn0004107',
    'FBgn0003525', 'FBgn0003124', 'FBgn0025564', 'FBgn0025948',
    'FBgn0027548', 'FBgn0005619', 'FBgn0011739', 'FBgn0002863',
    'FBgn0024822', 'FBgn0002610', 'FBgn0010309', 'FBgn0261823',
    'FBgn0036449',
]

FBGN_TO_SYMBOL = {v: k for k, v in FLYBASE_CELL_CYCLE_GENES.items()}

PHASE_ORDER  = ['g0/g1', 's', 'g2/m']
PHASE_COLORS = {'g0/g1': '#FF6B6B', 's': '#4ECDC4', 'g2/m': '#45B7D1'}




# ══════════════════════════════════════════════════════════════════════════════
# HELPERS
# ══════════════════════════════════════════════════════════════════════════════

def _get_leiden_col(adata):
    for col in ('leiden_ref', 'leiden'):
        if col in adata.obs.columns:
            return col
    raise KeyError(
        "No leiden column found. Expected 'leiden_ref' or 'leiden'. "
        f"Available: {list(adata.obs.columns)}"
    )


def _get_lognorm_matrix(adata):
    """
    Return a log-normalised dense matrix for marker scoring.
    Does NOT modify adata.X — operates on a copy.
    """
    X = adata.X
    if scipy.sparse.issparse(X):
        X = X.toarray().astype(np.float32)
    else:
        X = X.astype(np.float32).copy()

    mat_max = float(X.max())
    if mat_max > 20:
        print(f"  Data appears to be raw counts (max={mat_max:.1f}). "
              "Log-normalising for marker scoring (adata.X unchanged).")
        row_sums = X.sum(axis=1, keepdims=True)
        row_sums[row_sums == 0] = 1
        X = np.log1p(X / row_sums * 1e4)
    else:
        print(f"  Data appears log-normalised (max={mat_max:.3f}).")
    return X


def _find_marker_genes(adata, verbose=True):
    """
    Return lists of S and G2M marker gene IDs/symbols present in adata.var_names.
    Tries FBgn IDs first, then falls back to gene symbols.
    """
    s_present   = [g for g in S_GENES_FBGN   if g in adata.var_names]
    g2m_present = [g for g in G2M_GENES_FBGN if g in adata.var_names]

    # Fall back to gene symbols if no FBgn IDs found
    if not s_present and not g2m_present:
        for sym, fbgn in FLYBASE_CELL_CYCLE_GENES.items():
            if sym in adata.var_names:
                if fbgn in S_GENES_FBGN:
                    s_present.append(sym)
                elif fbgn in G2M_GENES_FBGN:
                    g2m_present.append(sym)

    if verbose:
        print(f"  S-phase markers found  : {len(s_present)}/{len(S_GENES_FBGN)}")
        if s_present:
            labels = [FBGN_TO_SYMBOL.get(g, g) for g in s_present[:10]]
            print("    " + ", ".join(labels) +
                  (f"  (+{len(s_present)-10} more)" if len(s_present) > 10 else ""))
        print(f"  G2/M markers found     : {len(g2m_present)}/{len(G2M_GENES_FBGN)}")
        if g2m_present:
            labels = [FBGN_TO_SYMBOL.get(g, g) for g in g2m_present[:10]]
            print("    " + ", ".join(labels) +
                  (f"  (+{len(g2m_present)-10} more)" if len(g2m_present) > 10 else ""))

    return s_present, g2m_present


# ══════════════════════════════════════════════════════════════════════════════
# SECTION 1 – CYCLUM PHASE ASSIGNMENT
# ══════════════════════════════════════════════════════════════════════════════

def assign_cell_cycle_stage_simple(pseudotime_flat):
    """Assign g0/g1, s, g2/m based on cyclum circular pseudotime."""
    print("\n[Cyclum] Assigning cell cycle stages from pseudotime...")

    if pseudotime_flat.max() <= 1:
        angles = pseudotime_flat * 2 * np.pi
    else:
        angles = ((pseudotime_flat - pseudotime_flat.min()) /
                  (pseudotime_flat.max() - pseudotime_flat.min())) * 2 * np.pi

    boundary1      = 2 * np.pi / 3
    boundary2      = 4 * np.pi / 3
    boundary_width = np.pi / 12

    phases = []
    for angle in angles:
        a = angle % (2 * np.pi)
        if a < boundary1:
            phases.append('g0/g1')
        elif a < boundary2:
            phases.append('s')
        else:
            phases.append('g2/m')

    # Light boundary smoothing
    n_cells = len(angles)
    if n_cells > 10:
        nn   = NearestNeighbors(n_neighbors=min(10, n_cells // 10))
        circ = np.column_stack([np.cos(angles), np.sin(angles)])
        nn.fit(circ)
        smoothed = phases.copy()
        changes  = 0
        for i, angle in enumerate(angles):
            a = angle % (2 * np.pi)
            near = (abs(a - boundary1) < boundary_width or
                    abs(a - boundary2) < boundary_width or
                    min(a, 2 * np.pi - a) < boundary_width)
            if near:
                _, indices = nn.kneighbors([circ[i]])
                nbr = [phases[j] for j in indices[0][1:]]
                counts = {}
                for p in nbr:
                    counts[p] = counts.get(p, 0) + 1
                if counts:
                    best = max(counts, key=counts.get)
                    if counts[best] > len(nbr) * 0.7 and best != phases[i]:
                        smoothed[i] = best
                        changes += 1
        phases = smoothed
        print(f"  Boundary smoothing changed {changes} assignments")

    # Confidence (distance from nearest phase boundary)
    confidence = np.ones(len(angles))
    for i, angle in enumerate(angles):
        a  = angle % (2 * np.pi)
        d1 = min(abs(a - boundary1), 2*np.pi - abs(a - boundary1))
        d2 = min(abs(a - boundary2), 2*np.pi - abs(a - boundary2))
        d0 = min(a, 2*np.pi - a)
        confidence[i] = min(1.0, min(d1, d2, d0) / (boundary_width * 2))

    print("\n  Cyclum phase distribution:")
    phase_s = pd.Series(phases)
    for phase in PHASE_ORDER:
        n = (phase_s == phase).sum()
        print(f"    {phase}: {n}  ({n / len(phases) * 100:.1f}%)")

    return phases, confidence, angles


# ══════════════════════════════════════════════════════════════════════════════
# SECTION 2 – DATA-DRIVEN VALIDATION: TOP DE GENES BETWEEN CYCLUM PHASES
# ══════════════════════════════════════════════════════════════════════════════

def validate_cyclum_phases(adata, output_dir, sample_name,
                            n_top_genes=5, n_umap_genes=6):
    """
    Visually validate cyclum phase assignments by identifying the genes most
    differentially expressed between cyclum phases and plotting their expression
    on UMAPs — no assumption about which specific genes should be present.

    Strategy
    --------
    1. Run sc.tl.rank_genes_groups (Wilcoxon) with cyclum_stage as groupby.
       This finds the top genes that distinguish each phase from the others.
    2. Select the top `n_top_genes` per phase (deduplicated), preferring genes
       that overlap with known FlyBase cell cycle markers where possible.
    3. Plot one UMAP panel per gene coloured by expression, with inset showing
       cyclum phase labels for reference.
    4. Also produce a heatmap of top DE genes sorted by cyclum phase.

    A good cyclum assignment should show clear spatial structure on the UMAP:
    top 's' genes high in one UMAP region, top 'g2/m' genes in another.

    Parameters
    ----------
    n_top_genes : int
        Top DE genes to extract per phase for heatmap / CSV (default 5).
    n_umap_genes : int
        Number of genes to show per phase on UMAP grid (default 6).
        Total UMAPs = n_umap_genes * 3 phases + 1 (phase label).

    Outputs
    -------
    {sample}_validation_de_genes.csv         top DE genes per phase
    {sample}_validation_umap_phase.pdf       UMAP coloured by cyclum phase
    {sample}_validation_umap_{phase}.pdf     UMAP grids per phase (expression)
    {sample}_validation_umap_pseudotime.pdf  UMAP coloured by pseudotime
    {sample}_validation_heatmap.pdf          DE gene heatmap sorted by phase
    """
    print("\n" + "=" * 60)
    print("DATA-DRIVEN VALIDATION: TOP DE GENES BETWEEN CYCLUM PHASES")
    print("=" * 60)

    if 'cyclum_stage' not in adata.obs.columns:
        print("  ERROR: 'cyclum_stage' not in adata.obs. Run cyclum first.")
        return None

    if 'X_umap' not in adata.obsm:
        print("  WARNING: No UMAP found in adata.obsm. "
              "UMAP plots will be skipped. Run sc.pp.neighbors + sc.tl.umap first.")

    os.makedirs(output_dir, exist_ok=True)

    # ── Ensure log-normalised layer exists for DE ─────────────────────────────
    # rank_genes_groups works on adata.X — normalise a copy if needed
    adata_de = adata.copy()
    X_check  = adata_de.X
    mat_max  = float(X_check.max() if not scipy.sparse.issparse(X_check)
                     else X_check.max())
    if mat_max > 20:
        print(f"  Raw counts detected (max={mat_max:.1f}). "
              "Normalising copy for DE (adata.X unchanged).")
        sc.pp.normalize_total(adata_de, target_sum=1e4)
        sc.pp.log1p(adata_de)
    else:
        print(f"  Data appears log-normalised (max={mat_max:.3f}).")

    # ── Known FlyBase cell cycle gene IDs/symbols present in dataset ──────────
    s_known, g2m_known = _find_marker_genes(adata, verbose=True)
    known_cc_genes = set(s_known + g2m_known)
    # Build a lookup: gene name → 'S marker' / 'G2M marker' / ''
    gene_type_label = {}
    for g in s_known:
        gene_type_label[g] = f"S ({FBGN_TO_SYMBOL.get(g, g)})"
    for g in g2m_known:
        gene_type_label[g] = f"G2M ({FBGN_TO_SYMBOL.get(g, g)})"

    # ── Differential expression: top genes per cyclum phase ───────────────────
    print("\n  Running Wilcoxon rank-sum DE (cyclum_stage groups)...")
    sc.tl.rank_genes_groups(
        adata_de,
        groupby='cyclum_stage',
        groups=PHASE_ORDER,
        reference='rest',
        method='wilcoxon',
        key_added='rank_genes_cyclum',
        pts=True,
    )

    # Extract top genes per phase
    de_rows   = []
    top_genes = {phase: [] for phase in PHASE_ORDER}

    for phase in PHASE_ORDER:
        if phase not in adata_de.obs['cyclum_stage'].unique():
            continue
        df = sc.get.rank_genes_groups_df(
            adata_de, group=phase, key='rank_genes_cyclum',
            pval_cutoff=None, log2fc_min=None,
        )
        df = df.sort_values('scores', ascending=False)

        # Annotate with known cell cycle gene info
        df['known_cc_gene'] = df['names'].isin(known_cc_genes)
        df['gene_label']    = df['names'].map(
            lambda g: gene_type_label.get(g, FBGN_TO_SYMBOL.get(g, g)))
        df['phase'] = phase
        de_rows.append(df.head(50))   # save top 50 per phase to CSV

        # For plotting: take top n_umap_genes, boosting known CC genes
        # by pulling them to the front if they rank in top 20
        top_all  = list(df['names'].head(20))
        top_cc   = [g for g in top_all if g in known_cc_genes]
        top_rest = [g for g in top_all if g not in known_cc_genes]
        # Interleave: put known CC genes first, fill rest to n_umap_genes
        prioritised = (top_cc + top_rest)[:n_umap_genes]
        top_genes[phase] = prioritised

        label_str = ", ".join([gene_type_label.get(g, FBGN_TO_SYMBOL.get(g, g))
                               for g in prioritised[:5]])
        print(f"  Top genes for '{phase}': {label_str}"
              + (f"  (+{len(prioritised)-5} more)" if len(prioritised) > 5 else ""))

    de_df = pd.concat(de_rows, ignore_index=True)
    de_df.to_csv(
        os.path.join(output_dir, f'{sample_name}_validation_de_genes.csv'),
        index=False)
    print(f"  DE results saved ({len(de_df)} rows)")

    # ── UMAP: cyclum phase labels + pseudotime ────────────────────────────────
    if 'X_umap' in adata.obsm:

        # 1. Phase label UMAP
        fig, ax = plt.subplots(figsize=(7, 6))
        sc.pl.umap(adata, color='cyclum_stage', ax=ax, show=False,
                   title=f'Cyclum phase — {sample_name}',
                   palette=PHASE_COLORS, frameon=False)
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir,
                                 f'{sample_name}_validation_umap_phase.pdf'),
                    dpi=300, bbox_inches='tight')
        plt.close()

        # 2. Pseudotime UMAP
        if 'cyclum_pseudotime' in adata.obs.columns:
            fig, ax = plt.subplots(figsize=(7, 6))
            sc.pl.umap(adata, color='cyclum_pseudotime', ax=ax, show=False,
                       title=f'Cyclum pseudotime — {sample_name}',
                       cmap='hsv', frameon=False)
            plt.tight_layout()
            plt.savefig(os.path.join(output_dir,
                                     f'{sample_name}_validation_umap_pseudotime.pdf'),
                        dpi=300, bbox_inches='tight')
            plt.close()

        # 3. Per-phase UMAP grids: top DE gene expression
        for phase in PHASE_ORDER:
            genes = top_genes.get(phase, [])
            if not genes:
                continue

            # Filter to genes actually in adata
            genes = [g for g in genes if g in adata.var_names]
            if not genes:
                continue

            n_genes = len(genes)
            # Layout: phase label in first panel, then gene expression panels
            n_cols  = min(4, n_genes + 1)
            n_rows  = int(np.ceil((n_genes + 1) / n_cols))
            fig, axes = plt.subplots(n_rows, n_cols,
                                      figsize=(5 * n_cols, 4.5 * n_rows))
            axes = np.array(axes).flatten()

            # First panel: phase label
            sc.pl.umap(adata, color='cyclum_stage', ax=axes[0], show=False,
                       title='Cyclum phase', palette=PHASE_COLORS,
                       frameon=False, legend_loc='on data', legend_fontsize=8)

            # Remaining panels: expression of top DE genes
            for i, gene in enumerate(genes):
                ax = axes[i + 1]
                gene_display = gene_type_label.get(gene, FBGN_TO_SYMBOL.get(gene, gene))
                # Flag if this is a known cell cycle gene
                cc_flag = ' ★' if gene in known_cc_genes else ''
                sc.pl.umap(adata, color=gene, ax=ax, show=False,
                           title=f'{gene_display}{cc_flag}',
                           cmap='viridis', frameon=False)

            # Hide unused axes
            for ax in axes[n_genes + 1:]:
                ax.set_visible(False)

            fig.suptitle(
                f"Top DE genes for cyclum phase '{phase}' — {sample_name}\n"
                f"(★ = known Drosophila cell cycle gene; "
                f"ranked by Wilcoxon score vs other phases)",
                fontsize=11, y=1.01,
            )
            plt.tight_layout()
            plt.savefig(os.path.join(output_dir,
                                     f'{sample_name}_validation_umap_{phase}.pdf'),
                        dpi=300, bbox_inches='tight')
            plt.close()
            print(f"  UMAP grid saved for phase '{phase}'")

    else:
        print("  Skipping UMAP plots (no X_umap).")

    # ── Heatmap: top DE genes sorted by cyclum phase ──────────────────────────
    # Collect unique top-n_top_genes per phase (in phase order)
    all_top = []
    for phase in PHASE_ORDER:
        df_phase = de_df[de_df['phase'] == phase].head(n_top_genes)
        all_top += list(df_phase['names'])
    heatmap_genes = list(dict.fromkeys(all_top))   # deduplicate, preserve order
    heatmap_genes = [g for g in heatmap_genes if g in adata.var_names]

    if heatmap_genes:
        # Sort cells: g0/g1 → s → g2/m, then by pseudotime
        sort_key = adata.obs['cyclum_stage'].map({'g0/g1': 0, 's': 1, 'g2/m': 2})
        if 'cyclum_pseudotime' in adata.obs.columns:
            sort_order = adata.obs.assign(_sk=sort_key).sort_values(
                ['_sk', 'cyclum_pseudotime']).index
        else:
            sort_order = adata.obs.assign(_sk=sort_key).sort_values('_sk').index

        # Expression matrix (log-normed)
        X_ln       = _get_lognorm_matrix(adata)
        var_list   = list(adata.var_names)
        gene_idx   = [var_list.index(g) for g in heatmap_genes]
        cell_idx   = [list(adata.obs_names).index(c) for c in sort_order]
        X_sub      = X_ln[np.ix_(cell_idx, gene_idx)]
        means      = X_sub.mean(axis=0)
        stds       = X_sub.std(axis=0) + 1e-10
        X_z        = ((X_sub - means) / stds).T       # genes × cells

        gene_labels = [gene_type_label.get(g, FBGN_TO_SYMBOL.get(g, g))
                       for g in heatmap_genes]

        # Assign gene-row colours by which phase they came from
        phase_of_gene = {}
        for phase in PHASE_ORDER:
            for g in de_df[de_df['phase'] == phase].head(n_top_genes)['names']:
                if g not in phase_of_gene:
                    phase_of_gene[g] = phase
        row_colors = pd.Series(
            [PHASE_COLORS.get(phase_of_gene.get(g, ''), '#CCCCCC')
             for g in heatmap_genes],
            index=gene_labels,
        )

        fig_h = max(6, len(heatmap_genes) * 0.4 + 2)
        fig, ax = plt.subplots(figsize=(12, fig_h))
        sns.heatmap(X_z, cmap='RdBu_r', center=0, vmin=-2, vmax=2,
                    xticklabels=False, yticklabels=gene_labels,
                    cbar_kws={'label': 'Z-score'}, ax=ax)

        # Colour gene labels by phase
        for ytick, gene in zip(ax.get_yticklabels(), heatmap_genes):
            phase = phase_of_gene.get(gene, '')
            ytick.set_color(PHASE_COLORS.get(phase, 'black'))
            if gene in known_cc_genes:
                ytick.set_fontweight('bold')

        ax.set_xlabel('Cells (sorted: g0/g1 → s → g2/m, then by pseudotime)')
        ax.set_ylabel(f'Top {n_top_genes} DE genes per phase')
        ax.set_title(
            f'Top DE genes by cyclum phase — {sample_name}\n'
            f'Colour = source phase  |  Bold = known Drosophila CC gene')

        # Phase colour strip
        phase_strip = np.array([[
            list(plt.cm.colors.to_rgb(PHASE_COLORS.get(p, '#CCCCCC')))
            for p in adata.obs.loc[sort_order, 'cyclum_stage']
        ]])
        ax2 = ax.inset_axes([0, -0.03, 1, 0.02])
        ax2.imshow(phase_strip, aspect='auto')
        ax2.axis('off')
        from matplotlib.patches import Patch
        ax2.legend(
            handles=[Patch(facecolor=c, label=p) for p, c in PHASE_COLORS.items()],
            loc='lower right', bbox_to_anchor=(1, -2), ncol=3,
            fontsize=8, frameon=True, title='Cyclum phase',
        )

        plt.tight_layout()
        plt.savefig(os.path.join(output_dir,
                                 f'{sample_name}_validation_heatmap.pdf'),
                    dpi=300, bbox_inches='tight')
        plt.close()
        print(f"  Heatmap saved ({len(heatmap_genes)} genes × {len(sort_order)} cells)")

    print(f"\n  Validation outputs saved to: {output_dir}")
    return de_df


# ══════════════════════════════════════════════════════════════════════════════
# SECTION 3 – CELL CYCLE GENE EXPRESSION ACROSS LEIDEN CLUSTERS
# ══════════════════════════════════════════════════════════════════════════════

def analyze_cc_genes_by_cluster(adata, output_dir, sample_name):
    """
    Test which FlyBase cell cycle genes are differentially expressed between
    Leiden clusters, and visualise their expression.

    Strategy
    --------
    - Restrict to the FlyBase S-phase and G2/M marker genes present in the dataset
    - Run Kruskal-Wallis per gene across clusters (non-parametric, no normality
      assumption), Bonferroni-corrected across the number of genes tested
    - Rank by H-statistic (effect size proxy) to find the most cluster-variable
      cell cycle genes
    - Plot: dot plot, heatmap, and UMAPs for the top significant genes

    Outputs
    -------
    {sample}_cc_cluster_stats.csv         per-gene KW results, ranked by H
    {sample}_cc_cluster_dotplot.pdf       dot plot: mean expression + % expressing
    {sample}_cc_cluster_heatmap.pdf       Z-score heatmap across clusters
    {sample}_cc_cluster_umap_{gene}.pdf   UMAP per top gene (top n_umap_genes)
    {sample}_cc_cluster_umap_grid.pdf     all top genes in one grid
    """
    print("\n" + "=" * 60)
    print("CELL CYCLE GENES: DIFFERENTIAL EXPRESSION ACROSS CLUSTERS")
    print("=" * 60)

    try:
        leiden_col = _get_leiden_col(adata)
    except KeyError as e:
        print(f"  ERROR: {e}")
        return None

    os.makedirs(output_dir, exist_ok=True)
    sc.settings.figdir = output_dir

    # ── Find which FlyBase CC genes are in the dataset ────────────────────────
    s_genes, g2m_genes = _find_marker_genes(adata, verbose=True)
    cc_genes = s_genes + g2m_genes

    if len(cc_genes) == 0:
        print("  ERROR: No FlyBase cell cycle genes found in dataset.")
        return None

    print(f"\n  Testing {len(cc_genes)} cell cycle genes across "
          f"{adata.obs[leiden_col].nunique()} clusters...")

    # ── Log-normalise for scoring ─────────────────────────────────────────────
    X_ln   = _get_lognorm_matrix(adata)
    var_list = list(adata.var_names)
    clusters = sorted(adata.obs[leiden_col].unique())

    # ── Kruskal-Wallis per gene across clusters ───────────────────────────────
    bonf_thr = 0.05 / len(cc_genes)
    rows = []
    for gene in cc_genes:
        gene_idx = var_list.index(gene)
        expr     = X_ln[:, gene_idx]
        groups   = [expr[adata.obs[leiden_col].values == c] for c in clusters]
        groups   = [g for g in groups if len(g) >= 3]
        if len(groups) < 2:
            continue
        h, p = kruskal(*groups)
        rows.append({
            'gene':             gene,
            'gene_symbol':      FBGN_TO_SYMBOL.get(gene, gene),
            'gene_type':        'S' if gene in s_genes else 'G2M',
            'KW_H':             h,
            'KW_p':             p,
            'bonf_significant': p < bonf_thr,
            'mean_expr':        float(expr.mean()),
            'pct_expressing':   float((expr > 0).mean() * 100),
        })

    stats_df = pd.DataFrame(rows).sort_values('KW_H', ascending=False)

    n_sig = stats_df['bonf_significant'].sum()
    print(f"\n  Bonferroni threshold: p < {bonf_thr:.2e}")
    print(f"  Significant genes: {n_sig}/{len(stats_df)}")
    print(f"\n  Top 10 most variable cell cycle genes across clusters:")
    print(stats_df[['gene_symbol', 'gene_type', 'KW_H', 'KW_p',
                     'bonf_significant']].head(10).to_string(index=False))

    stats_df.to_csv(
        os.path.join(output_dir, f'{sample_name}_cc_cluster_stats.csv'),
        index=False)

    # Use top significant genes for plots; fall back to top by H if none are sig
    plot_genes_df = (stats_df[stats_df['bonf_significant']]
                     if n_sig > 0 else stats_df.head(15))
    plot_genes = list(plot_genes_df['gene'].values)
    plot_symbols = list(plot_genes_df['gene_symbol'].values)

    if not plot_genes:
        print("  No genes to plot.")
        return stats_df

    # ── Dot plot: mean expression + % expressing per cluster ─────────────────
    # sc.pl.dotplot needs gene names in adata.var_names
    fig_title = (f"Cell cycle gene expression by cluster — {sample_name}\n"
                 f"({n_sig} Bonferroni-significant of {len(stats_df)} tested, "
                 f"ranked by Kruskal-Wallis H)")
    try:
        dp = sc.pl.dotplot(
            adata,
            var_names={
                'S-phase': [g for g in plot_genes if g in s_genes],
                'G2/M':    [g for g in plot_genes if g in g2m_genes],
            },
            groupby=leiden_col,
            gene_symbols=None,           # var_names are already IDs
            use_raw=False,
            show=False,
            return_fig=True,
            title=fig_title,
            var_group_rotation=0,
        )
        # Rename tick labels from FBgn IDs to symbols for readability
        ax_main = dp.get_axes()['mainplot_ax']
        ax_main.set_xticklabels(
            [FBGN_TO_SYMBOL.get(t.get_text(), t.get_text())
             for t in ax_main.get_xticklabels()],
            rotation=45, ha='right', fontsize=8,
        )
        dp.savefig(os.path.join(output_dir,
                                f'{sample_name}_cc_cluster_dotplot.pdf'),
                   bbox_inches='tight', dpi=300)
        plt.close()
        print("  Dot plot saved.")
    except Exception as e:
        print(f"  WARNING: dot plot failed ({e}). Skipping.")

    # ── Heatmap: mean Z-score per cluster for each CC gene ────────────────────
    # Build cluster × gene mean-expression matrix
    cluster_mean = pd.DataFrame(index=clusters, columns=plot_genes, dtype=float)
    for gene in plot_genes:
        gene_idx = var_list.index(gene)
        expr     = X_ln[:, gene_idx]
        for c in clusters:
            mask = adata.obs[leiden_col].values == c
            cluster_mean.loc[c, gene] = float(expr[mask].mean())

    # Z-score across clusters (per gene)
    cluster_mean_z = cluster_mean.apply(
        lambda col: (col - col.mean()) / (col.std() + 1e-10), axis=0)
    cluster_mean_z.columns = [FBGN_TO_SYMBOL.get(g, g) for g in plot_genes]

    # Annotate column headers with gene type
    col_colors = pd.Series(
        ['#4ECDC4' if g in s_genes else '#45B7D1' for g in plot_genes],
        index=cluster_mean_z.columns,
    )

    fig_w = max(8, len(plot_genes) * 0.55 + 2)
    fig_h = max(4, len(clusters) * 0.5 + 2)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    sns.heatmap(cluster_mean_z, cmap='RdBu_r', center=0, vmin=-2, vmax=2,
                annot=False, linewidths=0.4,
                cbar_kws={'label': 'Z-score (across clusters)'}, ax=ax)
    ax.set_xlabel('Cell cycle gene')
    ax.set_ylabel(f'Leiden cluster ({leiden_col})')
    ax.set_title(
        f'Mean CC gene expression by cluster (Z-scored) — {sample_name}\n'
        f'Blue ticks = S-phase  |  Teal ticks = G2/M')
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right', fontsize=8)

    # Colour x-axis tick labels by gene type
    for tick, gene in zip(ax.get_xticklabels(), plot_genes):
        tick.set_color('#4ECDC4' if gene in s_genes else '#45B7D1')
        if stats_df.loc[stats_df['gene'] == gene, 'bonf_significant'].values[0]:
            tick.set_fontweight('bold')

    plt.tight_layout()
    plt.savefig(os.path.join(output_dir,
                             f'{sample_name}_cc_cluster_heatmap.pdf'),
                dpi=300, bbox_inches='tight')
    plt.close()
    print("  Cluster heatmap saved.")

    # ── UMAP grid: top genes coloured by expression ───────────────────────────
    if 'X_umap' in adata.obsm:
        # All top genes in one grid
        n_genes = len(plot_genes)
        n_cols  = min(5, n_genes + 1)
        n_rows  = int(np.ceil((n_genes + 1) / n_cols))
        fig, axes = plt.subplots(n_rows, n_cols,
                                  figsize=(5 * n_cols, 4.5 * n_rows))
        axes = np.array(axes).flatten()

        # First panel: cluster labels
        sc.pl.umap(adata, color=leiden_col, ax=axes[0], show=False,
                   title=f'Clusters ({leiden_col})', frameon=False,
                   legend_loc='on data', legend_fontsize=8)

        for i, (gene, symbol) in enumerate(zip(plot_genes, plot_symbols)):
            ax  = axes[i + 1]
            h   = stats_df.loc[stats_df['gene'] == gene, 'KW_H'].values[0]
            sig = stats_df.loc[stats_df['gene'] == gene, 'bonf_significant'].values[0]
            gtype = 'S' if gene in s_genes else 'G2M'
            title = f'{symbol} ({gtype})\nH={h:.1f}{"*" if sig else ""}'
            sc.pl.umap(adata, color=gene, ax=ax, show=False,
                       title=title, cmap='viridis', frameon=False)

        for ax in axes[n_genes + 1:]:
            ax.set_visible(False)

        fig.suptitle(
            f'Cell cycle gene expression by cluster — {sample_name}\n'
            f'* = Bonferroni-significant (p < {bonf_thr:.1e})  |  '
            f'Ranked by Kruskal-Wallis H',
            fontsize=11, y=1.01,
        )
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir,
                                 f'{sample_name}_cc_cluster_umap_grid.pdf'),
                    dpi=300, bbox_inches='tight')
        plt.close()
        print(f"  UMAP grid saved ({n_genes} genes).")

        # Individual per-gene UMAPs
        umap_dir = os.path.join(output_dir, f'{sample_name}_umap_per_gene')
        os.makedirs(umap_dir, exist_ok=True)

        # Reference page: clusters + cyclum phase side by side
        fig, axes = plt.subplots(1, 2, figsize=(14, 6))
        sc.pl.umap(adata, color=leiden_col, ax=axes[0], show=False,
                   title=f'Leiden clusters ({leiden_col})', frameon=False,
                   legend_loc='on data')
        sc.pl.umap(adata, color='cyclum_stage', ax=axes[1], show=False,
                   title='Cyclum phase', frameon=False,
                   palette=PHASE_COLORS)
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir,
                                 f'{sample_name}_cc_cluster_umap_reference.pdf'),
                    dpi=300, bbox_inches='tight')
        plt.close()

        # One file per gene, named by gene symbol
        for gene, symbol in zip(plot_genes, plot_symbols):
            h     = stats_df.loc[stats_df['gene'] == gene, 'KW_H'].values[0]
            p     = stats_df.loc[stats_df['gene'] == gene, 'KW_p'].values[0]
            sig   = stats_df.loc[stats_df['gene'] == gene, 'bonf_significant'].values[0]
            gtype = 'S-phase' if gene in s_genes else 'G2/M'
            title = (f'{symbol}  [{gtype}]\n'
                     f'KW H={h:.2f}, p={p:.2e}'
                     f'{"  *Bonferroni-sig" if sig else ""}')

            fig, ax = plt.subplots(figsize=(7, 6))
            sc.pl.umap(adata, color=gene, ax=ax, show=False,
                       title=title, cmap='viridis', frameon=False)
            plt.tight_layout()

            safe_symbol = symbol.replace('/', '-').replace(' ', '_')
            plt.savefig(os.path.join(umap_dir,
                                     f'{sample_name}_umap_{safe_symbol}.pdf'),
                        dpi=300, bbox_inches='tight')
            plt.close()

        print(f"  Individual gene UMAPs saved to: {umap_dir}/")

    print(f"\n  Cell cycle cluster outputs saved to: {output_dir}")
    return stats_df


# ══════════════════════════════════════════════════════════════════════════════
# SECTION 4 – LEIDEN CLUSTER ~ CYCLUM PHASE ASSOCIATION
# ══════════════════════════════════════════════════════════════════════════════

def analyze_cluster_cellcycle_association(adata, output_dir, sample_name):
    """
    Test and visualise association between Leiden clusters and cyclum phase.
    S_score / G2M_score from the validation step are used as continuous measures.
    """
    from scipy.stats import chi2_contingency

    print("\n" + "=" * 60)
    print("LEIDEN CLUSTER ~ CYCLUM PHASE ASSOCIATION")
    print("=" * 60)

    try:
        leiden_col = _get_leiden_col(adata)
    except KeyError as e:
        print(f"  ERROR: {e}")
        return None

    if 'cyclum_stage' not in adata.obs.columns:
        print("  ERROR: 'cyclum_stage' not in adata.obs.")
        return None

    os.makedirs(output_dir, exist_ok=True)
    sc.settings.figdir = output_dir

    clusters      = sorted(adata.obs[leiden_col].unique())
    cmap          = plt.cm.get_cmap('tab20')
    leiden_colors = [cmap(i % 20) for i in range(len(clusters))]

    # ── Contingency table ─────────────────────────────────────────────────────
    contingency = pd.crosstab(adata.obs[leiden_col], adata.obs['cyclum_stage'])
    ordered_phases  = [p for p in PHASE_ORDER if p in contingency.columns]
    contingency     = contingency[ordered_phases]
    contingency_pct = contingency.div(contingency.sum(axis=1), axis=0) * 100

    # ── Chi-square ────────────────────────────────────────────────────────────
    chi2, p_value, dof, _ = chi2_contingency(contingency)
    n = contingency.sum().sum()
    cramers_v = np.sqrt(chi2 / (n * (min(contingency.shape) - 1)))
    if   cramers_v < 0.1: effect = "negligible"
    elif cramers_v < 0.3: effect = "weak"
    elif cramers_v < 0.5: effect = "moderate"
    else:                  effect = "strong"

    print(f"\n  Chi-square: chi2={chi2:.2f}, dof={dof}, p={p_value:.2e}")
    print(f"  Cramer's V: {cramers_v:.3f} ({effect} effect size)")
    print(f"  Clusters are "
          f"{'SIGNIFICANTLY' if p_value < 0.05 else 'NOT significantly'} "
          f"associated with cyclum phase (α=0.05)")

    print(f"\n  Phase distribution by cluster (%):")
    print(contingency_pct.round(1).to_string())

    # ── Kruskal-Wallis on continuous marker scores (only if present) ──────────
    kw_results = {}
    for score in ['S_score', 'G2M_score']:
        if score not in adata.obs.columns:
            continue
        groups = [adata.obs.loc[adata.obs[leiden_col] == c, score].dropna().values
                  for c in clusters]
        groups = [g for g in groups if len(g) >= 3]
        if len(groups) < 2:
            continue
        h, p = kruskal(*groups)
        kw_results[score] = (h, p)
        print(f"  Kruskal-Wallis {score}: H={h:.2f}, p={p:.2e} "
              f"({'sig' if p < 0.05 else 'ns'})")

    # ── Dominant phase per cluster ────────────────────────────────────────────
    dominant_phase = contingency_pct.idxmax(axis=1)
    max_pct        = contingency_pct.max(axis=1)
    print(f"\n  {'Cluster':<10} {'Dominant phase':<14} {'%':>6}  Status")
    print("  " + "-" * 50)
    for c in clusters:
        pct    = max_pct[c]
        status = ("STRONGLY ENRICHED" if pct > 50
                  else "ENRICHED"      if pct > 40
                  else "Mixed")
        print(f"  {str(c):<10} {dominant_phase[c]:<14} {pct:>6.1f}%  {status}")

    # ── FIGURES ───────────────────────────────────────────────────────────────

    # a) Phase-% heatmap
    fig, ax = plt.subplots(figsize=(10, max(5, len(clusters) * 0.5 + 2)))
    sns.heatmap(contingency_pct, annot=True, fmt='.1f', cmap='YlOrRd', ax=ax,
                linewidths=0.5, cbar_kws={'label': '% of cells in cluster'})
    ax.set_xlabel('Cyclum phase')
    ax.set_ylabel(f'Leiden cluster ({leiden_col})')
    ax.set_title(f"Cyclum phase distribution by cluster\n"
                 f"chi2={chi2:.2f}, p={p_value:.2e}, Cramer's V={cramers_v:.3f}")
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, f'{sample_name}_cluster_phase_heatmap.pdf'),
                dpi=300, bbox_inches='tight')
    plt.close()

    # b) Stacked bar
    fig, ax = plt.subplots(figsize=(max(10, len(clusters) * 0.9 + 4), 6))
    contingency_pct[ordered_phases].plot(
        kind='bar', stacked=True, ax=ax, width=0.8,
        color=[PHASE_COLORS[p] for p in ordered_phases])
    ax.set_xlabel(f'Leiden cluster ({leiden_col})', fontsize=12)
    ax.set_ylabel('% of cells', fontsize=12)
    ax.set_title(f'Cyclum phase composition by cluster\n'
                 f"chi2={chi2:.2f}, p={p_value:.2e}, Cramer's V={cramers_v:.3f}",
                 fontsize=13)
    ax.legend(title='Cyclum phase', bbox_to_anchor=(1.02, 1), loc='upper left')
    plt.xticks(rotation=0)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, f'{sample_name}_cluster_phase_stacked.pdf'),
                dpi=300, bbox_inches='tight')
    plt.close()

    # c) Violin + bar of marker scores by cluster
    if kw_results:
        n_scores = len(kw_results)
        fig, axes = plt.subplots(2, n_scores,
                                  figsize=(max(14, len(clusters) * 1.2) * n_scores // 2, 10))
        if n_scores == 1:
            axes = axes.reshape(2, 1)

        score_stats = adata.obs.groupby(leiden_col)[list(kw_results.keys())].agg(['mean', 'std'])

        for col_i, score in enumerate(kw_results):
            h_kw, p_kw = kw_results[score]
            # Violin
            sc.pl.violin(adata, score, groupby=leiden_col,
                         ax=axes[0, col_i], show=False, rotation=0)
            axes[0, col_i].set_title(
                f'{score} by cluster\nKW H={h_kw:.2f}, p={p_kw:.2e}')
            axes[0, col_i].axhline(0, color='k', linestyle='--', alpha=0.3)

            # Bar (mean ± SD)
            means = score_stats[score]['mean']
            stds  = score_stats[score]['std']
            axes[1, col_i].bar(range(len(clusters)), means, yerr=stds,
                                capsize=4, color=leiden_colors,
                                alpha=0.8, edgecolor='black')
            axes[1, col_i].set_xticks(range(len(clusters)))
            axes[1, col_i].set_xticklabels(clusters)
            axes[1, col_i].set_xlabel(f'Leiden cluster ({leiden_col})')
            axes[1, col_i].set_ylabel(f'Mean {score} ± SD')
            axes[1, col_i].set_title(f'Mean {score} per cluster')
            axes[1, col_i].axhline(0, color='k', linestyle='--', alpha=0.3)

        plt.suptitle(f'Marker scores by cluster — {sample_name}', fontsize=12)
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, f'{sample_name}_cluster_scores.pdf'),
                    dpi=300, bbox_inches='tight')
        plt.close()

    # d) UMAPs
    if 'X_umap' in adata.obsm:
        cols_to_plot = [leiden_col, 'cyclum_stage']
        if 'S_score' in adata.obs.columns:
            cols_to_plot += ['cyclum_pseudotime', 'S_score', 'G2M_score']

        fig, axes = plt.subplots(1, 2, figsize=(16, 6))
        sc.pl.umap(adata, color=leiden_col, ax=axes[0], show=False,
                   title=f'Leiden ({leiden_col})', frameon=False,
                   legend_loc='on data')
        sc.pl.umap(adata, color='cyclum_stage', ax=axes[1], show=False,
                   title='Cyclum phase', frameon=False,
                   palette=PHASE_COLORS)
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir,
                                 f'{sample_name}_umap_cluster_vs_phase.pdf'),
                    dpi=300, bbox_inches='tight')
        plt.close()

        if 'S_score' in adata.obs.columns:
            sc.pl.umap(adata,
                       color=['cyclum_pseudotime', 'S_score', 'G2M_score'],
                       save=f'_{sample_name}_continuous_scores.pdf',
                       cmap='viridis', show=False)

    # ── Save CSVs ──────────────────────────────────────────────────────────────
    contingency.to_csv(
        os.path.join(output_dir, f'{sample_name}_contingency_counts.csv'))
    contingency_pct.to_csv(
        os.path.join(output_dir, f'{sample_name}_contingency_pct.csv'))

    summary_df = pd.DataFrame({
        'Cluster':        clusters,
        'N_Cells':        [contingency.loc[c].sum() for c in clusters],
        'Dominant_Phase': [dominant_phase[c] for c in clusters],
        'Phase_Pct':      [max_pct[c] for c in clusters],
        **{f'Mean_{s}': [adata.obs.loc[adata.obs[leiden_col] == c, s].mean()
                          if s in adata.obs.columns else np.nan
                          for c in clusters]
           for s in ['S_score', 'G2M_score']},
    })
    summary_df.to_csv(
        os.path.join(output_dir, f'{sample_name}_cluster_summary.csv'), index=False)

    stats_row = {
        'chi2':      chi2,
        'dof':       dof,
        'p_value':   p_value,
        'cramers_v': cramers_v,
        'effect':    effect,
    }
    for s in ['S_score', 'G2M_score']:
        stats_row[f'kw_{s}_H'] = kw_results.get(s, (np.nan,))[0]
        stats_row[f'kw_{s}_p'] = kw_results.get(s, (np.nan, np.nan))[1]

    pd.DataFrame([stats_row]).to_csv(
        os.path.join(output_dir, f'{sample_name}_cluster_stats.csv'), index=False)

    print(f"\n  Cluster association outputs saved to: {output_dir}")

    return {
        'chi2':       chi2,
        'p_value':    p_value,
        'cramers_v':  cramers_v,
        'kw_results': kw_results,
        'summary':    summary_df,
    }


# ══════════════════════════════════════════════════════════════════════════════
# MAIN
# ══════════════════════════════════════════════════════════════════════════════

def main():
    parser = argparse.ArgumentParser(
        description='Cyclum cell cycle analysis with FlyBase marker validation '
                    'and Leiden cluster association testing',
    )
    parser.add_argument('--input',  '-i', required=False,
                        default='/private/groups/russelllab/jodie/scRNAseq/scripts/'
                                'snakemake_pipeline/results_kallisto_bustools/'
                                'filtered_h5ad/kallisto_JW18DOX-Ctrl-1_P.h5ad',
                        help='Input h5ad file')
    parser.add_argument('--output', '-o', required=False,
                        default='/private/groups/russelllab/jodie/scRNAseq/scripts/'
                                'snakemake_pipeline/results_kallisto_bustools/'
                                'filtered_h5ad/cyclum_JW18DOX-Ctrl-1_P',
                        help='Output directory')
    parser.add_argument('--sample', '-s', default='sample',
                        help='Prefix label for all output filenames')
    parser.add_argument('--epochs',      type=int,   default=800,
                        help='Cyclum training epochs (default: 800)')
    parser.add_argument('--rate',        type=float, default=2e-4,
                        help='Cyclum learning rate (default: 2e-4)')
    parser.add_argument('--n-top-genes',  type=int, default=5,
                        help='Top DE genes per phase saved to CSV / shown in '
                             'heatmap (default: 5)')
    parser.add_argument('--n-umap-genes', type=int, default=6,
                        help='Top DE genes per phase shown in UMAP grid '
                             '(default: 6)')
    parser.add_argument('--skip-cyclum', action='store_true',
                        help='Skip cyclum training; use cyclum_stage already in h5ad')
    parser.add_argument('--save-h5ad',   action='store_true',
                        help='Save annotated h5ad to output directory')

    args = parser.parse_args()
    os.makedirs(args.output, exist_ok=True)
    sc.settings.figdir = args.output

    # ── Load ──────────────────────────────────────────────────────────────────
    print(f"Loading {args.input}...")
    adata = sc.read_h5ad(args.input)
    print(f"  {adata.n_obs} cells × {adata.n_vars} genes")
    print(f"  obs columns: {list(adata.obs.columns)}")

    mtx = adata.X
    if scipy.sparse.issparse(mtx):
        mtx = mtx.toarray()

    # ── STEP 1: Cyclum ─────────────────────────────────────────────────────────
    if args.skip_cyclum and 'cyclum_stage' in adata.obs.columns:
        print("\n[Cyclum] Skipping training — using existing cyclum_stage.")
    else:
        print("\n[Cyclum] Training model...")
        model = cyclum.tuning.CyclumAutoTune(mtx)
        model.train(mtx, epochs=args.epochs, verbose=100, rate=args.rate)

        pseudotime      = model.predict_pseudotime(mtx)
        pseudotime_flat = pseudotime.flatten()
        print(f"  Pseudotime range: {pseudotime_flat.min():.3f} – {pseudotime_flat.max():.3f}")

        stages, confidence, angles = assign_cell_cycle_stage_simple(pseudotime_flat)

        adata.obs['cyclum_stage']      = stages
        adata.obs['cyclum_pseudotime'] = pseudotime_flat
        adata.obs['cyclum_confidence'] = confidence

        # Cyclum diagnostic plots
        color_map = {'g0/g1': 'red', 's': 'green', 'g2/m': 'blue'}
        fig = cyclum.illustration.plot_round_distr_color(
            pseudotime_flat, np.array(stages), color_map)
        plt.savefig(os.path.join(args.output, f'{args.sample}_cyclum_circular.pdf'),
                    dpi=300, bbox_inches='tight')
        plt.close()

        model.show_elbow()
        plt.savefig(os.path.join(args.output, f'{args.sample}_cyclum_elbow.pdf'),
                    dpi=300, bbox_inches='tight')
        plt.close()

        model.show_bar()
        plt.savefig(os.path.join(args.output, f'{args.sample}_cyclum_bar.pdf'),
                    dpi=300, bbox_inches='tight')
        plt.close()

        # Diagnostics panel
        fig, axes = plt.subplots(1, 3, figsize=(14, 4))
        axes[0].hist(confidence, bins=30, alpha=0.7, edgecolor='black')
        axes[0].set(xlabel='Confidence', ylabel='Cells',
                    title='Assignment confidence')
        axes[1].scatter(pseudotime_flat, confidence,
                        c=[color_map[s] for s in stages], alpha=0.5, s=8)
        axes[1].set(xlabel='Pseudotime', ylabel='Confidence',
                    title='Confidence vs Pseudotime')
        axes[2].hist(angles, bins=60, alpha=0.7, edgecolor='black')
        axes[2].axvline(2*np.pi/3, color='red',   linestyle='--',
                        alpha=0.7, label='G1/S')
        axes[2].axvline(4*np.pi/3, color='green', linestyle='--',
                        alpha=0.7, label='S/G2M')
        axes[2].set(xlabel='Angle (rad)', ylabel='Cells',
                    title='Circular distribution')
        axes[2].legend(fontsize=8)
        plt.tight_layout()
        plt.savefig(os.path.join(args.output,
                                 f'{args.sample}_cyclum_diagnostics.pdf'),
                    dpi=300, bbox_inches='tight')
        plt.close()

    # ── STEP 2: Validate cyclum phases with marker gene expression ────────────
    validate_cyclum_phases(adata, args.output, args.sample,
                           n_top_genes=args.n_top_genes,
                           n_umap_genes=args.n_umap_genes)

    # ── STEP 3: Which CC genes differ between clusters? ───────────────────────
    cc_cluster_results = analyze_cc_genes_by_cluster(
        adata, args.output, args.sample)

    # ── STEP 4: Leiden cluster ~ cyclum phase association ─────────────────────
    cluster_results = analyze_cluster_cellcycle_association(
        adata, args.output, args.sample)

    # ── STEP 5: Save h5ad ────────────────────────────────────────────────────
    if args.save_h5ad:
        out_h5ad = os.path.join(args.output, f'{args.sample}_cyclum_annotated.h5ad')
        adata.write_h5ad(out_h5ad)
        print(f"\nSaved annotated h5ad: {out_h5ad}")

    # ── Final summary ──────────────────────────────────────────────────────────
    print("\n" + "=" * 60)
    print("ANALYSIS COMPLETE — KEY FINDINGS")
    print("=" * 60)

    print("\nCyclum phase distribution:")
    for phase in PHASE_ORDER:
        n = (adata.obs['cyclum_stage'] == phase).sum()
        print(f"  {phase}: {n}  ({n / adata.n_obs * 100:.1f}%)")

    if 'S_score' in adata.obs.columns:
        print("\nMarker score means per cyclum phase (validation):")
        pm = adata.obs.groupby('cyclum_stage')[['S_score', 'G2M_score']].mean()
        print(pm.reindex([p for p in PHASE_ORDER if p in pm.index]).round(4))

    if cc_cluster_results is not None:
        n_sig = cc_cluster_results['bonf_significant'].sum()
        print(f"\nCell cycle genes significant across clusters: "
              f"{n_sig}/{len(cc_cluster_results)}")
        if n_sig > 0:
            top = cc_cluster_results[cc_cluster_results['bonf_significant']].head(5)
            print("  Top genes:")
            for _, r in top.iterrows():
                print(f"    {r['gene_symbol']} ({r['gene_type']})  "
                      f"H={r['KW_H']:.2f}, p={r['KW_p']:.2e}")


        print(f"\nCluster ~ cyclum phase:")
        print(f"  chi2={cluster_results['chi2']:.2f}, "
              f"p={cluster_results['p_value']:.2e}, "
              f"Cramer's V={cluster_results['cramers_v']:.3f}")

    print(f"\nAll outputs in: {args.output}")


if __name__ == '__main__':
    main()