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

# Expected enrichment direction for validation
EXPECTED_HIGH = {'S_score': 's', 'G2M_score': 'g2/m'}


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


def _mean_gene_set_score(X_lognorm, gene_indices):
    """Mean log-normalised expression across a gene set. Returns shape (n_cells,)."""
    if not gene_indices:
        return np.zeros(X_lognorm.shape[0])
    return X_lognorm[:, gene_indices].mean(axis=1)


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
# SECTION 2 – MARKER GENE VALIDATION OF CYCLUM PHASES
# ══════════════════════════════════════════════════════════════════════════════

def validate_cyclum_phases(adata, output_dir, sample_name):
    """
    Validate cyclum phase assignments using Drosophila cell cycle marker genes.

    Strategy
    --------
    Per-cell S_score and G2M_score are computed as the mean log-normalised
    expression of S-phase and G2/M-phase FlyBase marker genes respectively.
    These scores are stored in adata.obs but are NOT used to re-call phases —
    cyclum_stage remains the authoritative label.

    A good cyclum assignment should show:
      S_score   highest in cells assigned to 's'
      G2M_score highest in cells assigned to 'g2/m'
      Both      lowest  in cells assigned to 'g0/g1'

    Statistical tests (Kruskal-Wallis + pairwise Mann-Whitney U with Bonferroni
    correction) formally test whether scores differ significantly across phases.

    Outputs
    -------
    {sample}_validation_score_violin.pdf      score distributions per phase
    {sample}_validation_score_scatter.pdf     S vs G2M score scatter
    {sample}_validation_pseudotime_scores.pdf scores along pseudotime
    {sample}_validation_marker_heatmap.pdf    per-gene Z-score heatmap
    {sample}_validation_phase_means.csv       mean score per phase
    {sample}_validation_stats.csv             KW + pairwise test results
    """
    print("\n" + "=" * 60)
    print("MARKER GENE VALIDATION OF CYCLUM PHASES")
    print("=" * 60)

    if 'cyclum_stage' not in adata.obs.columns:
        print("  ERROR: 'cyclum_stage' not in adata.obs. Run cyclum first.")
        return None

    os.makedirs(output_dir, exist_ok=True)

    # ── Find marker genes ─────────────────────────────────────────────────────
    s_genes, g2m_genes = _find_marker_genes(adata)

    if len(s_genes) < 2 and len(g2m_genes) < 2:
        print("  ERROR: Fewer than 2 marker genes found. "
              "Check that var_names use FBgn IDs or matching gene symbols.")
        return None

    # ── Compute scores ────────────────────────────────────────────────────────
    X_ln    = _get_lognorm_matrix(adata)
    var_idx = list(adata.var_names)

    s_idx   = [var_idx.index(g) for g in s_genes   if g in var_idx]
    g2m_idx = [var_idx.index(g) for g in g2m_genes if g in var_idx]

    adata.obs['S_score']   = _mean_gene_set_score(X_ln, s_idx)
    adata.obs['G2M_score'] = _mean_gene_set_score(X_ln, g2m_idx)

    phases = adata.obs['cyclum_stage']

    # ── Mean scores per phase ─────────────────────────────────────────────────
    phase_means = adata.obs.groupby('cyclum_stage')[['S_score', 'G2M_score']].mean()
    # Reindex to canonical order
    phase_means = phase_means.reindex([p for p in PHASE_ORDER
                                        if p in phase_means.index])
    print("\n  Mean marker scores per cyclum phase:")
    print(phase_means.round(4).to_string())

    print("\n  Enrichment direction check (validation):")
    for score, expected in EXPECTED_HIGH.items():
        actual = phase_means[score].idxmax()
        ok = "PASS ✓" if actual == expected else "FAIL ✗ — check phase boundaries"
        print(f"    {score}: highest in '{actual}' (expected '{expected}') — {ok}")

    # ── Statistical tests ─────────────────────────────────────────────────────
    print("\n  Kruskal-Wallis (score ~ cyclum phase):")
    stats_rows = []
    for score in ['S_score', 'G2M_score']:
        phase_groups = [adata.obs.loc[phases == p, score].values
                        for p in PHASE_ORDER if p in phases.unique()]
        phase_groups = [g for g in phase_groups if len(g) >= 3]
        if len(phase_groups) < 2:
            continue
        h, p_kw = kruskal(*phase_groups)
        print(f"    {score}: H={h:.2f}, p={p_kw:.2e} "
              f"({'SIGNIFICANT' if p_kw < 0.05 else 'not significant'})")

        # Pairwise Mann-Whitney U with Bonferroni correction
        pairs    = list(combinations(PHASE_ORDER, 2))
        bonf_thr = 0.05 / len(pairs)
        for pa, pb in pairs:
            if pa not in phases.unique() or pb not in phases.unique():
                continue
            ga = adata.obs.loc[phases == pa, score].values
            gb = adata.obs.loc[phases == pb, score].values
            u, p_mw = mannwhitneyu(ga, gb, alternative='two-sided')
            stats_rows.append({
                'score':            score,
                'phase_A':          pa,
                'phase_B':          pb,
                'KW_H':             h,
                'KW_p':             p_kw,
                'MW_U':             u,
                'MW_p':             p_mw,
                'bonf_significant': p_mw < bonf_thr,
            })

    stats_df = pd.DataFrame(stats_rows)
    if not stats_df.empty:
        sig = stats_df[stats_df['bonf_significant']]
        print(f"\n  Bonferroni-significant pairwise differences "
              f"({sig.shape[0]}/{stats_df.shape[0]} pairs):")
        if sig.empty:
            print("    None")
        else:
            for _, row in sig.iterrows():
                print(f"    {row['score']}: {row['phase_A']} vs {row['phase_B']} "
                      f"U={row['MW_U']:.0f}, p={row['MW_p']:.2e} *")

    # ── Figure 1: violin plots of scores per cyclum phase ─────────────────────
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    for ax, score in zip(axes, ['S_score', 'G2M_score']):
        data_by_phase = [adata.obs.loc[phases == p, score].values
                         for p in PHASE_ORDER if p in phases.unique()]
        present_phases = [p for p in PHASE_ORDER if p in phases.unique()]
        vp = ax.violinplot(data_by_phase, positions=range(len(present_phases)),
                           showmedians=True, showextrema=False)
        for i, body in enumerate(vp['bodies']):
            body.set_facecolor(PHASE_COLORS[present_phases[i]])
            body.set_alpha(0.75)
        vp['cmedians'].set_color('black')
        ax.set_xticks(range(len(present_phases)))
        ax.set_xticklabels(present_phases)
        ax.set_xlabel('Cyclum phase')
        ax.set_ylabel(score)
        ax.set_title(f'{score} by Cyclum phase')

        # Annotate Bonferroni-significant pairs
        if not stats_df.empty:
            y_max   = adata.obs[score].max()
            y_range = adata.obs[score].max() - adata.obs[score].min()
            step    = y_range * 0.08
            for k_off, (_, row) in enumerate(
                stats_df[(stats_df['score'] == score) &
                          stats_df['bonf_significant']].iterrows()
            ):
                if row['phase_A'] not in present_phases or row['phase_B'] not in present_phases:
                    continue
                i = present_phases.index(row['phase_A'])
                j = present_phases.index(row['phase_B'])
                y = y_max + step * (k_off + 1)
                ax.plot([i, j], [y, y], color='black', linewidth=1)
                ax.text((i+j)/2, y + step*0.15, '*', ha='center',
                        fontsize=12, fontweight='bold')

    plt.suptitle(f'Marker gene score validation — {sample_name}\n'
                 f'(S_score = mean S-phase FlyBase markers; '
                 f'G2M_score = mean G2/M FlyBase markers)',
                 fontsize=10)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, f'{sample_name}_validation_score_violin.pdf'),
                dpi=300, bbox_inches='tight')
    plt.close()

    # ── Figure 2: S_score vs G2M_score scatter ────────────────────────────────
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    for phase in PHASE_ORDER:
        if phase not in phases.unique():
            continue
        mask = phases == phase
        axes[0].scatter(
            adata.obs.loc[mask, 'S_score'],
            adata.obs.loc[mask, 'G2M_score'],
            c=PHASE_COLORS[phase], label=phase,
            alpha=0.5, s=8, rasterized=True,
        )
    axes[0].axhline(0, color='k', linestyle='--', alpha=0.3)
    axes[0].axvline(0, color='k', linestyle='--', alpha=0.3)
    axes[0].set_xlabel('S score (FlyBase S-phase markers)')
    axes[0].set_ylabel('G2M score (FlyBase G2/M markers)')
    axes[0].set_title('Marker scores coloured by Cyclum phase')
    axes[0].legend(title='Cyclum phase', fontsize=9)

    if 'cyclum_pseudotime' in adata.obs.columns:
        pt      = adata.obs['cyclum_pseudotime'].values
        pt_norm = (pt - pt.min()) / (pt.max() - pt.min() + 1e-10)
        sc_pt   = axes[1].scatter(
            adata.obs['S_score'], adata.obs['G2M_score'],
            c=pt_norm, cmap='hsv', alpha=0.5, s=8, rasterized=True,
        )
        plt.colorbar(sc_pt, ax=axes[1], label='Cyclum pseudotime (normalised)')
    axes[1].axhline(0, color='k', linestyle='--', alpha=0.3)
    axes[1].axvline(0, color='k', linestyle='--', alpha=0.3)
    axes[1].set_xlabel('S score (FlyBase S-phase markers)')
    axes[1].set_ylabel('G2M score (FlyBase G2/M markers)')
    axes[1].set_title('Marker scores coloured by pseudotime')

    plt.suptitle(f'S score vs G2M score — {sample_name}', fontsize=12)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, f'{sample_name}_validation_score_scatter.pdf'),
                dpi=300, bbox_inches='tight')
    plt.close()

    # ── Figure 3: scores along pseudotime ────────────────────────────────────
    if 'cyclum_pseudotime' in adata.obs.columns:
        pt       = adata.obs['cyclum_pseudotime'].values
        sort_idx = np.argsort(pt)
        window   = max(1, len(pt) // 100)

        fig, axes = plt.subplots(2, 1, figsize=(14, 8), sharex=True)
        for ax, score, color in zip(axes,
                                    ['S_score', 'G2M_score'],
                                    ['#4ECDC4', '#45B7D1']):
            vals     = adata.obs[score].values[sort_idx]
            smoothed = pd.Series(vals).rolling(window, center=True,
                                               min_periods=1).mean().values
            ax.scatter(pt[sort_idx], vals, c=color, alpha=0.2, s=5, rasterized=True)
            ax.plot(pt[sort_idx], smoothed, color='black', linewidth=1.5,
                    label=f'Rolling mean (w={window})')
            ax.set_ylabel(score)
            ax.legend(fontsize=8)
            # Shade cyclum phase regions
            for phase, col in PHASE_COLORS.items():
                mask = adata.obs['cyclum_stage'].values == phase
                if mask.any():
                    pt_phase = pt[mask]
                    ax.axvspan(pt_phase.min(), pt_phase.max(),
                               alpha=0.08, color=col, label=phase)

        axes[1].set_xlabel('Cyclum pseudotime')
        axes[0].set_title(
            f'FlyBase marker scores along pseudotime — {sample_name}\n'
            f'Shaded regions = cyclum phase boundaries')
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir,
                                 f'{sample_name}_validation_pseudotime_scores.pdf'),
                    dpi=300, bbox_inches='tight')
        plt.close()

    # ── Figure 4: per-gene Z-score heatmap sorted by cyclum phase ─────────────
    all_markers = [g for g in (s_genes + g2m_genes) if g in var_idx]
    if all_markers:
        # Sort cells: g0/g1 → s → g2/m then by pseudotime within each phase
        sort_key = phases.map({'g0/g1': 0, 's': 1, 'g2/m': 2})
        if 'cyclum_pseudotime' in adata.obs.columns:
            sort_order = adata.obs.assign(_sk=sort_key).sort_values(
                ['_sk', 'cyclum_pseudotime']).index
        else:
            sort_order = adata.obs.assign(_sk=sort_key).sort_values('_sk').index

        cell_order_idx = [list(adata.obs_names).index(c) for c in sort_order]
        marker_col_idx = [var_idx.index(g) for g in all_markers]

        X_sub = X_ln[np.ix_(cell_order_idx, marker_col_idx)]
        means = X_sub.mean(axis=0)
        stds  = X_sub.std(axis=0) + 1e-10
        X_z   = ((X_sub - means) / stds).T   # genes × cells

        gene_labels   = [FBGN_TO_SYMBOL.get(g, g) for g in all_markers]
        n_s_present   = len([g for g in s_genes   if g in var_idx])
        n_g2m_present = len([g for g in g2m_genes if g in var_idx])

        fig_h = max(6, len(gene_labels) * 0.35 + 2)
        fig, ax = plt.subplots(figsize=(12, fig_h))

        sns.heatmap(X_z, cmap='RdBu_r', center=0, vmin=-2, vmax=2,
                    xticklabels=False, yticklabels=gene_labels,
                    cbar_kws={'label': 'Z-score'}, ax=ax)
        ax.set_xlabel('Cells (sorted: g0/g1 → s → g2/m, then by pseudotime)')
        ax.set_ylabel('Drosophila cell cycle marker genes')
        ax.set_title(f'Marker gene expression by cyclum phase — {sample_name}')

        # Dividing line between S and G2M gene blocks
        if 0 < n_s_present < len(gene_labels):
            ax.axhline(n_s_present, color='black', linewidth=1.5, linestyle='--')
            ax.text(X_z.shape[1] * 1.01, n_s_present / 2,
                    'S genes', va='center', fontsize=8, color='#4ECDC4',
                    fontweight='bold')
            ax.text(X_z.shape[1] * 1.01,
                    n_s_present + n_g2m_present / 2,
                    'G2M genes', va='center', fontsize=8, color='#45B7D1',
                    fontweight='bold')

        # Phase colour strip below heatmap
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
                                 f'{sample_name}_validation_marker_heatmap.pdf'),
                    dpi=300, bbox_inches='tight')
        plt.close()

    # ── Save CSVs ──────────────────────────────────────────────────────────────
    phase_means.to_csv(
        os.path.join(output_dir, f'{sample_name}_validation_phase_means.csv'))
    if not stats_df.empty:
        stats_df.to_csv(
            os.path.join(output_dir, f'{sample_name}_validation_stats.csv'),
            index=False)

    print(f"\n  Validation outputs saved to: {output_dir}")
    return stats_df


# ══════════════════════════════════════════════════════════════════════════════
# SECTION 3 – LEIDEN CLUSTER ~ CYCLUM PHASE ASSOCIATION
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

    # ── Kruskal-Wallis on continuous marker scores ────────────────────────────
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
    validate_cyclum_phases(adata, args.output, args.sample)

    # ── STEP 3: Leiden cluster ~ cyclum phase association ─────────────────────
    cluster_results = analyze_cluster_cellcycle_association(
        adata, args.output, args.sample)

    # ── STEP 4: Save h5ad ────────────────────────────────────────────────────
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

    if cluster_results:
        print(f"\nCluster ~ cyclum phase:")
        print(f"  chi2={cluster_results['chi2']:.2f}, "
              f"p={cluster_results['p_value']:.2e}, "
              f"Cramer's V={cluster_results['cramers_v']:.3f}")

    print(f"\nAll outputs in: {args.output}")


if __name__ == '__main__':
    main()