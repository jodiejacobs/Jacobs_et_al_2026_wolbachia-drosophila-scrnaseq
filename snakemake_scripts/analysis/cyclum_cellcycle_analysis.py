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
#
# Four-phase classification based on observed expression peaks in Drosophila
# scRNA-seq and known biology:
#
#   G0/G1  — arrest/exit markers; high when cells are quiescent or in G1
#   G1/S   — replication licensing factors transcribed in late G1 in
#             preparation for S phase; mRNA peaks BEFORE replication begins
#   S      — active replisome components and E2F targets expressed during
#             ongoing DNA synthesis; peak coincides with replication fork
#   G2/M   — mitotic entry and progression; peak at G2/M boundary
#
# FlyBase IDs verified against FlyBase r6 (2024).
FLYBASE_CELL_CYCLE_GENES = {
    # G0/G1 — quiescence and G1 arrest
    'dap':            'FBgn0010316',   # p21/p27 homolog; inhibits CycE/Cdk2
    'fzr':            'FBgn0262699',   # APC/C activator; degrades mitotic cyclins in G1
    'Rbf':            'FBgn0015799',   # Rb homolog; represses E2F targets in G1
    'Rbf2':           'FBgn0038390',   # Rb family member
    # G1/S — replication licensing; loaded in late G1, mRNA peaks before fork firing
    'Mcm2':           'FBgn0014861',
    'Mcm3':           'FBgn0284442',
    'Mcm5':           'FBgn0017577',
    'Mcm7':           'FBgn0020633',
    'Orc1':           'FBgn0286788',   # origin recognition complex
    'Orc2':           'FBgn0015270',
    'Orc6':           'FBgn0023180',
    'RPA1':           'FBgn0010173',   # ssDNA binding; pre-S loading
    'Pcna':           'FBgn0005655',   # sliding clamp; peaks in late G1
    'DNApol-delta':   'FBgn0019624',   # loaded at origins in G1
    'Dp':             'FBgn0011763',   # E2F partner; G1→S transcription factor
    'E2f2':           'FBgn0024371',   # repressive E2F; active at G1→S boundary
    'RnrS':           'FBgn0011704',   # ribonucleotide reductase small subunit
    'Rrp1':           'FBgn0004584',
    # S phase — active replisome; peak during ongoing DNA synthesis
    'E2f1':           'FBgn0011766',   # activating E2F; drives S-phase transcription
    'CycE':           'FBgn0010382',   # CycE/Cdk2 triggers G1→S; peaks mid-S
    'Cdk2':           'FBgn0004107',
    'Mcm6':           'FBgn0025815',   # MCM subunit that peaks in S (not late G1)
    'RnrL':           'FBgn0011703',   # large RnR subunit; peaks during replication
    'RPA2':           'FBgn0288834',
    'pol-alpha1':     'FBgn0259113',   # primase; active at replication forks
    # G2/M — mitotic entry and progression
    'CycA':           'FBgn0000404',
    'CycB':           'FBgn0000405',
    'CycB3':          'FBgn0015625',
    'Cdk1':           'FBgn0004106',
    'stg':            'FBgn0003525',   # Cdc25; activates Cdk1
    'polo':           'FBgn0003124',
    'aurA':           'FBgn0000147',
    'aurB':           'FBgn0024227',
    'Nek2':           'FBgn0029970',
    'Pbl':            'FBgn0003041',
    'Wee1':           'FBgn0011737',
    'myt':            'FBgn0040298',
    'BubR1':          'FBgn0263855',
    'Mad2':           'FBgn0035640',
    'Cdc20':          'FBgn0001086',
    'APC10':          'FBgn0034231',
}

# G0/G1: quiescence / G1 arrest markers
G0G1_GENES_FBGN = [
    'FBgn0010316',   # dap
    'FBgn0262699',   # fzr
    'FBgn0015799',   # Rbf
    'FBgn0038390',   # Rbf2
]

# G1/S: replication licensing — transcribed in late G1, peak before fork firing
G1S_GENES_FBGN = [
    'FBgn0014861',   # Mcm2
    'FBgn0284442',   # Mcm3
    'FBgn0017577',   # Mcm5
    'FBgn0020633',   # Mcm7
    'FBgn0286788',   # Orc1
    'FBgn0015270',   # Orc2
    'FBgn0023180',   # Orc6
    'FBgn0010173',   # RPA1
    'FBgn0005655',   # Pcna
    'FBgn0019624',   # DNApol-delta
    'FBgn0011763',   # Dp
    'FBgn0024371',   # E2f2
    'FBgn0011704',   # RnrS
    'FBgn0004584',   # Rrp1
]

# S phase: active replisome components; peak during ongoing DNA synthesis
S_GENES_FBGN = [
    'FBgn0011766',   # E2f1
    'FBgn0010382',   # CycE
    'FBgn0004107',   # Cdk2
    'FBgn0025815',   # Mcm6
    'FBgn0011703',   # RnrL
    'FBgn0288834',   # RPA2
    'FBgn0259113',   # pol-alpha1
]

# G2/M: mitotic entry and progression
G2M_GENES_FBGN = [
    'FBgn0000404',   # CycA
    'FBgn0000405',   # CycB
    'FBgn0015625',   # CycB3
    'FBgn0004106',   # Cdk1
    'FBgn0003525',   # stg
    'FBgn0003124',   # polo
    'FBgn0000147',   # aurA
    'FBgn0024227',   # aurB
    'FBgn0029970',   # Nek2
    'FBgn0003041',   # Pbl
    'FBgn0011737',   # Wee1
    'FBgn0040298',   # myt
    'FBgn0263855',   # BubR1
    'FBgn0035640',   # Mad2
    'FBgn0001086',   # Cdc20
    'FBgn0034231',   # APC10
]

# All CC genes in display order: G0/G1 → G1/S → S → G2/M
ALL_CC_GENES_FBGN = G0G1_GENES_FBGN + G1S_GENES_FBGN + S_GENES_FBGN + G2M_GENES_FBGN

FBGN_TO_SYMBOL = {v: k for k, v in FLYBASE_CELL_CYCLE_GENES.items()}

PHASE_ORDER  = ['g0/g1', 's', 'g2/m']
PHASE_COLORS = {'g0/g1': '#FF6B6B', 's': '#4ECDC4', 'g2/m': '#45B7D1'}

# Gene type labels and colours — four-phase classification
GENE_TYPE_COLORS = {
    'G0G1': '#FF6B6B',   # red    — quiescence / G1 arrest
    'G1S':  '#FFA500',   # orange — late G1 / replication licensing
    'S':    '#4ECDC4',   # teal   — active S phase
    'G2M':  '#45B7D1',   # blue   — G2/M
}


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



def _find_marker_genes(adata, verbose=True):
    """
    Return lists of G0G1, G1S, S, and G2M marker gene IDs present in dataset.
    Checks adata.raw first (full gene set), then adata.var_names (HVGs).
    Tries FBgn IDs first, then falls back to gene symbols.

    Four-phase classification:
      G0G1 — quiescence / G1 arrest markers
      G1S  — replication licensing factors; transcribed in late G1
      S    — active replisome components; peak during DNA synthesis
      G2M  — mitotic entry and progression
    """
    var_names = (list(adata.raw.var_names)
                 if adata.raw is not None
                 else list(adata.var_names))
    source = 'raw' if adata.raw is not None else 'X'

    g0g1_present = [g for g in G0G1_GENES_FBGN if g in var_names]
    g1s_present  = [g for g in G1S_GENES_FBGN  if g in var_names]
    s_present    = [g for g in S_GENES_FBGN    if g in var_names]
    g2m_present  = [g for g in G2M_GENES_FBGN  if g in var_names]

    # Fall back to gene symbols if no FBgn IDs found
    if not any([g0g1_present, g1s_present, s_present, g2m_present]):
        for sym, fbgn in FLYBASE_CELL_CYCLE_GENES.items():
            if sym in var_names:
                if fbgn in G0G1_GENES_FBGN:
                    g0g1_present.append(sym)
                elif fbgn in G1S_GENES_FBGN:
                    g1s_present.append(sym)
                elif fbgn in S_GENES_FBGN:
                    s_present.append(sym)
                elif fbgn in G2M_GENES_FBGN:
                    g2m_present.append(sym)

    if verbose:
        print(f"  Searching in adata.{source} ({len(var_names)} genes)")
        for label, present, total in [
            ("G0/G1 markers",   g0g1_present, G0G1_GENES_FBGN),
            ("G1/S markers",    g1s_present,  G1S_GENES_FBGN),
            ("S-phase markers", s_present,    S_GENES_FBGN),
            ("G2/M markers",    g2m_present,  G2M_GENES_FBGN),
        ]:
            print(f"  {label:<18}: {len(present)}/{len(total)}")
            if present:
                labels = [FBGN_TO_SYMBOL.get(g, g) for g in present[:10]]
                print("    " + ", ".join(labels) +
                      (f"  (+{len(present)-10} more)" if len(present) > 10 else ""))

    return g0g1_present, g1s_present, s_present, g2m_present


def _get_cc_expression(adata, cc_genes):
    """
    Extract log-normalised expression matrix for cc_genes only (cells × genes).
    Prefers adata.raw (full gene set) over adata.X (HVGs only), since cell
    cycle genes are frequently absent from HVG-filtered matrices.
    Returns (X_cc, genes_found, source_label).
    """
    def _extract_and_norm(X, var_names, genes):
        idx  = [list(var_names).index(g) for g in genes]
        Xsub = X[:, idx]
        if scipy.sparse.issparse(Xsub):
            Xsub = Xsub.toarray()
        Xsub = Xsub.astype(np.float32)
        if float(Xsub.max()) > 20:          # raw counts — normalise
            row_sums = Xsub.sum(axis=1, keepdims=True)
            row_sums[row_sums == 0] = 1
            Xsub = np.log1p(Xsub / row_sums * 1e4)
        return Xsub

    # Prefer .raw
    if adata.raw is not None:
        genes_found = [g for g in cc_genes if g in adata.raw.var_names]
        if genes_found:
            X_cc = _extract_and_norm(adata.raw.X, adata.raw.var_names, genes_found)
            print(f"  Using adata.raw  →  {len(genes_found)}/{len(cc_genes)} "
                  f"CC genes recovered")
            return X_cc, genes_found, 'raw'

    # Fall back to adata.X
    genes_found = [g for g in cc_genes if g in adata.var_names]
    if genes_found:
        n_hvg = len(genes_found)
        if n_hvg < len(cc_genes):
            print(f"  WARNING: adata.raw not found — using adata.X (HVGs only).")
            print(f"  Only {n_hvg}/{len(cc_genes)} CC genes present. Consider "
                  f"saving adata.raw before HVG filtering, or forcing CC genes "
                  f"into the HVG set.")
        X_cc = _extract_and_norm(adata.X, adata.var_names, genes_found)
        return X_cc, genes_found, 'X'

    return None, [], 'none'



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
    Validate cyclum phase assignments using Drosophila FlyBase cell cycle
    marker gene expression only.

    Strategy
    --------
    Rather than running DE across all genes, we extract expression of the
    FlyBase S-phase and G2/M marker genes directly (from adata.raw if
    available, to avoid HVG filtering removing them) and ask:
      - Do S-phase genes show higher expression in cyclum 's' cells?
      - Do G2/M genes show higher expression in cyclum 'g2/m' cells?

    Outputs
    -------
    {sample}_validation_umap_phase.pdf         UMAP coloured by cyclum phase
    {sample}_validation_umap_pseudotime.pdf    UMAP coloured by pseudotime
    {sample}_validation_umap_cc_genes.pdf      UMAP grid of all CC genes found
    {sample}_validation_umap_per_gene/         one PDF per CC gene
    {sample}_validation_heatmap.pdf            CC gene Z-score heatmap by phase
    {sample}_validation_cc_genes.csv           mean expression per phase per gene
    """
    print("\n" + "=" * 60)
    print("CELL CYCLE GENE VALIDATION OF CYCLUM PHASES")
    print("=" * 60)

    if 'cyclum_stage' not in adata.obs.columns:
        print("  ERROR: 'cyclum_stage' not in adata.obs. Run cyclum first.")
        return None

    os.makedirs(output_dir, exist_ok=True)

    # ── Find and extract CC gene expression ──────────────────────────────────
    g0g1_genes, g1s_genes, s_genes, g2m_genes = _find_marker_genes(adata, verbose=True)
    cc_genes  = g0g1_genes + g1s_genes + s_genes + g2m_genes

    if not cc_genes:
        print("  ERROR: No FlyBase cell cycle genes found. "
              "Check var_names use FBgn IDs or matching gene symbols.")
        return None

    X_cc, genes_found, source = _get_cc_expression(adata, cc_genes)

    if X_cc is None or not genes_found:
        print("  ERROR: Could not extract CC gene expression.")
        return None

    symbols_found = [FBGN_TO_SYMBOL.get(g, g) for g in genes_found]
    gene_type     = {g: ('G0G1' if g in g0g1_genes else 'G1S' if g in g1s_genes else 'S' if g in s_genes else 'G2M') for g in genes_found}

    # ── Mean expression per phase per gene ────────────────────────────────────
    phases = adata.obs['cyclum_stage'].values
    rows   = []
    for i, gene in enumerate(genes_found):
        for phase in PHASE_ORDER:
            mask = phases == phase
            rows.append({
                'gene':        gene,
                'gene_symbol': FBGN_TO_SYMBOL.get(gene, gene),
                'gene_type':   gene_type[gene],
                'phase':       phase,
                'mean_expr':   float(X_cc[mask, i].mean()),
                'pct_expr':    float((X_cc[mask, i] > 0).mean() * 100),
            })
    expr_df = pd.DataFrame(rows)

    print("\n  Mean expression per phase:")
    pivot = expr_df.pivot_table(
        index='gene_symbol', columns='phase',
        values='mean_expr', aggfunc='mean')
    pivot = pivot[[p for p in PHASE_ORDER if p in pivot.columns]]
    print(pivot.round(3).to_string())

    # Quick sanity check — which phase has highest mean for each gene?
    print("\n  Highest-expression phase per gene:")
    for _, row in pivot.iterrows():
        best = row.idxmax()
        gtype = expr_df.loc[expr_df['gene_symbol'] == row.name,
                             'gene_type'].values[0]
        expected = ('g0/g1' if gtype == 'G0G1' else
                    'g0/g1' if gtype == 'G1S'  else   # G1S peaks in g0/g1 by design
                    's'     if gtype == 'S'    else 'g2/m')
        flag = '✓' if best == expected else '✗'
        print(f"    {row.name:<20} {best}  {flag}  (expected {expected})")

    expr_df.to_csv(
        os.path.join(output_dir, f'{sample_name}_validation_cc_genes.csv'),
        index=False)

    # ── Heatmap: CC genes × cells sorted by cyclum phase ──────────────────────
    sort_key = pd.Series(phases, index=adata.obs_names).map(
        {'g0/g1': 0, 's': 1, 'g2/m': 2})
    if 'cyclum_pseudotime' in adata.obs.columns:
        sort_order = adata.obs.assign(_sk=sort_key).sort_values(
            ['_sk', 'cyclum_pseudotime']).index
    else:
        sort_order = adata.obs.assign(_sk=sort_key).sort_values('_sk').index

    cell_idx = [list(adata.obs_names).index(c) for c in sort_order]
    X_sorted = X_cc[cell_idx, :]
    means    = X_sorted.mean(axis=0)
    stds     = X_sorted.std(axis=0) + 1e-10
    X_z      = ((X_sorted - means) / stds).T   # genes × cells

    # Group S genes above G2M genes
    g0g1_idx_local = [i for i, g in enumerate(genes_found) if g in g0g1_genes]
    g1s_idx_local  = [i for i, g in enumerate(genes_found) if g in g1s_genes]
    s_idx_local    = [i for i, g in enumerate(genes_found) if g in s_genes]
    g2m_idx_local  = [i for i, g in enumerate(genes_found) if g in g2m_genes]
    row_order      = g0g1_idx_local + g1s_idx_local + s_idx_local + g2m_idx_local
    X_z_ordered   = X_z[row_order, :]
    labels_ordered = [symbols_found[i] for i in row_order]
    types_ordered  = [gene_type[genes_found[i]] for i in row_order]

    fig_h = max(5, len(genes_found) * 0.4 + 2)
    fig, ax = plt.subplots(figsize=(12, fig_h))
    sns.heatmap(X_z_ordered, cmap='RdBu_r', center=0, vmin=-2, vmax=2,
                xticklabels=False, yticklabels=labels_ordered,
                cbar_kws={'label': 'Z-score'}, ax=ax)

    for ytick, gtype in zip(ax.get_yticklabels(), types_ordered):
        ytick.set_color(GENE_TYPE_COLORS.get(gtype, 'black'))
        ytick.set_fontweight('bold')
        ytick.set_fontsize(9)

    # Dividing lines between gene type blocks
    n_g0g1 = len(g0g1_idx_local)
    n_g1s  = len(g1s_idx_local)
    n_s    = len(s_idx_local)
    for boundary in [n_g0g1, n_g0g1 + n_g1s, n_g0g1 + n_g1s + n_s]:
        if 0 < boundary < len(genes_found):
            ax.axhline(boundary, color='black', linewidth=1.5, linestyle='--')

    ax.set_xlabel('Cells (sorted: g0/g1 → s → g2/m, then by pseudotime)')
    ax.set_ylabel('FlyBase cell cycle genes')
    ax.set_title(
        f'Cell cycle gene expression by cyclum phase — {sample_name}\n'
        f'Red = G0/G1  |  Orange = G1/S  |  Teal = S  |  Blue = G2/M  |  '
        f'(from adata.{source})')

    # Phase colour strip
    phase_strip = np.array([[
        list(plt.cm.colors.to_rgb(PHASE_COLORS.get(p, '#CCCCCC')))
        for p in adata.obs.loc[sort_order, 'cyclum_stage']
    ]])
    ax2 = ax.inset_axes([0, -0.03, 1, 0.02])
    ax2.imshow(phase_strip, aspect='auto')
    ax2.axis('off')
    from matplotlib.patches import Patch as _Patch
    ax2.legend(
        handles=[_Patch(facecolor=c, label=p) for p, c in PHASE_COLORS.items()],
        loc='lower right', bbox_to_anchor=(1, -2), ncol=3,
        fontsize=8, frameon=True, title='Cyclum phase',
    )
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir,
                             f'{sample_name}_validation_heatmap.pdf'),
                dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  Heatmap saved ({len(genes_found)} genes).")

    # ── UMAPs ─────────────────────────────────────────────────────────────────
    if 'X_umap' not in adata.obsm:
        print("  Skipping UMAP plots (no X_umap in adata.obsm).")
        return expr_df

    # Phase + pseudotime reference UMAPs
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    sc.pl.umap(adata, color='cyclum_stage', ax=axes[0], show=False,
               title=f'Cyclum phase — {sample_name}',
               palette=PHASE_COLORS, frameon=False)
    if 'cyclum_pseudotime' in adata.obs.columns:
        sc.pl.umap(adata, color='cyclum_pseudotime', ax=axes[1], show=False,
                   title=f'Cyclum pseudotime — {sample_name}',
                   cmap='hsv', frameon=False)
    else:
        axes[1].set_visible(False)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir,
                             f'{sample_name}_validation_umap_phase.pdf'),
                dpi=300, bbox_inches='tight')
    plt.close()

    # CC genes need to be in adata.var_names for sc.pl.umap
    # If they came from .raw, temporarily add them to adata.obs for plotting
    genes_in_var  = [g for g in genes_found if g in adata.var_names]
    genes_via_obs = [g for g in genes_found if g not in adata.var_names]

    # Add raw-only genes as temporary obs columns
    for i, gene in enumerate(genes_found):
        if gene in genes_via_obs:
            col = f'_cc_{gene}'
            adata.obs[col] = X_cc[:, i]

    def _umap_color_key(gene):
        """Return the key to pass to sc.pl.umap color for this gene."""
        return gene if gene in genes_in_var else f'_cc_{gene}'

    # UMAP grid: all CC genes
    n_genes = len(genes_found)
    n_cols  = min(5, n_genes + 1)
    n_rows  = int(np.ceil((n_genes + 1) / n_cols))
    fig, axes = plt.subplots(n_rows, n_cols,
                              figsize=(5 * n_cols, 4.5 * n_rows))
    axes = np.array(axes).flatten()

    sc.pl.umap(adata, color='cyclum_stage', ax=axes[0], show=False,
               title='Cyclum phase', palette=PHASE_COLORS,
               frameon=False, legend_loc='on data', legend_fontsize=8)

    for i, (gene, symbol) in enumerate(zip(genes_found, symbols_found)):
        ax     = axes[i + 1]
        gtype  = gene_type[gene]
        color_key = _umap_color_key(gene)
        sc.pl.umap(adata, color=color_key, ax=ax, show=False,
                   title=f'{symbol} ({gtype})',
                   cmap='viridis', frameon=False)

    for ax in axes[n_genes + 1:]:
        ax.set_visible(False)

    fig.suptitle(
        f'FlyBase cell cycle gene expression — {sample_name}\n'
        f'Teal = S-phase  |  Blue = G2/M  |  source: adata.{source}',
        fontsize=11, y=1.01,
    )
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir,
                             f'{sample_name}_validation_umap_cc_genes.pdf'),
                dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  UMAP grid saved ({n_genes} CC genes).")

    # Individual per-gene UMAPs
    umap_dir = os.path.join(output_dir, f'{sample_name}_validation_umap_per_gene')
    os.makedirs(umap_dir, exist_ok=True)
    for gene, symbol in zip(genes_found, symbols_found):
        gtype     = gene_type[gene]
        color_key = _umap_color_key(gene)
        fig, ax   = plt.subplots(figsize=(7, 6))
        sc.pl.umap(adata, color=color_key, ax=ax, show=False,
                   title=f'{symbol}  [{gtype} marker]\n'
                         f'(adata.{source})',
                   cmap='viridis', frameon=False)
        plt.tight_layout()
        safe = symbol.replace('/', '-').replace(' ', '_')
        plt.savefig(os.path.join(umap_dir, f'{sample_name}_umap_{safe}.pdf'),
                    dpi=300, bbox_inches='tight')
        plt.close()
    print(f"  Individual gene UMAPs saved to: {umap_dir}/")

    # Clean up temporary obs columns
    for gene in genes_via_obs:
        col = f'_cc_{gene}'
        if col in adata.obs.columns:
            del adata.obs[col]

    print(f"\n  Validation outputs saved to: {output_dir}")
    return expr_df


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

    # ── Find CC genes and extract expression ──────────────────────────────────
    g0g1_genes, g1s_genes, s_genes, g2m_genes = _find_marker_genes(adata, verbose=True)
    cc_genes = s_genes + g2m_genes

    if not cc_genes:
        print("  ERROR: No FlyBase cell cycle genes found in dataset.")
        return None

    X_cc, genes_found, source = _get_cc_expression(adata, cc_genes)

    if X_cc is None or not genes_found:
        print("  ERROR: Could not extract CC gene expression.")
        return None

    g0g1_genes = [g for g in genes_found if g in g0g1_genes]
    g1s_genes  = [g for g in genes_found if g in g1s_genes]
    s_genes    = [g for g in genes_found if g in s_genes]
    g2m_genes  = [g for g in genes_found if g in g2m_genes]

    print(f"\n  Testing {len(genes_found)} cell cycle genes "
          f"(from adata.{source}) across "
          f"{adata.obs[leiden_col].nunique()} clusters...")

    # ── Kruskal-Wallis per gene across clusters ───────────────────────────────
    clusters  = sorted(adata.obs[leiden_col].unique())
    bonf_thr  = 0.05 / len(genes_found)
    leiden_vals = adata.obs[leiden_col].values
    rows = []
    for i, gene in enumerate(genes_found):
        expr   = X_cc[:, i]
        groups = [expr[leiden_vals == c] for c in clusters]
        groups = [g for g in groups if len(g) >= 3]
        if len(groups) < 2:
            continue
        h, p = kruskal(*groups)
        rows.append({
            'gene':             gene,
            'gene_symbol':      FBGN_TO_SYMBOL.get(gene, gene),
            'gene_type':        ('G0G1' if gene in g0g1_genes else 'G1S' if gene in g1s_genes else 'S' if gene in s_genes else 'G2M'),
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

    # ── Add raw-only genes as temporary obs columns for plotting ─────────────
    genes_in_var  = [g for g in genes_found if g in adata.var_names]
    genes_via_obs = [g for g in genes_found if g not in adata.var_names]
    gene_to_idx   = {g: i for i, g in enumerate(genes_found)}
    for gene in genes_via_obs:
        adata.obs[f'_cc_{gene}'] = X_cc[:, gene_to_idx[gene]]

    def _plot_key(gene):
        return gene if gene in genes_in_var else f'_cc_{gene}'

    # ── Dot plot ──────────────────────────────────────────────────────────────
    fig_title = (f"Cell cycle gene expression by cluster — {sample_name}\n"
                 f"({n_sig} Bonferroni-significant of {len(stats_df)} tested, "
                 f"ranked by Kruskal-Wallis H)")
    # dotplot only works with var_names keys; skip genes only in obs
    dp_g0g1_genes = [g for g in plot_genes if g in g0g1_genes and g in genes_in_var]
    dp_g1s_genes  = [g for g in plot_genes if g in g1s_genes  and g in genes_in_var]
    dp_s_genes    = [g for g in plot_genes if g in s_genes    and g in genes_in_var]
    dp_g2m_genes  = [g for g in plot_genes if g in g2m_genes  and g in genes_in_var]
    if dp_g0g1_genes or dp_g1s_genes or dp_s_genes or dp_g2m_genes:
        try:
            dp = sc.pl.dotplot(
                adata,
                var_names={k: v for k, v in
                           [('G0/G1', dp_g0g1_genes), ('G1/S', dp_g1s_genes), ('S-phase', dp_s_genes), ('G2/M', dp_g2m_genes)] if v},
                groupby=leiden_col,
                use_raw=False,
                show=False,
                return_fig=True,
                title=fig_title,
                var_group_rotation=0,
            )
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
    else:
        print("  Skipping dot plot (all CC genes sourced from .raw, "
              "not in adata.var_names).")

    # ── Clustermap: mean expression per cluster × CC gene ────────────────────
    # Design choices:
    #   - ALL CC genes shown (not just Bonferroni-sig), ordered G0G1→G1S→S→G2M
    #   - Raw mean log-normalised expression (not Z-scored) so colour reflects
    #     actual expression levels, not just rank within gene
    #   - Columns fixed in biological order; rows (clusters) dendrogrammed by
    #     Ward linkage so similar-CC-profile clusters group together
    #   - Dividing lines separate gene type blocks
    #   - Asterisk on x-tick = Bonferroni-significant KW test

    # Build full cluster × gene matrix (all CC genes, biological column order)
    gene_to_idx = {g: i for i, g in enumerate(genes_found)}
    col_order   = ([g for g in G0G1_GENES_FBGN + G1S_GENES_FBGN +
                    S_GENES_FBGN + G2M_GENES_FBGN
                    if g in genes_found])
    col_symbols = [FBGN_TO_SYMBOL.get(g, g) for g in col_order]

    cluster_mean = pd.DataFrame(index=clusters, columns=col_order, dtype=float)
    for gene in col_order:
        expr = X_cc[:, gene_to_idx[gene]]
        for c in clusters:
            cluster_mean.loc[c, gene] = float(expr[leiden_vals == c].mean())

    # Rename columns to gene symbols for display
    cluster_mean.columns = col_symbols

    # Row dendrogram only — cluster Leiden clusters by CC expression profile
    from matplotlib.patches import Patch
    fig_w = max(10, len(col_order) * 0.6 + 3)
    fig_h = max(5,  len(clusters)  * 0.5 + 3)

    g = sns.clustermap(
        cluster_mean,
        method='ward',
        metric='euclidean',
        row_cluster=True,       # cluster rows (Leiden clusters)
        col_cluster=False,      # keep biological gene order fixed
        cmap='viridis',         # sequential: shows absolute expression levels
        linewidths=0.3,
        figsize=(fig_w, fig_h),
        cbar_pos=(1.02, 0.35, 0.025, 0.3),
        cbar_kws={'label': 'Mean log-normalised expression'},
        dendrogram_ratio=(0.15, 0.0),
        yticklabels=True,
        xticklabels=True,
    )

    ax = g.ax_heatmap

    # X-tick labels: colour by gene type, bold + asterisk if Bonferroni-sig
    sig_genes = set(stats_df.loc[stats_df['bonf_significant'], 'gene'].values)
    for xtick, (sym, gene) in zip(ax.get_xticklabels(),
                                   zip(col_symbols, col_order)):
        gtype = ('G0G1' if gene in g0g1_genes else
                 'G1S'  if gene in g1s_genes  else
                 'S'    if gene in s_genes    else 'G2M')
        xtick.set_color(GENE_TYPE_COLORS[gtype])
        xtick.set_fontsize(8)
        xtick.set_rotation(45)
        xtick.set_ha('right')
        if gene in sig_genes:
            xtick.set_fontweight('bold')
            xtick.set_text(sym + ' *')

    # Vertical dividing lines between gene type blocks
    boundaries = []
    prev_type  = None
    for i, gene in enumerate(col_order):
        gtype = ('G0G1' if gene in g0g1_genes else
                 'G1S'  if gene in g1s_genes  else
                 'S'    if gene in s_genes    else 'G2M')
        if prev_type is not None and gtype != prev_type:
            boundaries.append(i)
        prev_type = gtype
    for b in boundaries:
        ax.axvline(b, color='white', linewidth=2.5, linestyle='-')

    # Gene type colour bar above x-axis
    type_strip_data = np.array([[
        matplotlib.colors.to_rgb(
            GENE_TYPE_COLORS['G0G1'] if g in g0g1_genes else
            GENE_TYPE_COLORS['G1S']  if g in g1s_genes  else
            GENE_TYPE_COLORS['S']    if g in s_genes    else
            GENE_TYPE_COLORS['G2M']
        ) for g in col_order
    ]])
    ax_strip = ax.inset_axes([0, 1.0, 1, 0.03])
    ax_strip.imshow(type_strip_data, aspect='auto', interpolation='none')
    ax_strip.set_axis_off()

    # Gene type block labels above strip
    block_starts = [0] + boundaries
    block_ends   = boundaries + [len(col_order)]
    block_labels = []
    for gene in [col_order[s] for s in block_starts]:
        block_labels.append(
            'G0/G1' if gene in g0g1_genes else
            'G1/S'  if gene in g1s_genes  else
            'S'     if gene in s_genes    else 'G2M'
        )
    ax_label = ax.inset_axes([0, 1.04, 1, 0.04])
    ax_label.set_xlim(0, len(col_order))
    ax_label.set_axis_off()
    for label, start, end in zip(block_labels, block_starts, block_ends):
        mid = (start + end) / 2
        gtype_key = ('G0G1' if label == 'G0/G1' else
                     'G1S'  if label == 'G1/S'  else
                     'S'    if label == 'S'      else 'G2M')
        ax_label.text(mid, 0.1, label, ha='center', va='bottom',
                      fontsize=9, fontweight='bold',
                      color=GENE_TYPE_COLORS[gtype_key])

    ax.set_xlabel('')
    ax.set_ylabel(f'Leiden cluster ({leiden_col})', fontsize=10)

    g.fig.suptitle(
        f'CC gene expression by cluster — {sample_name}\n'
        f'Mean log-normalised expression  |  '
        f'Clusters dendrogrammed by Ward linkage  |  * = Bonferroni-sig',
        fontsize=9, y=1.06,
    )

    # Legend
    legend_patches = [Patch(color=c, label=l) for l, c in [
        ('G0/G1',   GENE_TYPE_COLORS['G0G1']),
        ('G1/S',    GENE_TYPE_COLORS['G1S']),
        ('S-phase', GENE_TYPE_COLORS['S']),
        ('G2/M',    GENE_TYPE_COLORS['G2M']),
    ]]
    g.ax_heatmap.legend(
        handles=legend_patches, title='Gene type',
        loc='upper left', bbox_to_anchor=(1.08, 1.0),
        fontsize=8, frameon=True,
    )

    g.fig.savefig(os.path.join(output_dir,
                               f'{sample_name}_cc_cluster_heatmap.pdf'),
                  dpi=300, bbox_inches='tight')
    plt.close(g.fig)
    print("  Cluster heatmap saved (fixed column order).")

    # ── Second heatmap: fully clustered rows AND cols, Z-scored ──────────────
    # Z-scoring per gene removes absolute expression differences so the
    # dendrogram reflects which genes co-vary across clusters — useful for
    # seeing whether G2M genes all behave as a bloc vs. sub-groups.
    cluster_mean_z = cluster_mean.apply(
        lambda col: (col - col.mean()) / (col.std() + 1e-10), axis=0)

    symbol_to_gene = {FBGN_TO_SYMBOL.get(g, g): g for g in col_order}
    col_colors = pd.Series(
        {sym: (GENE_TYPE_COLORS['G0G1'] if symbol_to_gene[sym] in g0g1_genes
               else GENE_TYPE_COLORS['G1S'] if symbol_to_gene[sym] in g1s_genes
               else GENE_TYPE_COLORS['S']   if symbol_to_gene[sym] in s_genes
               else GENE_TYPE_COLORS['G2M'])
         for sym in cluster_mean_z.columns},
        name='Gene type',
    )

    g2 = sns.clustermap(
        cluster_mean_z,
        method='ward',
        metric='euclidean',
        row_cluster=True,
        col_cluster=True,
        col_colors=col_colors,
        cmap='RdBu_r',
        center=0, vmin=-2, vmax=2,
        linewidths=0.3,
        figsize=(fig_w, fig_h),
        cbar_pos=(1.02, 0.35, 0.025, 0.3),
        cbar_kws={'label': 'Z-score (across clusters)'},
        dendrogram_ratio=(0.15, 0.12),
        colors_ratio=0.03,
        yticklabels=True,
        xticklabels=True,
    )

    ax2 = g2.ax_heatmap
    reordered_genes = [col_order[i] for i in g2.dendrogram_col.reordered_ind]
    for xtick, gene in zip(ax2.get_xticklabels(), reordered_genes):
        sym   = FBGN_TO_SYMBOL.get(gene, gene)
        gtype = ('G0G1' if gene in g0g1_genes else
                 'G1S'  if gene in g1s_genes  else
                 'S'    if gene in s_genes    else 'G2M')
        xtick.set_color(GENE_TYPE_COLORS[gtype])
        xtick.set_fontsize(8)
        xtick.set_rotation(45)
        xtick.set_ha('right')
        if gene in sig_genes:
            xtick.set_fontweight('bold')
            xtick.set_text(sym + ' *')

    ax2.set_xlabel('')
    ax2.set_ylabel(f'Leiden cluster ({leiden_col})', fontsize=10)
    g2.fig.suptitle(
        f'CC gene expression by cluster (Z-scored, fully clustered) — {sample_name}\n'
        f'Both rows and cols dendrogrammed  |  * = Bonferroni-sig  |  Ward linkage',
        fontsize=9, y=1.04,
    )
    legend_patches2 = [Patch(color=c, label=l) for l, c in [
        ('G0/G1',   GENE_TYPE_COLORS['G0G1']),
        ('G1/S',    GENE_TYPE_COLORS['G1S']),
        ('S-phase', GENE_TYPE_COLORS['S']),
        ('G2/M',    GENE_TYPE_COLORS['G2M']),
    ]]
    g2.ax_heatmap.legend(
        handles=legend_patches2, title='Gene type',
        loc='upper left', bbox_to_anchor=(1.08, 1.0),
        fontsize=8, frameon=True,
    )
    g2.fig.savefig(os.path.join(output_dir,
                                f'{sample_name}_cc_cluster_heatmap_clustered.pdf'),
                   dpi=300, bbox_inches='tight')
    plt.close(g2.fig)
    print("  Cluster heatmap saved (fully clustered, Z-scored).")


    # ── UMAP grid: top genes coloured by expression ───────────────────────────
    if 'X_umap' in adata.obsm:
        # All top genes in one grid
        n_genes = len(plot_genes)
        n_cols  = min(5, n_genes + 1)
        n_rows  = int(np.ceil((n_genes + 1) / n_cols))
        fig, axes = plt.subplots(n_rows, n_cols,
                                  figsize=(5 * n_cols, 4.5 * n_rows))
        axes = np.array(axes).flatten()

        sc.pl.umap(adata, color=leiden_col, ax=axes[0], show=False,
                   title=f'Clusters ({leiden_col})', frameon=False,
                   legend_loc='on data', legend_fontsize=8)

        for i, (gene, symbol) in enumerate(zip(plot_genes, plot_symbols)):
            ax    = axes[i + 1]
            h     = stats_df.loc[stats_df['gene'] == gene, 'KW_H'].values[0]
            sig   = stats_df.loc[stats_df['gene'] == gene, 'bonf_significant'].values[0]
            gtype = ('G0G1' if gene in g0g1_genes else 'G1S' if gene in g1s_genes
                     else 'S' if gene in s_genes else 'G2M')
            title = f'{symbol} ({gtype})\nH={h:.1f}{"*" if sig else ""}'
            sc.pl.umap(adata, color=_plot_key(gene), ax=ax, show=False,
                       title=title, cmap='viridis', frameon=False)

        for ax in axes[n_genes + 1:]:
            ax.set_visible(False)

        fig.suptitle(
            f'Cell cycle gene expression by cluster — {sample_name}\n'
            f'* = Bonferroni-significant (p < {bonf_thr:.1e})  |  '
            f'Ranked by Kruskal-Wallis H  |  source: adata.{source}',
            fontsize=11, y=1.01,
        )
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir,
                                 f'{sample_name}_cc_cluster_umap_grid.pdf'),
                    dpi=300, bbox_inches='tight')
        plt.close()
        print(f"  UMAP grid saved ({n_genes} genes).")

        # Reference: clusters + cyclum phase side by side
        fig, axes = plt.subplots(1, 2, figsize=(14, 6))
        sc.pl.umap(adata, color=leiden_col, ax=axes[0], show=False,
                   title=f'Leiden clusters ({leiden_col})', frameon=False,
                   legend_loc='on data')
        sc.pl.umap(adata, color='cyclum_stage', ax=axes[1], show=False,
                   title='Cyclum phase', frameon=False, palette=PHASE_COLORS)
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir,
                                 f'{sample_name}_cc_cluster_umap_reference.pdf'),
                    dpi=300, bbox_inches='tight')
        plt.close()

        # Individual per-gene UMAPs
        umap_dir = os.path.join(output_dir, f'{sample_name}_umap_per_gene')
        os.makedirs(umap_dir, exist_ok=True)
        for gene, symbol in zip(plot_genes, plot_symbols):
            h     = stats_df.loc[stats_df['gene'] == gene, 'KW_H'].values[0]
            p     = stats_df.loc[stats_df['gene'] == gene, 'KW_p'].values[0]
            sig   = stats_df.loc[stats_df['gene'] == gene, 'bonf_significant'].values[0]
            gtype = ('G0/G1' if gene in g0g1_genes else 'G1/S' if gene in g1s_genes
                     else 'S-phase' if gene in s_genes else 'G2/M')
            title = (f'{symbol}  [{gtype}]\n'
                     f'KW H={h:.2f}, p={p:.2e}'
                     f'{"  *Bonferroni-sig" if sig else ""}')
            fig, ax = plt.subplots(figsize=(7, 6))
            sc.pl.umap(adata, color=_plot_key(gene), ax=ax, show=False,
                       title=title, cmap='viridis', frameon=False)
            plt.tight_layout()
            safe = symbol.replace('/', '-').replace(' ', '_')
            plt.savefig(os.path.join(umap_dir,
                                     f'{sample_name}_umap_{safe}.pdf'),
                        dpi=300, bbox_inches='tight')
            plt.close()
        print(f"  Individual gene UMAPs saved to: {umap_dir}/")

    # Clean up temporary obs columns
    for gene in genes_via_obs:
        col = f'_cc_{gene}'
        if col in adata.obs.columns:
            del adata.obs[col]

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
