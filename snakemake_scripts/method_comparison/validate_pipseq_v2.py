import os
import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from glob import glob
import argparse
import scipy.sparse
import re
import anndata as ad
import harmonypy
import warnings
import seaborn as sns
from scipy.stats import chi2_contingency, mannwhitneyu

def integrate(files, out_path, fig_dir, sample, batch_key, min_cells, min_genes, 
              calculate_titer=True, n_pcs=30, bio_condition=None, leiden_resolution=0.5):
    os.makedirs(fig_dir, exist_ok=True)

    adatas = []
    for file_path in files:
        adata = sc.read_h5ad(file_path)
        batch_name = os.path.splitext(os.path.basename(file_path))[0]
        adata.obs[batch_key] = batch_name
        adatas.append(adata)

    sc.settings.figdir = fig_dir
    combined = ad.concat(adatas, join='inner', merge='same', index_unique='-')

    # Metadata extraction
    combined.obs['cell_line']  = combined.obs[batch_key].str.extract(r'(JW18DOX|JW18wMel)')[0]
    combined.obs['treatment']  = combined.obs[batch_key].str.extract(r'-(Ctrl|SV)-')[0]
    combined.obs['timepoint']  = combined.obs[batch_key].str.extract(r'-(D\d+)-')[0]
    combined.obs['replicate']  = combined.obs[batch_key].str.extract(r'-(\d+)_')[0]  # fixed regex
    combined.obs['method']     = combined.obs[batch_key].str.extract(r'_(10x|pipseq)$')[0]
    combined.obs['bio_condition'] = combined.obs.apply(
        lambda row: f"{row['cell_line']}-{row['treatment']}-{row['timepoint']}"
                    if pd.notna(row['timepoint'])
                    else f"{row['cell_line']}-{row['treatment']}",
        axis=1
    )

    if bio_condition:
        print(f"Filtering to biological condition: {bio_condition}")
        combined = combined[combined.obs['bio_condition'] == bio_condition].copy()
        sample = f"{sample}_{bio_condition}"

    print(f"Combined data shape for {sample}: {combined.shape}")
    print(f"Methods: {combined.obs['method'].value_counts()}")

    # Remove bacterial genes
    bacteria_genes = ['GQX67_00940', 'GQX67_05945'] + [
        g for g in combined.var_names if g.startswith('16S_')]
    combined = combined[:, ~combined.var_names.isin(bacteria_genes)].copy()

    # --- Preprocessing (order matters) ---
    # Remove stored zeros, filter, then convert to dense
    if scipy.sparse.issparse(combined.X):
        combined.X.eliminate_zeros()
    sc.pp.filter_cells(combined, min_genes=min_genes)
    sc.pp.filter_cells(combined, min_counts=1)
    sc.pp.filter_genes(combined, min_cells=min_cells)
    sc.pp.filter_genes(combined, min_counts=1)

    if scipy.sparse.issparse(combined.X):
        combined.X = combined.X.toarray()
    combined.X = np.nan_to_num(combined.X.astype(np.float64), nan=0.0, posinf=0.0, neginf=0.0)

    # Drop zero-count cells/genes after dense conversion
    cell_sums = combined.X.sum(axis=1)
    gene_sums = combined.X.sum(axis=0)
    combined = combined[cell_sums > 0].copy()
    combined = combined[:, gene_sums > 0].copy()

    # Remove near-constant genes (cause NaN dispersion in HVG)
    gene_var  = combined.X.var(axis=0)
    gene_mean = combined.X.mean(axis=0)
    keep = gene_var > np.maximum(gene_mean * 1e-6, 1e-10)
    combined = combined[:, keep].copy()
    print(f"After filtering: {combined.n_obs} cells, {combined.n_vars} genes")

    # HVG on raw counts
    sc.pp.highly_variable_genes(combined, flavor='seurat', n_top_genes=2000,
                                 batch_key='method')

    # Normalize + log1p
    sc.pp.normalize_total(combined, target_sum=1e4)
    sc.pp.log1p(combined)
    combined.X = np.nan_to_num(combined.X, nan=0.0, posinf=0.0, neginf=0.0)

    # Store raw normalized counts for DE/visualization
    combined.raw = combined

    # Subset to HVGs, scale, PCA
    combined = combined[:, combined.var['highly_variable']].copy()
    sc.pp.scale(combined, max_value=10)
    sc.pp.pca(combined, n_comps=n_pcs)

    # Pre-correction UMAP
    combined_pre = combined.copy()
    sc.pp.neighbors(combined_pre, n_pcs=n_pcs)
    sc.tl.umap(combined_pre)
    sc.pl.umap(combined_pre, color=['method', 'bio_condition'],
               save=f'_{sample}_before_integration.pdf', ncols=2)
    del combined_pre

    # Harmony: correct for method + replicate simultaneously
    # Fills NaN in covariates with 'unknown' so Harmony doesn't silently skip cells
    for col in ['method', 'replicate']:
        combined.obs[col] = combined.obs[col].fillna('unknown').astype(str)

    print("Running Harmony batch correction (method + replicate) ...")
    ho = harmonypy.run_harmony(
        combined.obsm['X_pca'][:, :n_pcs],
        combined.obs,
        vars_use=['method', 'replicate'],
        max_iter_harmony=30,
        random_state=42,
    )
    combined.obsm['X_pca_harmony'] = ho.Z_corr.T

    # Build neighbor graph and UMAP on corrected embedding
    sc.pp.neighbors(combined, use_rep='X_pca_harmony', n_pcs=n_pcs)
    sc.tl.umap(combined)
    sc.tl.leiden(combined, resolution=leiden_resolution)

    print(f"Saving integrated object to {out_path}")
    combined.write(out_path)

    plot_ctrl_comparison(combined, fig_dir, sample)
    plot_full_dataset_analysis(combined, fig_dir, sample)

    print("\n" + "="*60)
    print("SUMMARY")
    print("="*60)
    print(f"Total cells: {combined.n_obs}")
    print(f"Total genes: {combined.n_vars}")
    print(f"Cells per method:\n{combined.obs['method'].value_counts()}")
    print(f"Clusters: {combined.obs['leiden'].nunique()}")
    print(f"Cells per cluster:\n{combined.obs['leiden'].value_counts().sort_index()}")


# --- plotting functions unchanged below ---

def plot_ctrl_comparison(combined, fig_dir, sample):
    print("\n" + "="*60)
    print("CTRL SAMPLES - METHOD COMPARISON")
    print("="*60)

    combined_ctrl = combined[combined.obs['treatment'] == 'Ctrl'].copy()
    if combined_ctrl.n_obs == 0:
        print("No Ctrl samples found, skipping comparison")
        return

    print(f"Ctrl samples: {combined_ctrl.n_obs} cells")
    print(f"Methods in Ctrl: {combined_ctrl.obs['method'].value_counts()}")

    clusters = sorted(combined_ctrl.obs['leiden'].unique())
    cmap = plt.cm.get_cmap('tab20')
    leiden_colors = [cmap(i % 20) for i in range(len(clusters))]

    sc.pl.umap(combined_ctrl, color='bio_condition',
               save=f'_{sample}_ctrl_by_bio_condition.pdf',
               title='Ctrl samples by biological condition')
    if 'wolbachia_titer' in combined_ctrl.obs.columns:
        sc.pl.umap(combined_ctrl, color='wolbachia_titer',
                   save=f'_{sample}_ctrl_by_titer.pdf',
                   title='Ctrl samples - Wolbachia titer',
                   vmax=np.percentile(combined_ctrl.obs['wolbachia_titer'].dropna(), 95))
    sc.pl.umap(combined_ctrl, color='leiden',
               save=f'_{sample}_ctrl_by_leiden.pdf',
               title='Ctrl samples - Leiden clusters', legend_loc='on data')
    sc.pl.umap(combined_ctrl, color='method',
               save=f'_{sample}_ctrl_by_method.pdf',
               title='Ctrl samples - Library prep method')

    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    sc.pl.umap(combined_ctrl[combined_ctrl.obs['method'] == '10x'],
               color='leiden', ax=axes[0], show=False,
               title='10x Genomics - Ctrl', frameon=False)
    sc.pl.umap(combined_ctrl[combined_ctrl.obs['method'] == 'pipseq'],
               color='leiden', ax=axes[1], show=False,
               title='PIPseq - Ctrl', frameon=False)
    plt.tight_layout()
    plt.savefig(os.path.join(fig_dir, f'umap_{sample}_ctrl_split_by_method.pdf'))
    plt.close()

    if 'phase' in combined_ctrl.obs.columns:
        sc.pl.umap(combined_ctrl, color='phase',
                   save=f'_{sample}_ctrl_by_cellcycle.pdf',
                   title='Ctrl samples - Cell cycle phase')
    elif 'cyclum_theta' in combined_ctrl.obs.columns:
        sc.pl.umap(combined_ctrl, color='cyclum_theta',
                   save=f'_{sample}_ctrl_by_cyclum_theta.pdf',
                   title='Ctrl samples - Cyclum theta', cmap='twilight')

    for cell_line, title in [('JW18wMel', 'wMel Ctrl'), ('JW18DOX', 'DOX Ctrl')]:
        if 'wolbachia_titer' not in combined_ctrl.obs.columns:
            break
        subset = combined_ctrl[combined_ctrl.obs['cell_line'] == cell_line]
        if subset.n_obs == 0:
            continue
        fig, ax = plt.subplots(figsize=(8, 6))
        sns.boxplot(data=subset.obs, x='method', y='wolbachia_titer',
                    ax=ax, palette=['#1f77b4', '#ff7f0e'])
        sns.stripplot(data=subset.obs, x='method', y='wolbachia_titer',
                      ax=ax, color='black', alpha=0.3, size=2)
        ax.set_title(f'{title} - Wolbachia titer comparison')
        t10 = subset.obs[subset.obs['method'] == '10x']['wolbachia_titer'].dropna()
        tpip = subset.obs[subset.obs['method'] == 'pipseq']['wolbachia_titer'].dropna()
        if len(t10) > 0 and len(tpip) > 0:
            _, p = mannwhitneyu(t10, tpip, alternative='two-sided')
            ax.text(0.5, 0.95, f'Mann-Whitney p = {p:.2e}',
                    transform=ax.transAxes, ha='center', va='top')
        plt.tight_layout()
        plt.savefig(os.path.join(fig_dir, f'titer_comparison_{sample}_{cell_line}_ctrl.pdf'))
        plt.close()

    if 'wolbachia_titer' in combined_ctrl.obs.columns:
        fig, ax = plt.subplots(figsize=(12, 6))
        plot_data = combined_ctrl.obs[['leiden', 'wolbachia_titer']].sort_values('leiden')
        clusters_sorted = sorted(plot_data['leiden'].unique())
        sns.stripplot(data=plot_data, x='leiden', y='wolbachia_titer',
                      color='black', alpha=0.3, size=2, ax=ax)
        bp = ax.boxplot([plot_data[plot_data['leiden'] == c]['wolbachia_titer'].values
                         for c in clusters_sorted],
                        positions=range(len(clusters_sorted)), widths=0.6,
                        patch_artist=True, medianprops=dict(color='black', linewidth=2))
        for patch, color in zip(bp['boxes'], leiden_colors):
            patch.set_facecolor(color); patch.set_alpha(0.7)
        ax.set_xlabel('Leiden Cluster'); ax.set_ylabel('Wolbachia Titer')
        ax.set_title('Ctrl samples - Wolbachia titer by cluster')
        plt.tight_layout()
        plt.savefig(os.path.join(fig_dir, f'boxplot_{sample}_ctrl_titer_by_cluster.pdf'))
        plt.close()

    genes_per_cluster = combined_ctrl.obs.groupby(['leiden', 'method'])['n_genes'].mean().unstack()
    fig, ax = plt.subplots(figsize=(10, 6))
    genes_per_cluster.plot(kind='bar', ax=ax, color=['#1f77b4', '#ff7f0e'])
    ax.set_xlabel('Leiden Cluster'); ax.set_ylabel('Mean genes per cell')
    ax.set_title('Ctrl samples - Genes per cluster by method')
    ax.legend(title='Method'); plt.xticks(rotation=0); plt.tight_layout()
    plt.savefig(os.path.join(fig_dir, f'genes_per_cluster_{sample}_ctrl.pdf'))
    plt.close()

    cluster_comp = pd.crosstab(combined_ctrl.obs['leiden'], combined_ctrl.obs['method'],
                                normalize='columns') * 100
    fig, ax = plt.subplots(figsize=(10, 6))
    cluster_comp.plot(kind='bar', ax=ax, color=['#1f77b4', '#ff7f0e'])
    ax.set_xlabel('Leiden Cluster'); ax.set_ylabel('Percentage of cells')
    ax.set_title('Ctrl samples - Cell distribution by method')
    ax.legend(title='Method'); plt.xticks(rotation=0); plt.tight_layout()
    plt.savefig(os.path.join(fig_dir, f'cell_percentage_per_cluster_{sample}_ctrl.pdf'))
    plt.close()

    contingency = pd.crosstab(combined_ctrl.obs['leiden'], combined_ctrl.obs['method'])
    chi2, p_value, dof, _ = chi2_contingency(contingency)
    print("\nCTRL SAMPLES - CLUSTER COMPOSITION BY METHOD")
    print(cluster_comp)
    print(f"\nChi-square test: χ² = {chi2:.2f}, p = {p_value:.2e}")
    print(f"Methods show {'SIGNIFICANT' if p_value < 0.05 else 'NO'} difference in cluster distribution")

    if 'wolbachia_titer' in combined_ctrl.obs.columns:
        fig, axes = plt.subplots(1, 2, figsize=(14, 5))
        vmax = np.percentile(combined_ctrl.obs['wolbachia_titer'].dropna(), 95)
        for ax, meth, title in zip(axes, ['10x', 'pipseq'],
                                   ['10x Genomics - Wolbachia titer', 'PIPseq - Wolbachia titer']):
            sc.pl.umap(combined_ctrl[combined_ctrl.obs['method'] == meth],
                       color='wolbachia_titer', ax=ax, show=False,
                       title=title, frameon=False, vmax=vmax)
        plt.tight_layout()
        plt.savefig(os.path.join(fig_dir, f'umap_{sample}_ctrl_titer_split.pdf'))
        plt.close()


def plot_full_dataset_analysis(combined, fig_dir, sample):
    print("\n" + "="*60)
    print("FULL DATASET ANALYSIS")
    print("="*60)

    clusters = sorted(combined.obs['leiden'].unique())
    cmap = plt.cm.get_cmap('tab20')
    leiden_colors = [cmap(i % 20) for i in range(len(clusters))]

    sc.pl.umap(combined, color='leiden',
               save=f'_{sample}_full_leiden_clusters.pdf',
               title='All samples - Leiden clusters', legend_loc='on data')

    if 'phase' in combined.obs.columns:
        sc.pl.umap(combined, color='phase',
                   save=f'_{sample}_full_by_cellcycle_phase.pdf',
                   title='All samples - Cell cycle phase')
        fig, ax = plt.subplots(figsize=(12, 6))
        pd.crosstab(combined.obs['leiden'], combined.obs['phase'],
                    normalize='index').mul(100).plot(kind='bar', stacked=True, ax=ax)
        ax.set_title('Cell cycle phase distribution by cluster')
        plt.tight_layout()
        plt.savefig(os.path.join(fig_dir, f'cellcycle_by_cluster_{sample}_full.pdf'))
        plt.close()

    if 'cyclum_theta' in combined.obs.columns:
        sc.pl.umap(combined, color='cyclum_theta',
                   save=f'_{sample}_full_by_cyclum_theta.pdf',
                   title='All samples - Cyclum theta', cmap='twilight')

    genes_per_cluster = combined.obs.groupby('leiden')['n_genes'].agg(['mean', 'std'])
    fig, ax = plt.subplots(figsize=(12, 6))
    ax.bar(range(len(genes_per_cluster)), genes_per_cluster['mean'],
           yerr=genes_per_cluster['std'], color=leiden_colors, alpha=0.7, capsize=5)
    ax.set_xticks(range(len(genes_per_cluster)))
    ax.set_xticklabels(genes_per_cluster.index)
    ax.set_xlabel('Leiden Cluster'); ax.set_ylabel('Mean genes per cell')
    ax.set_title('Gene count by cluster (all samples)')
    plt.tight_layout()
    plt.savefig(os.path.join(fig_dir, f'genes_per_cluster_{sample}_full.pdf'))
    plt.close()

    if 'wolbachia_titer' in combined.obs.columns:
        sc.pl.umap(combined, color='wolbachia_titer',
                   save=f'_{sample}_full_wolbachia_titer.pdf',
                   title='All samples - Wolbachia titer',
                   vmax=np.percentile(combined.obs['wolbachia_titer'].dropna(), 95))

        fig, ax = plt.subplots(figsize=(14, 6))
        plot_data = combined.obs[['leiden', 'wolbachia_titer']].sort_values('leiden')
        clusters_sorted = sorted(plot_data['leiden'].unique())
        sns.stripplot(data=plot_data, x='leiden', y='wolbachia_titer',
                      color='black', alpha=0.2, size=1, ax=ax)
        bp = ax.boxplot([plot_data[plot_data['leiden'] == c]['wolbachia_titer'].values
                         for c in clusters_sorted],
                        positions=range(len(clusters_sorted)), widths=0.6,
                        patch_artist=True, medianprops=dict(color='black', linewidth=2))
        for patch, color in zip(bp['boxes'], leiden_colors):
            patch.set_facecolor(color); patch.set_alpha(0.7)
        ax.set_xlabel('Leiden Cluster'); ax.set_ylabel('Wolbachia Titer')
        ax.set_title('Wolbachia titer by cluster (all samples)')
        plt.tight_layout()
        plt.savefig(os.path.join(fig_dir, f'boxplot_{sample}_full_titer_by_cluster.pdf'))
        plt.close()

        infected_per_cluster = combined.obs.groupby('leiden').apply(
            lambda x: (x['wolbachia_titer'] > 0).sum() / len(x) * 100)
        fig, ax = plt.subplots(figsize=(12, 6))
        ax.bar(range(len(infected_per_cluster)), infected_per_cluster.values,
               color=leiden_colors, alpha=0.7)
        ax.set_xticks(range(len(infected_per_cluster)))
        ax.set_xticklabels(infected_per_cluster.index)
        ax.set_ylabel('% Infected cells')
        ax.set_title('Percentage of Wolbachia-infected cells by cluster')
        plt.tight_layout()
        plt.savefig(os.path.join(fig_dir, f'infection_percentage_by_cluster_{sample}_full.pdf'))
        plt.close()

    sc.pl.umap(combined, color=['cell_line', 'treatment'],
               save=f'_{sample}_full_by_condition.pdf', ncols=2)

    fig, ax = plt.subplots(figsize=(14, 8))
    pd.crosstab(combined.obs['leiden'], combined.obs['bio_condition'],
                normalize='columns').mul(100).T.plot(
        kind='bar', stacked=True, ax=ax, color=leiden_colors, width=0.8)
    ax.set_xlabel('Biological Condition'); ax.set_ylabel('Percentage of cells')
    ax.set_title('Cluster composition by biological condition')
    ax.legend(title='Leiden Cluster', bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.xticks(rotation=45, ha='right'); plt.tight_layout()
    plt.savefig(os.path.join(fig_dir, f'cluster_composition_by_condition_{sample}_full.pdf'))
    plt.close()

    if 'wolbachia_titer' in combined.obs.columns and 'phase' in combined.obs.columns:
        fig, ax = plt.subplots(figsize=(10, 6))
        sns.boxplot(data=combined.obs, x='phase', y='wolbachia_titer', ax=ax)
        sns.stripplot(data=combined.obs, x='phase', y='wolbachia_titer',
                      ax=ax, color='black', alpha=0.2, size=1)
        ax.set_title('Wolbachia titer by cell cycle phase')
        plt.tight_layout()
        plt.savefig(os.path.join(fig_dir, f'titer_by_cellcycle_{sample}_full.pdf'))
        plt.close()

    cells_per_cluster = combined.obs['leiden'].value_counts().sort_index()
    fig, ax = plt.subplots(figsize=(12, 6))
    ax.bar(range(len(cells_per_cluster)), cells_per_cluster.values,
           color=leiden_colors, alpha=0.7)
    ax.set_xticks(range(len(cells_per_cluster)))
    ax.set_xticklabels(cells_per_cluster.index)
    ax.set_ylabel('Number of cells')
    ax.set_title('Cell count by cluster (all samples)')
    plt.tight_layout()
    plt.savefig(os.path.join(fig_dir, f'cells_per_cluster_{sample}_full.pdf'))
    plt.close()
    print(f"\nCells per cluster:\n{cells_per_cluster}")


def main():
    parser = argparse.ArgumentParser(description='Compare library prep methods with detailed analysis')
    parser.add_argument('--files', required=True, nargs='+')
    parser.add_argument('--sample', type=str, default='method_comparison')
    parser.add_argument('--bio_condition', type=str, default=None)
    parser.add_argument('--batch_key', type=str, default='batch')
    parser.add_argument('--min_cells', type=int, default=3)
    parser.add_argument('--min_genes', type=int, default=200)
    parser.add_argument('--out_path', type=str, default='integrated.h5ad')
    parser.add_argument('--fig_dir', type=str, default='figures')
    parser.add_argument('--n_pcs', type=int, default=30)
    parser.add_argument('--resolution', type=float, default=0.5)
    args = parser.parse_args()

    integrate(
        files=args.files,
        out_path=args.out_path,
        fig_dir=args.fig_dir,
        sample=args.sample,
        batch_key=args.batch_key,
        min_cells=args.min_cells,
        min_genes=args.min_genes,
        n_pcs=args.n_pcs,
        bio_condition=args.bio_condition,
        leiden_resolution=args.resolution,
    )

if __name__ == "__main__":
    main()