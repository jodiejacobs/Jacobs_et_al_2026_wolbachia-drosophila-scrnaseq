import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import chi2_contingency
from scipy.stats import chi2 as chi2_dist
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score, homogeneity_completeness_v_measure
import sys
import os

# ── Load ──────────────────────────────────────────────────────────────────────
adata = sc.read_h5ad(
    "/private/groups/russelllab/jodie/scRNAseq/Jacobs_et_al_2026_wolbachia-drosophila-scrnaseq/"
    "cell_culture_system/results/integrated/integrated_uninfected_with_cellcycle.h5ad"
)

leiden_col = "leiden"
phase_col  = "cyclum_stage"
pseudo_col = "cyclum_pseudotime"

outdir = "figures/cell_cycle_leiden"
os.makedirs(outdir, exist_ok=True)

obs = adata.obs[[leiden_col, phase_col, pseudo_col]].copy()
obs[leiden_col] = obs[leiden_col].astype(str)
obs[phase_col]  = obs[phase_col].astype(str)

# ── 1. Clustering metrics ─────────────────────────────────────────────────────
ari  = adjusted_rand_score(obs[phase_col], obs[leiden_col])
nmi  = normalized_mutual_info_score(obs[phase_col], obs[leiden_col])
hom, comp, vms = homogeneity_completeness_v_measure(obs[phase_col], obs[leiden_col])

# ── 2. Contingency table + Chi-squared ───────────────────────────────────────
ct = pd.crosstab(obs[leiden_col], obs[phase_col])
chi2, _, dof, _ = chi2_contingency(ct)
cramers_v = np.sqrt(chi2 / (ct.values.sum() * (min(ct.shape) - 1)))

# precise p-value
pval_precise = chi2_dist.sf(chi2, dof)
log10_pval   = chi2_dist.logsf(chi2, dof) / np.log(10)
if pval_precise == 0.0:
    pval_str = f"< {sys.float_info.min:.2e}  (log10(p) = {log10_pval:.1f})"
else:
    pval_str = f"{pval_precise:.6e}  (log10(p) = {log10_pval:.1f})"

# ── 3. Per-cluster summary ────────────────────────────────────────────────────
ct_norm = ct.div(ct.sum(axis=1), axis=0)
summary = ct_norm.idxmax(axis=1).rename("dominant_phase").to_frame()
summary["dominant_phase_frac"] = ct_norm.max(axis=1)
summary["n_cells"] = ct.sum(axis=1)
summary["n_phases_present"] = (ct > 0).sum(axis=1)
summary = summary.sort_index(key=lambda x: x.astype(int))

# ── 4. Print to stdout ────────────────────────────────────────────────────────
print("=== Clustering Agreement Metrics ===")
print(f"  Adjusted Rand Index (ARI):          {ari:.4f}   (1=perfect, 0=random)")
print(f"  Normalized Mutual Info (NMI):       {nmi:.4f}   (1=perfect, 0=none)")
print(f"  Homogeneity:                        {hom:.4f}   (each cluster = 1 phase?)")
print(f"  Completeness:                       {comp:.4f}  (each phase = 1 cluster?)")
print(f"  V-measure:                          {vms:.4f}   (harmonic mean of above)")
print(f"\n=== Chi-squared Test ===")
print(f"  Chi2={chi2:.1f}, dof={dof}, p={pval_str}")
print(f"  Cramer's V: {cramers_v:.4f}  (effect size; 1=perfect association)")
print(f"\n=== Per-cluster Summary ===")
print(summary.to_string())

# ── 5. Write stats to text file ───────────────────────────────────────────────
with open(f"{outdir}/leiden_phase_stats.txt", "w") as f:
    f.write("=== Clustering Agreement Metrics ===\n")
    f.write(f"  Adjusted Rand Index (ARI):          {ari:.4f}   (1=perfect, 0=random)\n")
    f.write(f"  Normalized Mutual Info (NMI):       {nmi:.4f}   (1=perfect, 0=none)\n")
    f.write(f"  Homogeneity:                        {hom:.4f}   (each cluster = 1 phase?)\n")
    f.write(f"  Completeness:                       {comp:.4f}  (each phase = 1 cluster?)\n")
    f.write(f"  V-measure:                          {vms:.4f}   (harmonic mean of above)\n")
    f.write(f"\n=== Chi-squared Test ===\n")
    f.write(f"  Chi2={chi2:.1f}, dof={dof}, p={pval_str}\n")
    f.write(f"  Cramer's V: {cramers_v:.4f}  (effect size; 1=perfect association)\n")
    f.write(f"\n=== Per-cluster Summary ===\n")
    f.write(summary.to_string())
    f.write("\n")

print(f"\nStats written to {outdir}/leiden_phase_stats.txt")

# ── 6. Heatmap: phase composition per cluster ─────────────────────────────────
fig, axes = plt.subplots(1, 2, figsize=(14, 5))

sns.heatmap(
    ct_norm, annot=True, fmt=".2f", cmap="Blues",
    linewidths=0.5, ax=axes[0], cbar_kws={"label": "Fraction of cluster"}
)
axes[0].set_title("Phase composition per Leiden cluster\n(row-normalized)")
axes[0].set_xlabel("Cyclum stage"); axes[0].set_ylabel("Leiden cluster")

sns.heatmap(
    ct, annot=True, fmt="d", cmap="Oranges",
    linewidths=0.5, ax=axes[1], cbar_kws={"label": "Cell count"}
)
axes[1].set_title("Raw cell counts per Leiden × Phase")
axes[1].set_xlabel("Cyclum stage"); axes[1].set_ylabel("Leiden cluster")

plt.tight_layout()
plt.savefig(f"{outdir}/leiden_phase_heatmap.pdf", bbox_inches="tight")
plt.close()

# ── 7. Cyclum pseudotime distribution per cluster ────────────────────────────
fig, ax = plt.subplots(figsize=(10, 4))
clusters = sorted(obs[leiden_col].unique(), key=lambda x: int(x))
data_by_cluster = [obs.loc[obs[leiden_col] == c, pseudo_col].values for c in clusters]

ax.violinplot(data_by_cluster, positions=range(len(clusters)), showmedians=True)
ax.set_xticks(range(len(clusters)))
ax.set_xticklabels(clusters)
ax.set_xlabel("Leiden cluster")
ax.set_ylabel("Cyclum pseudotime")
ax.set_title("Cyclum pseudotime distribution per Leiden cluster")
plt.tight_layout()
plt.savefig(f"{outdir}/leiden_pseudotime_violin.pdf", bbox_inches="tight")
plt.close()

summary.to_csv(f"{outdir}/leiden_phase_summary.csv")
print(f"Figures saved to {outdir}/")