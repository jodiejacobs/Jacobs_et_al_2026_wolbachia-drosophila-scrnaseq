import scanpy as sc
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np

adata = sc.read_h5ad("integrated.h5ad")

# Build a consistent color palette for leiden clusters
leiden_cats = adata.obs["leiden"].cat.categories.tolist()
n_clusters = len(leiden_cats)

# Use scanpy's default palette if available, otherwise generate one
if "leiden_colors" in adata.uns:
    cluster_colors = {c: col for c, col in zip(leiden_cats, adata.uns["leiden_colors"])}
else:
    cmap = plt.get_cmap("tab20", n_clusters)
    cluster_colors = {c: mcolors.to_hex(cmap(i)) for i, c in enumerate(leiden_cats)}

GREY = "#D3D3D3"
GREY_ALPHA = 0.3
POINT_SIZE = 8

umap_coords = adata.obsm["X_umap"]
samples = adata.obs["bio_condition"].cat.categories.tolist() \
    if hasattr(adata.obs["bio_condition"], "cat") \
    else adata.obs["bio_condition"].unique().tolist()

for sample in samples:
    fig, ax = plt.subplots(figsize=(6, 5))

    in_sample = adata.obs["bio_condition"] == sample
    out_sample = ~in_sample

    # Plot background (other samples) first
    ax.scatter(
        umap_coords[out_sample, 0],
        umap_coords[out_sample, 1],
        c=GREY,
        s=POINT_SIZE,
        alpha=GREY_ALPHA,
        linewidths=0,
        rasterized=True,
        label="_nolegend_",
    )

    # Plot foreground (focal sample) colored by leiden cluster
    for cluster in leiden_cats:
        mask = in_sample & (adata.obs["leiden"] == cluster)
        if mask.sum() == 0:
            continue
        ax.scatter(
            umap_coords[mask, 0],
            umap_coords[mask, 1],
            c=cluster_colors[cluster],
            s=POINT_SIZE,
            alpha=0.85,
            linewidths=0,
            rasterized=True,
            label=cluster,
        )

    ax.set_xlabel("UMAP 1", fontsize=10)
    ax.set_ylabel("UMAP 2", fontsize=10)
    ax.set_title(sample, fontsize=12)
    ax.tick_params(labelsize=8)
    ax.set_aspect("equal")

    # Legend — only show clusters present in this sample
    handles, labels = ax.get_legend_handles_labels()
    if handles:
        ax.legend(
            handles,
            labels,
            title="Leiden",
            bbox_to_anchor=(1.02, 1),
            loc="upper left",
            fontsize=7,
            title_fontsize=8,
            markerscale=1.5,
            frameon=False,
        )

    plt.tight_layout()
    safe_name = sample.replace("/", "_").replace(" ", "_")
    out_path = f"umap_leiden_{safe_name}_by_sample.pdf"
    fig.savefig(out_path, bbox_inches="tight", dpi=150)
    plt.close(fig)
    print(f"Saved: {out_path}")

print("Done.")

# --- Panel figure ---
import math

n_samples = len(samples)
ncols = min(4, n_samples)
nrows = math.ceil(n_samples / ncols)

fig, axes = plt.subplots(nrows, ncols, figsize=(ncols * 4, nrows * 3.5))
axes = np.array(axes).flatten()

for ax_idx, sample in enumerate(samples):
    ax = axes[ax_idx]

    in_sample = adata.obs["bio_condition"] == sample
    out_sample = ~in_sample

    ax.scatter(
        umap_coords[out_sample, 0],
        umap_coords[out_sample, 1],
        c=GREY,
        s=3,
        alpha=GREY_ALPHA,
        linewidths=0,
        rasterized=True,
    )

    for cluster in leiden_cats:
        mask = in_sample & (adata.obs["leiden"] == cluster)
        if mask.sum() == 0:
            continue
        ax.scatter(
            umap_coords[mask, 0],
            umap_coords[mask, 1],
            c=cluster_colors[cluster],
            s=3,
            alpha=0.85,
            linewidths=0,
            rasterized=True,
            label=cluster,
        )

    ax.set_title(sample, fontsize=9, pad=3)
    ax.set_xlabel("UMAP 1", fontsize=7)
    ax.set_ylabel("UMAP 2", fontsize=7)
    ax.tick_params(labelsize=6)
    ax.set_aspect("equal")

# Hide unused axes
for ax in axes[n_samples:]:
    ax.set_visible(False)

# Shared legend from last populated axis
handles, labels = axes[n_samples - 1].get_legend_handles_labels()
if handles:
    fig.legend(
        handles,
        labels,
        title="Leiden",
        loc="lower right",
        bbox_to_anchor=(1.0, 0.0),
        fontsize=7,
        title_fontsize=8,
        markerscale=2,
        frameon=False,
    )

fig.suptitle("UMAP by sample (Leiden clusters)", fontsize=12, y=1.01)
plt.tight_layout()
fig.savefig("umap_leiden_panel_by_sample.pdf", bbox_inches="tight", dpi=150)
plt.close(fig)
print("Saved: umap_leiden_panel_by_sample.pdf")
