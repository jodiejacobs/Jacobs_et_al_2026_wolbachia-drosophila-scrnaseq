#!/usr/bin/env python3
"""
run_velocity.py

RNA velocity (scVelo dynamical model) on kb nac (nascent/mature) counts,
aligned to the QC-filtered cells of results/integrated/integrated.h5ad.

Pseudotime-independent by design:
  - reads integrated.h5ad (upstream of SCEPTIC), never the SCEPTIC output
  - drops any obs column containing "pseudotime" before analysis
  - the kNN graph is built from X_pca_harmony (default) or from a PCA +
    Harmony embedding recomputed on the spliced layer (--neighbor-rep recompute)
  - the root/terminal states for latent time come from the velocity graph

Usage:
    python run_velocity.py \
        --integrated results/integrated/integrated.h5ad \
        --nac-files results/nac/*/adata_nac.h5ad \
        --platform pipseq \
        --out-h5ad results/velocity/velocity_pipseq.h5ad \
        --fig-dir results/velocity/pipseq \
        --n-jobs 16
"""
import argparse
import os

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
import scvelo as scv

SPLICED_KEYS = ["spliced", "mature"]
UNSPLICED_KEYS = ["unspliced", "nascent"]


def pick_layer(adata, candidates, path):
    for k in candidates:
        if k in adata.layers:
            return k
    raise KeyError(f"{path}: none of {candidates} in layers {list(adata.layers)}")


def load_nac(nac_files, obs):
    """Return spliced/unspliced CSR matrices in the row order of `obs`."""
    barcodes = obs.index.str.rsplit("-", n=1).str[0]
    per_sample, genes = {}, None

    for f in nac_files:
        sid = os.path.basename(os.path.dirname(f))
        mask = (obs["batch"] == sid).values
        if not mask.any():
            continue
        nac = ad.read_h5ad(f)
        s_key = pick_layer(nac, SPLICED_KEYS, f)
        u_key = pick_layer(nac, UNSPLICED_KEYS, f)
        idx = nac.obs_names.get_indexer(barcodes[mask])
        n_miss = int((idx < 0).sum())
        print(f"  {sid}: layers {s_key}/{u_key}, {mask.sum()} cells, {n_miss} unmatched barcodes")
        if n_miss / mask.sum() > 0.10:
            raise ValueError(f"{sid}: >10% of filtered barcodes missing from {f}; check barcode format")
        per_sample[sid] = (mask, idx, nac, s_key, u_key)
        genes = nac.var_names if genes is None else genes.union(nac.var_names)

    if not per_sample:
        raise ValueError("No nac file matched any obs['batch'] value; check sample IDs")
    genes = pd.Index(genes)
    n = obs.shape[0]
    blocks_s, blocks_u, rows = [], [], []

    for sid, (mask, idx, nac, s_key, u_key) in per_sample.items():
        col = genes.get_indexer(nac.var_names)
        keep = idx >= 0
        row_pos = np.where(mask)[0][keep]
        for key, store in ((s_key, blocks_s), (u_key, blocks_u)):
            m = sp.csr_matrix(nac.layers[key])[idx[keep]]
            # remap columns to the union gene index
            m = m.tocoo()
            store.append(sp.csr_matrix((m.data, (m.row, col[m.col])),
                                       shape=(len(row_pos), len(genes))))
        rows.append(row_pos)

    order = np.concatenate(rows)
    S_all = sp.vstack(blocks_s).tocsr()
    U_all = sp.vstack(blocks_u).tocsr()
    # place stacked rows back into obs order (unmatched cells stay all-zero)
    perm = sp.csr_matrix((np.ones(len(order)), (order, np.arange(len(order)))),
                         shape=(n, len(order)))
    S = (perm @ S_all).astype(np.float32).tocsr()
    U = (perm @ U_all).astype(np.float32).tocsr()
    matched = np.zeros(n, dtype=bool)
    matched[order] = True
    return S, U, genes, matched


def recompute_embedding(vdata, batch_vars, n_pcs=30):
    import harmonypy
    tmp = ad.AnnData(X=vdata.layers["spliced"].copy(), obs=vdata.obs.copy())
    tmp.var_names = vdata.var_names
    sc.pp.highly_variable_genes(tmp, flavor="seurat_v3", n_top_genes=2000)
    sc.pp.normalize_total(tmp)
    sc.pp.log1p(tmp)
    tmp = tmp[:, tmp.var["highly_variable"]].copy()
    sc.pp.scale(tmp, max_value=10)
    sc.pp.pca(tmp, n_comps=n_pcs)
    ho = harmonypy.run_harmony(tmp.obsm["X_pca"], tmp.obs, batch_vars, max_iter_harmony=20)
    Z = np.asarray(ho.Z_corr)
    return Z if Z.shape[0] == tmp.n_obs else Z.T  # harmonypy 2.x returns cells x PCs


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--integrated", required=True)
    p.add_argument("--nac-files", nargs="+", required=True)
    p.add_argument("--platform", required=True, help="pipseq, 10x, or all")
    p.add_argument("--out-h5ad", required=True)
    p.add_argument("--fig-dir", required=True)
    p.add_argument("--neighbor-rep", default="X_pca_harmony",
                   help="obsm key, or 'recompute' to build PCA+Harmony from the spliced layer")
    p.add_argument("--n-neighbors", type=int, default=30)
    p.add_argument("--n-top-genes", type=int, default=2000)
    p.add_argument("--min-shared-counts", type=int, default=20)
    p.add_argument("--cluster-key", default="leiden")
    p.add_argument("--n-jobs", type=int, default=8)
    args = p.parse_args()

    os.makedirs(args.fig_dir, exist_ok=True)
    scv.settings.figdir = args.fig_dir
    scv.settings.verbosity = 3

    print(f"Loading {args.integrated}")
    adata = sc.read_h5ad(args.integrated)
    if args.platform != "all":
        adata = adata[adata.obs["method"] == args.platform].copy()
    print(f"  {adata.n_obs} cells for platform={args.platform}")

    # trajectory pseudotimes only; cyclum (cell-cycle) pseudotime is kept
    pt_cols = [c for c in adata.obs.columns
               if "pseudotime" in c.lower() and not c.lower().startswith("cyclum")]
    if pt_cols:
        print(f"  Dropping pseudotime columns to keep velocity independent: {pt_cols}")
        adata.obs = adata.obs.drop(columns=pt_cols)

    print("Loading nac matrices")
    S, U, genes, matched = load_nac(args.nac_files, adata.obs)
    if not matched.all():
        print(f"  Dropping {(~matched).sum()} cells with no nac barcode match")

    vdata = ad.AnnData(X=S.copy(), obs=adata.obs.copy(),
                       var=pd.DataFrame(index=genes),
                       layers={"spliced": S, "unspliced": U})
    for key in ("X_pca_harmony", "X_umap"):
        if key in adata.obsm:
            vdata.obsm[key] = adata.obsm[key]
    vdata = vdata[matched].copy()

    # Wolbachia genes are intronless: nascent == mature, so velocity is undefined
    wmel = vdata.var_names.str.startswith("GQX67")
    print(f"  Removing {wmel.sum()} wMel (GQX67*) features")
    vdata = vdata[:, ~wmel].copy()

    # ── QC: unspliced fraction per sample, before any filtering ──────────────
    s_tot = np.asarray(vdata.layers["spliced"].sum(axis=1)).ravel()
    u_tot = np.asarray(vdata.layers["unspliced"].sum(axis=1)).ravel()
    vdata.obs["unspliced_fraction"] = u_tot / np.maximum(s_tot + u_tot, 1)
    qc = (vdata.obs.groupby("batch", observed=True)["unspliced_fraction"]
          .agg(["count", "median", "mean"]).rename(columns={"count": "n_cells"}))
    qc.to_csv(os.path.join(args.fig_dir, "unspliced_fraction_by_sample.csv"))
    print(qc.to_string())

    # ── Neighbour graph (no pseudotime input) ────────────────────────────────
    if args.neighbor_rep == "recompute":
        batch_vars = ["replicate"] if args.platform != "all" else ["method", "replicate"]
        vdata.obsm["X_harmony_spliced"] = recompute_embedding(vdata, batch_vars)
        rep = "X_harmony_spliced"
    else:
        rep = args.neighbor_rep
        if rep not in vdata.obsm:
            raise KeyError(f"{rep} not in obsm {list(vdata.obsm)}")
    print(f"  Neighbour graph from {rep}")

    # ── scVelo dynamical model ───────────────────────────────────────────────
    # scVelo >=0.3.4 dropped log1p / HVG from filter_and_normalize; do those in scanpy
    scv.pp.filter_and_normalize(vdata, min_shared_counts=args.min_shared_counts)
    sc.pp.log1p(vdata)
    if vdata.n_vars > args.n_top_genes:
        sc.pp.highly_variable_genes(vdata, n_top_genes=args.n_top_genes, flavor="seurat")
        vdata = vdata[:, vdata.var["highly_variable"]].copy()
    scv.pp.moments(vdata, n_neighbors=args.n_neighbors, use_rep=rep)
    scv.tl.recover_dynamics(vdata, n_jobs=args.n_jobs)
    scv.tl.velocity(vdata, mode="dynamical")
    scv.tl.velocity_graph(vdata, n_jobs=args.n_jobs)
    scv.tl.velocity_confidence(vdata)
    scv.tl.velocity_pseudotime(vdata)
    scv.tl.latent_time(vdata)

    n_vg = int(vdata.var["velocity_genes"].sum())
    print(f"  Velocity genes: {n_vg}")
    vdata.var[["fit_likelihood", "fit_alpha", "fit_beta", "fit_gamma", "velocity_genes"]] \
        .sort_values("fit_likelihood", ascending=False) \
        .to_csv(os.path.join(args.fig_dir, "velocity_gene_fits.csv"))

    # ── Figures ──────────────────────────────────────────────────────────────
    if "X_umap" in vdata.obsm:
        color = args.cluster_key if args.cluster_key in vdata.obs else None
        scv.pl.velocity_embedding_stream(vdata, basis="umap", color=color, show=False,
                                         save=f"stream_{args.platform}.pdf")
        scv.pl.scatter(vdata, basis="umap",
                       color=["latent_time", "velocity_pseudotime", "velocity_confidence",
                              "unspliced_fraction"],
                       show=False, save=f"latent_time_{args.platform}.pdf")

    vdata.write_h5ad(args.out_h5ad)
    print(f"Wrote {args.out_h5ad}")


if __name__ == "__main__":
    main()
