#!/usr/bin/env python3
"""
run_scribe.py

Velocity-coupled causal network inference with Scribe (Qiu et al. 2020,
Cell Systems), using Scribe.Scribe.causal_net_dynamics_coupling.

For each regulator x and target y, Scribe computes the restricted directed
information RDI(x -> y) = I(x_t ; y_{t+1} | y_t), where x_t and y_t are
kNN-smoothed spliced expression (layer "Ms") and y_{t+1} = Ms + velocity.

Gene selection uses velocity fits only (no pseudotime):
  targets    = top --n-targets velocity genes by scVelo fit_likelihood
  regulators = genes in --regulators (e.g. FlyTF) that are velocity genes,
               capped at --max-regulators by fit_likelihood

Modes
  --mode select : write the regulator and target lists
  --mode chunk  : RDI for one chunk of regulators (parallel Snakemake jobs)
  --mode merge  : combine chunks, add CLR-normalised scores, write edge list

Scribe-py is not pip-installable from its current packaging; put a clone
on PYTHONPATH (the Snakemake rule does this from config["scribe_repo"]).
"""
import argparse
import glob
import os

import numpy as np
import pandas as pd


def read_list(path):
    with open(path) as fh:
        return [l.strip() for l in fh if l.strip() and not l.startswith("#")]


def select_genes(args):
    import anndata as ad
    vdata = ad.read_h5ad(args.velocity, backed="r")
    var = vdata.var
    vg = var[var["velocity_genes"].astype(bool)].sort_values("fit_likelihood", ascending=False)
    exclude = set(read_list(args.exclude)) if args.exclude else set()
    vg = vg[~vg.index.isin(exclude)]
    print(f"Velocity genes after exclusions: {len(vg)}")

    targets = vg.index[: args.n_targets].tolist()
    tf = set(read_list(args.regulators))
    regs = [g for g in vg.index if g in tf][: args.max_regulators]
    print(f"Regulators: {len(regs)} (of {len(tf)} in list); targets: {len(targets)}")
    if not regs:
        raise ValueError("No regulators overlap velocity genes; check that the list uses the "
                         "same gene IDs as adata.var_names (e.g. FBgn IDs)")
    os.makedirs(args.outdir, exist_ok=True)
    pd.Series(regs).to_csv(os.path.join(args.outdir, "regulators.txt"), index=False, header=False)
    pd.Series(targets).to_csv(os.path.join(args.outdir, "targets.txt"), index=False, header=False)


def run_chunk(args):
    import anndata as ad
    from Scribe.Scribe import causal_net_dynamics_coupling

    regs = read_list(os.path.join(args.outdir, "regulators.txt"))
    targets = read_list(os.path.join(args.outdir, "targets.txt"))
    chunk = np.array_split(np.array(regs, dtype=object), args.n_chunks)[args.chunk].tolist()
    out = os.path.join(args.outdir, "chunks", f"rdi_chunk{args.chunk:03d}.csv")
    os.makedirs(os.path.dirname(out), exist_ok=True)
    if not chunk:
        pd.DataFrame(columns=targets).to_csv(out)
        return

    vdata = ad.read_h5ad(args.velocity)
    genes = sorted(set(chunk) | set(targets))
    vdata = vdata[:, genes].copy()
    if args.max_cells and vdata.n_obs > args.max_cells:
        rng = np.random.default_rng(args.seed)  # same seed -> same cells in every chunk
        keep = np.sort(rng.choice(vdata.n_obs, args.max_cells, replace=False))
        vdata = vdata[keep].copy()
    print(f"Chunk {args.chunk}: {len(chunk)} regulators x {len(targets)} targets, {vdata.n_obs} cells")

    causal_net_dynamics_coupling(vdata, TFs=chunk, Targets=targets,
                                 t0_key="Ms", t1_key="velocity",
                                 normalize=True, drop_zero_cells=True)
    vdata.uns["causal_net"]["RDI"].astype(float).to_csv(out)


def clr(mat):
    """Context likelihood of relatedness, z-scoring rows and columns (as in Scribe.CLR,
    with explicit row/column broadcasting)."""
    m = mat.to_numpy(dtype=float)
    zr = (m - m.mean(axis=1, keepdims=True)) / m.std(axis=1, ddof=1, keepdims=True)
    zc = (m - m.mean(axis=0, keepdims=True)) / m.std(axis=0, ddof=1, keepdims=True)
    zr, zc = np.clip(np.nan_to_num(zr), 0, None), np.clip(np.nan_to_num(zc), 0, None)
    return pd.DataFrame(np.sqrt((zr ** 2 + zc ** 2) / 2), index=mat.index, columns=mat.columns)


def merge(args):
    files = sorted(glob.glob(os.path.join(args.outdir, "chunks", "rdi_chunk*.csv")))
    rdi = pd.concat([pd.read_csv(f, index_col=0) for f in files])
    rdi = rdi.loc[~rdi.index.duplicated()]
    # self-pairs are skipped by Scribe; mask them so they do not enter CLR
    for g in rdi.index.intersection(rdi.columns):
        rdi.loc[g, g] = np.nan
    rdi_clr = clr(rdi.fillna(rdi.stack().median()))
    for g in rdi.index.intersection(rdi.columns):
        rdi_clr.loc[g, g] = np.nan

    rdi.to_csv(os.path.join(args.outdir, "rdi_matrix.csv"))
    rdi_clr.to_csv(os.path.join(args.outdir, "rdi_clr_matrix.csv"))

    edges = (rdi.stack().rename("rdi").to_frame()
             .join(rdi_clr.stack().rename("rdi_clr"))
             .reset_index().rename(columns={"level_0": "regulator", "level_1": "target"}))
    edges.columns = ["regulator", "target", "rdi", "rdi_clr"]
    edges = edges.sort_values("rdi_clr", ascending=False)
    edges.to_csv(os.path.join(args.outdir, "edges.csv"), index=False)
    print(f"{len(files)} chunks, {rdi.shape[0]} regulators x {rdi.shape[1]} targets")
    print(edges.head(20).to_string(index=False))


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--mode", choices=["select", "chunk", "merge"], required=True)
    p.add_argument("--velocity")
    p.add_argument("--outdir", required=True)
    p.add_argument("--regulators", help="file of regulator gene IDs, one per line")
    p.add_argument("--exclude", help="optional file of gene IDs to exclude (e.g. cell-cycle genes)")
    p.add_argument("--n-targets", type=int, default=300)
    p.add_argument("--max-regulators", type=int, default=100)
    p.add_argument("--chunk", type=int)
    p.add_argument("--n-chunks", type=int, default=10)
    p.add_argument("--max-cells", type=int, default=5000)
    p.add_argument("--seed", type=int, default=0)
    args = p.parse_args()
    {"select": select_genes, "chunk": run_chunk, "merge": merge}[args.mode](args)


if __name__ == "__main__":
    main()
