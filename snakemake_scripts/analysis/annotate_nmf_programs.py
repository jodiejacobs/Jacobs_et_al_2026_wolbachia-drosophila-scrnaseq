#!/usr/bin/env python3
"""
Annotate NMF programs with:
  1. Cell cycle score correlations (using existing scores in adata.obs)
  2. GSEA preranked enrichment against GO/KEGG gene sets
  3. FlyEnrichr API enrichment
  4. Summary heatmap of top annotations per program

Usage:
    python annotate_nmf_programs.py \
        --input data.h5ad \
        --program_dir /path/to/nmf_programs/ \
        --output_dir results/nmf_annotation \
        --mapping transcripts_to_genes.txt \
        --titer_var wolbachia_titer
"""

import argparse
import os
import glob
import time
import json

import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
import gseapy as gp
import requests

from scipy.stats import spearmanr, false_discovery_control


# ─────────────────────────────────────────────────────────────────────────────
# CLI
# ─────────────────────────────────────────────────────────────────────────────

def parse_args():
    parser = argparse.ArgumentParser(
        description="Annotate NMF programs via cell cycle correlations and GSEA",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument("--input",        required=True,  help="h5ad with NMF usage in obs")
    parser.add_argument("--program_dir",  required=True,  help="Dir containing Program_N_genes.csv files")
    parser.add_argument("--output_dir",   required=True,  help="Output directory")
    parser.add_argument("--mapping",      required=True,  help="transcripts_to_genes.txt (FBgn → symbol)")
    parser.add_argument("--titer_var",    default="wolbachia_titer",
                        help="Continuous titer column in obs (default: wolbachia_titer)")
    parser.add_argument("--cc_s_var",     default="S_score",
                        help="S phase score column (default: S_score)")
    parser.add_argument("--cc_g2m_var",   default="G2M_score",
                        help="G2M score column (default: G2M_score)")
    parser.add_argument("--cc_phase_var", default="phase",
                        help="Discrete phase column (default: phase)")
    parser.add_argument("--top_genes",    type=int, default=200,
                        help="Top N weighted genes per program for enrichment (default: 200)")
    parser.add_argument("--skip_gsea",    action="store_true",
                        help="Skip local GSEA (use FlyEnrichr only)")
    parser.add_argument("--skip_flyenrichr", action="store_true",
                        help="Skip FlyEnrichr API calls")
    return parser.parse_args()


# ─────────────────────────────────────────────────────────────────────────────
# Gene mapping
# ─────────────────────────────────────────────────────────────────────────────

def load_mapping(mapping_file):
    """Load FBgn → gene symbol from transcripts_to_genes.txt"""
    fbgn_to_sym = {}
    opener = __import__("gzip").open if mapping_file.endswith(".gz") else open
    mode = "rt" if mapping_file.endswith(".gz") else "r"
    with opener(mapping_file, mode) as fh:
        for line in fh:
            parts = line.strip().split("\t")
            if len(parts) >= 3 and parts[1].startswith("FBgn"):
                fbgn_to_sym[parts[1]] = parts[2]
    print(f"  Loaded {len(fbgn_to_sym):,} FBgn → symbol mappings")
    return fbgn_to_sym


# ─────────────────────────────────────────────────────────────────────────────
# Load program gene lists
# ─────────────────────────────────────────────────────────────────────────────

def load_program_genes(program_dir):
    """
    Load all Program_N_genes.csv files.
    Returns dict: program_name → DataFrame(gene, weight, rank)
    """
    csvs = sorted(glob.glob(os.path.join(program_dir, "Program_*_genes.csv")))
    programs = {}
    for path in csvs:
        name = os.path.basename(path).replace("_genes.csv", "")
        df = pd.read_csv(path)
        programs[name] = df
        print(f"  {name}: {len(df)} genes loaded")
    return programs


# ─────────────────────────────────────────────────────────────────────────────
# 1. Cell cycle correlations
# ─────────────────────────────────────────────────────────────────────────────

def cell_cycle_correlations(adata, output_dir, sample_name,
                             cc_s_var, cc_g2m_var, cc_phase_var, titer_var):
    """
    Correlate every Program_N usage score against:
      - S_score
      - G2M_score
      - wolbachia_titer (for reference)
    Also compute mean usage per cell cycle phase (boxplot).
    """
    print("\n" + "="*60)
    print("STEP 1: CELL CYCLE CORRELATIONS")
    print("="*60)

    program_cols = sorted(
        [c for c in adata.obs.columns if c.startswith("Program_")],
        key=lambda x: int(x.split("_")[1])
    )

    # Build list of continuous variables that actually exist
    cont_vars = {}
    for label, col in [("S_score", cc_s_var), ("G2M_score", cc_g2m_var),
                        ("titer", titer_var)]:
        if col in adata.obs.columns:
            cont_vars[label] = col
        else:
            print(f"  WARNING: '{col}' not found in obs, skipping")

    if not cont_vars:
        print("  No continuous variables found – skipping correlations")
        return None

    # ── Spearman correlations ────────────────────────────────────────────────
    rows = []
    for prog in program_cols:
        row = {"Program": prog}
        for label, col in cont_vars.items():
            valid = ~(adata.obs[prog].isna() | adata.obs[col].isna())
            r, p = spearmanr(adata.obs.loc[valid, col], adata.obs.loc[valid, prog])
            row[f"r_{label}"] = r
            row[f"p_{label}"] = p
        rows.append(row)

    corr_df = pd.DataFrame(rows)

    # FDR per variable
    for label in cont_vars:
        corr_df[f"FDR_{label}"] = false_discovery_control(corr_df[f"p_{label}"])

    corr_df.to_csv(os.path.join(output_dir, f"{sample_name}_cellcycle_correlations.csv"),
                   index=False)
    print(corr_df.to_string(index=False))

    # ── Heatmap of correlations ──────────────────────────────────────────────
    r_cols = [c for c in corr_df.columns if c.startswith("r_")]
    heat_data = corr_df.set_index("Program")[r_cols].rename(
        columns=lambda x: x.replace("r_", "")
    )

    fig, ax = plt.subplots(figsize=(max(4, len(r_cols) * 1.5), len(program_cols) * 0.5 + 1))
    sns.heatmap(
        heat_data, cmap="RdBu_r", center=0,
        annot=True, fmt=".2f", linewidths=0.5, ax=ax,
        cbar_kws={"label": "Spearman r"},
        vmin=-0.5, vmax=0.5
    )
    ax.set_title("NMF Program Correlations\n(S score / G2M score / Wolbachia titer)")
    ax.set_xlabel("")
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, f"{sample_name}_cellcycle_corr_heatmap.pdf"),
                dpi=300, bbox_inches="tight")
    plt.close()
    print(f"  Saved: {sample_name}_cellcycle_corr_heatmap.pdf")

    # ── Boxplots: usage per cell cycle phase ────────────────────────────────
    if cc_phase_var in adata.obs.columns:
        phases = sorted(adata.obs[cc_phase_var].dropna().unique())
        n_prog = len(program_cols)
        ncols = 4
        nrows = int(np.ceil(n_prog / ncols))

        fig, axes = plt.subplots(nrows, ncols,
                                 figsize=(ncols * 4, nrows * 3),
                                 sharey=False)
        axes = axes.flatten()

        palette = {"G1": "#4e9af1", "S": "#f5a623", "G2M": "#7ed321"}

        for i, prog in enumerate(program_cols):
            ax = axes[i]
            data_for_box = [
                adata.obs.loc[adata.obs[cc_phase_var] == ph, prog].dropna().values
                for ph in phases
            ]
            bp = ax.boxplot(data_for_box, labels=phases, patch_artist=True,
                            showfliers=False, widths=0.5)
            for patch, ph in zip(bp["boxes"], phases):
                patch.set_facecolor(palette.get(ph, "#aaa"))
                patch.set_alpha(0.7)
            ax.set_title(prog, fontsize=9)
            ax.set_ylabel("Usage", fontsize=7)
            ax.tick_params(labelsize=7)
            ax.grid(axis="y", alpha=0.3)

            # Annotate with r values
            for j, label in enumerate(["S_score", "G2M_score"]):
                if f"r_{label}" in corr_df.columns:
                    r_val = corr_df.loc[corr_df["Program"] == prog, f"r_{label}"].values
                    if len(r_val):
                        ax.text(0.98, 0.97 - j * 0.12,
                                f"{label.split('_')[0]} r={r_val[0]:.2f}",
                                transform=ax.transAxes, ha="right", va="top",
                                fontsize=6.5, color="#333")

        # Hide empty subplots
        for j in range(i + 1, len(axes)):
            axes[j].set_visible(False)

        legend_patches = [mpatches.Patch(color=v, alpha=0.7, label=k)
                          for k, v in palette.items() if k in phases]
        fig.legend(handles=legend_patches, loc="lower right", fontsize=9)

        plt.suptitle("NMF Program Usage by Cell Cycle Phase", fontsize=13, y=1.01)
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, f"{sample_name}_program_usage_by_phase.pdf"),
                    dpi=300, bbox_inches="tight")
        plt.close()
        print(f"  Saved: {sample_name}_program_usage_by_phase.pdf")

    return corr_df


# ─────────────────────────────────────────────────────────────────────────────
# 2. GSEA preranked (local, using gseapy)
# ─────────────────────────────────────────────────────────────────────────────

# Drosophila cell cycle gene sets (curated from GO:0007049 / FlyBase)
# These are gene symbols – used as a quick sanity check set
DROSOPHILA_CELLCYCLE_SETS = {
    "G1_S_transition": [
        "CycE", "dE2F1", "dDP", "rbf", "Cdc2", "CycD", "Cdk4",
        "stg", "twit", "pbl",
    ],
    "S_phase": [
        "PCNA", "RnrS", "RnrL", "mus209", "DNApol-alpha180", "DNApol-alpha50",
        "DNApol-delta", "DNApol-epsilon", "mcm2", "mcm3", "mcm5", "mcm7",
        "dup", "Orc1", "Orc2", "Orc5",
    ],
    "G2_M_transition": [
        "CycB", "Cdc2", "stg", "twine", "wee", "myt1",
    ],
    "M_phase_mitosis": [
        "aurA", "aurB", "polo", "BubR1", "Bub3", "mad2", "spd-2",
        "asp", "bw", "ncd", "Klp61F", "Klp67A", "pav", "tsr",
    ],
    "DNA_replication": [
        "Orc1", "Orc2", "Orc3", "Orc4", "Orc5", "Orc6",
        "dup", "mcm2", "mcm3", "mcm5", "mcm7",
    ],
}


def run_gsea_preranked(programs, fbgn_to_sym, output_dir, top_genes,
                       use_fly_sets=True):
    """
    For each program, rank genes by weight → run preranked GSEA against:
      - GO Biological Process (Fly)
      - KEGG Fly (if available)
      - Curated cell cycle gene sets

    Returns dict: program → results DataFrame
    """
    print("\n" + "="*60)
    print("STEP 2: GSEA PRERANKED ENRICHMENT")
    print("="*60)

    # Available Enrichr gene set libraries for Fly (gseapy)
    fly_libraries = [
        "GO_Biological_Process_2023",
        "KEGG_2021_Human",  # fallback – not ideal but broad coverage
    ]

    all_results = {}
    cc_summary = []   # track cell cycle hits per program

    for prog_name, prog_df in sorted(programs.items(),
                                      key=lambda x: int(x[0].split("_")[1])):
        print(f"\n  {prog_name} ({len(prog_df)} genes)")

        # Build ranked gene list: FBgn → symbol, weight as rank metric
        prog_df = prog_df.copy()

        # Separate FBgn genes and FBti (TE) elements
        fbgn_mask = prog_df["gene"].str.startswith("FBgn")
        fbti_mask = ~fbgn_mask

        n_te = fbti_mask.sum()
        if n_te > 0:
            print(f"    {n_te} transposable element loci (FBti) – excluded from GSEA")

        gene_df = prog_df[fbgn_mask].copy()
        gene_df["symbol"] = gene_df["gene"].map(fbgn_to_sym)
        gene_df = gene_df.dropna(subset=["symbol"])

        if len(gene_df) < 10:
            print(f"    Too few mapped genes ({len(gene_df)}), skipping GSEA")
            continue

        # Rank series: higher weight = higher rank
        rnk = gene_df.set_index("symbol")["weight"]
        # Remove duplicates (keep max weight)
        rnk = rnk.groupby(level=0).max().sort_values(ascending=False)

        print(f"    {len(rnk)} genes with symbols for GSEA")

        prog_results = {}

        # ── Curated cell cycle sets ──────────────────────────────────────────
        cc_hits = {}
        for cc_set_name, cc_genes in DROSOPHILA_CELLCYCLE_SETS.items():
            # Overlap with top genes
            top_set = set(rnk.head(top_genes).index)
            overlap = top_set & set(cc_genes)
            pct = len(overlap) / len(cc_genes) * 100 if cc_genes else 0
            cc_hits[cc_set_name] = {
                "overlap_n": len(overlap),
                "set_size": len(cc_genes),
                "pct_covered": pct,
                "genes": sorted(overlap),
            }

        cc_summary.append({
            "Program": prog_name,
            **{f"{k}_overlap": v["overlap_n"] for k, v in cc_hits.items()},
            **{f"{k}_pct": v["pct_covered"] for k, v in cc_hits.items()},
        })
        prog_results["curated_cc"] = cc_hits

        # ── gseapy preranked against GO ──────────────────────────────────────
        for lib in fly_libraries:
            try:
                pre_res = gp.prerank(
                    rnk=rnk,
                    gene_sets=lib,
                    threads=4,
                    min_size=10,
                    max_size=500,
                    permutation_num=100,   # keep fast; increase to 1000 for publication
                    outdir=None,
                    seed=42,
                    verbose=False,
                )
                res_df = pre_res.res2d
                if res_df is not None and len(res_df):
                    res_df = res_df.sort_values("FDR q-val")
                    prog_results[lib] = res_df
                    n_sig = (res_df["FDR q-val"] < 0.25).sum()
                    print(f"    {lib}: {n_sig} terms FDR<0.25")

                    # Save per-program per-library
                    out_path = os.path.join(
                        output_dir,
                        f"{prog_name}_gsea_{lib.replace(' ', '_')}.csv"
                    )
                    res_df.to_csv(out_path, index=False)

            except Exception as e:
                print(f"    WARNING: GSEA failed for {lib}: {e}")

        all_results[prog_name] = prog_results

    # ── Cell cycle summary table ─────────────────────────────────────────────
    if cc_summary:
        cc_df = pd.DataFrame(cc_summary)
        cc_df.to_csv(os.path.join(output_dir, "program_cellcycle_overlap.csv"),
                     index=False)
        print("\n  Cell cycle overlap summary (top genes vs curated sets):")
        pct_cols = [c for c in cc_df.columns if c.endswith("_pct")]
        print(cc_df[["Program"] + pct_cols].to_string(index=False))

    return all_results, cc_df if cc_summary else None


# ─────────────────────────────────────────────────────────────────────────────
# 3. FlyEnrichr API
# ─────────────────────────────────────────────────────────────────────────────

def flyenrichr_query(gene_symbols, library, description="query"):
    """Submit to FlyEnrichr and return results DataFrame."""
    SUBMIT_URL = "https://maayanlab.cloud/FlyEnrichr/addList"
    ENRICH_URL = "https://maayanlab.cloud/FlyEnrichr/enrich"

    payload = {
        "list": (None, "\n".join(gene_symbols)),
        "description": (None, description),
    }
    try:
        r = requests.post(SUBMIT_URL, files=payload, timeout=30)
        r.raise_for_status()
        uid = r.json()["userListId"]

        r2 = requests.get(f"{ENRICH_URL}?userListId={uid}&backgroundType={library}",
                          timeout=30)
        r2.raise_for_status()
        data = r2.json()

        if library not in data:
            return None

        rows = []
        for entry in data[library]:
            rows.append({
                "term":           entry[1],
                "p_value":        entry[2],
                "z_score":        entry[3],
                "combined_score": entry[4],
                "genes":          entry[5],
                "adj_p_value":    entry[6],
            })
        return pd.DataFrame(rows)

    except Exception as e:
        print(f"      FlyEnrichr error ({library}): {e}")
        return None


def run_flyenrichr(programs, fbgn_to_sym, output_dir, top_genes):
    """
    Run FlyEnrichr for each program's top genes.
    Focuses on cell cycle–relevant libraries.
    """
    print("\n" + "="*60)
    print("STEP 3: FLYENRICHR ANNOTATION")
    print("="*60)

    libraries = [
        "GO_Biological_Process_2018",
        "GO_Molecular_Function_2018",
        "KEGG_2019",
    ]

    all_rows = []

    for prog_name, prog_df in sorted(programs.items(),
                                      key=lambda x: int(x[0].split("_")[1])):
        print(f"\n  {prog_name}")

        # FBgn only, convert to symbols
        gene_df = prog_df[prog_df["gene"].str.startswith("FBgn")].copy()
        gene_df["symbol"] = gene_df["gene"].map(fbgn_to_sym)
        gene_df = gene_df.dropna(subset=["symbol"]).head(top_genes)
        symbols = gene_df["symbol"].tolist()

        if len(symbols) < 5:
            print(f"    Too few symbols ({len(symbols)}), skipping")
            continue

        for lib in libraries:
            print(f"    {lib} ...", end=" ", flush=True)
            res = flyenrichr_query(symbols, lib, description=f"{prog_name}_{lib}")
            if res is not None and len(res):
                res["program"] = prog_name
                res["library"] = lib
                all_rows.append(res)
                n_sig = (res["adj_p_value"] < 0.05).sum()
                print(f"{n_sig} significant terms")

                # Save per-program
                out = os.path.join(output_dir,
                                   f"{prog_name}_flyenrichr_{lib.split('_')[0]}.csv")
                res.to_csv(out, index=False)
            else:
                print("no results")
            time.sleep(0.3)   # be polite to the API

    if not all_rows:
        print("  No FlyEnrichr results obtained")
        return None

    combined = pd.concat(all_rows, ignore_index=True)
    combined.to_csv(os.path.join(output_dir, "all_flyenrichr_results.csv"), index=False)
    print(f"\n  Saved: all_flyenrichr_results.csv ({len(combined)} rows)")
    return combined


# ─────────────────────────────────────────────────────────────────────────────
# 4. Summary visualisation
# ─────────────────────────────────────────────────────────────────────────────

CELL_CYCLE_KEYWORDS = [
    "cell cycle", "mitosis", "mitotic", "meiosis", "meiotic",
    "dna replication", "dna repair", "s phase", "g1", "g2",
    "spindle", "chromosome segregation", "kinetochore",
    "cyclin", "cdk", "checkpoint", "rb", "e2f",
]


def is_cc_term(term):
    t = term.lower()
    return any(kw in t for kw in CELL_CYCLE_KEYWORDS)


def build_summary(corr_df, flyenrichr_df, cc_overlap_df, output_dir, sample_name):
    """
    Build and plot a combined annotation summary per program.
    """
    print("\n" + "="*60)
    print("STEP 4: ANNOTATION SUMMARY")
    print("="*60)

    if corr_df is None:
        print("  No correlation data available")
        return

    programs = corr_df["Program"].tolist()

    # ── Flag cell cycle programs ─────────────────────────────────────────────
    cc_flags = {}

    # Flag from FlyEnrichr GO terms
    if flyenrichr_df is not None:
        fly_cc = flyenrichr_df[
            (flyenrichr_df["adj_p_value"] < 0.05) &
            (flyenrichr_df["term"].apply(is_cc_term))
        ]
        for prog in programs:
            n_cc_terms = len(fly_cc[fly_cc["program"] == prog])
            cc_flags[prog] = cc_flags.get(prog, 0) + n_cc_terms

    # Flag from curated overlap
    if cc_overlap_df is not None:
        pct_cols = [c for c in cc_overlap_df.columns if c.endswith("_pct")]
        for _, row in cc_overlap_df.iterrows():
            if row[pct_cols].max() > 20:   # >20% overlap with any CC set
                cc_flags[row["Program"]] = cc_flags.get(row["Program"], 0) + 5

    # ── Scatter: S_score r vs G2M_score r, coloured by titer r ──────────────
    r_cols = {c.replace("r_", ""): c for c in corr_df.columns if c.startswith("r_")}

    if "S_score" in r_cols and "G2M_score" in r_cols:
        fig, ax = plt.subplots(figsize=(7, 6))

        s_r   = corr_df[r_cols["S_score"]].values
        g2m_r = corr_df[r_cols["G2M_score"]].values
        titer_r = (corr_df[r_cols.get("titer", r_cols["S_score"])].values
                   if "titer" in r_cols else np.zeros(len(corr_df)))

        sc_plot = ax.scatter(s_r, g2m_r, c=titer_r, cmap="RdBu_r",
                             s=120, vmin=-0.4, vmax=0.4, edgecolors="k",
                             linewidths=0.5, zorder=3)
        plt.colorbar(sc_plot, ax=ax, label="Spearman r (titer)")

        for i, prog in enumerate(corr_df["Program"]):
            prog_num = prog.split("_")[1]
            color = "red" if cc_flags.get(prog, 0) > 0 else "black"
            ax.annotate(prog_num, (s_r[i], g2m_r[i]),
                        fontsize=7, ha="center", va="bottom",
                        color=color, fontweight="bold" if color == "red" else "normal")

        ax.axhline(0, color="gray", lw=0.8, ls="--")
        ax.axvline(0, color="gray", lw=0.8, ls="--")
        ax.set_xlabel("Spearman r (S score)", fontsize=12)
        ax.set_ylabel("Spearman r (G2M score)", fontsize=12)
        ax.set_title("Program correlations with cell cycle scores\n"
                     "red labels = likely cell cycle programs", fontsize=11)
        ax.grid(alpha=0.2)

        red_patch = mpatches.Patch(color="red", label="Likely cell cycle program")
        ax.legend(handles=[red_patch], fontsize=9)

        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, f"{sample_name}_cc_scatter.pdf"),
                    dpi=300, bbox_inches="tight")
        plt.close()
        print(f"  Saved: {sample_name}_cc_scatter.pdf")

    # ── Print annotation table ───────────────────────────────────────────────
    print("\n  PROGRAM ANNOTATION SUMMARY")
    print("  " + "-"*70)

    header = f"{'Program':<14}"
    for col in [c for c in corr_df.columns if c.startswith("r_")]:
        header += f"  {col:<12}"
    header += "  CC_flag  Interpretation"
    print("  " + header)
    print("  " + "-"*70)

    for _, row in corr_df.iterrows():
        prog = row["Program"]
        line = f"{prog:<14}"
        for col in [c for c in corr_df.columns if c.startswith("r_")]:
            line += f"  {row[col]:>+.3f}      "

        cc_score = cc_flags.get(prog, 0)
        flag = "  *** CC  " if cc_score > 4 else ("  *  CC  " if cc_score > 0 else "         ")

        # Quick interpretation
        interp = ""
        if "r_S_score" in corr_df.columns and "r_G2M_score" in corr_df.columns:
            s  = row["r_S_score"]
            g  = row["r_G2M_score"]
            ti = row.get("r_titer", 0)
            if abs(s) > 0.3 or abs(g) > 0.3:
                phase = "S-phase" if abs(s) > abs(g) else "G2M"
                dirn  = "enriched" if max(s, g) > 0 else "depleted"
                interp = f"{phase} {dirn}"
            if abs(ti) > 0.2:
                interp += f" | titer {'↑' if ti > 0 else '↓'}"

        print(f"  {line}{flag}  {interp}")

    # ── Top cell cycle GO terms per program ──────────────────────────────────
    if flyenrichr_df is not None:
        print("\n  TOP CELL CYCLE GO TERMS PER PROGRAM (FlyEnrichr, adj_p<0.05)")
        print("  " + "-"*70)
        fly_cc_sig = flyenrichr_df[
            (flyenrichr_df["adj_p_value"] < 0.05) &
            (flyenrichr_df["term"].apply(is_cc_term))
        ].sort_values("adj_p_value")

        for prog in programs:
            subset = fly_cc_sig[fly_cc_sig["program"] == prog].head(3)
            if len(subset):
                print(f"\n  {prog}:")
                for _, r in subset.iterrows():
                    print(f"    [{r['library'].split('_')[0]:<4}] "
                          f"{r['term'][:60]:<60} p={r['adj_p_value']:.2e}")


# ─────────────────────────────────────────────────────────────────────────────
# Main
# ─────────────────────────────────────────────────────────────────────────────

def main():
    args = parse_args()
    os.makedirs(args.output_dir, exist_ok=True)
    sample_name = os.path.basename(args.output_dir).rstrip("/") or "nmf"

    # ── Load ─────────────────────────────────────────────────────────────────
    print("="*60)
    print("LOADING DATA")
    print("="*60)
    adata = sc.read_h5ad(args.input)
    print(f"  {adata.n_obs:,} cells × {adata.n_vars:,} genes")

    fbgn_to_sym = load_mapping(args.mapping)
    programs    = load_program_genes(args.program_dir)
    print(f"  {len(programs)} programs loaded")

    # ── Step 1: Cell cycle correlations ──────────────────────────────────────
    corr_df = cell_cycle_correlations(
        adata, args.output_dir, sample_name,
        args.cc_s_var, args.cc_g2m_var, args.cc_phase_var, args.titer_var,
    )

    # ── Step 2: GSEA preranked ────────────────────────────────────────────────
    cc_overlap_df = None
    if not args.skip_gsea:
        _, cc_overlap_df = run_gsea_preranked(
            programs, fbgn_to_sym, args.output_dir, args.top_genes
        )

    # ── Step 3: FlyEnrichr ───────────────────────────────────────────────────
    flyenrichr_df = None
    if not args.skip_flyenrichr:
        flyenrichr_df = run_flyenrichr(
            programs, fbgn_to_sym, args.output_dir, args.top_genes
        )

    # ── Step 4: Summary ──────────────────────────────────────────────────────
    build_summary(corr_df, flyenrichr_df, cc_overlap_df,
                  args.output_dir, sample_name)

    print("\n" + "="*60)
    print("DONE")
    print("="*60)
    print(f"\nOutputs in: {args.output_dir}/")
    print("  program_cellcycle_overlap.csv   – curated CC gene set overlap")
    print("  *_cellcycle_corr_heatmap.pdf    – S/G2M/titer correlation heatmap")
    print("  *_program_usage_by_phase.pdf    – usage boxplots per phase")
    print("  *_cc_scatter.pdf                – S vs G2M correlation scatter")
    print("  all_flyenrichr_results.csv      – all FlyEnrichr hits")
    print("  Program_N_gsea_*.csv            – per-program GSEA results")


if __name__ == "__main__":
    main()
