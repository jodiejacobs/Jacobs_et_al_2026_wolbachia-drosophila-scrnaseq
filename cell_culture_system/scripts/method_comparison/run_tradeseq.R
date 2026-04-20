#!/usr/bin/env Rscript
# run_tradeseq.R
# ==============
# Standalone tradeSeq GAM analysis for single-trajectory pseudotime data.
# Reads counts CSV + pseudotime CSV written by pseudotime_gene_importance.py.
#
# Usage:
#   Rscript run_tradeseq.R \
#       --counts   results/.../tradeseq_inputs/counts_genesXcells.csv \
#       --pt       results/.../tradeseq_inputs/pseudotime.csv \
#       --outdir   results/.../wolbachia_infection \
#       --nknots   6 \
#       --nworkers 16
#
# Outputs (written to --outdir):
#   tradeseq_sce.rds                  : fitted SingleCellExperiment (re-use for tests)
#   tradeseq_association.csv          : associationTest results (non-constant along PT)
#   tradeseq_startvsend.csv           : startVsEndTest results
#   tradeseq_smooth_predictions.csv   : predictSmooth for top 300 sig genes

suppressPackageStartupMessages({
    library(tradeSeq)
    library(BiocParallel)
    library(Matrix)
})

# ─────────────────────────────────────────────────────────────────────────────
# Argument parsing
# ─────────────────────────────────────────────────────────────────────────────

args <- commandArgs(trailingOnly = TRUE)

parse_arg <- function(flag, default = NULL) {
    idx <- which(args == flag)
    if (length(idx) == 0) {
        if (is.null(default)) stop(sprintf("Required argument %s not provided", flag))
        return(default)
    }
    return(args[idx + 1])
}

counts_file <- parse_arg("--counts")
pt_file     <- parse_arg("--pt")
out_dir     <- parse_arg("--outdir")
n_knots     <- as.integer(parse_arg("--nknots",   "6"))
n_workers   <- as.integer(parse_arg("--nworkers", "8"))

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

cat("=== tradeSeq standalone ===\n")
cat(sprintf("  counts   : %s\n", counts_file))
cat(sprintf("  pt       : %s\n", pt_file))
cat(sprintf("  outdir   : %s\n", out_dir))
cat(sprintf("  nknots   : %d\n", n_knots))
cat(sprintf("  nworkers : %d\n", n_workers))

# ─────────────────────────────────────────────────────────────────────────────
# Load data
# ─────────────────────────────────────────────────────────────────────────────

cat("\n[1/5] Loading counts ...\n")
counts <- read.csv(counts_file, row.names = 1, check.names = FALSE)
counts <- as.matrix(counts)
cat(sprintf("  Shape: %d genes x %d cells\n", nrow(counts), ncol(counts)))
cat(sprintf("  Count range: %.1f – %.1f\n", min(counts), max(counts)))

# Round to integers (in case of expm1 recovery producing near-integers)
counts <- round(counts)
storage.mode(counts) <- "integer"

cat("\n[2/5] Loading pseudotime ...\n")
pt_df  <- read.csv(pt_file, row.names = 1)
pt_vec <- setNames(pt_df[[1]], rownames(pt_df))
cat(sprintf("  Cells with pseudotime: %d\n", length(pt_vec)))
cat(sprintf("  Pseudotime range: %.3f – %.3f\n", min(pt_vec), max(pt_vec)))

# ─────────────────────────────────────────────────────────────────────────────
# Align cells
# ─────────────────────────────────────────────────────────────────────────────

shared <- intersect(colnames(counts), names(pt_vec))
cat(sprintf("  Shared cells: %d\n", length(shared)))
if (length(shared) == 0) stop("No shared cells between counts and pseudotime!")

counts <- counts[, shared, drop = FALSE]
pt_vec <- pt_vec[shared]

# Remove genes with all-zero counts (would cause fitGAM to fail)
gene_sums <- rowSums(counts)
n_zero <- sum(gene_sums == 0)
if (n_zero > 0) {
    cat(sprintf("  Removing %d all-zero genes\n", n_zero))
    counts <- counts[gene_sums > 0, , drop = FALSE]
}

# Remove genes expressed in fewer than 5 cells
n_cells_expr <- rowSums(counts > 0)
n_sparse <- sum(n_cells_expr < 5)
if (n_sparse > 0) {
    cat(sprintf("  Removing %d genes expressed in <5 cells\n", n_sparse))
    counts <- counts[n_cells_expr >= 5, , drop = FALSE]
}

cat(sprintf("  Final: %d genes x %d cells\n", nrow(counts), ncol(counts)))

pt_mat <- matrix(pt_vec, ncol = 1,
                 dimnames = list(shared, "pseudotime"))
wt_mat <- matrix(1, nrow = length(shared), ncol = 1,
                 dimnames = list(shared, "w1"))

# ─────────────────────────────────────────────────────────────────────────────
# Fit GAMs
# ─────────────────────────────────────────────────────────────────────────────

cat(sprintf("\n[3/5] Fitting GAMs (%d genes, %d knots, %d workers) ...\n",
            nrow(counts), n_knots, n_workers))
cat("  This is the slow step — expect 1-4 hrs depending on gene count\n")

set.seed(42)
BPPARAM <- MulticoreParam(workers   = n_workers,
                          progressbar = TRUE)

sce <- tryCatch(
    fitGAM(counts      = counts,
           pseudotime  = pt_mat,
           cellWeights = wt_mat,
           nknots      = n_knots,
           parallel    = TRUE,
           BPPARAM     = BPPARAM,
           verbose     = FALSE),
    error = function(e) {
        cat(sprintf("  ERROR in fitGAM: %s\n", conditionMessage(e)))
        quit(status = 1)
    }
)

sce_path <- file.path(out_dir, "tradeseq_sce.rds")
saveRDS(sce, sce_path)
cat(sprintf("  Saved SCE: %s\n", sce_path))

# ─────────────────────────────────────────────────────────────────────────────
# Association test — is expression non-constant along pseudotime?
# ─────────────────────────────────────────────────────────────────────────────

cat("\n[4/5] Running tests ...\n")
cat("  associationTest ...\n")
assoc      <- tryCatch(
    as.data.frame(associationTest(sce, lineages = FALSE)),
    error = function(e) { cat(sprintf("  ERROR: %s\n", conditionMessage(e))); NULL }
)

if (!is.null(assoc)) {
    assoc$gene <- rownames(assoc)
    assoc$padj <- p.adjust(assoc$pvalue, method = "BH")
    assoc      <- assoc[order(assoc$waldStat, decreasing = TRUE), ]
    write.csv(assoc, file.path(out_dir, "tradeseq_association.csv"))
    n_sig <- sum(assoc$padj < 0.05, na.rm = TRUE)
    cat(sprintf("  Significant genes (padj<0.05): %d / %d\n", n_sig, nrow(assoc)))
}

cat("  startVsEndTest ...\n")
sve <- tryCatch(
    as.data.frame(startVsEndTest(sce)),
    error = function(e) { cat(sprintf("  ERROR: %s\n", conditionMessage(e))); NULL }
)

if (!is.null(sve)) {
    sve$gene <- rownames(sve)
    sve$padj <- p.adjust(sve$pvalue, method = "BH")
    sve      <- sve[order(sve$waldStat, decreasing = TRUE), ]
    write.csv(sve, file.path(out_dir, "tradeseq_startvsend.csv"))
    cat(sprintf("  Significant genes (padj<0.05): %d / %d\n",
                sum(sve$padj < 0.05, na.rm = TRUE), nrow(sve)))
}

# ─────────────────────────────────────────────────────────────────────────────
# Smooth predictions for top dynamic genes
# ─────────────────────────────────────────────────────────────────────────────

cat("\n[5/5] Smooth predictions ...\n")
if (!is.null(assoc)) {
    top_genes <- head(assoc$gene[assoc$padj < 0.05], 300)
    if (length(top_genes) == 0) {
        top_genes <- head(assoc$gene, 50)
        cat("  No sig genes at padj<0.05 — using top 50 by waldStat\n")
    }
    cat(sprintf("  Predicting smooth curves for %d genes ...\n", length(top_genes)))
    yhat <- tryCatch(
        predictSmooth(sce, gene = top_genes, nPoints = 100, tidy = FALSE),
        error = function(e) { cat(sprintf("  ERROR: %s\n", conditionMessage(e))); NULL }
    )
    if (!is.null(yhat)) {
        write.csv(yhat, file.path(out_dir, "tradeseq_smooth_predictions.csv"))
        cat(sprintf("  Saved smooth predictions\n"))
    }
}

cat("\n=== tradeSeq complete ===\n")
cat(sprintf("Outputs written to: %s\n", out_dir))
