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
#   Rscript scripts/method_comparison/run_tradeseq.R \
#     --counts  results/pseudotime_genes/wolbachia_infection/tradeseq_inputs/counts_genesXcells.csv  \
#     --pt      results/pseudotime_genes/wolbachia_infection/tradeseq_inputs/pseudotime.csv    \
#     --outdir  results/pseudotime_genes/wolbachia_infection/ \
#     --nknots  6 \
#     --nworkers 16 \
#     --genes results/pseudotime_genes/wolbachia_infection/tradeseq_inputs/custom_genes.csv
#
# Outputs (written to --outdir):
#   tradeseq_sce.rds                    : fitted SCE (re-use without re-fitting)
#   tradeseq_association.csv            : associationTest results
#   tradeseq_startvsend.csv             : startVsEndTest results
#   tradeseq_smooth_predictions.csv     : predictSmooth for top 300 sig genes
#   tradeseq_association_volcano.pdf    : waldStat vs -log10(padj)
#   tradeseq_smooth_heatmap.pdf         : top 100 genes heatmap ordered by peak pseudotime
#   tradeseq_heatmap_genes_ordered.csv  : heatmap genes with peak stage + stats
#   tradeseq_smooth_curves_top9.pdf     : individual GAM curves for top 9 genes
#   tradeseq_smooth_curves_manual.pdf   : GAM curves for --genes list (if provided)
#   tradeseq_startvsend_barplot.pdf     : top 30 start-vs-end genes by direction

suppressPackageStartupMessages({
    library(tradeSeq)
    library(SingleCellExperiment)
    library(BiocParallel)
    library(Matrix)
    library(ggplot2)
    library(patchwork)
})

has_ggrepel <- requireNamespace("ggrepel", quietly = TRUE)

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

counts_file  <- parse_arg("--counts")
pt_file      <- parse_arg("--pt")
out_dir      <- parse_arg("--outdir")
n_knots      <- as.integer(parse_arg("--nknots",   "6"))
n_workers    <- as.integer(parse_arg("--nworkers", "8"))
genes_file   <- parse_arg("--genes", default = NULL)

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

cat("=== tradeSeq standalone ===\n")
cat(sprintf("  counts   : %s\n", counts_file))
cat(sprintf("  pt       : %s\n", pt_file))
cat(sprintf("  outdir   : %s\n", out_dir))
cat(sprintf("  nknots   : %d\n", n_knots))
cat(sprintf("  nworkers : %d\n", n_workers))
cat(sprintf("  ggrepel  : %s\n", if (has_ggrepel) "yes" else "no (gene labels disabled)"))
cat(sprintf("  genes    : %s\n", if (!is.null(genes_file)) genes_file else "none"))

# ─────────────────────────────────────────────────────────────────────────────
# Plotting functions
# ─────────────────────────────────────────────────────────────────────────────

plot_volcano <- function(assoc, out_dir) {
    cat("  Volcano plot ...\n")
    assoc$neg_log10p <- -log10(assoc$padj + 1e-300)

    p <- ggplot(assoc, aes(x = waldStat, y = neg_log10p, colour = sig)) +
        geom_point(size = 0.8, alpha = 0.6) +
        scale_colour_manual(values = c("FALSE" = "#bdc3c7", "TRUE" = "#c0392b"),
                            labels = c("ns", "padj < 0.05")) +
        geom_hline(yintercept = -log10(0.05), linetype = "dashed",
                   colour = "grey40", linewidth = 0.5) +
        labs(x      = "Wald statistic",
             y      = expression(-log[10](p[adj])),
             title  = "Association test: genes dynamic along pseudotime",
             colour = NULL) +
        theme_bw(base_size = 11) +
        theme(legend.position = "top")

    if (has_ggrepel) {
        top20 <- head(assoc[assoc$sig, ], 20)
        if (nrow(top20) > 0)
            p <- p + ggrepel::geom_text_repel(
                data = top20, aes(label = gene),
                size = 2.5, max.overlaps = 20, colour = "black")
    }

    out <- file.path(out_dir, "tradeseq_association_volcano.pdf")
    ggsave(out, p, width = 8, height = 6)
    cat(sprintf("    Saved: %s\n", out))
}


plot_smooth_heatmap <- function(assoc, yhat, out_dir) {
    cat("  Smooth heatmap ...\n")
    top_genes <- head(assoc$gene[assoc$sig], 100)
    top_genes <- top_genes[top_genes %in% rownames(yhat)]

    if (length(top_genes) < 5) {
        cat(sprintf("  Only %d sig genes in yhat — using top 50 by waldStat\n",
                    length(top_genes)))
        top_genes <- head(assoc$gene, 50)
        top_genes <- top_genes[top_genes %in% rownames(yhat)]
    }

    if (length(top_genes) < 5) {
        cat("  Too few genes for heatmap — skipping\n")
        return(invisible(NULL))
    }

    mat   <- as.matrix(yhat[top_genes, , drop = FALSE])
    mat_z <- t(scale(t(mat)))
    mat_z[is.nan(mat_z)] <- 0

    # Order genes by peak pseudotime (column of max smoothed expression)
    peak_idx   <- apply(mat_z, 1, which.max)
    gene_order <- order(peak_idx)
    mat_z      <- mat_z[gene_order, , drop = FALSE]
    peak_idx   <- peak_idx[gene_order]

    df_long <- data.frame(
        gene     = rep(rownames(mat_z), times = ncol(mat_z)),
        pt_point = rep(seq_len(ncol(mat_z)), each = nrow(mat_z)),
        expr     = as.vector(mat_z)
    )
    df_long$gene <- factor(df_long$gene, levels = rev(rownames(mat_z)))

    p <- ggplot(df_long, aes(x = pt_point, y = gene, fill = expr)) +
        geom_tile() +
        scale_fill_viridis_c(option = "viridis", name = "Scaled\nexpr") +
        scale_x_continuous(expand = c(0, 0), name = "Pseudotime ->") +
        theme_minimal(base_size = 8) +
        theme(axis.text.y  = element_text(
                  size = ifelse(length(top_genes) > 60, 4, 6)),
              axis.text.x  = element_blank(),
              axis.ticks.x = element_blank(),
              panel.grid   = element_blank()) +
        labs(title = sprintf("Top %d dynamic genes ordered by peak pseudotime",
                             length(top_genes)),
             y = NULL)

    h   <- max(5, length(top_genes) * 0.12 + 1)
    out <- file.path(out_dir, "tradeseq_smooth_heatmap.pdf")
    ggsave(out, p, width = 10, height = h, limitsize = FALSE)
    cat(sprintf("    Saved: %s\n", out))

    # Export genes in heatmap order (top to bottom = early to late peak)
    n_pts     <- ncol(mat_z)
    peak_frac <- peak_idx / n_pts
    stage     <- cut(peak_frac,
                     breaks = c(0, 1/3, 2/3, 1),
                     labels = c("early", "middle", "late"),
                     include.lowest = TRUE)

    heatmap_df <- assoc[match(rownames(mat_z), assoc$gene),
                        c("gene", "waldStat", "pvalue", "padj", "sig")]
    heatmap_df$heatmap_rank  <- seq_len(nrow(heatmap_df))
    heatmap_df$peak_pt_index <- peak_idx
    heatmap_df$peak_pt_frac  <- round(peak_frac, 4)
    heatmap_df$stage         <- as.character(stage)
    heatmap_df <- heatmap_df[, c("heatmap_rank", "gene", "stage",
                                  "peak_pt_index", "peak_pt_frac",
                                  "waldStat", "pvalue", "padj", "sig")]
    out_csv <- file.path(out_dir, "tradeseq_heatmap_genes_ordered.csv")
    write.csv(heatmap_df, out_csv, row.names = FALSE)
    cat(sprintf("    Saved: %s\n", out_csv))
}


plot_smooth_curves <- function(assoc, sce, out_dir, manual_genes = character(0)) {
    cat("  Smooth curves for top 9 genes ...\n")
    top9 <- head(assoc$gene[assoc$sig], 9)
    if (length(top9) == 0) {
        top9 <- head(assoc$gene, 9)
        cat("  No sig genes — using top 9 by waldStat\n")
    }
    if (length(top9) == 0) {
        cat("  No genes available — skipping smooth curves\n")
        return(invisible(NULL))
    }

    sce_counts <- as.matrix(assay(sce, "counts"))

    make_curve_grid <- function(genes, label) {
        plots <- list()
        for (g in genes) {
            if (!g %in% rownames(sce)) {
                cat(sprintf("    WARNING: %s not in SCE — skipping\n", g))
                next
            }
            tryCatch({
                p <- plotSmoothers(sce, counts = sce_counts, gene = g) +
                    labs(title = g, x = "Pseudotime",
                         y = "log-normalised expression") +
                    theme_bw(base_size = 9) +
                    theme(legend.position = "none")
                plots[[g]] <- p
            }, error = function(e) {
                cat(sprintf("    WARNING: %s failed: %s\n", g, conditionMessage(e)))
            })
        }
        if (length(plots) == 0) {
            cat(sprintf("  All %s genes failed — skipping\n", label))
            return(invisible(NULL))
        }
        combined <- wrap_plots(plots, ncol = 3)
        out <- file.path(out_dir, sprintf("tradeseq_smooth_curves_%s.pdf", label))
        ggsave(out, combined,
               width  = min(3, length(plots)) * 5,
               height = ceiling(length(plots) / 3) * 5)
        cat(sprintf("    Saved: %s\n", out))
    }

    make_curve_grid(top9, "top9")

    if (length(manual_genes) > 0) {
        cat(sprintf("  Smooth curves for %d manual genes ...\n", length(manual_genes)))
        make_curve_grid(manual_genes, "manual")
    }
}


plot_startvsend <- function(sve, out_dir) {
    cat("  Start vs end barplot ...\n")
    top_sve <- head(sve[sve$sig, ], 30)

    if (nrow(top_sve) == 0) {
        cat("  No significant start-vs-end genes — skipping\n")
        return(invisible(NULL))
    }

    lfc_col <- grep("^logFC", colnames(top_sve), value = TRUE)[1]
    top_sve$direction <- if (!is.na(lfc_col)) {
        ifelse(top_sve[[lfc_col]] > 0, "Up at end", "Down at end")
    } else {
        "Dynamic"
    }
    top_sve$gene <- factor(top_sve$gene, levels = rev(top_sve$gene))

    p <- ggplot(top_sve, aes(x = waldStat, y = gene, fill = direction)) +
        geom_col(alpha = 0.8) +
        scale_fill_manual(values = c("Up at end"   = "#c0392b",
                                     "Down at end" = "#2980b9",
                                     "Dynamic"     = "#888888")) +
        labs(x     = "Wald statistic",
             y     = NULL,
             title = "Start vs End test - top dynamic genes",
             fill  = NULL) +
        theme_bw(base_size = 10) +
        theme(legend.position = "top")

    out <- file.path(out_dir, "tradeseq_startvsend_barplot.pdf")
    ggsave(out, p, width = 7, height = max(4, nrow(top_sve) * 0.25 + 2))
    cat(sprintf("    Saved: %s\n", out))
}


plot_custom_gene_curves <- function(genes_df, sce, assoc, out_dir) {
    cat("  Custom gene smooth curves ...\n")

    genes_out  <- file.path(out_dir, "custom_gene_curves")
    dir.create(genes_out, showWarnings = FALSE, recursive = TRUE)

    sce_counts <- as.matrix(assay(sce, "counts"))
    available  <- rownames(sce_counts)
    plots      <- list()

    for (i in seq_len(nrow(genes_df))) {
        gene_name  <- as.character(genes_df$Gene[i])
        flybase_id <- as.character(genes_df$FlyBaseId[i])

        idx <- match(flybase_id, available)
        if (is.na(idx)) idx <- match(tolower(gene_name), tolower(available))

        if (is.na(idx)) {
            cat(sprintf("    SKIP: '%s' (%s) not found in SCE\n", gene_name, flybase_id))
            next
        }

        matched_gene <- available[idx]

        subtitle <- NULL
        if (!is.null(assoc) && matched_gene %in% assoc$gene) {
            ar       <- assoc[assoc$gene == matched_gene, ]
            padj_val <- ar$padj
            padj_str <- if (is.na(padj_val)) "NA" else if (padj_val < 0.001) sprintf("%.2e", padj_val) else sprintf("%.4f", padj_val)
            subtitle <- sprintf("waldStat = %.2f  |  padj = %s  |  sig = %s",
                                ar$waldStat, padj_str, ifelse(isTRUE(ar$sig), "YES", "NO"))
        }

        tryCatch({
            p <- plotSmoothers(sce, counts = sce_counts, gene = matched_gene) +
                labs(title    = sprintf("%s (%s)", gene_name, flybase_id),
                     subtitle = subtitle,
                     x = "Pseudotime", y = "log-normalised expression") +
                theme_bw(base_size = 9) +
                theme(legend.position = "none")
            plots[[gene_name]] <- p
        }, error = function(e) {
            cat(sprintf("    WARNING: %s failed: %s\n", gene_name, conditionMessage(e)))
        })
    }

    if (length(plots) == 0) {
        cat("  All custom genes failed or not found — skipping\n")
        return(invisible(NULL))
    }

    combined <- wrap_plots(plots, ncol = 3)
    out <- file.path(out_dir, "tradeseq_smooth_curves_custom.pdf")
    ggsave(out, combined,
           width  = min(3, length(plots)) * 5,
           height = ceiling(length(plots) / 3) * 5)
    cat(sprintf("    Saved: %s\n", out))
}

# ─────────────────────────────────────────────────────────────────────────────
# Load data
# ─────────────────────────────────────────────────────────────────────────────

cat("\n[1/6] Loading counts ...\n")
counts <- as.matrix(read.csv(counts_file, row.names = 1, check.names = FALSE))
cat(sprintf("  Shape: %d genes x %d cells\n", nrow(counts), ncol(counts)))
cat(sprintf("  Count range: %.1f - %.1f\n", min(counts), max(counts)))
counts <- round(counts)
storage.mode(counts) <- "integer"

cat("\n[2/6] Loading pseudotime ...\n")
pt_df  <- read.csv(pt_file, row.names = 1)
pt_vec <- setNames(pt_df[[1]], rownames(pt_df))
cat(sprintf("  Cells: %d  Range: %.3f - %.3f\n",
            length(pt_vec), min(pt_vec), max(pt_vec)))

# ─────────────────────────────────────────────────────────────────────────────
# Align and filter
# ─────────────────────────────────────────────────────────────────────────────

shared <- intersect(colnames(counts), names(pt_vec))
cat(sprintf("  Shared cells: %d\n", length(shared)))
if (length(shared) == 0) stop("No shared cells between counts and pseudotime!")

counts <- counts[, shared, drop = FALSE]
pt_vec <- pt_vec[shared]

gene_sums <- rowSums(counts)
n_zero    <- sum(gene_sums == 0)
if (n_zero > 0) {
    cat(sprintf("  Removing %d all-zero genes\n", n_zero))
    counts <- counts[gene_sums > 0, , drop = FALSE]
}

n_cells_expr <- rowSums(counts > 0)
n_sparse     <- sum(n_cells_expr < 5)
if (n_sparse > 0) {
    cat(sprintf("  Removing %d genes expressed in <5 cells\n", n_sparse))
    counts <- counts[n_cells_expr >= 5, , drop = FALSE]
}
cat(sprintf("  Final: %d genes x %d cells\n", nrow(counts), ncol(counts)))

pt_mat <- matrix(pt_vec, ncol = 1, dimnames = list(shared, "pseudotime"))
wt_mat <- matrix(1, nrow = length(shared), ncol = 1, dimnames = list(shared, "w1"))

# ─────────────────────────────────────────────────────────────────────────────
# Fit GAMs
# ─────────────────────────────────────────────────────────────────────────────

cat(sprintf("\n[3/6] Fitting GAMs (%d genes, %d knots, %d workers) ...\n",
            nrow(counts), n_knots, n_workers))
cat("  This is the slow step — expect 1-4 hrs depending on gene count\n")

set.seed(42)
BPPARAM <- MulticoreParam(workers = n_workers, progressbar = TRUE)

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

saveRDS(sce, file.path(out_dir, "tradeseq_sce.rds"))
cat(sprintf("  Saved: %s\n", file.path(out_dir, "tradeseq_sce.rds")))

# ─────────────────────────────────────────────────────────────────────────────
# Statistical tests
# ─────────────────────────────────────────────────────────────────────────────

cat("\n[4/6] Running tests ...\n")

cat("  associationTest ...\n")
assoc <- tryCatch(
    as.data.frame(associationTest(sce, lineages = FALSE)),
    error = function(e) { cat(sprintf("  ERROR: %s\n", conditionMessage(e))); NULL }
)
if (!is.null(assoc)) {
    assoc$gene <- rownames(assoc)
    assoc$padj <- p.adjust(assoc$pvalue, method = "BH")
    assoc$sig  <- assoc$padj < 0.05
    assoc      <- assoc[order(assoc$waldStat, decreasing = TRUE), ]
    write.csv(assoc, file.path(out_dir, "tradeseq_association.csv"))
    cat(sprintf("  Significant genes (padj<0.05): %d / %d\n",
                sum(assoc$sig, na.rm = TRUE), nrow(assoc)))
}

cat("  startVsEndTest ...\n")
sve <- tryCatch(
    as.data.frame(startVsEndTest(sce)),
    error = function(e) { cat(sprintf("  ERROR: %s\n", conditionMessage(e))); NULL }
)
if (!is.null(sve)) {
    sve$gene <- rownames(sve)
    sve$padj <- p.adjust(sve$pvalue, method = "BH")
    sve$sig  <- sve$padj < 0.05
    sve      <- sve[order(sve$waldStat, decreasing = TRUE), ]
    write.csv(sve, file.path(out_dir, "tradeseq_startvsend.csv"))
    cat(sprintf("  Significant genes (padj<0.05): %d / %d\n",
                sum(sve$sig, na.rm = TRUE), nrow(sve)))
}

# ─────────────────────────────────────────────────────────────────────────────
# Smooth predictions
# ─────────────────────────────────────────────────────────────────────────────

cat("\n[5/6] Smooth predictions ...\n")
yhat <- NULL
if (!is.null(assoc)) {
    top_genes <- head(assoc$gene[assoc$sig], 300)
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
        cat("  Saved smooth predictions\n")
    }
}

# ─────────────────────────────────────────────────────────────────────────────
# Plots
# ─────────────────────────────────────────────────────────────────────────────

cat("\n[6/6] Generating plots ...\n")

if (!is.null(assoc))                    plot_volcano(assoc, out_dir)
if (!is.null(assoc) && !is.null(yhat)) plot_smooth_heatmap(assoc, yhat, out_dir)
if (!is.null(assoc))                    plot_smooth_curves(assoc, sce, out_dir)
if (!is.null(sve))                      plot_startvsend(sve, out_dir)
if (!is.null(genes_file) && !is.null(assoc)) {
    cat("  Loading custom gene list ...\n")
    genes_df <- tryCatch(
        read.csv(genes_file, stringsAsFactors = FALSE, strip.white = TRUE),
        error = function(e) {
            cat(sprintf("  WARNING: Could not read genes file: %s\n",
                        conditionMessage(e)))
            NULL
        }
    )
    if (!is.null(genes_df) && all(c("Gene", "FlyBaseId") %in% colnames(genes_df))) {
        plot_custom_gene_curves(genes_df, sce, assoc, out_dir)
    } else {
        cat("  WARNING: genes file must have columns 'Gene' and 'FlyBaseId'\n")
    }
}

cat("\n=== tradeSeq complete ===\n")
cat(sprintf("Outputs written to: %s\n", out_dir))