#!/usr/bin/env Rscript
# =============================================================================
# proteomics_multiomics.R
#
# Multi-omics integration for CPTAC-style proteogenomics data. Three canonical
# analyses, each optional (runs only if the required input files are provided):
#
#   1. mRNA-PROTEIN correlation: per-gene correlation of RNA vs protein across
#      samples. Identifies post-transcriptionally regulated genes.
#   2. KINASE ACTIVITY from phosphoproteome: infers per-sample kinase activity
#      by scoring phosphosite changes against kinase-substrate sets, then tests
#      tumor vs normal. (Uses a simple substrate-set enrichment approach.)
#   3. CNV-PROTEIN cis effects: per-gene correlation of copy-number vs its OWN
#      protein abundance. Shows how genomic amplification propagates to protein.
#
# All inputs are LinkedOmics .cct (genes/sites x samples). Tumor and normal are
# provided as separate files (as downloaded) and merged with _Tumor/_Normal tags.
#
# Config keys (in addition to the standard ones):
#   rna_tumor / rna_normal        : RNAseq .cct files (log2 FPKM)
#   prot_tumor / prot_normal      : proteome .cct files (already have these)
#   phospho_tumor / phospho_normal: phosphoproteome site-level .cct files
#   cnv_tumor                     : SCNV gene-level .cct (tumor; log ratios)
#   kinase_gmt                    : optional GMT of kinase->substrate sites;
#                                   if absent, kinase step is skipped
#
# Usage:
#   Rscript proteomics_multiomics.R config.yaml
#
# Author: <your name>   |   License: MIT
# =============================================================================

suppressPackageStartupMessages({
  library(yaml); library(dplyr); library(ggplot2)
  library(ComplexHeatmap); library(circlize); library(matrixStats)
})
have_fgsea <- requireNamespace("fgsea", quietly = TRUE)
if (have_fgsea) suppressPackageStartupMessages(library(fgsea))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("Usage: Rscript proteomics_multiomics.R config.yaml")
cfg <- yaml::read_yaml(args[1])
if (is.null(cfg$outdir)) cfg$outdir <- "results/proteomics"
if (is.null(cfg$seed)) cfg$seed <- 42
set.seed(cfg$seed)
outdir <- file.path(cfg$outdir, "multiomics")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
log_msg <- function(...) cat(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), ..., "\n")

# ---- helpers ----
read_cct <- function(path) if (!is.null(path) && file.exists(path))
  as.matrix(read.delim(path, row.names = 1, check.names = FALSE)) else NULL

# Merge tumor + normal matrices with suffixes; return list(mat, type)
merge_tn <- function(tfile, nfile) {
  t <- read_cct(tfile); n <- read_cct(nfile)
  if (is.null(t) || is.null(n)) return(NULL)
  g <- intersect(rownames(t), rownames(n)); t <- t[g, , drop = FALSE]; n <- n[g, , drop = FALSE]
  colnames(t) <- paste0(colnames(t), "_Tumor"); colnames(n) <- paste0(colnames(n), "_Normal")
  list(mat = cbind(t, n), type = c(rep("Tumor", ncol(t)), rep("Normal", ncol(n))))
}

# =============================================================================
# 1. mRNA-PROTEIN CORRELATION
# =============================================================================
mrna_protein <- function(cfg) {
  rna  <- merge_tn(cfg$rna_tumor,  cfg$rna_normal)
  prot <- merge_tn(cfg$prot_tumor, cfg$prot_normal)
  if (is.null(rna) || is.null(prot)) { log_msg("RNA or protein files missing; skipping mRNA-protein"); return(invisible()) }
  log_msg("mRNA-protein correlation")

  # Shared genes and shared samples
  genes   <- intersect(rownames(rna$mat), rownames(prot$mat))
  samples <- intersect(colnames(rna$mat), colnames(prot$mat))
  R <- rna$mat[genes, samples, drop = FALSE]
  P <- prot$mat[genes, samples, drop = FALSE]
  log_msg(sprintf("  %d shared genes, %d shared samples", length(genes), length(samples)))

  # Per-gene Spearman correlation across samples
  cors <- sapply(seq_along(genes), function(i) {
    x <- R[i, ]; y <- P[i, ]; ok <- is.finite(x) & is.finite(y)
    if (sum(ok) < 10) return(NA_real_)
    suppressWarnings(cor(x[ok], y[ok], method = "spearman"))
  })
  df <- data.frame(gene = genes, rna_protein_cor = cors)
  df <- df[!is.na(df$rna_protein_cor), ]
  write.csv(df[order(-df$rna_protein_cor), ], file.path(outdir, "mrna_protein_correlation.csv"), row.names = FALSE)

  # Distribution plot (CPTAC-style: most genes positive, long tail)
  p <- ggplot(df, aes(rna_protein_cor)) +
    geom_histogram(bins = 60, fill = "steelblue", color = "white") +
    geom_vline(xintercept = median(df$rna_protein_cor), linetype = "dashed", color = "firebrick") +
    labs(title = "mRNA-protein correlation per gene",
         subtitle = sprintf("median Spearman = %.2f", median(df$rna_protein_cor)),
         x = "Spearman correlation (RNA vs protein)", y = "genes") + theme_bw()
  ggsave(file.path(outdir, "mrna_protein_correlation_hist.pdf"), p, width = 7, height = 5)
  log_msg(sprintf("  median mRNA-protein correlation = %.2f", median(df$rna_protein_cor)))
}

# =============================================================================
# 2. KINASE ACTIVITY from phosphoproteome
# =============================================================================
# Rank phosphosites by tumor-vs-normal change, then enrich kinase substrate sets
# (from a GMT) to estimate which kinases are more/less active in tumor.
kinase_activity <- function(cfg) {
  ph <- merge_tn(cfg$phospho_tumor, cfg$phospho_normal)
  if (is.null(ph)) { log_msg("Phospho files missing; skipping kinase activity"); return(invisible()) }
  if (is.null(cfg$kinase_gmt) || !file.exists(cfg$kinase_gmt)) {
    log_msg("No kinase_gmt provided; writing phosphosite tumor-vs-normal stats only")
  }
  log_msg("Kinase activity from phosphoproteome")

  type <- factor(ph$type, levels = c("Normal","Tumor"))
  M <- ph$mat
  # Per-site tumor-vs-normal t statistic (simple, robust for a ranking)
  tstat <- apply(M, 1, function(v) {
    tv <- v[type == "Tumor"]; nv <- v[type == "Normal"]
    tv <- tv[is.finite(tv)]; nv <- nv[is.finite(nv)]
    if (length(tv) < 3 || length(nv) < 3) return(NA_real_)
    tryCatch(t.test(tv, nv)$statistic, error = function(e) NA_real_)
  })
  site_df <- data.frame(site = rownames(M), t_tumor_vs_normal = tstat)
  site_df <- site_df[is.finite(site_df$t_tumor_vs_normal), ]
  write.csv(site_df[order(-site_df$t_tumor_vs_normal), ],
            file.path(outdir, "phosphosite_tumor_vs_normal.csv"), row.names = FALSE)

  # If a kinase-substrate GMT is given, enrich to get kinase activity scores
  if (!is.null(cfg$kinase_gmt) && file.exists(cfg$kinase_gmt) && have_fgsea) {
    read_gmt <- function(path) {
      lines <- readLines(path, warn = FALSE); s <- list()
      for (ln in lines) { p <- strsplit(ln, "\t")[[1]]; if (length(p) < 3) next
        g <- p[-c(1,2)]; s[[p[1]]] <- g[nzchar(g)] }
      s
    }
    ksets <- read_gmt(cfg$kinase_gmt)
    ranks <- setNames(site_df$t_tumor_vs_normal, site_df$site)
    res <- tryCatch(fgsea(ksets, ranks, minSize = 3, maxSize = 500), error = function(e) NULL)
    if (!is.null(res) && nrow(res) > 0) {
      out <- as.data.frame(res); out$leadingEdge <- sapply(out$leadingEdge, paste, collapse = ";")
      out <- out[order(-out$NES), ]
      write.csv(out, file.path(outdir, "kinase_activity.csv"), row.names = FALSE)
      top <- rbind(head(out, 15), tail(out, 15))
      top$pathway <- factor(top$pathway, levels = top$pathway[order(top$NES)])
      p <- ggplot(top, aes(NES, pathway, fill = padj < 0.1)) + geom_col() +
        scale_fill_manual(values = c("grey70","firebrick"), name = "padj<0.1") +
        labs(title = "Kinase activity (tumor vs normal)", x = "NES", y = NULL) +
        theme_bw() + theme(axis.text.y = element_text(size = 6))
      ggsave(file.path(outdir, "kinase_activity_barplot.pdf"), p, width = 8, height = 7)
      log_msg(sprintf("  kinase activity: %d kinases scored", nrow(out)))
    }
  } else {
    log_msg("  (install fgsea + provide kinase_gmt to get kinase activity scores)")
  }
}

# =============================================================================
# 3. CNV-PROTEIN cis effects
# =============================================================================
# Per-gene correlation of copy-number (tumor) vs its own protein abundance.
cnv_protein_cis <- function(cfg) {
  cnv  <- read_cct(cfg$cnv_tumor)
  prot <- read_cct(cfg$prot_tumor)
  if (is.null(cnv) || is.null(prot)) { log_msg("CNV or tumor protein missing; skipping CNV-protein"); return(invisible()) }
  log_msg("CNV-protein cis effects (tumor)")

  genes   <- intersect(rownames(cnv), rownames(prot))
  samples <- intersect(colnames(cnv), colnames(prot))
  C <- cnv[genes, samples, drop = FALSE]
  P <- prot[genes, samples, drop = FALSE]
  log_msg(sprintf("  %d shared genes, %d shared tumor samples", length(genes), length(samples)))

  cors <- sapply(seq_along(genes), function(i) {
    x <- C[i, ]; y <- P[i, ]; ok <- is.finite(x) & is.finite(y)
    if (sum(ok) < 10) return(NA_real_)
    suppressWarnings(cor(x[ok], y[ok], method = "spearman"))
  })
  df <- data.frame(gene = genes, cnv_protein_cis_cor = cors)
  df <- df[!is.na(df$cnv_protein_cis_cor), ]
  write.csv(df[order(-df$cnv_protein_cis_cor), ], file.path(outdir, "cnv_protein_cis.csv"), row.names = FALSE)

  p <- ggplot(df, aes(cnv_protein_cis_cor)) +
    geom_histogram(bins = 60, fill = "darkorange", color = "white") +
    geom_vline(xintercept = 0, linetype = "dashed") +
    geom_vline(xintercept = median(df$cnv_protein_cis_cor), linetype = "dashed", color = "firebrick") +
    labs(title = "CNV-protein cis correlation per gene",
         subtitle = sprintf("median = %.2f; %d%% of genes positive",
                            median(df$cnv_protein_cis_cor),
                            round(100*mean(df$cnv_protein_cis_cor > 0))),
         x = "Spearman (copy number vs own protein)", y = "genes") + theme_bw()
  ggsave(file.path(outdir, "cnv_protein_cis_hist.pdf"), p, width = 7, height = 5)
  log_msg(sprintf("  median CNV-protein cis = %.2f", median(df$cnv_protein_cis_cor)))
}

# ---- MAIN ----
main <- function(cfg) {
  mrna_protein(cfg)
  kinase_activity(cfg)
  cnv_protein_cis(cfg)
  writeLines(capture.output(sessionInfo()), file.path(outdir, "sessionInfo_multiomics.txt"))
  log_msg("Multi-omics analysis complete. Outputs in:", outdir)
}

if (sys.nframe() == 0) main(cfg)
