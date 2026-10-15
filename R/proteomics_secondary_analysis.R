#!/usr/bin/env Rscript
# =============================================================================
# proteomics_secondary_analysis.R
#
# Standalone SECONDARY / exploratory analyses for a proteomics abundance matrix,
# complementing the main differential-abundance pipeline. Modules:
#
#   1. Sample correlation + hierarchical clustering (QC / global structure)
#   2. Unsupervised sample clustering (k-means) -> candidate subtypes
#   3. ssGSEA per-sample pathway scores (GSVA) + pathway x sample heatmap
#   4. Co-expression modules (WGCNA-lite: correlation -> hierarchical modules)
#
# Reuses the SAME config as the main pipeline (abundance, metadata, sample_col,
# condition, organism, etc.). Reads the already-prepared abundance matrix.
#
# Usage:
#   Rscript proteomics_secondary_analysis.R config.yaml
#
# Author: <your name>   |   License: MIT
# =============================================================================

suppressPackageStartupMessages({
  library(yaml); library(dplyr); library(ggplot2)
  library(ComplexHeatmap); library(circlize); library(matrixStats)
})
have_gsva   <- requireNamespace("GSVA", quietly = TRUE)
have_msig   <- requireNamespace("msigdbr", quietly = TRUE)
have_impute <- requireNamespace("impute", quietly = TRUE)
if (have_msig) suppressPackageStartupMessages(library(msigdbr))

# ---- CONFIG ----
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("Usage: Rscript proteomics_secondary_analysis.R config.yaml")
cfg <- yaml::read_yaml(args[1])
defs <- list(organism = "human", outdir = "results/proteomics", min_valid = 0.5,
             n_clusters = 3, n_top_var = 2000, n_modules = 8, seed = 42)
for (k in names(defs)) if (is.null(cfg[[k]])) cfg[[k]] <- defs[[k]]
set.seed(cfg$seed)
outdir <- file.path(cfg$outdir, "secondary")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
log_msg <- function(...) cat(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), ..., "\n")
msig_species <- if (cfg$organism == "mouse") "Mus musculus" else "Homo sapiens"

# ---- LOAD + light prep (filter, impute to complete matrix for these analyses) ----
load_matrix <- function(cfg) {
  log_msg("Loading abundance + metadata")
  mat <- as.matrix(read.delim(cfg$abundance, row.names = 1, check.names = FALSE))
  sep <- if (grepl("\\.tsv$", cfg$metadata, ignore.case = TRUE)) "\t" else ","
  meta <- read.delim(cfg$metadata, sep = sep, stringsAsFactors = FALSE, check.names = FALSE)
  rownames(meta) <- meta[[cfg$sample_col]]
  common <- intersect(colnames(mat), rownames(meta))
  mat <- mat[, common, drop = FALSE]; meta <- meta[common, , drop = FALSE]
  if (!isTRUE(cfg$is_log)) { mat[mat == 0] <- NA; mat <- log2(mat) }
  # keep well-quantified proteins
  mat <- mat[rowMeans(!is.na(mat)) >= cfg$min_valid, , drop = FALSE]
  # impute to complete (needed for correlation/clustering/ssGSEA)
  if (anyNA(mat)) {
    if (have_impute) { log_msg("  KNN imputation"); mat <- impute::impute.knn(mat)$data }
    else { log_msg("  dropping NA rows (impute not available)"); mat <- mat[rowSums(is.na(mat)) == 0, , drop = FALSE] }
  }
  list(mat = mat, meta = meta)
}

# ---- 1. SAMPLE CORRELATION + HIERARCHICAL CLUSTERING ----
sample_correlation <- function(mat, meta, cfg) {
  log_msg("Sample correlation + hierarchical clustering")
  cm <- cor(mat, method = "pearson")
  ann <- HeatmapAnnotation(group = meta[[cfg$condition]])
  pdf(file.path(outdir, "sample_correlation_heatmap.pdf"), width = 9, height = 8)
  draw(Heatmap(cm, name = "Pearson", top_annotation = ann,
               show_row_names = ncol(cm) <= 40, show_column_names = ncol(cm) <= 40,
               col = colorRamp2(c(min(cm), median(cm), 1), c("navy","white","firebrick")),
               column_title = "Sample-sample correlation (hierarchical)"))
  dev.off()
  write.csv(cm, file.path(outdir, "sample_correlation_matrix.csv"))
}

# ---- 2. UNSUPERVISED SAMPLE CLUSTERING (k-means on top-variable proteins) ----
sample_clustering <- function(mat, meta, cfg) {
  log_msg(sprintf("Unsupervised clustering (k = %d) on top %d variable proteins", cfg$n_clusters, cfg$n_top_var))
  vars <- rowVars(mat)
  top <- mat[order(vars, decreasing = TRUE)[1:min(cfg$n_top_var, nrow(mat))], , drop = FALSE]
  km <- kmeans(t(top), centers = cfg$n_clusters, nstart = 25)
  cl <- factor(km$cluster)

  # PCA colored by both cluster and known group, to see if clusters track biology
  p <- prcomp(t(top), scale. = TRUE)
  pv <- round(100 * (p$sdev^2 / sum(p$sdev^2))[1:2])
  df <- data.frame(PC1 = p$x[,1], PC2 = p$x[,2], cluster = cl, group = meta[[cfg$condition]])
  g1 <- ggplot(df, aes(PC1, PC2, color = cluster)) + geom_point(size = 2, alpha = 0.8) +
    labs(title = "Unsupervised clusters", x = sprintf("PC1 (%d%%)", pv[1]), y = sprintf("PC2 (%d%%)", pv[2])) + theme_bw()
  g2 <- ggplot(df, aes(PC1, PC2, color = group)) + geom_point(size = 2, alpha = 0.8) +
    labs(title = "Known group", x = sprintf("PC1 (%d%%)", pv[1]), y = sprintf("PC2 (%d%%)", pv[2])) + theme_bw()
  ggsave(file.path(outdir, "clusters_pca.pdf"), g1, width = 6, height = 5)
  ggsave(file.path(outdir, "groups_pca.pdf"), g2, width = 6, height = 5)

  # Cross-tab clusters vs known group (do subtypes align with tumor/normal?)
  ct <- table(cluster = cl, group = meta[[cfg$condition]])
  write.csv(as.data.frame.matrix(ct), file.path(outdir, "cluster_vs_group.csv"))
  write.csv(data.frame(sample = colnames(mat), cluster = cl), file.path(outdir, "sample_clusters.csv"), row.names = FALSE)
  log_msg("  cluster vs group cross-tab written")
}

# ---- 3. ssGSEA PER-SAMPLE PATHWAY SCORES (GSVA) ----
run_ssgsea <- function(mat, meta, cfg) {
  if (!have_gsva || !have_msig) { log_msg("GSVA/msigdbr not installed; skipping ssGSEA"); return(invisible()) }
  log_msg("ssGSEA per-sample pathway scores (Hallmark)")
  h <- msigdbr(species = msig_species, category = "H")
  gs <- split(h$gene_symbol, h$gs_name)
  gs <- lapply(gs, function(g) intersect(g, rownames(mat)))
  gs <- gs[sapply(gs, length) >= 5]
  if (length(gs) == 0) { log_msg("  no gene sets overlap the matrix; skipping"); return(invisible()) }

  es <- tryCatch({
    par <- GSVA::ssgseaParam(mat, gs)          # new GSVA API
    GSVA::gsva(par)
  }, error = function(e) {
    tryCatch(GSVA::gsva(mat, gs, method = "ssgsea"), error = function(e2) NULL)  # old API fallback
  })
  if (is.null(es)) { log_msg("  ssGSEA failed; skipping"); return(invisible()) }

  write.csv(es, file.path(outdir, "ssGSEA_hallmark_scores.csv"))
  m <- t(scale(t(es)))
  ann <- HeatmapAnnotation(group = meta[[cfg$condition]])
  pdf(file.path(outdir, "ssGSEA_hallmark_heatmap.pdf"), width = 11, height = 10)
  draw(Heatmap(m, name = "z", top_annotation = ann, show_column_names = ncol(m) <= 50,
               row_names_gp = gpar(fontsize = 6),
               col = colorRamp2(c(-2,0,2), c("navy","white","firebrick")),
               column_title = "ssGSEA Hallmark pathway activity"))
  dev.off()
  log_msg(sprintf("  ssGSEA done (%d pathways x %d samples)", nrow(es), ncol(es)))
}

# ---- 4. CO-EXPRESSION MODULES (WGCNA-lite) ----
# Correlation-based module detection without the full WGCNA dependency:
# correlate top-variable proteins, cluster them hierarchically, cut into modules,
# then relate each module's eigenprotein (PC1) to the sample groups.
coexpression_modules <- function(mat, meta, cfg) {
  log_msg(sprintf("Co-expression modules (top %d proteins, %d modules)", cfg$n_top_var, cfg$n_modules))
  vars <- rowVars(mat)
  top <- mat[order(vars, decreasing = TRUE)[1:min(cfg$n_top_var, nrow(mat))], , drop = FALSE]

  # Protein-protein correlation -> distance -> hierarchical clustering
  pc <- cor(t(top))
  d <- as.dist(1 - pc)
  hc <- hclust(d, method = "average")
  modules <- cutree(hc, k = cfg$n_modules)
  write.csv(data.frame(protein = rownames(top), module = modules),
            file.path(outdir, "coexpression_modules.csv"), row.names = FALSE)

  # Module eigenprotein = PC1 of each module; correlate with group
  grp_num <- as.numeric(factor(meta[[cfg$condition]]))
  eig <- sapply(sort(unique(modules)), function(m) {
    sub <- top[modules == m, , drop = FALSE]
    if (nrow(sub) < 2) return(rep(NA, ncol(sub)))
    prcomp(t(sub), scale. = TRUE)$x[, 1]
  })
  colnames(eig) <- paste0("module", sort(unique(modules)))
  eig_cor <- apply(eig, 2, function(e) if (all(is.na(e))) NA else cor(e, grp_num))
  write.csv(data.frame(module = names(eig_cor), group_correlation = eig_cor),
            file.path(outdir, "module_group_correlation.csv"), row.names = FALSE)

  # Heatmap of module eigenproteins across samples
  m <- t(scale(eig))
  ann <- HeatmapAnnotation(group = meta[[cfg$condition]])
  pdf(file.path(outdir, "module_eigenprotein_heatmap.pdf"), width = 10, height = 5)
  draw(Heatmap(m, name = "z", top_annotation = ann, show_column_names = ncol(m) <= 50,
               col = colorRamp2(c(-2,0,2), c("navy","white","firebrick")),
               column_title = "Co-expression module eigenproteins"))
  dev.off()

  mod_sizes <- table(modules)
  log_msg(sprintf("  %d modules (sizes: %s)", length(mod_sizes), paste(mod_sizes, collapse = ", ")))
}

# ---- MAIN ----
main <- function(cfg) {
  d <- load_matrix(cfg)
  sample_correlation(d$mat, d$meta, cfg)
  sample_clustering(d$mat, d$meta, cfg)
  run_ssgsea(d$mat, d$meta, cfg)
  coexpression_modules(d$mat, d$meta, cfg)
  writeLines(capture.output(sessionInfo()), file.path(outdir, "sessionInfo_secondary.txt"))
  log_msg("Secondary analysis complete. Outputs in:", outdir)
}

if (sys.nframe() == 0) main(cfg)
