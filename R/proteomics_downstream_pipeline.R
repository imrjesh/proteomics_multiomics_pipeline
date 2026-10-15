#!/usr/bin/env Rscript
# =============================================================================
# proteomics_downstream_pipeline.R
#
# Industry-standard bulk PROTEOMICS downstream analysis from a protein-abundance
# matrix. Two input branches (set `platform` in the config):
#
#   platform: "TMT"  -> isobaric-labeling data (e.g. CPTAC LUAD).
#                       Input is typically already log2-ratio / normalized
#                       abundance (protein x sample). Handles batch/plex.
#   platform: "DIA"  -> label-free DIA data (e.g. DIA-NN report from PRIDE).
#                       Input is linear intensities; pipeline log2-transforms
#                       and normalizes.
#
# SCOPE:
#   Starts FROM a protein-group abundance matrix (the output of the spectral
#   search: FragPipe/DIA-NN/Spectronaut for DIA, or the CPTAC/PDC processed
#   matrix for TMT). Raw-file spectral search is an UPSTREAM step, documented
#   in the README, and not repeated here (mirrors starting RNA-seq from counts).
#
# Downstream (shared across both branches):
#   QC (missingness, intensity distributions, PCA, correlation) -> filtering ->
#   log2 + normalization -> missing-value imputation -> differential abundance
#   (limma) -> volcano/heatmap -> enrichment (ORA + GSEA: GO/Reactome/Hallmark).
#
# Usage:
#   Rscript proteomics_downstream_pipeline.R config.yaml
#
# Author: <your name>   |   License: MIT
# =============================================================================

suppressPackageStartupMessages({
  library(yaml); library(limma); library(dplyr); library(tibble); library(tidyr)
  library(ggplot2); library(ComplexHeatmap); library(circlize); library(matrixStats)
  library(fgsea); library(msigdbr)
})
have_cp      <- requireNamespace("clusterProfiler", quietly = TRUE)
have_impute  <- requireNamespace("impute", quietly = TRUE)         # KNN imputation
have_sva     <- requireNamespace("sva", quietly = TRUE)            # ComBat batch correction
suppressPackageStartupMessages({ if (have_cp) library(clusterProfiler) })

# -----------------------------------------------------------------------------
# 0. CONFIG
# -----------------------------------------------------------------------------
# config.yaml keys:
#   platform:    "TMT" or "DIA"
#   abundance:   protein x sample matrix (proteins as rownames/first column)
#   metadata:    sample sheet (.csv/.tsv); rows = samples
#   sample_col:  metadata column with sample IDs (match abundance columns)
#   gene_col:    (optional) column in abundance mapping protein -> gene SYMBOL;
#                if absent, rownames are assumed to be gene symbols
#   organism:    "human" or "mouse"
#   outdir:      output directory
#   condition:   metadata column defining groups (e.g. tumor/normal)
#   contrasts:   list of {name, numerator, denominator} on `condition`
#   batch_col:   (optional) metadata column for batch/plex (TMT) -> ComBat
#   is_log:      TRUE if abundance is already log-scale (typical for CPTAC TMT);
#                FALSE for linear DIA intensities (pipeline will log2 them)
#   min_valid:   min fraction of samples a protein must be quantified in (default 0.5)
#   impute:      "knn" | "mindet" | "none" (default "knn")
#   normalize:   "median" | "quantile" | "none" (default "median")
#   padj:        adjusted-p cutoff (default 0.05)
#   lfc:         |log2FC| cutoff (default 1)
#   seed:        random seed (default 42)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("Usage: Rscript proteomics_downstream_pipeline.R config.yaml")
cfg <- yaml::read_yaml(args[1])
defs <- list(platform = "TMT", organism = "human", outdir = "results/proteomics",
             is_log = TRUE, min_valid = 0.5, impute = "knn", normalize = "median",
             padj = 0.05, lfc = 1, seed = 42, gene_col = NULL, batch_col = NULL)
for (k in names(defs)) if (is.null(cfg[[k]])) cfg[[k]] <- defs[[k]]
set.seed(cfg$seed)
for (s in c("qc","de","enrichment")) dir.create(file.path(cfg$outdir, s), showWarnings = FALSE, recursive = TRUE)
log_msg <- function(...) cat(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), ..., "\n")

msig_species <- if (cfg$organism == "mouse") "Mus musculus" else "Homo sapiens"
org_pkg <- if (cfg$organism == "mouse") "org.Mm.eg.db" else "org.Hs.eg.db"
if (have_cp) suppressPackageStartupMessages(library(org_pkg, character.only = TRUE))

# -----------------------------------------------------------------------------
# 1. LOAD  (branch-aware)
# -----------------------------------------------------------------------------
load_data <- function(cfg) {
  log_msg("Loading abundance matrix + metadata (platform = ", cfg$platform, ")")
  mat <- as.matrix(read.delim(cfg$abundance, row.names = 1, check.names = FALSE))
  sep <- if (grepl("\\.tsv$", cfg$metadata, ignore.case = TRUE)) "\t" else ","
  meta <- read.delim(cfg$metadata, sep = sep, stringsAsFactors = FALSE, check.names = FALSE)
  rownames(meta) <- meta[[cfg$sample_col]]
  common <- intersect(colnames(mat), rownames(meta))
  if (length(common) < 4) stop("Fewer than 4 shared samples between abundance and metadata.")
  mat <- mat[, common, drop = FALSE]; meta <- meta[common, , drop = FALSE]

  # DIA branch: linear intensities -> log2. TMT branch: usually already log.
  if (cfg$platform == "DIA" && !isTRUE(cfg$is_log)) {
    log_msg("  DIA: log2-transforming linear intensities")
    mat[mat == 0] <- NA
    mat <- log2(mat)
  } else if (!isTRUE(cfg$is_log)) {
    mat[mat == 0] <- NA; mat <- log2(mat)
  }
  log_msg(sprintf("  %d proteins x %d samples", nrow(mat), ncol(mat)))
  list(mat = mat, meta = meta)
}

# -----------------------------------------------------------------------------
# 2. QC  (missingness, intensity distributions, PCA, correlation)
# -----------------------------------------------------------------------------
run_qc <- function(mat, meta, cfg) {
  log_msg("QC")
  # Missingness per sample
  miss <- colMeans(is.na(mat))
  pdf(file.path(cfg$outdir, "qc", "missingness_per_sample.pdf"), width = 8, height = 4)
  barplot(miss, las = 2, ylab = "fraction missing", main = "Missing values per sample",
          names.arg = rep("", length(miss)))
  dev.off()
  # Intensity boxplots (sampled columns if many)
  idx <- if (ncol(mat) > 60) sort(sample(ncol(mat), 60)) else seq_len(ncol(mat))
  pdf(file.path(cfg$outdir, "qc", "intensity_boxplots.pdf"), width = 10, height = 5)
  boxplot(mat[, idx], las = 2, outline = FALSE, names = rep("", length(idx)),
          main = "Log2 intensity distribution (subset)", ylab = "log2 abundance")
  dev.off()
}

pca_plot <- function(mat, meta, cfg, tag = "postnorm") {
  # PCA needs complete data; use proteins present in all samples
  complete <- mat[rowSums(is.na(mat)) == 0, , drop = FALSE]
  if (nrow(complete) < 10) return(invisible())
  p <- prcomp(t(complete), scale. = TRUE)
  pv <- round(100 * (p$sdev^2 / sum(p$sdev^2))[1:2])
  df <- data.frame(PC1 = p$x[,1], PC2 = p$x[,2], group = meta[[cfg$condition]])
  if (!is.null(cfg$batch_col)) df$batch <- factor(meta[[cfg$batch_col]])
  g <- ggplot(df, aes(PC1, PC2, color = group)) + geom_point(size = 2, alpha = 0.8) +
    labs(x = sprintf("PC1 (%d%%)", pv[1]), y = sprintf("PC2 (%d%%)", pv[2]),
         title = paste0("PCA (", tag, ")")) + theme_bw()
  ggsave(file.path(cfg$outdir, "qc", paste0("pca_", tag, ".pdf")), g, width = 6, height = 5)
}

# -----------------------------------------------------------------------------
# 3. FILTER  (valid-value threshold) + NORMALIZE + IMPUTE
# -----------------------------------------------------------------------------
filter_normalize_impute <- function(mat, meta, cfg) {
  # Keep proteins quantified in >= min_valid fraction of samples
  keep <- rowMeans(!is.na(mat)) >= cfg$min_valid
  log_msg(sprintf("Filter: kept %d / %d proteins (>= %.0f%% valid)", sum(keep), nrow(mat), 100*cfg$min_valid))
  mat <- mat[keep, , drop = FALSE]

  # Normalization (on log scale)
  if (cfg$normalize == "median") {
    med <- apply(mat, 2, median, na.rm = TRUE)
    mat <- sweep(mat, 2, med - mean(med), "-")
  } else if (cfg$normalize == "quantile") {
    mat <- limma::normalizeQuantiles(mat)
  }

  # Optional batch correction (common for TMT plexes)
  if (!is.null(cfg$batch_col) && have_sva && cfg$batch_col %in% colnames(meta)) {
    log_msg("  ComBat batch correction on ", cfg$batch_col)
    # ComBat needs complete data; impute first for the correction, then re-NA is overkill
    tmp <- mat; tmp[is.na(tmp)] <- min(mat, na.rm = TRUE)
    mod <- model.matrix(~ factor(meta[[cfg$condition]]))
    mat <- sva::ComBat(dat = tmp, batch = factor(meta[[cfg$batch_col]]), mod = mod)
  }

  # Missing-value imputation
  if (cfg$impute == "knn" && have_impute) {
    log_msg("  KNN imputation")
    mat <- impute::impute.knn(mat)$data
  } else if (cfg$impute == "mindet") {
    # Minimal-value (left-censored) imputation: MNAR assumption, common in DIA/LFQ
    log_msg("  MinDet (left-censored) imputation")
    for (j in seq_len(ncol(mat))) {
      v <- mat[, j]; q <- quantile(v, 0.01, na.rm = TRUE)
      v[is.na(v)] <- q; mat[, j] <- v
    }
  } else {
    # none: drop proteins with any NA so limma can run
    mat <- mat[rowSums(is.na(mat)) == 0, , drop = FALSE]
    log_msg(sprintf("  No imputation; %d complete proteins retained", nrow(mat)))
  }
  mat
}

# -----------------------------------------------------------------------------
# 4. DIFFERENTIAL ABUNDANCE (limma)  — per contrast
# -----------------------------------------------------------------------------
run_limma <- function(mat, meta, ct, cfg) {
  grp <- factor(meta[[cfg$condition]])
  design <- model.matrix(~ 0 + grp); colnames(design) <- levels(grp)
  fit <- lmFit(mat, design)
  contr <- makeContrasts(contrasts = paste0(ct$numerator, "-", ct$denominator), levels = design)
  fit2 <- eBayes(contrasts.fit(fit, contr))
  tt <- topTable(fit2, number = Inf, sort.by = "P")
  tt$protein <- rownames(tt)
  tt$significant <- !is.na(tt$adj.P.Val) & tt$adj.P.Val < cfg$padj & abs(tt$logFC) >= cfg$lfc
  tt
}

plot_de <- function(tt, ct, cfg) {
  cdir <- file.path(cfg$outdir, "de", ct$name); dir.create(cdir, showWarnings = FALSE, recursive = TRUE)
  write.csv(tt, file.path(cdir, "differential_abundance.csv"), row.names = FALSE)
  # Volcano
  p <- ggplot(tt, aes(logFC, -log10(P.Value), color = significant)) +
    geom_point(size = 0.7, alpha = 0.6) +
    scale_color_manual(values = c("grey70", "firebrick")) +
    geom_vline(xintercept = c(-cfg$lfc, cfg$lfc), linetype = "dashed") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
    labs(title = paste0(ct$name, " (limma)"), x = "log2 fold-change", y = "-log10 p") +
    theme_bw() + theme(legend.position = "none")
  ggsave(file.path(cdir, "volcano.pdf"), p, width = 6, height = 5)
  list(dir = cdir, sig = tt$protein[which(tt$significant)])
}

de_heatmap <- function(mat, meta, sig, ct, cdir, cfg) {
  sig <- intersect(sig, rownames(mat)); if (length(sig) < 2) return(invisible())
  m <- t(scale(t(mat[sig, , drop = FALSE])))
  ann <- HeatmapAnnotation(group = meta[[cfg$condition]])
  pdf(file.path(cdir, "DA_heatmap.pdf"), width = 8, height = 9)
  draw(Heatmap(m, name = "z", top_annotation = ann, show_row_names = length(sig) <= 60,
               show_column_names = ncol(m) <= 50,
               col = colorRamp2(c(-2,0,2), c("navy","white","firebrick")),
               column_title = paste0(ct$name, " differentially abundant proteins")))
  dev.off()
}

# -----------------------------------------------------------------------------
# 5. ENRICHMENT  (fgsea GSEA + optional clusterProfiler ORA)
# -----------------------------------------------------------------------------
get_genesets <- function(cfg) {
  h  <- msigdbr(species = msig_species, category = "H")
  re <- msigdbr(species = msig_species, category = "C2", subcategory = "CP:REACTOME")
  go <- msigdbr(species = msig_species, category = "C5", subcategory = "GO:BP")
  to <- function(d) split(d$gene_symbol, d$gs_name)
  list(Hallmark = to(h), Reactome = to(re), GO_BP = to(go))
}

run_enrichment <- function(tt, genesets, ct, cdir, cfg) {
  # Rank by signed -log10 p (needs gene symbols as names)
  r <- tt %>% dplyr::filter(!is.na(P.Value)) %>%
    dplyr::mutate(stat = sign(logFC) * -log10(P.Value)) %>%
    dplyr::arrange(dplyr::desc(stat))
  ranks <- setNames(r$stat, r$protein)
  write.table(data.frame(r$protein, r$stat), file.path(cdir, paste0(ct$name, ".rnk")),
              sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  for (db in names(genesets)) {
    res <- tryCatch(fgsea(genesets[[db]], ranks, minSize = 10, maxSize = 500), error = function(e) NULL)
    if (!is.null(res) && nrow(res) > 0) {
      out <- as.data.frame(res); out$leadingEdge <- sapply(out$leadingEdge, paste, collapse = ";")
      write.csv(out[order(out$padj), ], file.path(cdir, paste0("fgsea_", db, ".csv")), row.names = FALSE)
    }
  }
  if (have_cp) {
    sig <- tt$protein[which(tt$significant)]
    if (length(sig) >= 10) {
      ego <- tryCatch(enrichGO(sig, OrgDb = get(org_pkg), keyType = "SYMBOL", ont = "BP",
                               qvalueCutoff = 0.05, readable = TRUE), error = function(e) NULL)
      if (!is.null(ego) && nrow(as.data.frame(ego)) > 0) {
        write.csv(as.data.frame(ego), file.path(cdir, "ORA_GO_BP.csv"), row.names = FALSE)
        ggsave(file.path(cdir, "ORA_GO_dotplot.pdf"), dotplot(ego, showCategory = 20), width = 8, height = 7)
      }
    }
  }
}

# -----------------------------------------------------------------------------
# MAIN
# -----------------------------------------------------------------------------
main <- function(cfg) {
  d <- load_data(cfg)
  run_qc(d$mat, d$meta, cfg)
  pca_plot(d$mat, d$meta, cfg, tag = "prenorm")
  mat <- filter_normalize_impute(d$mat, d$meta, cfg)
  pca_plot(mat, d$meta, cfg, tag = "postnorm")

  if (is.null(cfg$contrasts)) stop("No contrasts in config.")
  genesets <- get_genesets(cfg)
  for (ct in cfg$contrasts) {
    log_msg("==== Contrast: ", ct$name, " ====")
    keep <- d$meta[[cfg$condition]] %in% c(ct$numerator, ct$denominator)
    tt <- run_limma(mat[, keep, drop = FALSE], d$meta[keep, , drop = FALSE], ct, cfg)
    res <- plot_de(tt, ct, cfg)
    de_heatmap(mat[, keep, drop = FALSE], d$meta[keep, , drop = FALSE], res$sig, ct, res$dir, cfg)
    run_enrichment(tt, genesets, ct, res$dir, cfg)
    log_msg(sprintf("  %s: %d differentially abundant proteins", ct$name, length(res$sig)))
  }
  writeLines(capture.output(sessionInfo()), file.path(cfg$outdir, "sessionInfo.txt"))
  log_msg("Pipeline complete. Outputs in:", cfg$outdir)
}

if (sys.nframe() == 0) main(cfg)
