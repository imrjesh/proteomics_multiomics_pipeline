#!/usr/bin/env Rscript
# =============================================================================
# prepare_cptac.R
# Reshapes LinkedOmics CPTAC-LUAD proteome files into the abundance + metadata
# tables the proteomics pipeline expects.
#
# Inputs (LinkedOmics .cct, genes x samples, TMT log2 ratios):
#   HS_CPTAC_LUAD_proteome_ratio_NArm_TUMOR.cct   (110 tumor samples)
#   HS_CPTAC_LUAD_proteome_ratio_NArm_NORMAL.cct  (101 normal samples)
#
# Outputs (into data/):
#   cptac_luad_abundance.tsv   genes x (tumor+normal) samples, columns renamed
#   cptac_luad_metadata.tsv    sample, sample_type (Tumor/Normal), patient
#
# Note: tumor and normal files share patient IDs (matched samples), so columns
# are renamed with a _Tumor / _Normal suffix to keep them distinct after merge.
#
# Usage: Rscript prepare_cptac.R
# =============================================================================

tumor_file  <- "HS_CPTAC_LUAD_proteome_ratio_NArm_TUMOR.cct"
normal_file <- "HS_CPTAC_LUAD_proteome_ratio_NArm_NORMAL.cct"
dir.create("data", showWarnings = FALSE, recursive = TRUE)

read_cct <- function(path) {
  m <- read.delim(path, row.names = 1, check.names = FALSE)
  as.matrix(m)   # genes x samples
}

cat("Reading tumor and normal proteome matrices...\n")
tum <- read_cct(tumor_file)
nor <- read_cct(normal_file)

cat(sprintf("  tumor:  %d genes x %d samples\n", nrow(tum), ncol(tum)))
cat(sprintf("  normal: %d genes x %d samples\n", nrow(nor), ncol(nor)))

# Keep genes present in BOTH, in the same order
common_genes <- intersect(rownames(tum), rownames(nor))
tum <- tum[common_genes, , drop = FALSE]
nor <- nor[common_genes, , drop = FALSE]
cat(sprintf("  common genes: %d\n", length(common_genes)))

# Rename columns with suffix so matched patient IDs don't collide on merge
colnames(tum) <- paste0(colnames(tum), "_Tumor")
colnames(nor) <- paste0(colnames(nor), "_Normal")

# Merge side by side (same gene order guaranteed above)
mat <- cbind(tum, nor)
cat(sprintf("  merged: %d genes x %d samples\n", nrow(mat), ncol(mat)))

# Build metadata
meta <- data.frame(
  sample      = colnames(mat),
  sample_type = c(rep("Tumor", ncol(tum)), rep("Normal", ncol(nor))),
  patient     = sub("_(Tumor|Normal)$", "", colnames(mat)),
  stringsAsFactors = FALSE
)

# Write outputs (protein as first column, matching pipeline's read.delim)
write.table(data.frame(protein = rownames(mat), mat, check.names = FALSE),
            "data/cptac_luad_abundance.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(meta, "data/cptac_luad_metadata.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

cat("\nWrote:\n")
cat("  data/cptac_luad_abundance.tsv (", nrow(mat), "genes x", ncol(mat), "samples )\n")
cat("  data/cptac_luad_metadata.tsv (", nrow(meta), "samples:",
    sum(meta$sample_type == "Tumor"), "Tumor,", sum(meta$sample_type == "Normal"), "Normal )\n")
