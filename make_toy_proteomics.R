#!/usr/bin/env Rscript
# =============================================================================
# make_toy_proteomics.R
# Generates a synthetic proteomics dataset so the downstream pipeline runs
# end-to-end with no external data. Mimics a TMT-style log2 abundance matrix
# (tumor vs normal) with realistic missing values (MNAR: low-abundance proteins
# more likely to be missing). Uses real human gene symbols so enrichment works.
#
# Produces (into data/):
#   toy_protein_abundance.tsv   proteins x samples (gene symbols as rownames)
#   toy_metadata.tsv            sample sheet: sample, sample_type, plex
#
# Usage: Rscript make_toy_proteomics.R [n_proteins] [n_per_group]
#   default: 3000 proteins, 50 tumor + 50 normal
# =============================================================================
set.seed(42)
args <- commandArgs(trailingOnly = TRUE)
n_prot <- ifelse(length(args) >= 1, as.integer(args[1]), 3000L)
per_grp <- ifelse(length(args) >= 2, as.integer(args[2]), 50L)
dir.create("data", showWarnings = FALSE, recursive = TRUE)

# Sample sheet: tumor/normal, split across two TMT "plexes" (batch effect)
samples <- c(sprintf("T%03d", 1:per_grp), sprintf("N%03d", 1:per_grp))
grp <- c(rep("Tumor", per_grp), rep("Normal", per_grp))
plex <- rep(c("plex1", "plex2"), length.out = length(samples))
meta <- data.frame(sample = samples, sample_type = grp, plex = plex, stringsAsFactors = FALSE)

# Gene symbols from org.Hs.eg.db if available, else synthetic
genes <- tryCatch({ suppressPackageStartupMessages(library(org.Hs.eg.db))
  sample(keys(org.Hs.eg.db, "SYMBOL"), n_prot) },
  error = function(e) sprintf("PROT%04d", seq_len(n_prot)))

# Baseline log2 abundance per protein
base <- rnorm(n_prot, mean = 20, sd = 2)
mat <- matrix(rep(base, length(samples)), nrow = n_prot)
colnames(mat) <- samples; rownames(mat) <- genes

# True differential proteins (10%): up or down in tumor
n_diff <- round(0.10 * n_prot)
up <- sample(n_prot, n_diff %/% 2); dn <- sample(setdiff(seq_len(n_prot), up), n_diff %/% 2)
fc_up <- runif(length(up), 1, 3); fc_dn <- runif(length(dn), 1, 3)

for (j in seq_along(samples)) {
  mat[, j] <- base + rnorm(n_prot, 0, 0.5)                  # biological + technical noise
  if (grp[j] == "Tumor") { mat[up, j] <- mat[up, j] + fc_up; mat[dn, j] <- mat[dn, j] - fc_dn }
  if (plex[j] == "plex2") mat[, j] <- mat[, j] + 0.6        # batch shift (ComBat should remove)
}

# MNAR missingness: lower-abundance proteins more likely missing (~15% overall)
prob_missing <- plogis(-(rowMeans(mat) - mean(mat)) * 0.8) * 0.3
for (i in seq_len(n_prot)) {
  miss <- runif(length(samples)) < prob_missing[i]
  mat[i, miss] <- NA
}

write.table(data.frame(protein = rownames(mat), mat, check.names = FALSE),
            "data/toy_protein_abundance.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(meta, "data/toy_metadata.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

cat(sprintf("Wrote toy proteomics data to data/\n  %d proteins x %d samples (%d Tumor, %d Normal)\n",
            n_prot, length(samples), per_grp, per_grp))
cat(sprintf("  true differential: %d up, %d down in Tumor; ~%.0f%% missing (MNAR)\n",
            length(up), length(dn), 100*mean(is.na(mat))))
