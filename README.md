# proteomics_multiomics_pipeline
Multi-omics proteogenomics pipeline (CPTAC LUAD): limma differential abundance, QC, subtype clustering, ssGSEA, co-expression modules, mRNA-protein &amp; CNV-protein integration
# Multi-Omics Proteogenomics Pipeline (CPTAC LUAD)

A config-driven pipeline for downstream analysis of cancer proteomics data,
with multi-omics integration. Built and tested on CPTAC lung adenocarcinoma
(LUAD) TMT proteomics, with branches for DIA label-free data.

> Scope: starts from protein/omics abundance matrices (the output of spectral
> search tools like FragPipe/DIA-NN, or processed CPTAC matrices). Raw-file
> spectral search is an upstream step, documented below, not run here.

## Analyses

**Primary (differential abundance)**
- QC (missingness, intensity distributions, PCA, correlation)
- Filtering, normalization, missing-value imputation (KNN / MinDet)
- Differential abundance with limma (multi-contrast)
- Volcano, heatmap, ranked files, GSEA (Hallmark/Reactome/GO) + ORA

**Secondary (exploratory)**
- Sample correlation + hierarchical clustering
- Unsupervised sample clustering (candidate subtypes)
- ssGSEA per-sample pathway activity
- Co-expression modules (WGCNA-lite)

**Multi-omics integration**
- mRNA-protein correlation (post-transcriptional regulation)
- CNV-protein cis effects (genomic amplification to protein)
- Phosphosite tumor-vs-normal; kinase activity (optional, needs substrate sets)

## Platform branches

| Branch | Input | Set in config |
|--------|-------|---------------|
| TMT | log2-ratio matrix (e.g. CPTAC) | `platform: TMT`, `is_log: true` |
| DIA | linear intensities (DIA-NN report) | `platform: DIA`, `is_log: false` |

## Quickstart (toy data — clone and run)

```bash
# Apple Silicon: prefix with CONDA_SUBDIR=osx-64
mamba env create -f environment.yml
conda activate proteomics-downstream

Rscript make_toy_proteomics.R
Rscript R/proteomics_downstream_pipeline.R config/config_toy.yaml
```

## Running on CPTAC LUAD (real data)

1. Download from LinkedOmics (https://linkedomics.org/data_download/CPTAC-LUAD/):
   - Proteome Tumor + Normal (TMT log ratio, gene level)
   - RNAseq Tumor + Normal (log2 FPKM)   [for multi-omics]
   - Phosphoproteome Tumor + Normal (site level)   [for multi-omics]
   - SCNV gene-level (log ratios)   [for multi-omics]
2. Place the `.cct` files in the project root (or `data/`).
3. Reshape proteome into pipeline inputs:
```bash
   Rscript prepare_cptac.R
```
4. Run the three stages:
```bash
   Rscript R/proteomics_downstream_pipeline.R config/config_cptac_luad.yaml
   Rscript R/proteomics_secondary_analysis.R  config/config_cptac_luad.yaml
   Rscript R/proteomics_multiomics.R          config/config_cptac_luad.yaml
```

See `METHODS.md` for details of each step and `INSTALL.md` for full setup.

## Upstream (primary search — documented, not run here)

Raw instrument files (.raw/.d) are processed with spectral-search tools
(FragPipe/MSFragger or DIA-NN for DIA; CPTAC Common Data Analysis Pipeline for
the TMT matrices used here) to produce the abundance matrices this pipeline
consumes. This mirrors starting an RNA-seq pipeline from a count matrix.

## Data note

CPTAC data is public (LinkedOmics / PDC). Large `.cct` files are not committed;
download them via the links above. A small synthetic dataset is included so the
pipeline runs immediately after cloning.

## License
MIT
