## =====================================================================
## Reference list of packages the pipeline uses, grouped by source.
## The AUTHORITATIVE environment is renv (renv.lock) — prefer
## `renv::restore()`. This file is a human-readable manifest and a
## fallback installer, extracted from the analysis scripts.
## =====================================================================

## --- Pipeline machinery ----------------------------------------------
# install.packages(c("targets", "callr", "here", "renv"))

## --- CRAN ------------------------------------------------------------
cran <- c(
  "dplyr","tidyr","tibble","readr","readxl","ggplot2","reshape2",
  "Rtsne","umap","circlize","RColorBrewer","pheatmap","viridis",
  "patchwork","gridExtra","cowplot","ggpubr","lubridate","forestplot",
  "survival","survminer","timeROC","rms","Hmisc","pec",
  "glmnet","randomForestSRC","gbm","plsRcox","survivalsvm","superpc",
  "foreach","doParallel","msigdbr","devtools","remotes","GOplot"
)
# install.packages(cran)

## --- Bioconductor ----------------------------------------------------
bioc <- c(
  "TCGAbiolinks","SummarizedExperiment","biomaRt","maftools","EDASeq",
  "DESeq2","org.Hs.eg.db","AnnotationDbi","edgeR","minfi",
  "IlluminaHumanMethylation450kanno.ilmn12.hg19","GEOquery","sva","limma",
  "Biobase","ComplexHeatmap","dorothea","viper","GSVA","illuminaHumanv3.db",
  "hgu133plus2.db","clusterProfiler","enrichplot","EnhancedVolcano",
  "depmap","ExperimentHub"
)
# if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager")
# BiocManager::install(bioc)

## --- GitHub ----------------------------------------------------------
# remotes::install_github("xlucpu/MOVICS")
# remotes::install_github("IOBR/IOBR")
# remotes::install_github("danro9685/CIMLR", ref = "R")   # used by MOVICS
# also: CoxBoost, MAFDash  (CRAN archive / GitHub, version-dependent)

## --- Local vendored tarballs (see data/README.md) --------------------
# install.packages("pRRophetic_0.5.tar.gz", repos = NULL, type = "source")
# install.packages("oncoPredict_1.2.tar.gz", repos = NULL, type = "source")
# install.packages("IMvigor210CoreBiologies_1.0.1.tar.gz", repos = NULL, type = "source")
