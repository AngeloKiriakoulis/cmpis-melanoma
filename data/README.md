# Data manifest

Large and regenerable files are **not committed** to git. This lists what they are,
where they go, and how to obtain them. All data are public.

## Downloaded automatically by the pipeline
- **TCGA-SKCM multi-omics** (RNA-seq, methylation 450K, somatic mutation, miRNA,
  clinical): `TCGA_Preprocessing.R` via `TCGAbiolinks` into `GDCdata/`, cached to
  `skcm.tcga.RData`.
- **GEO bulk cohorts** (GSE15605, GSE22153, GSE46517, GSE54467, GSE19234, GSE98394,
  GSE65904, GSE53118): `GEO.R` via `GEOquery`; batch-corrected output
  `scaled_expression_data.csv`.
- **GEO single-cell datasets** (GSE72056 Tirosh, GSE115978 Jerby-Arnon, GSE120575
  Sade-Feldman): Stage 13, into `data/sc/` (~1 GB).
- **Checkpoint-inhibitor cohorts** (Liu 2019, Riaz 2017, Hugo 2016, Gide 2019 from the
  CRI iAtlas harmonised studies; Van Allen 2015): Stage 15, from the public cBioPortal
  API (https://www.cbioportal.org), cached in `data/icb/`.
- **GSE78220 / GSE91061** (Stage 8 only): via `GEOquery`.

## Manual downloads required
- **`Expression_Public_25Q3_subsetted.csv`** (~526 MB): DepMap 25Q3 cell-line
  expression, subset to the genes used. Source: https://depmap.org/portal/download/
- **`Drug_sensitivity_AUC_(Sanger_GDSC2)_subsetted.csv`**: GDSC2 drug AUC (committed;
  original source https://www.cancerrxgene.org / DepMap).
- **`GSE98394_expression.txt`** (~17 MB): validation expression matrix.

## Vendored R packages (local tarballs, not committed)
Not on CRAN/Bioconductor in the required versions; install after `renv::restore()`:
- `pRRophetic_0.5.tar.gz`
- `oncoPredict_1.2.tar.gz`
- `IMvigor210CoreBiologies_1.0.0.tar.gz` / `1.0.1.tar.gz`
- `DESeq_1.12.1.tar.gz`

```r
install.packages("pRRophetic_0.5.tar.gz", repos = NULL, type = "source")
# repeat for the others
```

## Regenerable intermediates (not committed)
`*.RData` stage caches (`TCGA_Preprocessing.RData` ~5.8 GB, `TCGA_Clustering.RData`
~1.9 GB, `umap.RData`, `mutation&deseq.RData`, `Cox.RData`, `GEO.RData`,
`Validation_PCDRGS.RData`, `SKCM_ML_Algorithms.RData`), `de.csv`, `expr_est_raw.txt`,
`TCGA_SKCM_gene_expression_counts.csv`, and the generated MOVICS output directories
(`1-1-SKCM_MPCDS-*`, `1-GSEA-result`, most of `skcm_analysis_results`).
