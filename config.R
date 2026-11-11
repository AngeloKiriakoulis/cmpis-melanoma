## =====================================================================
## Central configuration for the CMPIS melanoma pipeline.
## Single source of truth for paths, filenames, parameters and the seed.
## Sourced at the start of every stage (see R/run_stage.R) so settings
## and determinism apply uniformly.
## =====================================================================

suppressWarnings(suppressMessages({
  if (requireNamespace("here", quietly = TRUE)) library(here)
}))

## Reproducibility -----------------------------------------------------
SEED <- 20240101L
set.seed(SEED)
options(stringsAsFactors = FALSE)

cfg <- list(
  seed = SEED,

  ## Cache behaviour ---------------------------------------------------
  ## TRUE  = a stage is SKIPPED when its declared output files already
  ##         exist (adopts the caches already in your folder, so the
  ##         first run does not re-download ~6 GB or re-run finished work).
  ## FALSE = always re-run every stage from scratch.
  ## To force a single stage to rebuild, delete its output file(s).
  reuse_cache = TRUE,

  ## curated inputs (tracked in git) -----------------------------------
  pcd_gene_set  = "PCD GENES.xlsx",   # full MSigDB PCD gene set
  pcdrg_genes   = "pcdrg_genes.txt",  # final 8 CMPIS genes

  ## stage caches / key intermediate outputs ---------------------------
  tcga_rdata    = "skcm.tcga.RData",              # Stage 1a output
  scaled_expr   = "scaled_expression_data.csv",   # Stage 1b output
  tcga_expr     = "tcga_skcm_expression.txt",     # Stage 2 outputs
  tcga_clinical = "tcga_skcm_clinical.txt",
  meta_cohort   = "meta_cohort.txt",
  geo_cohorts   = c("gse19234_cohort.txt", "gse98394_cohort.txt",
                    "gse65904_cohort.txt", "gse53118_cohort.txt"),

  ## external data (manual download — see data/README.md) --------------
  depmap_expr   = "Expression_Public_25Q3_subsetted.csv",
  gdsc_drug     = "Drug_sensitivity_AUC_(Sanger_GDSC2)_subsetted.csv",

  ## machine-learning stage --------------------------------------------
  ml_output_dir = "skcm_analysis_results",

  ## sensitivity analysis (Stage 6) --------------
  sens_output   = file.path("results", "pcd_sensitivity.rds"),

  ## classifier benchmark (Stage 7) -------------
  bench_output  = file.path("results", "classifier_benchmark.rds"),

  ## anti-PD-1 validation (Stage 8) --------------
  icb_output    = file.path("results", "icb_validation.rds"),

  ## mechanistic network figure (Stage 9) --------
  net_output    = file.path("results", "mechanistic_network.rds"),

  ## CMPIS weights (Stage 10) ------------------------------------
  ## "corrected" = final weights (default; used throughout the paper)
  ## "original"  = alternative ensemble weights, applied inverted (sensitivity)
  cmpis_version = "corrected",
  cmpis_output  = file.path("results", "cmpis_reweight.rds"),

  ## CMPIS clinical panels (Stage 11)
  clin_output   = file.path("results", "cmpis_clinical.rds"),

  ## single-cell validation (Stage 13: Fig 7)
  sc_output     = file.path("results", "singlecell.rds"),

  ## parameters --------------------------------------------------------
  n_clusters    = 2L                  # k selected by CPI + gap statistic
)
