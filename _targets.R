## =====================================================================
## {targets} pipeline for the CMPIS melanoma study.
##
##   Run:      source("run.R")   or   targets::tar_make()
##   Inspect:  targets::tar_visnetwork()   /   targets::tar_manifest()
##
## Each stage wraps an existing analysis script (see R/run_stage.R) and is
## gated on its declared input/output files: an unchanged stage is
## SKIPPED, and the heavy Stage 1 downloads (multi-GB TCGA/GEO) run only
## when their output cache is missing. Referencing an upstream target in a
## stage's body is what creates the dependency edge.
##
## Supplementary scripts UMAP&Clinical.R and Mutation&DESeq2.R are not in
## the core DAG; add them as extra targets if you want them in the run.
## =====================================================================

library(targets)

source("config.R")
source("R/run_stage.R")

tar_option_set(
  packages = character(0),   # each stage sources its own library() calls
  format   = "rds",
  error    = "abort"
)

list(
  ## tracked curated input -------------------------------------------------
  tar_target(pcd_gene_set, cfg$pcd_gene_set, format = "file"),

  ## Stage 1a — TCGA-SKCM multi-omics download + preprocess ----------------
  tar_target(
    tcga_preprocess,
    run_stage("TCGA_Preprocessing.R", outputs = cfg$tcga_rdata),
    format = "file"
  ),

  ## Stage 1b — GEO cohorts + ComBat batch correction ---------------------
  tar_target(
    geo_scaled,
    run_stage("GEO.R", outputs = cfg$scaled_expr),
    format = "file"
  ),

  ## Stage 2 — MOVICS subtyping (CS1/CS2) + characterisation --------------
  tar_target(
    subtypes,
    {
      c(tcga_preprocess, geo_scaled, pcd_gene_set)   # dependency edges
      run_stage(
        "TCGA_Clustering.R",
        outputs = c(cfg$tcga_expr, cfg$tcga_clinical,
                    cfg$meta_cohort, cfg$pcdrg_genes)
      )
    },
    format = "file"
  ),

  ## Stage 4 — CMPIS machine-learning signature --------------------------
  tar_target(
    cmpis_model,
    {
      c(subtypes)                                    # dependency edge
      run_stage(
        "SKCM_ML_Algorithms.R",
        outputs = file.path(cfg$ml_output_dir, "final_model_info.rds")
      )
    },
    format = "file"
  ),

  ## Stage 5 — validation, immunotherapy, pharmacogenomics ---------------
  tar_target(
    validation,
    {
      c(tcga_preprocess, geo_scaled, cmpis_model)    # dependency edges
      run_stage(
        "Validation_PCDRGS.R",
        outputs = "Plot_D_MPCDS_by_response.pdf"
      )
    },
    format = "file"
  ),

  ## Stage 6 — PCD gene-set sensitivity analysis ----
  ## Per-pathway Supplementary Table 1 + sensitivity of the CS1/CS2
  ## two-cluster solution to PCD pathway choice. Depends on the subtypes
  ## (reference CS1/CS2 labels) and the curated PCD gene set; reuses the
  ## existing skcm.tcga / TCGA_Clustering caches, no downloads.
  tar_target(
    pcd_sensitivity,
    {
      c(subtypes, pcd_gene_set)                      # dependency edges
      run_stage(
        "R/06_PCD_Sensitivity.R",
        outputs = cfg$sens_output
      )
    },
    format = "file"
  ),

  ## Stage 7 — benchmark CS1/CS2 (+MPCDS) vs established classifiers ------
  ## Depends on the subtypes (CS1/CS2 reference) and
  ## the CMPIS model (MPCDS coefficients). Pure downstream analysis on the
  ## existing TCGA cache — no downloads.
  tar_target(
    classifier_benchmark,
    {
      c(subtypes, cmpis_model, cmpis_reweight)       # dependency edges
      run_stage(
        "R/07_Classifier_Benchmark.R",
        outputs = cfg$bench_output
      )
    },
    format = "file"
  ),

  ## Stage 8 — anti-PD-1 cohort validation ---------
  ## Applies MPCDS to real anti-PD-1 melanoma cohorts (GSE78220, GSE91061)
  ## with RECIST response. Depends on the CMPIS model (MPCDS coefficients).
  tar_target(
    icb_validation,
    {
      c(cmpis_model, cmpis_reweight)
      run_stage(
        "R/08_ICB_Validation.R",
        outputs = cfg$icb_output
      )
    },
    format = "file"
  ),

  ## Stage 9 — mechanistic network figure ----------
  tar_target(
    mechanistic_network,
    {
      c(subtypes, cmpis_model)
      run_stage(
        "R/09_Mechanistic_Network.R",
        outputs = cfg$net_output
      )
    },
    format = "file"
  ),

  ## Stage 10 — corrected CMPIS weights + re-validation ------------------
  tar_target(
    cmpis_reweight,
    {
      c(subtypes, cmpis_model)
      run_stage("R/10_CMPIS_Reweight.R", outputs = cfg$cmpis_output)
    },
    format = "file"
  ),

  ## Stage 11 — CMPIS manuscript panels (Figs S3-S6) on corrected score --
  tar_target(
    cmpis_clinical,
    {
      c(subtypes, cmpis_reweight)
      run_stage("R/11_CMPIS_Clinical.R", outputs = cfg$clin_output)
    },
    format = "file"
  )
)
