## =====================================================================
## R/helpers.R — convenience helpers for module-by-module runs.
## Sourced by run.R. You can also source this file alone to get the
## helpers without running the pipeline:  source("R/helpers.R")
## =====================================================================

## Stage (target) names, for reference:
##   tcga_preprocess  Stage 1a  (TCGA_Preprocessing.R)
##   geo_scaled       Stage 1b  (GEO.R)
##   subtypes         Stage 2   (TCGA_Clustering.R)
##   cmpis_model      Stage 4   (SKCM_ML_Algorithms.R)
##   validation       Stage 5   (Validation_PCDRGS.R)

## rerun(): re-run ONE stage (and anything downstream of it), reusing all
## other cached results. Use this after editing that stage's script — it
## invalidates only that target, then rebuilds what depends on it.
##   rerun("validation")
rerun <- function(stage) {
  stopifnot(is.character(stage), length(stage) == 1L)
  targets::tar_invalidate(tidyselect::any_of(stage))
  targets::tar_make()
}

## status(): per-stage status of the last run (built / skipped / errored).
status <- function() targets::tar_progress()

## graph(): open the dependency graph + what is currently up to date.
graph <- function() targets::tar_visnetwork()
