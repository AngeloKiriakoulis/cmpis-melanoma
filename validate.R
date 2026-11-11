#!/usr/bin/env Rscript
## =====================================================================
## Pre-flight check: confirms the pipeline is wired correctly WITHOUT
## running any stage (no downloads, no compute). Run this first.
##   source("validate.R")
## =====================================================================

ok <- TRUE
for (p in c("targets", "callr", "here")) {
  if (!requireNamespace(p, quietly = TRUE)) {
    message("MISSING package: ", p, "  ->  renv::install('", p, "')"); ok <- FALSE
  }
}
if (!ok) stop("Install the missing pipeline packages above, then re-run validate.R")

source("config.R")
message("Config loaded. reuse_cache = ", cfg$reuse_cache, "; seed = ", cfg$seed)

## Print the stage graph (parses _targets.R; does not run anything)
print(targets::tar_manifest(fields = c("name", "command")))

message("\nWhich stage outputs already exist (these will be SKIPPED on run):")
checks <- c(
  tcga_rdata    = cfg$tcga_rdata,
  scaled_expr   = cfg$scaled_expr,
  tcga_expr     = cfg$tcga_expr,
  tcga_clinical = cfg$tcga_clinical,
  meta_cohort   = cfg$meta_cohort,
  pcdrg_genes   = cfg$pcdrg_genes,
  ml_model      = file.path(cfg$ml_output_dir, "final_model_info.rds")
)
for (nm in names(checks)) {
  message(sprintf("  [%s] %s", if (file.exists(checks[[nm]])) "present" else "MISSING", checks[[nm]]))
}

message("\nWiring OK. When ready, run:  source('run.R')")
