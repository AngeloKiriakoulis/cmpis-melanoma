#!/usr/bin/env Rscript
## =====================================================================
## One-command entrypoint for the CMPIS pipeline.
##   Rscript run.R          (from a shell, at the project root)
##   source("run.R")        (from an R session)
## =====================================================================

## On a fresh checkout, restore the pinned environment first:
if (!requireNamespace("renv", quietly = TRUE)) install.packages("renv")
## renv::restore()   # <- uncomment on first setup to install pinned packages

## Pipeline machinery (add to renv if missing: renv::install(c("targets","callr","here")))
for (p in c("targets", "callr", "here")) {
  if (!requireNamespace(p, quietly = TRUE)) {
    stop("Missing package '", p, "'. Install it with: renv::install('", p, "')",
         call. = FALSE)
  }
}

## Load convenience helpers: rerun("<stage>"), status(), graph()
source("R/helpers.R")

targets::tar_make()

## After a run:
##   targets::tar_visnetwork()   # see the DAG + what is up to date
##   targets::tar_progress()     # per-stage status
