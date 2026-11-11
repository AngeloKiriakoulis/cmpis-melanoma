## =====================================================================
## run_stage(): execute a legacy analysis script as a pipeline stage.
##
## Each of the original scripts is a large, interactive analysis file.
## Rather than rewriting them, we run each in a CLEAN child R process
## (via callr) from the project root, with a few environment-mutating and
## interactive calls neutralised so the stage runs unattended:
##
##   * setwd()             -> no-op  (keep the project root as the wd; this
##                            also fixes scripts that setwd() to an old
##                            absolute path)
##   * install.packages()  -> no-op  (renv provides packages; the pipeline
##   * BiocManager::install-> no-op   must not install mid-run)
##   * View()              -> no-op  (non-interactive)
##
## Cache behaviour: if `skip_if_exists` is TRUE (default, driven by
## cfg$reuse_cache) and every declared output already exists, the stage is
## SKIPPED and its outputs returned as-is. This lets the first run adopt
## the caches already in your folder instead of recomputing them. Delete a
## stage's output file(s) to force it to rebuild.
##
## Returning the output paths lets {targets} track them with format="file".
## =====================================================================

run_stage <- function(script, outputs = character(0), config = "config.R",
                      skip_if_exists = if (exists("cfg") && !is.null(cfg$reuse_cache))
                                         isTRUE(cfg$reuse_cache) else TRUE) {

  if (isTRUE(skip_if_exists) && length(outputs) > 0 && all(file.exists(outputs))) {
    message("--- skipping '", script, "': cached output(s) already present (",
            paste(outputs, collapse = ", "), ").")
    return(outputs)
  }

  callr::r(
    function(script, config) {
      if (file.exists(config)) source(config, local = FALSE)

      assign("setwd",            function(...) invisible(NULL), envir = globalenv())
      assign("install.packages", function(...) invisible(NULL), envir = globalenv())
      assign("View",             function(...) invisible(NULL), envir = globalenv())
      if (requireNamespace("BiocManager", quietly = TRUE)) {
        try(assignInNamespace("install", function(...) invisible(NULL),
                              ns = "BiocManager"), silent = TRUE)
      }

      message(">>> running stage script: ", script)
      source(script, local = FALSE, echo = TRUE, max.deparse.length = Inf)
      invisible(TRUE)
    },
    args    = list(script = script, config = config),
    show    = TRUE,
    spinner = FALSE
  )

  missing <- outputs[!file.exists(outputs)]
  if (length(missing)) {
    stop("Stage '", script, "' did not produce expected output(s): ",
         paste(missing, collapse = ", "), call. = FALSE)
  }
  outputs
}
