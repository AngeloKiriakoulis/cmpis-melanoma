## =====================================================================
## R/cmpis_weights.R — single source of truth for the CMPIS/MPCDS score.
##
## Final weights (Stage 10, "corrected" in the code): Cox coefficient x
##   normalised RSF importance; sign from Cox, magnitude modulated by RSF;
##   used without inversion.
## Alternative ensemble weights ("original" in the code, sensitivity only):
##   (Cox_std + RSFimportance_std)/2, applied with risk_direction = "invert".
##
## load_cmpis()  -> list(coef, sign, version)
## score_cmpis() -> per-sample score, z-scaled, higher = higher risk
## =====================================================================

load_cmpis <- function(prefer = if (exists("cfg") && !is.null(cfg$cmpis_version))
                                   cfg$cmpis_version else "corrected") {
  corr <- file.path("skcm_analysis_results", "final_model_info_corrected.rds")
  orig <- file.path("skcm_analysis_results", "final_model_info.rds")
  if (identical(prefer, "corrected") && file.exists(corr)) {
    m <- readRDS(corr)
    return(list(coef = m$coefficients, sign = +1, version = "corrected"))
  }
  if (file.exists(orig)) {
    f <- readRDS(orig); co <- f$coefficients
    if (is.null(co)) co <- f$ensemble_coefficients
    return(list(coef = co, sign = -1, version = "original (inverted)"))
  }
  NULL
}

score_cmpis <- function(expr, cm = load_cmpis()) {       # expr: genes x samples
  if (is.null(cm)) return(NULL)
  g <- intersect(names(cm$coef), rownames(expr)); if (length(g) < 3) return(NULL)
  m <- as.matrix(expr[g, , drop = FALSE]); storage.mode(m) <- "double"
  if (max(m, na.rm = TRUE) > 50) m <- log2(m + 1)
  z <- t(scale(t(m))); z[is.na(z)] <- 0                   # per-gene z within cohort
  s <- cm$sign * colSums(z * cm$coef[g])
  setNames(as.numeric(scale(s)), colnames(expr))
}
