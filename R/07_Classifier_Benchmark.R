## =====================================================================
## Stage 7 — Benchmark CS1/CS2 (and MPCDS) vs established melanoma /
## immune classifiers
## ---------------------------------------------------------------------
## Compares the PCD subtypes and the CMPIS (earlier name: MPCDS) with
## established signatures: prognostic value beyond them, and association
## with immune state.
##
## Comparators computed on TCGA-SKCM mRNA:
##   * IFN-gamma signature (Ayers et al. 2017, 6-gene)
##   * T-cell-inflamed GEP (Ayers et al. 2017, 18-gene)
##   * ESTIMATE ImmuneScore / StromalScore   (needs 'estimate' pkg; gated)
##   * Tsoi differentiation axis (melanocytic vs undifferentiated poles)
##   * MITF program (MITF-high vs MITF-low)
##   * TIDE score            (read from results/tide*.csv if provided; gated)
##
## PCD variables: CS1/CS2 (published subtypes) and MPCDS (8-gene score).
##
## Analyses: univariable Cox + C-index; multivariable Cox (does CS / MPCDS
## stay independent?); nested-model likelihood-ratio tests; C-index
## comparison (base vs base+CS vs base+MPCDS); decision-curve net benefit
## (gated); signature-correlation heatmap.
##
## Outputs (results/):
##   Supplementary_Table_3_benchmark_univariable.csv
##   Supplementary_Table_3_benchmark_multivariable.csv
##   Supplementary_Table_3_incremental_value.csv
##   Figure_S_benchmark.pdf
##   classifier_benchmark.rds   (declared stage output)
##
## Reuses caches — no downloads:
##   skcm.tcga.RData, TCGA_Clustering.RData (cmoic.skcm),
##   skcm_analysis_results/final_model_info.rds
##
## Run:  source("config.R"); source("R/07_Classifier_Benchmark.R")
##   or: source("R/helpers.R"); rerun("classifier_benchmark")
## =====================================================================

suppressWarnings(suppressMessages({
  library(survival)
  have_survcomp <- requireNamespace("survcomp", quietly = TRUE)
  have_gsva     <- requireNamespace("GSVA",     quietly = TRUE)
  have_estimate <- requireNamespace("estimate", quietly = TRUE)
  have_dcurves  <- requireNamespace("dcurves",  quietly = TRUE)
  have_pheatmap <- requireNamespace("pheatmap", quietly = TRUE)
}))
set.seed(if (exists("cfg") && !is.null(cfg$seed)) cfg$seed else 20240101L)
results_dir <- "results"; if (!dir.exists(results_dir)) dir.create(results_dir, recursive = TRUE)
msg <- function(...) message(sprintf(...))

## ---------------------------------------------------------------------
## 1. Load expression, clinical, CS labels, MPCDS coefficients
## ---------------------------------------------------------------------
stopifnot(file.exists("skcm.tcga.RData")); load("skcm.tcga.RData")   # -> skcm.tcga
expr <- as.matrix(skcm.tcga$mRNA.expr); storage.mode(expr) <- "double"
if (max(expr, na.rm = TRUE) > 50) expr <- log2(expr + 1)             # -> log space
clin <- skcm.tcga$clin.info
clin <- clin[!is.na(clin$futime) & clin$futime > 0, , drop = FALSE]

## reference CS1/CS2
CS <- NULL
if (file.exists("TCGA_Clustering.RData")) {
  e <- new.env(); tryCatch(load("TCGA_Clustering.RData", envir = e), error = function(x) NULL)
  if (exists("cmoic.skcm", envir = e)) {
    cr <- get("cmoic.skcm", envir = e)$clust.res
    CS <- setNames(factor(ifelse(cr$clust == 1, "CS1", "CS2"), c("CS1","CS2")), cr$samID)
  }
}
if (is.null(CS)) stop("CS1/CS2 labels not found in TCGA_Clustering.RData — run Stage 2 first.")

## MPCDS (8-gene weighted score, inverted — matches Validation_PCDRGS.R)
pcdrg_genes <- c("CFLAR","CYBB","CYLD","GYS1","IRF1","NFKBIA","PIM2","VCAM1")
source("R/cmpis_weights.R")
cmpis <- load_cmpis(); mpcds <- score_cmpis(expr, cmpis)
if (!is.null(cmpis)) msg("MPCDS weights: %s", cmpis$version)
if (is.null(mpcds)) msg("NOTE: MPCDS unavailable (final_model_info coefficients missing) — CS1/CS2 remains the primary PCD variable.")

## ---------------------------------------------------------------------
## 2. Comparator signature scores
## ---------------------------------------------------------------------
zmean <- function(genes) {
  g <- intersect(genes, rownames(expr)); if (!length(g)) return(setNames(rep(NA_real_, ncol(expr)), colnames(expr)))
  m <- expr[g, , drop = FALSE]; m <- m[apply(m,1,function(x) sd(x,na.rm=TRUE)>0), , drop = FALSE]
  setNames(colMeans(t(scale(t(m))), na.rm = TRUE), colnames(expr))
}
sig <- list()
sig$IFNg <- zmean(c("IFNG","STAT1","IDO1","CXCL9","CXCL10","HLA-DRA"))
sig$GEP  <- zmean(c("CCL5","CD27","CD274","CD276","CD8A","CMKLR1","CXCL9","CXCR6",
                    "HLA-DQA1","HLA-DRB1","HLA-E","IDO1","LAG3","NKG7","PDCD1LG2",
                    "PSMB10","STAT1","TIGIT"))
mel  <- c("MITF","MLANA","PMEL","TYR","TYRP1","DCT","SOX10","GPR143")
undf <- c("AXL","WNT5A","EGFR","NGFR","NT5E","FN1","TGFBI","SERPINE1")
sig$Tsoi_diff <- zmean(mel) - zmean(undf)          # + = differentiated (melanocytic)
sig$MITF <- zmean(c("MITF","MLANA","PMEL","TYR","TYRP1","DCT","GPR143"))

## ESTIMATE — ssGSEA (GSVA) on the ESTIMATE immune/stromal signatures.
## Uses the gene sets shipped with the 'estimate' package but bypasses its
## estimateScore()/filterCommonGenes(), which hit a 'common_genes' lazy-data
## bug in v1.0.13 (that also breaks IOBR's estimate method).
est_ok <- FALSE
if (!requireNamespace("estimate", quietly = TRUE))
  try(utils::install.packages("estimate", repos = "http://R-Forge.R-project.org"), silent = TRUE)
if (requireNamespace("estimate", quietly = TRUE) && requireNamespace("GSVA", quietly = TRUE)) {
  tryCatch({
    e2 <- new.env(); utils::data("SI_geneset", package = "estimate", envir = e2)
    sidf <- as.data.frame(get("SI_geneset", envir = e2), stringsAsFactors = FALSE)
    extract <- function(pat) {                    # robust to row/col/first-col layouts
      cand <- character(0); rn <- rownames(sidf)
      if (!is.null(rn) && any(grepl(pat, rn, ignore.case = TRUE)))
        cand <- as.character(unlist(sidf[grep(pat, rn, ignore.case = TRUE)[1], ]))
      else if (any(grepl(pat, colnames(sidf), ignore.case = TRUE)))
        cand <- as.character(sidf[[grep(pat, colnames(sidf), ignore.case = TRUE)[1]]])
      else if (any(grepl(pat, sidf[[1]], ignore.case = TRUE)))
        cand <- as.character(unlist(sidf[grep(pat, sidf[[1]], ignore.case = TRUE)[1], -1]))
      cand <- trimws(cand); cand <- cand[nzchar(cand) & toupper(cand) != "NA"]
      intersect(unique(cand), rownames(expr))
    }
    gs <- list(ESTIMATE_Stromal = extract("Stromal"), ESTIMATE_Immune = extract("Immune"))
    if (min(lengths(gs)) < 10) stop(sprintf("SI_geneset matched too few genes (stromal=%d, immune=%d)",
                                            length(gs$ESTIMATE_Stromal), length(gs$ESTIMATE_Immune)))
    eu <- expr[, !duplicated(colnames(expr)), drop = FALSE]
    sc <- GSVA::gsva(GSVA::ssgseaParam(exprData = eu, geneSets = gs))
    for (nm in rownames(sc)) sig[[nm]] <- setNames(as.numeric(sc[nm, ]), colnames(eu))[colnames(expr)]
    est_ok <- TRUE
    msg("ESTIMATE (ssGSEA on estimate signatures: %d stromal, %d immune genes) computed.",
        length(gs$ESTIMATE_Stromal), length(gs$ESTIMATE_Immune))
  }, error = function(x) msg("ESTIMATE ssGSEA failed (%s).", conditionMessage(x)))
}
if (!est_ok) msg("NOTE: ESTIMATE unavailable — panel proceeds without it.")

## TIDE (gated on a provided file)
tide_file <- list.files(results_dir, pattern = "^tide.*\\.csv$", full.names = TRUE, ignore.case = TRUE)
if (length(tide_file)) {
  td <- tryCatch(read.csv(tide_file[1], check.names = FALSE), error = function(x) NULL)
  if (!is.null(td)) {
    idc <- names(td)[1]; sc <- grep("TIDE", names(td), ignore.case = TRUE, value = TRUE)[1]
    if (!is.na(sc)) { v <- setNames(td[[sc]], td[[idc]]); sig$TIDE <- v[colnames(expr)] }
    msg("TIDE read from %s", tide_file[1])
  }
} else {
  ti <- expr - apply(expr, 1, median, na.rm = TRUE)          # log2, gene-median-centered
  utils::write.table(data.frame(Gene = rownames(ti), ti, check.names = FALSE),
                     file.path(results_dir, "tide_input.txt"), sep = "\t",
                     quote = FALSE, row.names = FALSE)
  msg("NOTE: no results/tide*.csv — wrote results/tide_input.txt. Upload to tide.dfci.harvard.edu (Cancer=Melanoma, Previous ICB=No), save result as results/tide.csv, re-run.")
}

## ---------------------------------------------------------------------
## 3. Assemble modelling frame
## ---------------------------------------------------------------------
samples <- Reduce(intersect, list(colnames(expr), names(CS), rownames(clin)))
msg("Benchmark cohort: %d samples", length(samples))
cl <- clin[samples, , drop = FALSE]
## detect survival columns robustly (names vary across cached clin.info)
tcol <- intersect(c("futime","time","OS.time","OS_time","days"), names(cl))[1]
scol <- intersect(c("fustat","status","OS","event","vital_status"), names(cl))[1]
if (is.na(tcol) || is.na(scol))
  stop("Could not find survival columns in clin.info (have: ",
       paste(names(cl), collapse = ", "), ")")
dat <- data.frame(time   = suppressWarnings(as.numeric(cl[[tcol]])),
                  status = suppressWarnings(as.integer(cl[[scol]])),
                  row.names = samples)
dat$CS <- CS[samples]
if (!is.null(mpcds)) dat$CMPIS <- as.numeric(scale(mpcds[samples]))
## age is optional
acol <- intersect(c("age","Age","age_at_diagnosis"), names(cl))[1]
if (!is.na(acol)) {
  dat$age <- suppressWarnings(as.numeric(cl[[acol]]))
  if (all(is.na(dat$age))) dat$age <- NULL
  else dat$age[is.na(dat$age)] <- median(dat$age, na.rm = TRUE)
}
sig_names <- names(sig)
for (s in sig_names) dat[[s]] <- as.numeric(scale(sig[[s]][samples]))
## keep only rows with valid survival
dat <- dat[!is.na(dat$time) & dat$time > 0 & !is.na(dat$status), , drop = FALSE]
msg("Modelling on %d samples with valid survival", nrow(dat))

pcd_vars <- c("CS", if (!is.null(mpcds)) "CMPIS")
covars   <- if ("age" %in% names(dat)) "age" else character(0)
comparators <- sig_names[sapply(sig_names, function(s) sum(!is.na(dat[[s]])) > 20)]

cidx <- function(lp, d) {
  if (have_survcomp) return(survcomp::concordance.index(
    x = lp, surv.time = d$time, surv.event = d$status, method = "noether")$c.index)
  as.numeric(summary(coxph(Surv(time,status)~lp, cbind(d, lp=lp)))$concordance[1])
}

## ---------------------------------------------------------------------
## 4. Univariable Cox + C-index
## ---------------------------------------------------------------------
uni <- do.call(rbind, lapply(c(pcd_vars, comparators), function(v) {
  f <- as.formula(sprintf("Surv(time,status) ~ %s", v))
  cx <- tryCatch(coxph(f, dat), error = function(x) NULL); if (is.null(cx)) return(NULL)
  sm <- summary(cx); lp <- predict(cx, type = "lp")
  data.frame(variable = v, HR = round(sm$coef[1,2],3),
             lower = round(sm$conf.int[1,3],3), upper = round(sm$conf.int[1,4],3),
             p = signif(sm$coef[1,5],3), C_index = round(cidx(lp, dat),3),
             stringsAsFactors = FALSE)
}))
utils::write.csv(uni, file.path(results_dir,"Supplementary_Table_3_benchmark_univariable.csv"), row.names = FALSE)
print(uni)

## ---------------------------------------------------------------------
## 5. Multivariable Cox (does CS / MPCDS stay independent?)
## ---------------------------------------------------------------------
multi_tbl <- list()
fit_multi <- function(pcd) {
  terms <- c(pcd, comparators, covars)
  f <- as.formula(sprintf("Surv(time,status) ~ %s", paste(terms, collapse = " + ")))
  cx <- coxph(f, dat); sm <- summary(cx)
  data.frame(model = pcd, term = rownames(sm$coef),
             HR = round(sm$coef[,2],3),
             lower = round(sm$conf.int[,3],3), upper = round(sm$conf.int[,4],3),
             p = signif(sm$coef[,5],3), row.names = NULL, stringsAsFactors = FALSE)
}
for (pcd in pcd_vars) multi_tbl[[pcd]] <- tryCatch(fit_multi(pcd), error = function(x) NULL)
multi <- do.call(rbind, multi_tbl)
utils::write.csv(multi, file.path(results_dir,"Supplementary_Table_3_benchmark_multivariable.csv"), row.names = FALSE)
print(multi)

## ---------------------------------------------------------------------
## 6. Incremental value: nested-model LR test + C-index comparison
## ---------------------------------------------------------------------
base_terms <- c(comparators, covars)
base_f <- as.formula(sprintf("Surv(time,status) ~ %s", paste(base_terms, collapse = " + ")))
base_cx <- coxph(base_f, dat, x = TRUE)
inc_rows <- list()
for (pcd in pcd_vars) {
  full_f <- as.formula(sprintf("Surv(time,status) ~ %s", paste(c(base_terms, pcd), collapse = " + ")))
  full_cx <- coxph(full_f, dat, x = TRUE)
  lr <- anova(base_cx, full_cx, test = "LRT")
  lr_p <- lr$`Pr(>|Chi|)`[2]; if (is.null(lr_p)) lr_p <- lr[2, ncol(lr)]
  c_base <- cidx(predict(base_cx, type="lp"), dat)
  c_full <- cidx(predict(full_cx, type="lp"), dat)
  comp_p <- NA_real_
  if (have_survcomp) comp_p <- tryCatch(survcomp::cindex.comp(
      survcomp::concordance.index(predict(full_cx,type="lp"), dat$time, dat$status, method="noether"),
      survcomp::concordance.index(predict(base_cx,type="lp"), dat$time, dat$status, method="noether")
    )$p.value, error = function(x) NA_real_)
  inc_rows[[pcd]] <- data.frame(
    added = pcd, LR_p = signif(lr_p,3),
    C_base = round(c_base,3), C_full = round(c_full,3), dC = round(c_full-c_base,3),
    Cindex_comp_p = signif(comp_p,3), stringsAsFactors = FALSE)
}
inc <- do.call(rbind, inc_rows)
utils::write.csv(inc, file.path(results_dir,"Supplementary_Table_3_incremental_value.csv"), row.names = FALSE)
print(inc)

## correlation of PCD vars with comparators (association)
cor_mat <- sapply(comparators, function(s)
  sapply(pcd_vars, function(p) {
    x <- if (p == "CS") as.numeric(dat$CS) else dat[[p]]
    suppressWarnings(cor(x, dat[[s]], method = "spearman", use = "pairwise.complete.obs"))
  }))
cor_mat <- rbind(cor_mat); rownames(cor_mat) <- pcd_vars

## ---------------------------------------------------------------------
## 7. Decision-curve analysis (gated) — net benefit at 5 years
## ---------------------------------------------------------------------
dca_ok <- FALSE
if (have_dcurves) {
  dca_res <- tryCatch({
    horizon <- 1825
    d2 <- dat; d2$surv_obj <- with(d2, Surv(time, status))
    mk <- function(cx) 1 - summary(survfit(cx, newdata = d2), times = horizon)$surv[1, ]
    d2$base   <- mk(base_cx)
    d2$plusCS <- mk(coxph(as.formula(sprintf("Surv(time,status)~%s", paste(c(base_terms,"CS"),collapse="+"))), dat))
    dcurves::dca(surv_obj ~ base + plusCS, data = d2, time = horizon)
  }, error = function(x) { msg("DCA failed: %s", conditionMessage(x)); NULL })
  dca_ok <- !is.null(dca_res)
} else msg("NOTE: 'dcurves' not installed — decision-curve panel skipped.")

## ---------------------------------------------------------------------
## 8. Figure
## ---------------------------------------------------------------------
figf <- file.path(results_dir, "Figure_S_benchmark.pdf")
pdf(figf, width = 11, height = 9)
op <- par(mfrow = c(2,2), mar = c(7,5,3,1))
## A: univariable C-index
o <- order(uni$C_index)
barplot(uni$C_index[o], names.arg = uni$variable[o], las = 2, horiz = FALSE,
        ylim = c(0.4, max(0.7, max(uni$C_index)+.03)), xpd = FALSE, col = "#1F4E79",
        ylab = "C-index", main = "A. Univariable prognostic C-index")
abline(h = 0.5, lty = 2)
## B: multivariable HR (CS model) forest
cm <- multi[multi$model == "CS", ]
if (nrow(cm)) {
  yy <- seq_len(nrow(cm)); par(mar = c(4,8,3,1))
  plot(cm$HR, yy, xlim = range(c(cm$lower, cm$upper, 1), na.rm=TRUE), pch = 19,
       yaxt = "n", ylab = "", xlab = "HR (95% CI)", main = "B. Multivariable Cox (CS model)")
  segments(cm$lower, yy, cm$upper, yy); abline(v = 1, lty = 2)
  axis(2, yy, cm$term, las = 2, cex.axis = .8)
}
## C: incremental C-index
par(mar = c(7,5,3,1))
if (exists("inc")) {
  vals <- c(inc$C_base[1], inc$C_full); nm <- c("base", paste0("+", inc$added))
  bc <- barplot(vals, names.arg = nm, las = 2, ylim = c(0.4, max(vals)+.03), xpd = FALSE,
                col = c("grey70", rep("#C55A11", nrow(inc))), ylab = "C-index",
                main = "C. Incremental prognostic value")
  text(bc, vals, round(vals,3), pos = 3, cex = .8)
}
## D: correlation heatmap PCD vs comparators
par(mar = c(8,6,3,2))
image(seq_len(ncol(cor_mat)), seq_len(nrow(cor_mat)), t(cor_mat[nrow(cor_mat):1,,drop=FALSE]),
      axes = FALSE, zlim = c(-1,1), col = hcl.colors(21,"Blue-Red 3"),
      xlab = "", ylab = "", main = "D. Spearman corr: PCD vs signatures")
axis(1, seq_len(ncol(cor_mat)), colnames(cor_mat), las = 2, cex.axis = .8)
axis(2, seq_len(nrow(cor_mat)), rev(rownames(cor_mat)), las = 2, cex.axis = .8)
par(op); dev.off()
msg("Wrote %s", figf)

## ---------------------------------------------------------------------
## 9. Declared stage output
## ---------------------------------------------------------------------
out <- list(univariable = uni, multivariable = multi, incremental = inc,
            correlation = cor_mat, comparators = comparators, pcd_vars = pcd_vars,
            n_samples = length(samples), dca = dca_ok,
            figure = figf, generated = Sys.time())
saveRDS(out, file.path(results_dir, "classifier_benchmark.rds"))
msg("Stage 7 complete. Outputs in %s/", results_dir)
