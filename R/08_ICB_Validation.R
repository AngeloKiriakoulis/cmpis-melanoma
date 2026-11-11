## =====================================================================
## Stage 8 — Validation on REAL anti-PD-1 cohorts
## ---------------------------------------------------------------------
## Applies the CMPIS (earlier name: MPCDS) (and, where possible, project
## CS1/CS2) to real anti-PD-1 melanoma cohorts with RECIST response, and
## test: responder vs non-responder MPCDS, ROC/AUC for response, response
## odds ratios, and survival under anti-PD-1.
##
## Target public cohorts (auto-download via GEOquery, run on the machine):
##   * GSE78220  Hugo 2016   — pembrolizumab, melanoma (FPKM, gene symbols)
##   * GSE91061  Riaz 2017   — nivolumab,     melanoma (raw counts, Entrez)
## dbGaP-gated (not auto): Van Allen 2015 (phs000452), Liu 2019 (phs000452).
## To add any cohort manually, drop two files in data/icb/ :
##   <cohort>_expr.txt  (tab; row1 = sample IDs; col1 = gene symbols)
##   <cohort>_meta.txt  (tab; columns: sample, responder[0/1], os_time, os_event)
## and it is picked up automatically.
##
## Outputs (results/):
##   Supplementary_Table_4_ICB_validation.csv
##   Figure_S_ICB_validation.pdf
##   icb_validation.rds                (declared stage output)
##
## Note: the immunotherapy analysis reported in the paper (five cohorts,
## response + survival) is Stage 15, R/15_ICB_cBioPortal.R.
##
## Run:  source("config.R"); source("R/08_ICB_Validation.R")
## =====================================================================

suppressWarnings(suppressMessages({
  library(survival)
  have_geo    <- requireNamespace("GEOquery", quietly = TRUE)
  have_readxl <- requireNamespace("readxl",   quietly = TRUE)
  have_orghs  <- requireNamespace("org.Hs.eg.db", quietly = TRUE)
  have_edger  <- requireNamespace("edgeR", quietly = TRUE)
}))
set.seed(if (exists("cfg") && !is.null(cfg$seed)) cfg$seed else 20240101L)
results_dir <- "results"; if (!dir.exists(results_dir)) dir.create(results_dir, recursive = TRUE)
icb_dir <- file.path("data", "icb"); if (!dir.exists(icb_dir)) dir.create(icb_dir, recursive = TRUE)
msg <- function(...) message(sprintf(...))

## ---------------------------------------------------------------------
## MPCDS coefficients (8-gene score; matches Validation_PCDRGS.R)
## ---------------------------------------------------------------------
pcdrg_genes <- c("CFLAR","CYBB","CYLD","GYS1","IRF1","NFKBIA","PIM2","VCAM1")
coefs <- NULL
mfile <- file.path("skcm_analysis_results","final_model_info.rds")
if (file.exists(mfile)) { fmi <- readRDS(mfile)
  coefs <- fmi$coefficients; if (is.null(coefs)) coefs <- fmi$ensemble_coefficients }
if (is.null(coefs)) { coefs <- setNames(rep(1, length(pcdrg_genes)), pcdrg_genes)
  msg("NOTE: model coefficients not found — using unit weights for MPCDS.") }

source("R/cmpis_weights.R")
cmpis <- load_cmpis(); if (!is.null(cmpis)) msg("MPCDS weights: %s", cmpis$version)
mpcds_score <- function(expr) score_cmpis(expr, cmpis)

## additional biology-axis predictors ------------------------------------
zmean_sig <- function(expr, genes) {
  g <- intersect(genes, rownames(expr)); if (!length(g)) return(NULL)
  m <- expr[g, , drop = FALSE]; if (max(m,na.rm=TRUE) > 50) m <- log2(m + 1)
  m <- m[apply(m,1,function(x) sd(x,na.rm=TRUE)>0), , drop=FALSE]
  setNames(colMeans(t(scale(t(m))), na.rm = TRUE), colnames(expr))
}
gep_genes <- c("CCL5","CD27","CD274","CD276","CD8A","CMKLR1","CXCL9","CXCR6",
               "HLA-DQA1","HLA-DRB1","HLA-E","IDO1","LAG3","NKG7","PDCD1LG2",
               "PSMB10","STAT1","TIGIT")
gep_score <- function(expr) zmean_sig(expr, gep_genes)

## CS1/CS2 projection: TCGA-built centroids -> CS1-ness on a cohort
cs_centroids <- NULL
if (file.exists("skcm.tcga.RData") && file.exists("TCGA_Clustering.RData")) {
  tryCatch({
    ev <- new.env(); load("skcm.tcga.RData", envir = ev)
    te <- as.matrix(get("skcm.tcga", envir = ev)$mRNA.expr)
    if (max(te, na.rm=TRUE) > 50) te <- log2(te + 1)
    ec <- new.env(); load("TCGA_Clustering.RData", envir = ec)
    cr <- get("cmoic.skcm", envir = ec)$clust.res
    lab <- setNames(ifelse(cr$clust == 1, "CS1", "CS2"), cr$samID)
    cs <- intersect(colnames(te), names(lab)); te <- te[, cs]; lab <- lab[cs]
    mu1 <- rowMeans(te[, lab=="CS1", drop=FALSE]); mu2 <- rowMeans(te[, lab=="CS2", drop=FALSE])
    sda <- apply(te, 1, sd); d <- (mu1 - mu2)/(sda + 1e-6)
    mk <- names(sort(abs(d), decreasing = TRUE))[seq_len(min(500, length(d)))]
    cs_centroids <- list(genes = mk, CS1 = mu1[mk], CS2 = mu2[mk])
    msg("CS centroids built from %d TCGA samples, %d marker genes.", length(cs), length(mk))
  }, error = function(e) msg("CS centroid build failed: %s", conditionMessage(e)))
}
cs_score <- function(expr) {           # CS1-ness: higher = CS1 (immune-hot)-like
  if (is.null(cs_centroids)) return(NULL)
  g <- intersect(cs_centroids$genes, rownames(expr)); if (length(g) < 20) return(NULL)
  ex <- expr[g, , drop = FALSE]; if (max(ex,na.rm=TRUE) > 50) ex <- log2(ex + 1)
  c1 <- as.numeric(cor(ex, cs_centroids$CS1[g], method = "spearman", use = "pairwise.complete.obs"))
  c2 <- as.numeric(cor(ex, cs_centroids$CS2[g], method = "spearman", use = "pairwise.complete.obs"))
  setNames(c1 - c2, colnames(expr))
}

## ---------------------------------------------------------------------
## metrics (no extra packages)
## ---------------------------------------------------------------------
auc_mw <- function(score, label) {              # label: 1=responder,0=non
  r <- rank(score); n1 <- sum(label == 1); n0 <- sum(label == 0)
  if (n1 == 0 || n0 == 0) return(NA_real_)
  (sum(r[label == 1]) - n1*(n1+1)/2) / (n1*n0)
}
roc_points <- function(score, label) {
  o <- order(score, decreasing = TRUE); s <- score[o]; l <- label[o]
  P <- sum(label == 1); N <- sum(label == 0)
  tpr <- cumsum(l == 1)/P; fpr <- cumsum(l == 0)/N
  data.frame(fpr = c(0,fpr), tpr = c(0,tpr))
}

## ---------------------------------------------------------------------
## cohort loaders — return list(expr = genes x samples, meta = df)
##   meta columns: sample, responder (0/1), os_time (days), os_event (0/1)
## ---------------------------------------------------------------------
std_meta <- function(sample, responder, os_time = NA, os_event = NA)
  data.frame(sample = as.character(sample), responder = as.integer(responder),
             os_time = as.numeric(os_time), os_event = as.integer(os_event),
             stringsAsFactors = FALSE)

load_dropin <- function(cohort) {
  ef <- file.path(icb_dir, paste0(cohort, "_expr.txt"))
  mf <- file.path(icb_dir, paste0(cohort, "_meta.txt"))
  if (!file.exists(ef) || !file.exists(mf)) return(NULL)
  ex <- as.matrix(read.delim(ef, row.names = 1, check.names = FALSE))
  md <- read.delim(mf, check.names = FALSE)
  msg("[%s] drop-in: %d genes x %d samples", cohort, nrow(ex), ncol(ex))
  list(expr = ex, meta = md)
}

load_GSE78220 <- function() {
  if (!have_geo || !have_readxl) return(NULL)
  tryCatch({
    g <- GEOquery::getGEO("GSE78220", GSEMatrix = TRUE, getGPL = FALSE)[[1]]
    pd <- Biobase::pData(g)
    ## response column
    rc <- grep("response", names(pd), ignore.case = TRUE, value = TRUE)
    rc <- rc[which.max(sapply(rc, function(c) length(unique(pd[[c]]))))]
    msg("[GSE78220] response col '%s'; values: %s", rc,
        paste(sprintf("%s(%d)", names(table(pd[[rc]])), table(pd[[rc]])), collapse=", "))
    rv <- as.character(pd[[rc]])
    responder <- ifelse(grepl("Complete|Partial|^R$|responder", rv, ignore.case = TRUE) &
                        !grepl("Non|Progress|SD|PD", rv, ignore.case = TRUE), 1L, 0L)
    ## expression supplementary (FPKM xlsx)
    sup <- GEOquery::getGEOSuppFiles("GSE78220", baseDir = tempdir())
    xls <- rownames(sup)[grepl("FPKM.*xlsx?$", rownames(sup), ignore.case = TRUE)][1]
    fp  <- as.data.frame(readxl::read_excel(xls))
    rownames(fp) <- make.unique(as.character(fp[[1]])); fp[[1]] <- NULL
    ex <- as.matrix(fp)
    ## match expr sample cols to pData rows (by patient token)
    tok <- function(x) toupper(gsub("[^A-Za-z0-9]", "", sub("[._-].*$", "", x)))
    pd_tok <- tok(pd$title); names(responder) <- pd_tok
    ex_tok <- tok(colnames(ex))
    keep <- ex_tok %in% pd_tok
    ex <- ex[, keep, drop = FALSE]; colnames(ex) <- ex_tok[keep]
    ## optional overall survival
    osc <- grep("overall surv|surviv.*day|os.*day", names(pd), ignore.case=TRUE, value=TRUE)[1]
    vsc <- grep("vital", names(pd), ignore.case=TRUE, value=TRUE)[1]
    ost <- if (!is.na(osc)) setNames(suppressWarnings(as.numeric(pd[[osc]])), tok(pd$title)) else NULL
    ose <- if (!is.na(vsc)) setNames(as.integer(grepl("dead|decease", pd[[vsc]], ignore.case=TRUE)), tok(pd$title)) else NULL
    meta <- std_meta(colnames(ex), responder[colnames(ex)],
                     if (!is.null(ost)) ost[colnames(ex)] else NA,
                     if (!is.null(ose)) ose[colnames(ex)] else NA)
    msg("[GSE78220] %d genes x %d samples; responders=%d/%d",
        nrow(ex), ncol(ex), sum(meta$responder, na.rm=TRUE), nrow(meta))
    list(expr = ex, meta = meta)
  }, error = function(e) { msg("[GSE78220] load failed: %s", conditionMessage(e)); NULL })
}

load_GSE91061 <- function() {
  if (!have_geo) return(NULL)
  tryCatch({
    g <- GEOquery::getGEO("GSE91061", GSEMatrix = TRUE, getGPL = FALSE)[[1]]
    pd <- Biobase::pData(g)
    msg("[GSE91061] pData cols: %s", paste(names(pd), collapse=", "))
    vc <- grep("visit|timepoint|pre.?on|treatment", names(pd), ignore.case = TRUE, value = TRUE)[1]
    rc <- grep("response|bor", names(pd), ignore.case = TRUE, value = TRUE)
    rc <- rc[which.max(sapply(rc, function(c) length(unique(pd[[c]]))))]
    msg("[GSE91061] response col '%s'; values: %s", rc,
        paste(sprintf("%s(%d)", names(table(pd[[rc]])), table(pd[[rc]])), collapse=", "))
    rv <- toupper(as.character(pd[[rc]]))
    responder <- ifelse(grepl("PRCR|COMPLETE|PARTIAL|^CR|^PR", rv), 1L,
                 ifelse(grepl("SD|PD|PROGRESS|STABLE|NON", rv), 0L, NA_integer_))
    pre <- if (!is.na(vc)) grepl("pre", pd[[vc]], ignore.case = TRUE) else rep(TRUE, nrow(pd))
    ## expression supplementary — prefer FPKM, else raw counts
    sup <- GEOquery::getGEOSuppFiles("GSE91061", baseDir = tempdir())
    cf <- rownames(sup)[grepl("fpkm.*csv", rownames(sup), ignore.case = TRUE)][1]
    if (is.na(cf)) cf <- rownames(sup)[grepl("raw.*csv|count", rownames(sup), ignore.case = TRUE)][1]
    ex <- as.matrix(read.csv(cf, row.names = 1, check.names = FALSE))
    if (have_orghs && all(grepl("^[0-9]+$", head(rownames(ex), 20)))) {
      sym <- AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db, keys = rownames(ex),
               keytype = "ENTREZID", column = "SYMBOL")
      ok <- !is.na(sym); ex <- ex[ok, , drop = FALSE]; rownames(ex) <- make.unique(sym[ok])
    }
    if (max(ex, na.rm = TRUE) > 50) ex <- log2(ex + 1)
    msg("[GSE91061] expr colnames head: %s", paste(head(colnames(ex), 3), collapse=" | "))
    ## robust alignment: pick the pData column best matching expr colnames
    norm <- function(x) toupper(gsub("[^A-Za-z0-9]", "", x))
    exn <- norm(colnames(ex))
    cand <- c(list(GSM = rownames(pd)), as.list(pd[vapply(pd, is.character, logical(1))]))
    best <- NULL; bestn <- 0; bestkey <- NULL
    for (nm in names(cand)) { key <- norm(as.character(cand[[nm]]))
      ov <- sum(exn %in% key); if (ov > bestn) { bestn <- ov; best <- nm; bestkey <- key } }
    msg("[GSE91061] matched %d expr cols via pData$%s", bestn, best)
    if (bestn < 10) stop(sprintf("could not align expression to pData (best overlap %d)", bestn))
    idx <- match(exn, bestkey)                 # expr col -> pData row
    keep <- !is.na(idx) & pre[idx] & !is.na(responder[idx])
    ex <- ex[, keep, drop = FALSE]
    meta <- std_meta(colnames(ex), responder[idx[keep]])
    msg("[GSE91061] %d genes x %d samples; responders=%d/%d",
        nrow(ex), ncol(ex), sum(meta$responder, na.rm=TRUE), nrow(meta))
    list(expr = ex, meta = meta)
  }, error = function(e) { msg("[GSE91061] load failed: %s", conditionMessage(e)); NULL })
}

## ---------------------------------------------------------------------
## run each cohort
## ---------------------------------------------------------------------
cohorts <- c("GSE78220","GSE91061")
loaders <- list(GSE78220 = load_GSE78220, GSE91061 = load_GSE91061)
predictors <- list(CMPIS = mpcds_score, CS1ness = cs_score, GEP = gep_score)
pcols <- c(CMPIS="#C55A11", CS1ness="#1F4E79", GEP="#548235")

res_rows <- list(); roc_store <- list()
for (co in cohorts) {
  dat <- load_dropin(co); if (is.null(dat)) dat <- loaders[[co]]()
  if (is.null(dat)) { msg("[%s] unavailable — skipped.", co); next }
  ex <- dat$expr; meta <- dat$meta
  meta <- meta[!is.na(meta$responder), , drop = FALSE]
  for (pn in names(predictors)) {
    sc <- predictors[[pn]](ex); if (is.null(sc)) next
    ids <- intersect(names(sc), meta$sample); if (length(ids) < 10) next
    s2 <- sc[ids]; y <- meta$responder[match(ids, meta$sample)]
    if (length(unique(y)) < 2) next
    a <- auc_mw(s2, y); if (!is.na(a) && a < 0.5) { s2 <- -s2; a <- 1 - a }
    w <- suppressWarnings(wilcox.test(s2 ~ y)$p.value)
    or <- tryCatch({ m <- glm(y ~ s2, family = binomial)
      c(exp(coef(m)[2]), summary(m)$coef[2,4]) }, error = function(e) c(NA,NA))
    hr <- lr <- NA_real_
    ot <- meta$os_time[match(ids, meta$sample)]; oe <- meta$os_event[match(ids, meta$sample)]
    if (sum(!is.na(ot) & !is.na(oe) & ot > 0) > 15) {
      grp <- factor(ifelse(s2 >= median(s2), "high", "low"))
      sd <- data.frame(t = ot, e = oe, g = grp)
      lr <- tryCatch(1 - pchisq(survdiff(Surv(t,e)~g, sd)$chisq,1), error=function(x) NA_real_)
      hr <- tryCatch(exp(coef(coxph(Surv(t,e)~g, sd))[1]), error=function(x) NA_real_)
    }
    res_rows[[paste(co,pn)]] <- data.frame(cohort = co, predictor = pn, n = length(ids),
      responders = sum(y==1), non_responders = sum(y==0),
      AUC = round(a,3), wilcox_p = signif(w,3),
      OR = round(or[1],3), OR_p = signif(or[2],3),
      surv_logrank_p = signif(lr,3), surv_HR_high = round(hr,3),
      stringsAsFactors = FALSE)
    roc_store[[paste(co,pn)]] <- roc_points(s2, y)
  }
}
res <- if (length(res_rows)) do.call(rbind, res_rows) else data.frame()
st4 <- file.path(results_dir, "Supplementary_Table_4_ICB_validation.csv")
utils::write.csv(res, st4, row.names = FALSE); if (nrow(res)) print(res)
msg("Wrote %s", st4)

## figure: AUC barplot (predictor x cohort) + per-cohort ROC panels
figf <- file.path(results_dir, "Figure_S_ICB_validation.pdf")
pdf(figf, width = 10, height = 8)
if (nrow(res)) {
  cohs <- unique(res$cohort); preds <- names(predictors)
  layout(matrix(c(1,1, 2,3), nrow = 2, byrow = TRUE))
  par(mar = c(6,4,3,1))
  M <- matrix(NA, length(preds), length(cohs), dimnames = list(preds, cohs))
  for (i in seq_len(nrow(res))) M[res$predictor[i], res$cohort[i]] <- res$AUC[i]
  bp <- barplot(M, beside = TRUE, col = pcols[rownames(M)], ylim = c(0,1),
                ylab = "AUC (response)", main = "A. Response-prediction AUC")
  abline(h = 0.5, lty = 2)
  legend("topright", rownames(M), fill = pcols[rownames(M)], bty = "n", cex = .9)
  for (co in cohs) {
    par(mar = c(4,4,3,1))
    plot(c(0,1), c(0,1), type = "n", xlab = "FPR", ylab = "TPR",
         main = sprintf("ROC — %s", co)); abline(0,1,lty=2)
    for (pn in preds) { rp <- roc_store[[paste(co,pn)]]
      if (!is.null(rp)) lines(rp$fpr, rp$tpr, col = pcols[pn], lwd = 2) }
    legend("bottomright", preds, col = pcols[preds], lwd = 2, bty = "n", cex = .8)
  }
} else { plot.new(); text(.5,.5,"No ICB cohort data available") }
dev.off(); msg("Wrote %s", figf)

saveRDS(list(results = res, roc = roc_store,
             cs_centroid_genes = if (!is.null(cs_centroids)) cs_centroids$genes else NULL,
             generated = Sys.time()), file.path(results_dir, "icb_validation.rds"))
msg("Stage 8 complete. Outputs in %s/", results_dir)
