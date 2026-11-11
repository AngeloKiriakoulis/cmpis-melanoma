## =====================================================================
## Stage 15 — Real immunotherapy validation in independent ICB cohorts
## ---------------------------------------------------------------------
## Applies CS1/CS2 (nearest-centroid to the TCGA subtypes), CS1-likeness
## and the CMPIS to melanoma immune-checkpoint cohorts with real clinical
## response and survival, pulled from cBioPortal (iAtlas-harmonised
## expression + clinical data; public API, no dbGaP needed):
##
##   Liu 2019   (Van Allen lab, DFCI)   anti-PD-1        mel_iatlas_liu_2019
##   Riaz 2017  (CheckMate-038)          nivolumab        mel_iatlas_riaz_nivolumab_2017
##   Hugo 2016  (UCLA)                   pembrolizumab    mel_iatlas_hugo_ucla_2016
##   Gide 2019  (MIA)                    anti-PD-1 mono   mel_iatlas_gide_2019
##   Van Allen 2015 (DFCI)               ipilimumab       skcm_dfci_2015   (anti-CTLA-4, secondary)
##
## Pre-specified directions (no post-hoc flipping):
##   CMPIS     higher = higher risk  -> expected LOWER in responders
##   CS1-ness  higher = CS1-like     -> expected HIGHER in responders
##   GEP       T-cell-inflamed GEP (Ayers 2017), positive-control reference
## Per cohort: responder vs non-responder (Wilcoxon), AUC + 95% CI (DeLong),
## logistic OR per SD, CS1 vs CS2 response (Fisher OR), OS/PFS by median
## split (log-rank) and Cox HR per SD. Across anti-PD-1 cohorts:
## inverse-variance fixed + DerSimonian-Laird random-effects meta-analysis,
## and a pooled analysis with cohort-stratified logistic / Cox models.
##
## Outputs (results/):
##   Supplementary_Table_11_ICB_cohorts.csv      per cohort x predictor
##   Supplementary_Table_11b_ICB_meta.csv        meta-analysis + pooled
##   Supplementary_Table_11c_ICB_CS_response.csv CS1 vs CS2 response rates
##   ICB_cohort_scores.csv                       per-patient scores
##   Figure_S9_ICB_response.pdf / Figure_S9b_ICB_survival.pdf
##   icb_cbioportal.rds
## Cache: data/icb/cbio_<study>.rds (delete to re-download)
##
## Run:  source("config.R"); source("R/15_ICB_cBioPortal.R")
## =====================================================================

suppressWarnings(suppressMessages({
  library(survival); library(httr); library(jsonlite)
  have_proc  <- requireNamespace("pROC", quietly = TRUE)
  have_orghs <- requireNamespace("org.Hs.eg.db", quietly = TRUE)
}))
set.seed(if (exists("cfg") && !is.null(cfg$seed)) cfg$seed else 20240101L)
results_dir <- "results"; if (!dir.exists(results_dir)) dir.create(results_dir, recursive = TRUE)
icb_dir <- file.path("data", "icb"); if (!dir.exists(icb_dir)) dir.create(icb_dir, recursive = TRUE)
msg <- function(...) message(sprintf(...))
CBIO <- "https://www.cbioportal.org/api"

## ---------------------------------------------------------------------
## scores
## ---------------------------------------------------------------------
source("R/cmpis_weights.R")
cmpis <- if (exists("cmpis_override")) cmpis_override else load_cmpis()
if (is.null(cmpis)) stop("CMPIS weights not found (skcm_analysis_results/final_model_info_corrected.rds)")
msg("CMPIS weights: %s", cmpis$version)

gep_genes <- c("CCL5","CD27","CD274","CD276","CD8A","CMKLR1","CXCL9","CXCR6",
               "HLA-DQA1","HLA-DRB1","HLA-E","IDO1","LAG3","NKG7","PDCD1LG2",
               "PSMB10","STAT1","TIGIT")
logx <- function(m) { m <- as.matrix(m); storage.mode(m) <- "double"
  m[m < 0 & !is.na(m)] <- 0; if (max(m, na.rm = TRUE) > 50) m <- log2(m + 1); m }
gep_score <- function(ex) { g <- intersect(gep_genes, rownames(ex)); if (length(g) < 10) return(NULL)
  m <- logx(ex[g, , drop = FALSE]); m <- m[apply(m, 1, sd, na.rm = TRUE) > 0, , drop = FALSE]
  setNames(colMeans(t(scale(t(m))), na.rm = TRUE), colnames(ex)) }

## CS1/CS2 centroids from TCGA (same construction as Stage 8)
cs_file <- file.path(results_dir, "cs_centroids.rds")
cs_centroids <- if (exists("cs_override")) cs_override else NULL
if (is.null(cs_centroids) && file.exists("skcm.tcga.RData") && file.exists("TCGA_Clustering.RData")) {
  ev <- new.env(); load("skcm.tcga.RData", envir = ev)
  te <- as.matrix(get("skcm.tcga", envir = ev)$mRNA.expr)
  if (max(te, na.rm = TRUE) > 50) te <- log2(te + 1)
  ec <- new.env(); load("TCGA_Clustering.RData", envir = ec)
  cr <- get("cmoic.skcm", envir = ec)$clust.res
  lab <- setNames(ifelse(cr$clust == 1, "CS1", "CS2"), cr$samID)
  cs <- intersect(colnames(te), names(lab)); te <- te[, cs]; lab <- lab[cs]
  mu1 <- rowMeans(te[, lab == "CS1", drop = FALSE]); mu2 <- rowMeans(te[, lab == "CS2", drop = FALSE])
  d <- (mu1 - mu2) / (apply(te, 1, sd) + 1e-6)
  mk <- names(sort(abs(d), decreasing = TRUE))[seq_len(min(500, length(d)))]
  cs_centroids <- list(genes = mk, CS1 = mu1[mk], CS2 = mu2[mk])
  saveRDS(cs_centroids, cs_file)
  msg("CS centroids: %d TCGA samples, %d marker genes", length(cs), length(mk))
  tcga_check <- te[mk, ]
}
if (is.null(cs_centroids) && file.exists(cs_file)) cs_centroids <- readRDS(cs_file)
if (is.null(cs_centroids)) stop("CS centroids unavailable (need skcm.tcga.RData + TCGA_Clustering.RData)")
cs_score <- function(ex) {
  ## Nearest-centroid on gene-centred data (as for PAM50): expression is centred per gene within the
  ## cohort and correlated with the TCGA CS1-minus-CS2 centroid difference. This removes platform /
  ## sample-type offsets between TCGA and the ICB cohorts; > 0 = relatively CS1-like.
  g <- intersect(cs_centroids$genes, rownames(ex)); if (length(g) < 50) return(NULL)
  e <- logx(ex[g, , drop = FALSE]); e <- e - rowMeans(e, na.rm = TRUE)
  d <- cs_centroids$CS1[g] - cs_centroids$CS2[g]
  setNames(as.numeric(cor(e, d, method = "spearman", use = "pairwise.complete.obs")), colnames(ex))
}
if (exists("tcga_check")) {             # sanity check: does the gene-centred classifier recover TCGA CS1/CS2?
  pr <- ifelse(cs_score(tcga_check) > 0, "CS1", "CS2")
  msg("CS classifier on TCGA (in-sample): %.1f%% agreement with CS1/CS2; kappa-free table CS1->CS1 %d, CS2->CS2 %d",
      100 * mean(pr == lab), sum(pr == "CS1" & lab == "CS1"), sum(pr == "CS2" & lab == "CS2"))
  rm(tcga_check)
}
cs_score_abs <- function(ex) {           # absolute (uncentred) version, reported for transparency
  g <- intersect(cs_centroids$genes, rownames(ex)); if (length(g) < 50) return(NULL)
  e <- logx(ex[g, , drop = FALSE])
  c1 <- cor(e, cs_centroids$CS1[g], method = "spearman", use = "pairwise.complete.obs")
  c2 <- cor(e, cs_centroids$CS2[g], method = "spearman", use = "pairwise.complete.obs")
  setNames(as.numeric(c1 - c2), colnames(ex))
}

## ---------------------------------------------------------------------
## cBioPortal API
## ---------------------------------------------------------------------
api_get <- function(path, query = list()) {
  for (k in 1:4) {
    r <- tryCatch(GET(paste0(CBIO, path), query = query, timeout(120),
                      add_headers(Accept = "application/json")), error = function(e) NULL)
    if (!is.null(r) && status_code(r) == 200)
      return(fromJSON(content(r, "text", encoding = "UTF-8"), simplifyVector = TRUE))
    if (!is.null(r)) msg("  GET %s -> HTTP %s: %s", path, status_code(r), substr(content(r, "text", encoding = "UTF-8"), 1, 200))
    if (!is.null(r) && status_code(r) %in% c(400, 404)) break
    Sys.sleep(2 * k)
  }
  stop("cBioPortal GET failed: ", path)
}
api_post <- function(path, body, query = list()) {
  for (k in 1:4) {
    r <- tryCatch(POST(paste0(CBIO, path), query = query, timeout(300), encode = "raw",
                       body = toJSON(body, auto_unbox = TRUE),
                       add_headers(`Content-Type` = "application/json", Accept = "application/json")),
                  error = function(e) NULL)
    if (!is.null(r) && status_code(r) == 200)
      return(fromJSON(content(r, "text", encoding = "UTF-8"), simplifyVector = TRUE))
    if (!is.null(r)) msg("  POST %s -> HTTP %s: %s", path, status_code(r), substr(content(r, "text", encoding = "UTF-8"), 1, 200))
    if (!is.null(r) && status_code(r) %in% c(400, 404)) break
    Sys.sleep(2 * k)
  }
  stop("cBioPortal POST failed: ", path)
}
clin_wide <- function(study, type) {
  d <- api_get(sprintf("/studies/%s/clinical-data", study),
               list(clinicalDataType = type, projection = "SUMMARY", pageSize = "100000"))
  if (!length(d) || !nrow(d)) return(NULL)
  key <- if (type == "SAMPLE") "sampleId" else "patientId"
  ids <- unique(d[[key]]); att <- unique(d$clinicalAttributeId)
  w <- data.frame(id = ids, stringsAsFactors = FALSE)
  for (a in att) { s <- d[d$clinicalAttributeId == a, ]; w[[a]] <- s$value[match(ids, s[[key]])] }
  if (type == "SAMPLE") w$patientId <- d$patientId[match(ids, d$sampleId)]
  w
}
sym2entrez <- function(sym) {
  if (!have_orghs) stop("org.Hs.eg.db needed for gene-ID mapping")
  m <- suppressMessages(AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db, sym, "ENTREZID", "SYMBOL"))
  m[!is.na(m)]
}
fetch_expr <- function(study, profile, genes) {
  ez <- sym2entrez(genes); out <- list()
  for (ch in split(ez, ceiling(seq_along(ez) / 150))) {
    d <- api_post(sprintf("/molecular-profiles/%s/molecular-data/fetch", profile),
                  list(entrezGeneIds = I(as.integer(ch)), sampleListId = paste0(study, "_all")),
                  list(projection = "SUMMARY"))
    if (length(d) && nrow(d)) out[[length(out) + 1]] <- d[, c("entrezGeneId", "sampleId", "patientId", "value")]
  }
  d <- do.call(rbind, out); if (is.null(d)) return(NULL)
  d$gene <- names(ez)[match(as.character(d$entrezGeneId), ez)]
  smp <- unique(d$sampleId); gn <- unique(d$gene)
  m <- matrix(NA_real_, length(gn), length(smp), dimnames = list(gn, smp))
  m[cbind(match(d$gene, gn), match(d$sampleId, smp))] <- as.numeric(d$value)
  list(expr = m, sample2patient = setNames(d$patientId[match(smp, d$sampleId)], smp))
}
get_study <- function(study, profile) {
  cf <- file.path(icb_dir, paste0("cbio_", study, ".rds"))
  if (file.exists(cf)) { msg("[%s] cached", study); return(readRDS(cf)) }
  msg("[%s] downloading from cBioPortal ...", study)
  genes <- unique(c(names(cmpis$coef), gep_genes, cs_centroids$genes))
  ef <- file.path(icb_dir, paste0("cbio_", study, "_expr.rds"))
  ex <- if (file.exists(ef)) readRDS(ef) else { e <- fetch_expr(study, profile, genes); saveRDS(e, ef); e }
  x <- list(expr = ex,
            sample = clin_wide(study, "SAMPLE"), patient = clin_wide(study, "PATIENT"),
            fetched = Sys.time())
  saveRDS(x, cf); x
}

## ---------------------------------------------------------------------
## cohort definitions -> data.frame(patient, sample, responder, os_t, os_e, pfs_t, pfs_e)
## ---------------------------------------------------------------------
surv_ev <- function(s) suppressWarnings(as.integer(sub(":.*", "", s)))
num <- function(x) suppressWarnings(as.numeric(x))
lgl <- function(x) { x <- toupper(as.character(x)); ifelse(x %in% c("TRUE","YES","1"), 1L,
                     ifelse(x %in% c("FALSE","NO","0"), 0L, NA_integer_)) }
cohorts <- list(
  Liu2019   = list(study = "mel_iatlas_liu_2019", profile = "mel_iatlas_liu_2019_rna_seq_mrna",
                   label = "Liu 2019 (Van Allen lab)", drug = "anti-PD-1", antiPD1 = TRUE),
  Riaz2017  = list(study = "mel_iatlas_riaz_nivolumab_2017", profile = "mel_iatlas_riaz_nivolumab_2017_rna_seq_mrna",
                   label = "Riaz 2017", drug = "nivolumab", antiPD1 = TRUE),
  Hugo2016  = list(study = "mel_iatlas_hugo_ucla_2016", profile = "mel_iatlas_hugo_ucla_2016_rna_seq_mrna",
                   label = "Hugo 2016", drug = "pembrolizumab", antiPD1 = TRUE),
  Gide2019  = list(study = "mel_iatlas_gide_2019", profile = "mel_iatlas_gide_2019_rna_seq_mrna",
                   label = "Gide 2019", drug = "anti-PD-1 (monotherapy)", antiPD1 = TRUE),
  VanAllen2015 = list(study = "skcm_dfci_2015", profile = "skcm_dfci_2015_rna_seq_mrna",
                   label = "Van Allen 2015 (anti-CTLA-4)", drug = "ipilimumab (anti-CTLA-4)", antiPD1 = FALSE))
if (exists("cohort_subset")) cohorts <- cohorts[cohort_subset]

build_meta <- function(nm, x) {
  P <- x$patient; S <- x$sample; s2p <- x$expr$sample2patient
  smp <- colnames(x$expr$expr)
  df <- data.frame(sample = smp, patient = unname(s2p[smp]), stringsAsFactors = FALSE)
  if (!is.null(S)) {
    si <- match(df$sample, S$id)
    if ("SAMPLE_TREATMENT" %in% names(S)) df$treat <- S$SAMPLE_TREATMENT[si]
    if ("RESPONSE" %in% names(S))         df$recist_s <- S$RESPONSE[si]
  }
  pi <- match(df$patient, P$id); pc <- function(a) if (a %in% names(P)) P[[a]][pi] else rep(NA, nrow(df))
  if (nm == "VanAllen2015") {
    br <- toupper(pc("DURABLE_CLINICAL_BENEFIT"))
    df$responder <- ifelse(br %in% c("CR","PR"), 1L, ifelse(br %in% c("SD","PD"), 0L, NA_integer_))
    df$pfs_t <- num(pc("DFS_MONTHS")); df$pfs_e <- surv_ev(pc("DFS_STATUS"))
  } else {
    df$responder <- lgl(pc("RESPONDER"))
    if (all(is.na(df$responder)) && !is.null(df$recist_s))
      df$responder <- ifelse(grepl("Complete|Partial", df$recist_s), 1L,
                      ifelse(grepl("Stable|Progress", df$recist_s), 0L, NA_integer_))
    df$pfs_t <- num(pc("PFS_MONTHS")); df$pfs_e <- surv_ev(pc("PFS_STATUS"))
    df$target <- pc("ICI_TARGET"); df$rx <- pc("ICI_RX")
  }
  df$os_t <- num(pc("OS_MONTHS")); df$os_e <- surv_ev(pc("OS_STATUS"))
  ## pre-treatment only; one sample per patient
  if (!is.null(df$treat) && any(!is.na(df$treat))) df <- df[is.na(df$treat) | grepl("^pre", df$treat, ignore.case = TRUE), ]
  if (nm == "Gide2019" && !is.null(df$target)) df <- df[!is.na(df$target) & toupper(df$target) %in% c("PD1","PD-1"), ]
  df <- df[!duplicated(df$patient), ]
  df
}

## ---------------------------------------------------------------------
## statistics
## ---------------------------------------------------------------------
auc_ci <- function(score, y) {     # score oriented: higher = predicted responder
  if (have_proc) {
    r <- pROC::roc(y, score, levels = c(0, 1), direction = "<", quiet = TRUE)
    ci <- as.numeric(pROC::ci.auc(r, method = "delong"))
    return(c(auc = ci[2], lo = ci[1], hi = ci[3], se = (ci[3] - ci[1]) / (2 * 1.96)))
  }
  rk <- rank(score); n1 <- sum(y == 1); n0 <- sum(y == 0)
  a <- (sum(rk[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
  q1 <- a / (2 - a); q2 <- 2 * a^2 / (1 + a)               # Hanley-McNeil
  se <- sqrt((a * (1 - a) + (n1 - 1) * (q1 - a^2) + (n0 - 1) * (q2 - a^2)) / (n1 * n0))
  c(auc = a, lo = max(0, a - 1.96 * se), hi = min(1, a + 1.96 * se), se = se)
}
roc_pts <- function(score, y) { o <- order(score, decreasing = TRUE); l <- y[o]
  data.frame(fpr = c(0, cumsum(l == 0) / sum(y == 0)), tpr = c(0, cumsum(l == 1) / sum(y == 1))) }
coef_ci <- function(fit, i = 1) { b <- coef(fit)[i]; se <- sqrt(diag(vcov(fit)))[i]
  c(est = unname(exp(b)), lo = unname(exp(b - 1.96 * se)), hi = unname(exp(b + 1.96 * se)),
    p = unname(2 * pnorm(-abs(b / se))), logb = unname(b), se = unname(se)) }
cox_sd <- function(t, e, x) { ok <- !is.na(t) & !is.na(e) & !is.na(x) & t > 0
  if (sum(ok) < 15 || sum(e[ok]) < 5) return(NULL)
  coef_ci(coxph(Surv(t[ok], e[ok]) ~ as.numeric(scale(x[ok])))) }
km_p <- function(t, e, g) { ok <- !is.na(t) & !is.na(e) & !is.na(g) & t > 0
  if (sum(ok) < 15 || length(unique(g[ok])) < 2) return(NA_real_)
  1 - pchisq(survdiff(Surv(t[ok], e[ok]) ~ g[ok])$chisq, 1) }
meta_iv <- function(b, se) {
  ok <- is.finite(b) & is.finite(se) & se > 0; b <- b[ok]; se <- se[ok]; k <- length(b)
  if (k < 2) return(NULL)
  w <- 1 / se^2; bf <- sum(w * b) / sum(w); sef <- sqrt(1 / sum(w))
  Q <- sum(w * (b - bf)^2); tau2 <- max(0, (Q - (k - 1)) / (sum(w) - sum(w^2) / sum(w)))
  wr <- 1 / (se^2 + tau2); br <- sum(wr * b) / sum(wr); ser <- sqrt(1 / sum(wr))
  data.frame(k = k, fixed = bf, fixed_se = sef, fixed_p = 2 * pnorm(-abs(bf / sef)),
             random = br, random_se = ser, random_p = 2 * pnorm(-abs(br / ser)),
             Q = Q, Q_p = 1 - pchisq(Q, k - 1), I2 = if (Q > 0) max(0, (Q - (k - 1)) / Q) else 0, tau2 = tau2)
}

## ---------------------------------------------------------------------
## run
## ---------------------------------------------------------------------
pred_dir <- c(CMPIS = -1, CS1ness = 1, GEP = 1)   # sign that makes "higher = responder"
rows <- list(); csrows <- list(); rocs <- list(); allsc <- list()
for (nm in names(cohorts)) {
  co <- cohorts[[nm]]
  x <- tryCatch(get_study(co$study, co$profile), error = function(e) { msg("[%s] FAILED: %s", nm, conditionMessage(e)); NULL })
  if (is.null(x) || is.null(x$expr)) next
  md <- build_meta(nm, x); ex <- x$expr$expr[, md$sample, drop = FALSE]
  ex <- ex[rowSums(!is.na(ex)) > 0.8 * ncol(ex), , drop = FALSE]
  sc <- data.frame(cohort = nm, md[, c("patient","sample","responder","os_t","os_e","pfs_t","pfs_e")],
                   CMPIS = unname(score_cmpis(ex, cmpis)[md$sample]),
                   CS1ness = unname(cs_score(ex)[md$sample]),
                   GEP = unname(gep_score(ex)[md$sample]), stringsAsFactors = FALSE)
  sc$CS <- ifelse(sc$CS1ness > 0, "CS1", "CS2")
  sc$CS_absolute <- ifelse(cs_score_abs(ex)[md$sample] > 0, "CS1", "CS2")
  allsc[[nm]] <- sc
  y <- sc$responder; ny <- !is.na(y)
  msg("[%s] %s: %d pre-treatment patients, %d with response (%d R / %d NR); OS %d; PFS %d; CMPIS genes %d/8",
      nm, co$drug, nrow(sc), sum(ny), sum(y == 1, na.rm = TRUE), sum(y == 0, na.rm = TRUE),
      sum(!is.na(sc$os_t)), sum(!is.na(sc$pfs_t)), length(intersect(names(cmpis$coef), rownames(ex))))
  for (pn in names(pred_dir)) {
    s <- sc[[pn]]; ok <- ny & !is.na(s)
    r <- data.frame(cohort = nm, label = co$label, treatment = co$drug, predictor = pn,
                    n = sum(ok), responders = sum(y[ok] == 1), non_responders = sum(y[ok] == 0))
    if (sum(ok) >= 10 && length(unique(y[ok])) == 2) {
      so <- pred_dir[[pn]] * s[ok]; a <- auc_ci(so, y[ok])
      r$median_R <- median(s[ok][y[ok] == 1]); r$median_NR <- median(s[ok][y[ok] == 0])
      r$wilcox_p <- suppressWarnings(wilcox.test(s[ok] ~ y[ok])$p.value)
      r$AUC <- a[["auc"]]; r$AUC_lo <- a[["lo"]]; r$AUC_hi <- a[["hi"]]; r$AUC_se <- a[["se"]]
      o <- coef_ci(glm(y[ok] ~ as.numeric(scale(s[ok])), family = binomial), 2)
      r$OR_perSD <- o[["est"]]; r$OR_lo <- o[["lo"]]; r$OR_hi <- o[["hi"]]; r$OR_p <- o[["p"]]
      r$logOR <- o[["logb"]]; r$logOR_se <- o[["se"]]
      rocs[[paste(nm, pn)]] <- roc_pts(so, y[ok])
    }
    for (ep in c("os", "pfs")) {
      t <- sc[[paste0(ep, "_t")]]; e <- sc[[paste0(ep, "_e")]]
      h <- cox_sd(t, e, s); g <- ifelse(s >= median(s, na.rm = TRUE), "high", "low")
      r[[paste0(toupper(ep), "_n")]] <- sum(!is.na(t) & !is.na(e) & !is.na(s))
      r[[paste0(toupper(ep), "_events")]] <- sum(e[!is.na(t) & !is.na(s)], na.rm = TRUE)
      r[[paste0(toupper(ep), "_HR_perSD")]] <- if (is.null(h)) NA else h[["est"]]
      r[[paste0(toupper(ep), "_HR_lo")]] <- if (is.null(h)) NA else h[["lo"]]
      r[[paste0(toupper(ep), "_HR_hi")]] <- if (is.null(h)) NA else h[["hi"]]
      r[[paste0(toupper(ep), "_HR_p")]]  <- if (is.null(h)) NA else h[["p"]]
      r[[paste0(toupper(ep), "_logHR")]] <- if (is.null(h)) NA else h[["logb"]]
      r[[paste0(toupper(ep), "_logHR_se")]] <- if (is.null(h)) NA else h[["se"]]
      r[[paste0(toupper(ep), "_logrank_p_median")]] <- km_p(t, e, g)
    }
    rows[[paste(nm, pn)]] <- r
  }
  ## CS1 vs CS2 class
  ok <- ny & !is.na(sc$CS)
  tb <- table(factor(sc$CS[ok], c("CS1","CS2")), factor(y[ok], 1:0))
  ft <- tryCatch(fisher.test(tb), error = function(e) NULL)
  csrows[[nm]] <- data.frame(cohort = nm, label = co$label, treatment = co$drug,
    CS1_n = sum(tb["CS1", ]), CS1_resp = tb["CS1", "1"], CS1_rate = tb["CS1", "1"] / max(1, sum(tb["CS1", ])),
    CS2_n = sum(tb["CS2", ]), CS2_resp = tb["CS2", "1"], CS2_rate = tb["CS2", "1"] / max(1, sum(tb["CS2", ])),
    OR_CS1_vs_CS2 = if (is.null(ft)) NA else unname(ft$estimate),
    OR_lo = if (is.null(ft)) NA else ft$conf.int[1], OR_hi = if (is.null(ft)) NA else ft$conf.int[2],
    fisher_p = if (is.null(ft)) NA else ft$p.value,
    OS_logrank_p_CS = km_p(sc$os_t, sc$os_e, sc$CS), PFS_logrank_p_CS = km_p(sc$pfs_t, sc$pfs_e, sc$CS))
}
if (!length(rows)) stop("No cohort could be loaded — check the internet connection to www.cbioportal.org")
bind <- function(l) { cols <- unique(unlist(lapply(l, names)))
  do.call(rbind, lapply(l, function(d) { for (c in setdiff(cols, names(d))) d[[c]] <- NA; d[, cols] })) }
res <- bind(rows); csres <- bind(csrows); scores <- bind(allsc)
apd1 <- names(cohorts)[vapply(cohorts, `[[`, logical(1), "antiPD1")]

## meta-analysis + pooled (anti-PD-1 cohorts only)
mrows <- list()
for (pn in names(pred_dir)) {
  d <- res[res$predictor == pn & res$cohort %in% apd1, ]
  for (spec in list(c("OR", "logOR", "logOR_se"), c("OS_HR", "OS_logHR", "OS_logHR_se"),
                    c("PFS_HR", "PFS_logHR", "PFS_logHR_se"), c("AUC", "AUC", "AUC_se"))) {
    m <- meta_iv(d[[spec[2]]], d[[spec[3]]]); if (is.null(m)) next
    tr <- if (spec[1] == "AUC") identity else exp
    mrows[[paste(pn, spec[1])]] <- data.frame(predictor = pn, measure = paste0(spec[1], if (spec[1] != "AUC") "_perSD" else ""),
      k = m$k, fixed = tr(m$fixed), fixed_lo = tr(m$fixed - 1.96 * m$fixed_se), fixed_hi = tr(m$fixed + 1.96 * m$fixed_se),
      fixed_p = if (spec[1] == "AUC") 2 * pnorm(-abs((m$fixed - 0.5) / m$fixed_se)) else m$fixed_p,
      random = tr(m$random), random_lo = tr(m$random - 1.96 * m$random_se), random_hi = tr(m$random + 1.96 * m$random_se),
      random_p = if (spec[1] == "AUC") 2 * pnorm(-abs((m$random - 0.5) / m$random_se)) else m$random_p,
      I2 = m$I2, Q_p = m$Q_p)
  }
  ## pooled, cohort-stratified (scores z-scaled within cohort)
  p <- scores[scores$cohort %in% apd1, ]
  p$z <- ave(p[[pn]], p$cohort, FUN = function(v) as.numeric(scale(v)))
  q <- p[!is.na(p$responder) & !is.na(p$z), ]
  if (nrow(q) > 30) {
    f <- glm(responder ~ z + factor(cohort), family = binomial, data = q); o <- coef_ci(f, 2)
    mrows[[paste(pn, "pooledOR")]] <- data.frame(predictor = pn, measure = "OR_perSD_pooled_stratified",
      k = length(unique(q$cohort)), fixed = o[["est"]], fixed_lo = o[["lo"]], fixed_hi = o[["hi"]], fixed_p = o[["p"]])
  }
  for (ep in c("os", "pfs")) {
    q <- p[!is.na(p[[paste0(ep, "_t")]]) & !is.na(p[[paste0(ep, "_e")]]) & !is.na(p$z) & p[[paste0(ep, "_t")]] > 0, ]
    if (nrow(q) > 30) {
      f <- coxph(Surv(q[[paste0(ep, "_t")]], q[[paste0(ep, "_e")]]) ~ q$z + strata(q$cohort)); o <- coef_ci(f)
      mrows[[paste(pn, ep, "pooled")]] <- data.frame(predictor = pn, measure = paste0(toupper(ep), "_HR_perSD_pooled_stratified"),
        k = length(unique(q$cohort)), fixed = o[["est"]], fixed_lo = o[["lo"]], fixed_hi = o[["hi"]], fixed_p = o[["p"]])
    }
  }
}
## CMPIS adjusted for the T-cell-inflamed GEP (pooled, cohort-stratified)
p <- scores[scores$cohort %in% apd1, ]
p$zc <- ave(p$CMPIS, p$cohort, FUN = function(v) as.numeric(scale(v)))
p$zg <- ave(p$GEP, p$cohort, FUN = function(v) as.numeric(scale(v)))
for (ep in c("os", "pfs")) {
  q <- p[!is.na(p[[paste0(ep, "_t")]]) & !is.na(p[[paste0(ep, "_e")]]) & p[[paste0(ep, "_t")]] > 0 & !is.na(p$zc) & !is.na(p$zg), ]
  if (nrow(q) > 30) { f <- coxph(Surv(q[[paste0(ep, "_t")]], q[[paste0(ep, "_e")]]) ~ q$zc + q$zg + strata(q$cohort))
    for (i in 1:2) { o <- coef_ci(f, i)
      mrows[[paste("adj", ep, i)]] <- data.frame(predictor = c("CMPIS", "GEP")[i],
        measure = paste0(toupper(ep), "_HR_perSD_pooled_stratified_CMPIS+GEP"), k = length(unique(q$cohort)),
        fixed = o[["est"]], fixed_lo = o[["lo"]], fixed_hi = o[["hi"]], fixed_p = o[["p"]]) } }
}
q <- p[!is.na(p$responder) & !is.na(p$zc) & !is.na(p$zg), ]
if (nrow(q) > 30) { f <- glm(responder ~ zc + zg + factor(cohort), family = binomial, data = q)
  for (i in 2:3) { o <- coef_ci(f, i)
    mrows[[paste("adjOR", i)]] <- data.frame(predictor = c("CMPIS", "GEP")[i - 1], measure = "OR_perSD_pooled_stratified_CMPIS+GEP",
      k = length(unique(q$cohort)), fixed = o[["est"]], fixed_lo = o[["lo"]], fixed_hi = o[["hi"]], fixed_p = o[["p"]]) } }

## CS1 vs CS2 response, pooled Mantel-Haenszel
p <- scores[scores$cohort %in% apd1 & !is.na(scores$responder) & !is.na(scores$CS), ]
if (nrow(p) > 30) {
  mh <- tryCatch(mantelhaen.test(table(factor(p$CS, c("CS1","CS2")), factor(p$responder, 1:0), p$cohort), exact = FALSE),
                 error = function(e) NULL)
  if (!is.null(mh)) mrows[["CS MH"]] <- data.frame(predictor = "CS1_vs_CS2", measure = "response_OR_MantelHaenszel",
    k = length(unique(p$cohort)), fixed = unname(mh$estimate), fixed_lo = mh$conf.int[1], fixed_hi = mh$conf.int[2], fixed_p = mh$p.value)
  for (ep in c("os", "pfs")) {
    q <- p[!is.na(p[[paste0(ep, "_t")]]) & !is.na(p[[paste0(ep, "_e")]]) & p[[paste0(ep, "_t")]] > 0, ]
    if (nrow(q) > 30) { f <- coxph(Surv(q[[paste0(ep, "_t")]], q[[paste0(ep, "_e")]]) ~ I(q$CS == "CS2") + strata(q$cohort)); o <- coef_ci(f)
      mrows[[paste("CS", ep)]] <- data.frame(predictor = "CS2_vs_CS1", measure = paste0(toupper(ep), "_HR_pooled_stratified"),
        k = length(unique(q$cohort)), fixed = o[["est"]], fixed_lo = o[["lo"]], fixed_hi = o[["hi"]], fixed_p = o[["p"]]) }
  }
}
meta <- bind(mrows)

rnd <- function(d) { for (c in names(d)) if (is.numeric(d[[c]])) d[[c]] <- signif(d[[c]], 3); d }
write.csv(rnd(res),   file.path(results_dir, "Supplementary_Table_11_ICB_cohorts.csv"), row.names = FALSE)
write.csv(rnd(meta),  file.path(results_dir, "Supplementary_Table_11b_ICB_meta.csv"), row.names = FALSE)
write.csv(rnd(csres), file.path(results_dir, "Supplementary_Table_11c_ICB_CS_response.csv"), row.names = FALSE)
write.csv(scores,     file.path(results_dir, "ICB_cohort_scores.csv"), row.names = FALSE)

## ---------------------------------------------------------------------
## figures
## ---------------------------------------------------------------------
pcol <- c(CMPIS = "#C55A11", CS1ness = "#1F4E79", GEP = "#7F7F7F")
plab <- c(CMPIS = "CMPIS (low = responder)", CS1ness = "CS1-likeness", GEP = "T-cell GEP (reference)")
cl <- vapply(cohorts, `[[`, "", "label"); cn <- names(allsc); K <- length(cn)
forest <- function(d, est, lo, hi, xlab, main, mt = NULL) {
  lab <- c(cl[d$cohort], if (!is.null(mt)) "Random-effects (anti-PD-1)")
  e <- c(d[[est]], mt[1]); l <- c(d[[lo]], mt[2]); h <- c(d[[hi]], mt[3]); n <- length(e); yy <- n:1
  ok <- is.finite(e)
  xr <- range(c(l[ok], h[ok], 1), na.rm = TRUE); xr <- c(max(0.05, xr[1]), min(20, xr[2]))
  par(mar = c(4, 11, 3, 1))
  plot(NA, xlim = xr, ylim = c(0.5, n + 0.5), log = "x", yaxt = "n", xlab = xlab, ylab = "", main = main, cex.main = 1)
  abline(v = 1, lty = 2, col = "grey50")
  segments(pmax(l, xr[1]), yy, pmin(h, xr[2]), yy, lwd = 2, col = ifelse(seq_len(n) == n & !is.null(mt), "black", "grey30"))
  points(e, yy, pch = ifelse(seq_len(n) == n & !is.null(mt), 18, 15), cex = ifelse(seq_len(n) == n & !is.null(mt), 2, 1.3))
  axis(2, at = yy, labels = lab, las = 1, cex.axis = .8)
}
mt_of <- function(pn, m) { r <- meta[meta$predictor == pn & meta$measure == m, ]
  if (!nrow(r)) NULL else c(r$random, r$random_lo, r$random_hi) }

pdf(file.path(results_dir, "Figure_S9_ICB_response.pdf"), width = 13, height = 11)
layout(matrix(c(seq_len(K), rep(0, 5 - K), K + seq_len(K), rep(0, 5 - K), 2*K + c(1, 1, 2, 2, 3)), 3, 5, byrow = TRUE))
for (nm in cn) {                                    # A: CMPIS by response
  sc <- allsc[[nm]]; ok <- !is.na(sc$responder)
  par(mar = c(4, 4, 3, 1))
  r <- res[res$cohort == nm & res$predictor == "CMPIS", ]
  boxplot(CMPIS ~ factor(responder, 0:1, c("NR", "R")), data = sc[ok, ], col = c("#F4B183", "#9DC3E6"),
          outline = FALSE, ylab = "CMPIS (z)", xlab = "", main = sprintf("%s\nWilcoxon p = %s", cl[nm], signif(r$wilcox_p, 2)), cex.main = .9)
  stripchart(CMPIS ~ factor(responder, 0:1), data = sc[ok, ], vertical = TRUE, add = TRUE, method = "jitter", pch = 16, cex = .6, col = "#00000088")
}
for (nm in cn) {                                    # B: ROC
  par(mar = c(4, 4, 3, 1))
  plot(c(0, 1), c(0, 1), type = "n", xlab = "1 - specificity", ylab = "Sensitivity", main = cl[nm], cex.main = .9); abline(0, 1, lty = 3)
  lg <- c()
  for (pn in names(pred_dir)) { rp <- rocs[[paste(nm, pn)]]; if (is.null(rp)) next
    lines(rp$fpr, rp$tpr, col = pcol[pn], lwd = 2)
    a <- res[res$cohort == nm & res$predictor == pn, ]; lg <- c(lg, setNames(sprintf("%s %.2f", pn, a$AUC), pn)) }
  legend("bottomright", lg, col = pcol[names(lg)], lwd = 2, bty = "n", cex = .75, title = "AUC")
}
d <- res[res$predictor == "CMPIS", ]
forest(d, "OR_perSD", "OR_lo", "OR_hi", "OR for response per SD CMPIS", "C. CMPIS: odds of response", mt_of("CMPIS", "OR_perSD"))
d <- res[res$predictor == "CS1ness", ]
forest(d, "OR_perSD", "OR_lo", "OR_hi", "OR for response per SD CS1-likeness", "D. CS1-likeness: odds of response", mt_of("CS1ness", "OR_perSD"))
par(mar = c(6, 4, 3, 1))                             # E: CS1 vs CS2 response rate
M <- rbind(CS1 = csres$CS1_rate, CS2 = csres$CS2_rate) * 100; colnames(M) <- sub(" \\(.*", "", csres$label)
bp <- barplot(M, beside = TRUE, col = c("#1F4E79", "#C55A11"), ylim = c(0, 100), las = 2, cex.names = .75,
              ylab = "Response rate (%)", main = "E. Response by subtype")
text(colMeans(bp), apply(M, 2, max, na.rm = TRUE) + 6, sprintf("p=%s", signif(csres$fisher_p, 2)), cex = .7)
legend("topright", c("CS1", "CS2"), fill = c("#1F4E79", "#C55A11"), bty = "n", cex = .8)
dev.off()

pdf(file.path(results_dir, "Figure_S9b_ICB_survival.pdf"), width = 13, height = 9)
layout(matrix(c(seq_len(K), rep(0, 5 - K), K + c(1, 1, 2, 2, 3)), 2, 5, byrow = TRUE))
for (nm in cn) {
  sc <- allsc[[nm]]; useP <- all(is.na(sc$os_t)); t <- if (useP) sc$pfs_t else sc$os_t; e <- if (useP) sc$pfs_e else sc$os_e
  ok <- !is.na(t) & !is.na(e) & !is.na(sc$CMPIS) & t > 0
  par(mar = c(4, 4, 3, 1))
  if (sum(ok) < 15) { plot.new(); title(sprintf("%s\n(no survival data)", cl[nm]), cex.main = .9); next }
  g <- factor(ifelse(sc$CMPIS[ok] >= median(sc$CMPIS[ok]), "High", "Low"), c("Low", "High"))
  f <- survfit(Surv(t[ok], e[ok]) ~ g)
  plot(f, col = c("#1F4E79", "#C55A11"), lwd = 2, mark.time = TRUE, xlab = "Months", ylab = if (useP) "PFS" else "Overall survival",
       main = sprintf("%s\nlog-rank p = %s", cl[nm], signif(km_p(t[ok], e[ok], g), 2)), cex.main = .9)
  legend("topright", c("CMPIS Low", "CMPIS High"), col = c("#1F4E79", "#C55A11"), lwd = 2, bty = "n", cex = .8)
}
d <- res[res$predictor == "CMPIS" & is.finite(res$OS_HR_perSD), ]
if (nrow(d)) forest(d, "OS_HR_perSD", "OS_HR_lo", "OS_HR_hi", "HR for death per SD CMPIS", "F. CMPIS and overall survival", mt_of("CMPIS", "OS_HR_perSD"))
d <- res[res$predictor == "CS1ness" & is.finite(res$OS_HR_perSD), ]
if (nrow(d)) forest(d, "OS_HR_perSD", "OS_HR_lo", "OS_HR_hi", "HR for death per SD CS1-likeness", "G. CS1-likeness and overall survival", mt_of("CS1ness", "OS_HR_perSD"))
par(mar = c(1, 1, 3, 1)); plot.new(); title("Notes", cex.main = .9)
text(0, .9, paste("KM: CMPIS split at the cohort median.", "Forest: Cox HR per SD (95% CI);",
  "diamond = random-effects pooled", "estimate over anti-PD-1 cohorts.", "Van Allen 2015 = anti-CTLA-4,", "shown but not pooled.", sep = "\n"),
  adj = c(0, 1), cex = .85)
dev.off()

saveRDS(list(results = res, meta = meta, cs = csres, scores = scores, roc = rocs, cohorts = cohorts,
             cmpis = cmpis$version, generated = Sys.time()), file.path(results_dir, "icb_cbioportal.rds"))

options(width = 200)
cat("\n=== Per cohort (AUC oriented a priori; OR/HR per SD of the native score) ===\n")
print(rnd(res[, c("cohort","predictor","n","responders","AUC","AUC_lo","AUC_hi","wilcox_p","OR_perSD","OR_p",
                  "OS_n","OS_HR_perSD","OS_HR_p","OS_logrank_p_median","PFS_n","PFS_HR_perSD","PFS_HR_p")]), row.names = FALSE)
cat("\n=== CS1 vs CS2 (gene-centred nearest centroid) ===\n"); print(rnd(csres[, -2]), row.names = FALSE)
cat("\nAbsolute (uncentred) assignment, CS1 / CS2 counts per cohort:\n"); print(table(scores$cohort, scores$CS_absolute))
cat("\n=== Meta-analysis / pooled (anti-PD-1 cohorts) ===\n"); print(rnd(meta), row.names = FALSE)
msg("Stage 15 complete. Outputs in %s/", results_dir)
