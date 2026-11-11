## =====================================================================
## Stage 11 — CMPIS clinical panels: survival, clinical models, immune
##            correlates and drug sensitivity
## ---------------------------------------------------------------------
## Headless version of the CMPIS analyses of Validation_PCDRGS.R.
## Score = R/cmpis_weights.R. Groups = median split of the CMPIS
## (High = higher risk).
##
##  S3  Kaplan-Meier, CMPIS High vs Low: TCGA + GSE19234, GSE98394
##      (selection cohorts) + GSE65904, GSE53118 (independent)
##  S4  TCGA: time-dependent ROC (1/3/5/10 y), uni/multivariable Cox
##      (CMPIS + age + sex + stage), nomogram + calibration (rms)
##  S5  TCGA immune correlates by group: TMB, antigen-presentation score
##      immune-evasion proxy (IOBR exhaustion/exclusion ssGSEA),
##      real TIDE if results/tide*.csv exists, CS subtype overlap
##  S6  Drug sensitivity: GDSC2 AUC (DepMap) predicted into TCGA by
##      ComBat-homogenised kernel ridge (LOO-tuned), High vs Low
##
## Tables: Supplementary_Table_8..11 (results/)
## Run: source("config.R"); source("R/11_CMPIS_Clinical.R")
## =====================================================================

suppressMessages(suppressWarnings({ library(survival) }))
set.seed(if (exists("cfg") && !is.null(cfg$seed)) cfg$seed else 20240101L)
results_dir <- "results"; if (!dir.exists(results_dir)) dir.create(results_dir, recursive = TRUE)
msg  <- function(...) message(sprintf(...))
has  <- function(p) requireNamespace(p, quietly = TRUE)
fmtp <- function(p) ifelse(is.na(p), "NA", ifelse(p < 1e-4, formatC(p, format = "e", digits = 1), formatC(p, format = "g", digits = 2)))
source("R/cmpis_weights.R")
cmpis <- load_cmpis(); if (is.null(cmpis)) stop("CMPIS weights not found — run Stage 10 first.")
msg("CMPIS weights in use: %s", cmpis$version)
out <- list(version = cmpis$version)

## ---------------------------------------------------------------------
## 0. TCGA data
## ---------------------------------------------------------------------
stopifnot(file.exists("skcm.tcga.RData")); load("skcm.tcga.RData")
expr <- as.matrix(skcm.tcga$mRNA.expr); storage.mode(expr) <- "double"
if (max(expr, na.rm = TRUE) > 50) expr <- log2(expr + 1)
expr <- expr[, !duplicated(colnames(expr)), drop = FALSE]
clin <- as.data.frame(skcm.tcga$clin.info)
msg("clin.info columns: %s", paste(names(clin), collapse = ", "))
pick <- function(cands) { x <- intersect(cands, names(clin)); if (length(x)) x[1] else NA }
tcol <- pick(c("futime","time","OS.time","OS_time")); scol <- pick(c("fustat","status","OS","event"))
acol <- pick(c("Age","age","age_at_diagnosis","age_at_initial_pathologic_diagnosis"))
gcol <- pick(c("Gender","gender","sex","Sex"))
stcol <- pick(c("Stage","stage","ajcc_pathologic_stage","pathologic_stage","tumor_stage","ajcc_pathologic_tumor_stage"))
if (is.na(tcol) || is.na(scol)) stop("survival columns not found in clin.info")
sm <- intersect(colnames(expr), rownames(clin))
cmp <- score_cmpis(expr[, sm], cmpis)
d <- data.frame(row.names = sm, time = as.numeric(clin[sm, tcol]), status = as.integer(clin[sm, scol]),
                CMPIS = cmp[sm])
if (!is.na(acol)) d$Age <- suppressWarnings(as.numeric(clin[sm, acol]))
if (!is.na(acol) && median(d$Age, na.rm = TRUE) > 200) d$Age <- d$Age / 365.25     # days -> years
if (!is.na(gcol)) d$Gender <- factor(tolower(as.character(clin[sm, gcol])))
if (!is.na(stcol)) {                                          # collapse to I / II / III / IV
  s <- toupper(as.character(clin[sm, stcol])); s <- sub("^STAGE\\s*", "", s)
  s <- ifelse(grepl("^IV", s), "IV", ifelse(grepl("^III", s), "III", ifelse(grepl("^II", s), "II",
       ifelse(grepl("^I", s) | grepl("^0", s), "I", NA))))
  d$Stage <- factor(s, levels = c("I","II","III","IV"))
}
d <- d[!is.na(d$time) & d$time > 0 & !is.na(d$status) & !is.na(d$CMPIS), , drop = FALSE]
d$years <- d$time / 365.25
d$Group <- factor(ifelse(d$CMPIS > median(d$CMPIS), "High", "Low"), levels = c("Low","High"))
msg("TCGA: %d samples with survival (%d events); age=%s sex=%s stage=%s", nrow(d), sum(d$status),
    !is.null(d$Age), !is.null(d$Gender), !is.null(d$Stage))

## ---------------------------------------------------------------------
## 1. S3 — Kaplan-Meier by CMPIS group (TCGA + GEO)
## ---------------------------------------------------------------------
cohorts <- list(`TCGA-SKCM (training)` = d[, c("years","status","CMPIS")])
geo <- c(GSE19234 = "selection", GSE98394 = "selection", GSE65904 = "independent", GSE53118 = "independent")
for (g in names(geo)) {
  f <- sprintf("%s_cohort.txt", tolower(g)); if (!file.exists(f)) { msg("missing %s", f); next }
  x <- read.delim(f, check.names = FALSE); x <- x[!is.na(x$time) & x$time > 0 & !is.na(x$status), ]
  genes <- setdiff(names(x), c("time","status"))
  sc <- score_cmpis(t(as.matrix(x[, genes])), cmpis)
  unit_div <- if (max(x$time) < 400) 12 else 365.25            # months vs days
  cohorts[[sprintf("%s (%s)", g, geo[[g]])]] <- data.frame(years = x$time / unit_div, status = x$status, CMPIS = sc)
}
km_tab <- do.call(rbind, lapply(names(cohorts), function(nm) {
  z <- cohorts[[nm]]; z$Group <- factor(ifelse(z$CMPIS > median(z$CMPIS), "High", "Low"), c("Low","High"))
  lr <- survdiff(Surv(years, status) ~ Group, z); p <- 1 - pchisq(lr$chisq, 1)
  cg <- summary(coxph(Surv(years, status) ~ Group, z))
  cc <- summary(coxph(Surv(years, status) ~ scale(CMPIS), z))
  mm <- tryCatch(summary(survfit(Surv(years, status) ~ Group, z))$table[, "median"], error = function(e) c(NA, NA))
  data.frame(cohort = nm, n = nrow(z), events = sum(z$status), logrank_p = signif(p, 3),
             median_OS_Low_y = round(unname(mm[1]), 2), median_OS_High_y = round(unname(mm[2]), 2),
             HR_high_vs_low = round(cg$conf.int[1,1], 3), HR_lower = round(cg$conf.int[1,3], 3),
             HR_upper = round(cg$conf.int[1,4], 3), HR_per_SD = round(cc$conf.int[1,1], 3),
             HR_per_SD_p = signif(cc$coefficients[1,5], 3), C_index = round(cc$concordance[1], 3))
}))
print(km_tab); out$km <- km_tab
pdf(file.path(results_dir, "Figure_S3_CMPIS_KM.pdf"), width = 12, height = 7.5)
op <- par(mfrow = c(2, 3), mar = c(4.2, 4.2, 3, 1))
for (nm in names(cohorts)) {
  z <- cohorts[[nm]]; z$Group <- factor(ifelse(z$CMPIS > median(z$CMPIS), "High", "Low"), c("Low","High"))
  fit <- survfit(Surv(years, status) ~ Group, z); r <- km_tab[km_tab$cohort == nm, ]
  plot(fit, col = c("#1F4E79","#B4232A"), lwd = 2, mark.time = TRUE, xlab = "Years", ylab = "Overall survival",
       main = nm, cex.main = 1)
  legend("topright", c(sprintf("CMPIS-Low (n=%d)", sum(z$Group == "Low")), sprintf("CMPIS-High (n=%d)", sum(z$Group == "High"))),
         col = c("#1F4E79","#B4232A"), lwd = 2, bty = "n", cex = 0.85)
  legend("bottomleft", c(sprintf("log-rank p = %s", fmtp(r$logrank_p)), sprintf("HR (High vs Low) = %.2f", r$HR_high_vs_low),
         sprintf("C-index = %.3f", r$C_index)), bty = "n", cex = 0.85)
}
par(op); dev.off(); msg("Wrote Figure_S3_CMPIS_KM.pdf")

## ---------------------------------------------------------------------
## 2. S4 — timeROC, Cox, nomogram, calibration (TCGA)
## ---------------------------------------------------------------------
tp <- c(1, 3, 5, 10)
roc_tab <- NULL
if (has("timeROC")) tryCatch({
  mk <- list(CMPIS = d$CMPIS)
  if (!is.null(d$Age))   mk$Age   <- d$Age
  if (!is.null(d$Stage)) mk$Stage <- as.numeric(d$Stage)
  roc_tab <- do.call(rbind, lapply(names(mk), function(k) {
    ok <- !is.na(mk[[k]])
    r <- timeROC::timeROC(T = d$years[ok], delta = d$status[ok], marker = mk[[k]][ok], cause = 1,
                          times = tp, iid = TRUE, weighting = "marginal")
    ci <- tryCatch(confint(r)$CI_AUC, error = function(e) matrix(NA, length(tp), 2))
    data.frame(marker = k, years = tp, AUC = round(r$AUC, 3), lower = round(ci[,1]/100, 3), upper = round(ci[,2]/100, 3))
  }))
  print(roc_tab)
}, error = function(e) msg("timeROC failed: %s", conditionMessage(e)))
out$timeROC <- roc_tab

cov <- intersect(c("Age","Gender","Stage"), names(d))
dm <- droplevels(d[complete.cases(d[, c("CMPIS", cov)]), , drop = FALSE]); dm$CMPIS_SD <- as.numeric(scale(dm$CMPIS))
uni <- do.call(rbind, lapply(c("CMPIS_SD", cov), function(v) {
  s <- summary(coxph(as.formula(paste("Surv(years, status) ~", v)), dm))
  data.frame(model = "univariable", term = rownames(s$conf.int), HR = round(s$conf.int[,1], 3),
             lower = round(s$conf.int[,3], 3), upper = round(s$conf.int[,4], 3), p = signif(s$coefficients[,5], 3))
}))
mfit <- coxph(as.formula(paste("Surv(years, status) ~", paste(c("CMPIS_SD", cov), collapse = " + "))), dm)
sm_ <- summary(mfit)
mv <- data.frame(model = "multivariable", term = rownames(sm_$conf.int), HR = round(sm_$conf.int[,1], 3),
                 lower = round(sm_$conf.int[,3], 3), upper = round(sm_$conf.int[,4], 3), p = signif(sm_$coefficients[,5], 3))
cox_tab <- rbind(uni, mv); rownames(cox_tab) <- NULL
base_fit <- if (length(cov)) coxph(as.formula(paste("Surv(years, status) ~", paste(cov, collapse = " + "))), dm) else NULL
lr_p <- if (!is.null(base_fit)) signif(anova(base_fit, mfit)[2, "Pr(>|Chi|)"], 3) else NA
cidx <- c(clinical = if (!is.null(base_fit)) round(base_fit$concordance[["concordance"]], 3) else NA,
          clinical_plus_CMPIS = round(mfit$concordance[["concordance"]], 3))
msg("Multivariable n=%d; LR p (add CMPIS to %s) = %s; C %s -> %s", nrow(dm), paste(cov, collapse = "+"),
    fmtp(lr_p), cidx[1], cidx[2])
print(cox_tab); out$cox <- cox_tab; out$lr_p <- lr_p; out$cindex <- cidx

pdf(file.path(results_dir, "Figure_S4_CMPIS_clinical.pdf"), width = 12, height = 9)
op <- par(mfrow = c(1, 2), mar = c(6.5, 4.5, 3, 1))
## (a) time-dependent AUC
if (!is.null(roc_tab)) {
  cols <- c(CMPIS = "#B4232A", Age = "#777777", Stage = "#1F4E79")
  plot(NA, xlim = range(tp), ylim = c(0.3, 0.9), xlab = "Years", ylab = "Time-dependent AUC",
       main = "(a) Time-dependent ROC (TCGA)", xaxt = "n"); axis(1, tp); abline(h = 0.5, lty = 3)
  for (k in unique(roc_tab$marker)) { r <- roc_tab[roc_tab$marker == k, ]
    lines(r$years, r$AUC, type = "b", pch = 19, lwd = 2, col = cols[k])
    arrows(r$years, r$lower, r$years, r$upper, angle = 90, code = 3, length = .03, col = cols[k]) }
  legend("topleft", names(cols)[names(cols) %in% roc_tab$marker], col = cols[names(cols) %in% roc_tab$marker],
         lwd = 2, pch = 19, bty = "n")
} else plot.new()
mtext("Markers oriented a priori (higher = higher risk); AUC < 0.5 = marker runs the opposite way in this cohort (age).",
      side = 1, line = 5, cex = 0.6)
## (b) forest plot (multivariable)
par(mar = c(4.5, 8.5, 3, 8))
fp <- mv; k <- nrow(fp)
plot(NA, xlim = range(log(c(fp$lower, fp$upper, 0.5, 2)), na.rm = TRUE), ylim = c(0.5, k + 0.5), yaxt = "n",
     xaxt = "n", xlab = "Hazard ratio (log scale)", ylab = "", main = sprintf("(b) Multivariable Cox (n = %d)", nrow(dm)))
at <- c(0.25, 0.5, 1, 2, 4); axis(1, log(at), at); abline(v = 0, lty = 2)
segments(log(fp$lower), k:1, log(fp$upper), k:1, lwd = 2); points(log(fp$HR), k:1, pch = 15, cex = 1.3)
axis(2, k:1, sub("CMPIS_SD", "CMPIS (per SD)", fp$term), las = 1, cex.axis = 0.8, tick = FALSE)
text(par("usr")[2], k:1, sprintf("%.2f  p=%s", fp$HR, fmtp(fp$p)), adj = -0.05, cex = 0.75, xpd = TRUE)
## (c) nomogram + (d) calibration
nomo_ok <- FALSE
if (has("rms")) tryCatch({
  suppressMessages(library(rms))
  dn <- dm; dn$CMPIS <- dn$CMPIS_SD; assign("dd", rms::datadist(dn), envir = .GlobalEnv); options(datadist = "dd")
  f <- rms::cph(as.formula(paste("Surv(years, status) ~", paste(c("CMPIS", cov), collapse = " + "))),
                data = dn, x = TRUE, y = TRUE, surv = TRUE, time.inc = 3)
  S <- rms::Survival(f)
  nom <- rms::nomogram(f, fun = list(function(x) S(1, x), function(x) S(3, x), function(x) S(5, x)),
                       funlabel = c("1-y OS", "3-y OS", "5-y OS"), fun.at = c(.9, .8, .6, .4, .2), lp = FALSE)
  par(op); op <- par(mfrow = c(1, 1), mar = c(2, 2, 3, 1))
  plot(nom, cex.axis = 0.7, cex.var = 0.8); title("(c) Nomogram (CMPIS per SD)", line = 1)
  par(mfrow = c(1, 1), mar = c(4.5, 4.5, 3, 1))
  cols <- c("#1F4E79","#B4232A","#2E7D32"); first <- TRUE
  for (i in seq_along(c(1, 3, 5))) { u <- c(1, 3, 5)[i]
    fu <- rms::cph(formula(f), data = dn, x = TRUE, y = TRUE, surv = TRUE, time.inc = u)
    cal <- rms::calibrate(fu, u = u, cmethod = "KM", method = "boot", B = 200, m = floor(nrow(dn)/4))
    pr <- cal[, "mean.predicted"]; ob <- cal[, "KM"]
    if (first) { plot(pr, ob, xlim = c(0,1), ylim = c(0,1), type = "b", pch = 19, col = cols[i], lwd = 2,
                      xlab = "Nomogram-predicted OS", ylab = "Observed OS (KM)", main = "(d) Calibration (bootstrap, B = 200)")
                 abline(0, 1, lty = 2, col = "grey50"); first <- FALSE
    } else lines(pr, ob, type = "b", pch = 19, col = cols[i], lwd = 2)
  }
  legend("topleft", c("1 year","3 years","5 years"), col = cols, lwd = 2, pch = 19, bty = "n")
  nomo_ok <- TRUE
}, error = function(e) msg("nomogram/calibration failed: %s", conditionMessage(e)))
par(op); dev.off(); msg("Wrote Figure_S4_CMPIS_clinical.pdf")

## ---------------------------------------------------------------------
## 3. S5 — immune / genomic correlates by CMPIS group (TCGA)
## ---------------------------------------------------------------------
bio <- data.frame(row.names = rownames(d), Group = d$Group, CMPIS = d$CMPIS)
id12 <- function(x) substr(x, 1, 12)
## TMB (mutations / 30 Mb, as in the original script) from skcm.tcga$maf if present
mafn <- grep("maf", names(skcm.tcga), ignore.case = TRUE, value = TRUE)
if (length(mafn)) tryCatch({
  m <- skcm.tcga[[mafn[1]]]; if (methods::is(m, "MAF")) m <- m@data
  m <- as.data.frame(m); tb <- table(id12(m$Tumor_Sample_Barcode)) / 30
  bio$TMB <- as.numeric(tb[id12(rownames(bio))])
}, error = function(e) msg("TMB failed: %s", conditionMessage(e)))
if (is.null(bio$TMB)) msg("NOTE: no MAF in skcm.tcga — TMB skipped.")
## ssGSEA scores
ss <- function(gs) {
  gs <- lapply(gs, intersect, rownames(expr)); gs <- gs[lengths(gs) >= 3]
  if (!length(gs) || !has("GSVA")) return(NULL)
  sc <- GSVA::gsva(GSVA::ssgseaParam(exprData = expr[, rownames(d)], geneSets = gs)); sc
}
ap <- tryCatch(ss(list(AP = c("HLA-A","HLA-B","HLA-C","HLA-DRA","HLA-DRB1","HLA-DQA1","HLA-DQB1",
                              "TAP1","TAP2","PSMB8","PSMB9","NLRC5","CIITA"))), error = function(e) NULL)
if (!is.null(ap)) bio$Antigen_presentation <- as.numeric(ap[1, rownames(bio)])
if (has("IOBR")) tryCatch({
  e <- new.env()
  for (k in c("sig_group","signature_collection")) try(utils::data(list = k, package = "IOBR", envir = e), silent = TRUE)
  if (is.null(e$sig_group)) e$sig_group <- getExportedValue("IOBR", "sig_group")
  if (is.null(e$signature_collection)) e$signature_collection <- getExportedValue("IOBR", "signature_collection")
  sc_all <- e$signature_collection
  ex <- ss(sc_all[intersect(e$sig_group$immu_exhaustion, names(sc_all))])
  xc <- ss(sc_all[intersect(e$sig_group$immu_exclusion,  names(sc_all))])
  z <- scale(t(rbind(ex, xc)))
  bio$Immune_evasion_proxy <- as.numeric(rowMeans(z, na.rm = TRUE)[rownames(bio)])
}, error = function(e) msg("IOBR evasion proxy failed: %s", conditionMessage(e)))
tide_file <- list.files(results_dir, pattern = "^tide.*\\.csv$", full.names = TRUE, ignore.case = TRUE)
if (length(tide_file)) tryCatch({
  td <- read.csv(tide_file[1], check.names = FALSE); idc <- names(td)[1]
  tc <- intersect(c("TIDE","Dysfunction","Exclusion"), names(td))
  for (k in tc) bio[[paste0("TIDE_", k)]] <- as.numeric(td[[k]][match(id12(rownames(bio)), id12(td[[idc]]))])
}, error = function(e) msg("TIDE file unreadable: %s", conditionMessage(e)))
## CS overlap
if (file.exists("TCGA_Clustering.RData")) tryCatch({
  e <- new.env(); load("TCGA_Clustering.RData", envir = e); cr <- get("cmoic.skcm", envir = e)$clust.res
  CS <- setNames(ifelse(cr$clust == 1, "CS1", "CS2"), cr$samID); bio$CS <- CS[rownames(bio)]
  tb <- table(bio$Group, bio$CS); print(tb); out$cs_overlap <- tb
  msg("CMPIS-High vs CS2: Fisher p = %s", fmtp(fisher.test(tb)$p.value))
}, error = function(e) msg("CS overlap failed: %s", conditionMessage(e)))
vars <- setdiff(names(bio), c("Group","CMPIS","CS"))
bio_tab <- do.call(rbind, lapply(vars, function(v) {
  x <- bio[[v]]; ok <- !is.na(x)
  data.frame(variable = v, n = sum(ok), median_Low = signif(median(x[ok & bio$Group == "Low"]), 3),
             median_High = signif(median(x[ok & bio$Group == "High"]), 3),
             wilcox_p = signif(wilcox.test(x[ok] ~ bio$Group[ok])$p.value, 3),
             spearman_rho_vs_CMPIS = round(cor(x[ok], bio$CMPIS[ok], method = "spearman"), 3))
}))
print(bio_tab); out$biomarkers <- bio_tab
pdf(file.path(results_dir, "Figure_S5_CMPIS_immune.pdf"), width = 3.2 * max(1, length(vars)), height = 4.6)
op <- par(mfrow = c(1, max(1, length(vars))), mar = c(3, 4.2, 3.5, 1))
for (v in vars) { r <- bio_tab[bio_tab$variable == v, ]
  boxplot(bio[[v]] ~ bio$Group, col = c("#EAF1F8","#FBEBEC"), border = c("#1F4E79","#B4232A"), outline = FALSE,
          ylab = gsub("_", " ", v), xlab = "", names = c("CMPIS-Low","CMPIS-High"),
          main = sprintf("%s\nWilcoxon p = %s", gsub("_", " ", v), fmtp(r$wilcox_p)), cex.main = 0.9)
  stripchart(bio[[v]] ~ bio$Group, vertical = TRUE, method = "jitter", add = TRUE, pch = 16, cex = 0.35,
             col = adjustcolor(c("#1F4E79","#B4232A"), 0.4)) }
par(op); dev.off(); msg("Wrote Figure_S5_CMPIS_immune.pdf")

## ---------------------------------------------------------------------
## 4. S6 — drug sensitivity (GDSC2 AUC; lower AUC = more sensitive)
## ---------------------------------------------------------------------
drug_tab <- NULL
ef <- "Expression_Public_25Q3_subsetted.csv"; df <- "Drug_sensitivity_AUC_(Sanger_GDSC2)_subsetted.csv"
panel <- c("BEZ235","PICTILISIB","GSK-2141795","ALPELISIB","TASELISIB","AFURESERTIB","IPATASERTIB",
           "VORINOSTAT","ENTINOSTAT","(+)-JQ-1","PACLITAXEL","DOCETAXEL","CIS-DDP","FLUOROURACIL",
           "DABRAFENIB","TRAMETINIB","PALBOCICLIB","RIBOCICLIB","NILOTINIB","DASATINIB","BORTEZOMIB")
if (file.exists(ef) && file.exists(df) && has("data.table")) tryCatch({
  msg("Drug model: reading DepMap expression + GDSC2 AUC ...")
  E  <- data.table::fread(ef, data.table = FALSE, check.names = FALSE)
  rownames(E) <- E[[1]]; E <- as.matrix(E[, -(1:7)]); storage.mode(E) <- "double"      # cell lines x genes
  A  <- data.table::fread(df, data.table = FALSE, check.names = FALSE)
  rownames(A) <- A[[1]]; A <- A[, -(1:7), drop = FALSE]
  base <- toupper(sub(" \\(GDSC2:[0-9]+\\)$", "", names(A)))
  keep <- which(base %in% panel)
  keep <- keep[order(-colSums(!is.na(A[, keep, drop = FALSE])))]; keep <- keep[!duplicated(base[keep])]  # best-covered duplicate
  cl <- intersect(rownames(E), rownames(A)); E <- E[cl, ]; A <- A[cl, keep, drop = FALSE]
  genes <- intersect(colnames(E), rownames(expr))
  tr <- t(E[, genes]); te <- expr[genes, rownames(d)]
  keepg <- apply(tr, 1, sd) > 0 & apply(te, 1, sd) > 0; tr <- tr[keepg, ]; te <- te[keepg, ]
  msg("Drug model: %d cell lines, %d TCGA samples, %d shared genes, %d drugs", ncol(tr), ncol(te), nrow(tr), ncol(A))
  comb <- cbind(tr, te); batch <- rep(c("cells","tcga"), c(ncol(tr), ncol(te)))
  if (has("sva")) comb <- suppressMessages(sva::ComBat(comb, batch = batch)) else {
    msg("sva missing — per-batch gene centring instead of ComBat")
    comb[, batch == "cells"] <- comb[, batch == "cells"] - rowMeans(comb[, batch == "cells"])
    comb[, batch == "tcga"]  <- comb[, batch == "tcga"]  - rowMeans(comb[, batch == "tcga"]) }
  v <- apply(comb, 1, var); comb <- comb[v >= quantile(v, 0.2), ]          # drop 20% lowest-variance genes (oncoPredict default)
  X <- t((comb - rowMeans(comb)) / apply(comb, 1, sd))
  Xtr <- X[batch == "cells", ]; Xte <- X[batch == "tcga", ]
  Kall <- tcrossprod(Xtr) / ncol(Xtr); Kte_all <- tcrossprod(Xte, Xtr) / ncol(Xtr)
  lam <- 10^seq(-3, 3, length.out = 31)
  pred <- matrix(NA, nrow(Xte), ncol(A), dimnames = list(rownames(Xte), base[keep]))
  qual <- data.frame(drug = base[keep], n_cells = NA, lambda = NA, LOO_r = NA, cell_sd = NA)
  for (j in seq_len(ncol(A))) {
    y <- as.numeric(A[, j]); ok <- !is.na(y); yc <- y[ok] - mean(y[ok])
    ei <- eigen(Kall[ok, ok], symmetric = TRUE); U <- ei$vectors; dv <- pmax(ei$values, 0); Uy <- crossprod(U, yc)
    loo <- sapply(lam, function(l) { h <- dv / (dv + l); yh <- U %*% (h * Uy)
      Hii <- rowSums(U^2 * rep(h, each = nrow(U))); mean(((yc - yh) / (1 - Hii))^2) })
    l <- lam[which.min(loo)]; h <- dv / (dv + l); yh <- U %*% (h * Uy); Hii <- rowSums(U^2 * rep(h, each = nrow(U)))
    loo_pred <- yc - (yc - yh) / (1 - Hii)
    alpha <- U %*% (Uy / (dv + l))
    pred[, j] <- as.numeric(Kte_all[, ok] %*% alpha) + mean(y[ok])
    qual[j, 2:5] <- c(sum(ok), l, round(cor(loo_pred, yc), 3), sd(y[ok]))
  }
  ## Three groupings, so we can see WHERE any drug signal comes from:
  ##   the CMPIS, the alternative ensemble weighting (sensitivity), and CS2 vs CS1
  cell_sd <- setNames(qual$cell_sd, qual$drug)
  gr <- list(CMPIS_corrected = d[rownames(pred), "Group"] == "High")
  org <- tryCatch(score_cmpis(expr[, rownames(pred)], load_cmpis("original")), error = function(e) NULL)
  if (!is.null(org)) gr$CMPIS_original <- org[rownames(pred)] > median(org)
  if (!is.null(bio$CS)) { cs <- bio[rownames(pred), "CS"]; gr$CS2_vs_CS1 <- ifelse(is.na(cs), NA, cs == "CS2") }
  one <- function(isHigh, tag) {
    okk <- !is.na(isHigh)
    do.call(rbind, lapply(colnames(pred), function(k) {
      x <- pred[okk, k]; g <- isHigh[okk]
      data.frame(drug = k, grouping = tag,
                 median_AUC_Low = signif(median(x[!g]), 4), median_AUC_High = signif(median(x[g]), 4),
                 delta_SD = round((median(x[g]) - median(x[!g])) / cell_sd[[k]], 3),
                 more_sensitive = ifelse(median(x[g]) < median(x[!g]), "High", "Low"),
                 wilcox_p = signif(wilcox.test(x ~ g)$p.value, 3))
    }))
  }
  drug_long <- do.call(rbind, lapply(names(gr), function(n) one(gr[[n]], n)))
  drug_long <- merge(drug_long, qual, by = "drug")
  ## a drug model is only interpretable if it predicts held-out cell lines (LOO r >= 0.3)
  ## and its ridge penalty is not pinned at the top of the grid (= intercept-only model)
  drug_long$reliable <- drug_long$LOO_r >= 0.3 & drug_long$lambda < max(lam) * 0.99
  drug_long$FDR <- NA_real_
  for (n in unique(drug_long$grouping)) {
    i <- drug_long$grouping == n & drug_long$reliable
    drug_long$FDR[i] <- signif(p.adjust(drug_long$wilcox_p[i], "BH"), 3)
  }
  drug_long <- drug_long[order(drug_long$grouping, drug_long$wilcox_p), ]; rownames(drug_long) <- NULL
  drug_tab <- drug_long; out$drugs <- drug_long
  print(drug_long[drug_long$reliable, c("grouping","drug","delta_SD","more_sensitive","wilcox_p","FDR","LOO_r")])
  msg("Drug models NOT predictive (LOO r < 0.3 or lambda at bound; excluded from inference): %s",
      paste(unique(drug_long$drug[!drug_long$reliable]), collapse = ", "))
  ## figure: page 1 = corrected-CMPIS boxplots (unreliable models flagged); page 2 = effect sizes by grouping
  pdf(file.path(results_dir, "Figure_S6_CMPIS_drugs.pdf"), width = 14, height = 8)
  dd1 <- drug_long[drug_long$grouping == "CMPIS_corrected", ]
  gf  <- factor(gr$CMPIS_corrected, c(FALSE, TRUE), c("Low", "High"))
  op <- par(mfrow = c(ceiling(nrow(dd1) / 7), 7), mar = c(2.5, 3.5, 3.8, 0.5))
  for (k in dd1$drug) {
    r <- dd1[dd1$drug == k, ]
    boxplot(pred[, k] ~ gf, col = c("#EAF1F8", "#FBEBEC"), border = c("#1F4E79", "#B4232A"),
            outline = FALSE, ylab = "", xlab = "",
            main = sprintf("%s\n%s | LOO r %.2f | %.2f SD", k,
                           if (r$reliable) paste("FDR", fmtp(r$FDR)) else "MODEL NOT PREDICTIVE", r$LOO_r, r$delta_SD),
            cex.main = 0.7, col.main = if (r$reliable) "black" else "#B4232A")
  }
  par(op)
  rel <- unique(drug_long$drug[drug_long$reliable])
  if (length(rel) > 0) {
    M <- sapply(names(gr), function(n) { z <- drug_long[drug_long$grouping == n, ]; setNames(z$delta_SD, z$drug)[rel] })
    M <- matrix(M, nrow = length(rel), dimnames = list(rel, names(gr)))
    M <- M[order(M[, 1]), , drop = FALSE]
    cl3 <- c("#B4232A", "#999999", "#1F4E79")[seq_len(ncol(M))]
    op <- par(mar = c(4.5, 9, 3, 1))
    barplot(t(M), beside = TRUE, horiz = TRUE, las = 1, col = cl3,
            xlab = "Predicted AUC shift, High - Low (in SD of cell-line AUC; negative = High more sensitive)",
            main = "Predicted drug-sensitivity shift (drugs with a predictive model only)")
    abline(v = 0); legend("bottomright", colnames(M), fill = cl3, bty = "n", cex = 0.85)
    par(op)
  }
  dev.off(); msg("Wrote Figure_S6_CMPIS_drugs.pdf (predicted GDSC2 AUC; lower = more sensitive)")
}, error = function(e) msg("Drug sensitivity failed: %s", conditionMessage(e))) else
  msg("NOTE: DepMap/GDSC2 files or data.table missing — drug panel skipped.")

## ---------------------------------------------------------------------
## 5. tables
## ---------------------------------------------------------------------
w <- function(x, f) if (!is.null(x)) { utils::write.csv(x, file.path(results_dir, f), row.names = FALSE); msg("Wrote %s", f) }
w(km_tab,   "Supplementary_Table_8_CMPIS_KM.csv")
w(rbind(cox_tab, data.frame(model = "LR test", term = "add CMPIS to clinical", HR = NA, lower = cidx[1], upper = cidx[2], p = lr_p)),
  "Supplementary_Table_9_CMPIS_Cox.csv")
w(roc_tab,  "Supplementary_Table_9b_CMPIS_timeROC.csv")
w(bio_tab,  "Supplementary_Table_10_CMPIS_immune.csv")
w(drug_tab, "Supplementary_Table_11_CMPIS_drugs.csv")
saveRDS(out, file.path(results_dir, "cmpis_clinical.rds"))
msg("Stage 11 complete.")
