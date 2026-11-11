## =====================================================================
## Stage 12 — Cross-check of the drug-sensitivity panel with the
##            oncoPredict / pRRophetic ridge pipeline (Fig S7)
## ---------------------------------------------------------------------
## Stage 11 predicts GDSC2 AUC with a kernel-ridge implementation. This
## stage re-predicts a 12-drug panel with the SAME recipe as the original
## analysis (calcPhenotype_FIXED taken verbatim from Validation_PCDRGS.R:
## ComBat-homogenise -> drop 20% low-variance genes -> ridge::linearRidge),
## then compares High-vs-Low for the three groupings (CMPIS, alternative
## ensemble weighting, CS2 vs CS1) and checks agreement with Stage 11.
##
## Output: results/Supplementary_Table_12_drug_oncoPredict.csv,
##         results/Figure_S7_oncoPredict_check.pdf, results/drug_oncopredict.rds
## Run:    source("config.R"); source("R/12_Drug_oncoPredict_Check.R")   (~15-30 min)
## =====================================================================

suppressMessages(suppressWarnings(library(survival)))
set.seed(if (exists("cfg") && !is.null(cfg$seed)) cfg$seed else 20240101L)
results_dir <- "results"; if (!dir.exists(results_dir)) dir.create(results_dir, recursive = TRUE)
msg  <- function(...) message(sprintf(...))
has  <- function(p) requireNamespace(p, quietly = TRUE)
fmtp <- function(p) ifelse(is.na(p), "NA", ifelse(p < 1e-4, formatC(p, format = "e", digits = 1), formatC(p, format = "g", digits = 2)))
source("R/cmpis_weights.R")

panel <- c("TRAMETINIB","DABRAFENIB","VORINOSTAT","BORTEZOMIB","BEZ235","PICTILISIB",
           "GSK-2141795","(+)-JQ-1","PACLITAXEL","DOCETAXEL","DASATINIB","IPATASERTIB")

## ---------------------------------------------------------------------
## 0. the original calcPhenotype_FIXED, extracted verbatim from the old script
## ---------------------------------------------------------------------
if (!has("ridge")) {
  msg("Installing 'ridge' (required by the original ridge pipeline) ...")
  try(install.packages("ridge", repos = "https://cloud.r-project.org"), silent = TRUE)
}
if (!has("ridge")) stop("Package 'ridge' is not available. Run install.packages('ridge') and re-run.")
getfun <- function(nm) {
  for (pk in c("oncoPredict", "pRRophetic")) if (has(pk)) {
    f <- tryCatch(utils::getFromNamespace(nm, pk), error = function(e) NULL); if (!is.null(f)) return(f)
  }
  NULL
}
for (nm in c("homogenizeData", "doVariableSelection")) {
  f <- getfun(nm); if (is.null(f)) stop(sprintf("'%s' not found in oncoPredict/pRRophetic - install oncoPredict.", nm))
  assign(nm, f, envir = .GlobalEnv)
}
src <- readLines("Validation_PCDRGS.R", warn = FALSE)
a <- grep("^calcPhenotype_FIXED <- function", src)[1]
if (is.na(a)) stop("calcPhenotype_FIXED not found in Validation_PCDRGS.R")
b <- a - 1 + which(trimws(src[a:length(src)], "right") == "}")[1]
eval(parse(text = src[a:b]), envir = .GlobalEnv)
msg("Loaded calcPhenotype_FIXED (Validation_PCDRGS.R lines %d-%d)", a, b)

## ---------------------------------------------------------------------
## 1. TCGA tumours + the three groupings
## ---------------------------------------------------------------------
stopifnot(file.exists("skcm.tcga.RData")); load("skcm.tcga.RData")
expr <- as.matrix(skcm.tcga$mRNA.expr); storage.mode(expr) <- "double"
if (max(expr, na.rm = TRUE) > 50) expr <- log2(expr + 1)
expr <- expr[, !duplicated(colnames(expr)), drop = FALSE]
clin <- as.data.frame(skcm.tcga$clin.info)
tcol <- intersect(c("futime","time","OS.time"), names(clin))[1]; scol <- intersect(c("fustat","status","OS"), names(clin))[1]
sm <- intersect(colnames(expr), rownames(clin))
ok <- !is.na(clin[sm, tcol]) & clin[sm, tcol] > 0 & !is.na(clin[sm, scol]); sm <- sm[ok]     # same 303 as Stage 11
cm_cor <- score_cmpis(expr[, sm], load_cmpis("corrected"))
cm_org <- score_cmpis(expr[, sm], load_cmpis("original"))
gr <- list(CMPIS_corrected = cm_cor[sm] > median(cm_cor))
if (!is.null(cm_org)) gr$CMPIS_original <- cm_org[sm] > median(cm_org)
if (file.exists("TCGA_Clustering.RData")) tryCatch({
  e <- new.env(); load("TCGA_Clustering.RData", envir = e); cr <- get("cmoic.skcm", envir = e)$clust.res
  cs <- setNames(ifelse(cr$clust == 1, "CS1", "CS2"), cr$samID)[sm]; gr$CS2_vs_CS1 <- ifelse(is.na(cs), NA, cs == "CS2")
}, error = function(e) msg("CS labels unavailable: %s", conditionMessage(e)))
msg("Groupings: %s (n = %d tumours)", paste(names(gr), collapse = ", "), length(sm))

## ---------------------------------------------------------------------
## 2. GDSC2 training data (identical selection to Stage 11)
## ---------------------------------------------------------------------
ef <- "Expression_Public_25Q3_subsetted.csv"; df <- "Drug_sensitivity_AUC_(Sanger_GDSC2)_subsetted.csv"
stopifnot(file.exists(ef), file.exists(df), has("data.table"))
E <- data.table::fread(ef, data.table = FALSE, check.names = FALSE); rownames(E) <- E[[1]]
E <- as.matrix(E[, -(1:7)]); storage.mode(E) <- "double"
A <- data.table::fread(df, data.table = FALSE, check.names = FALSE); rownames(A) <- A[[1]]; A <- A[, -(1:7), drop = FALSE]
base <- toupper(sub(" \\(GDSC2:[0-9]+\\)$", "", names(A)))
keep <- which(base %in% panel); keep <- keep[order(-colSums(!is.na(A[, keep, drop = FALSE])))]; keep <- keep[!duplicated(base[keep])]
cl <- intersect(rownames(E), rownames(A)); E <- E[cl, ]; A <- A[cl, keep, drop = FALSE]; colnames(A) <- base[keep]
genes <- intersect(colnames(E), rownames(expr))
tr <- t(E[, genes]); te <- expr[genes, sm]
kg <- apply(tr, 1, sd) > 0 & apply(te, 1, sd) > 0; tr <- tr[kg, ]; te <- te[kg, ]
msg("Training: %d cell lines x %d genes; test: %d tumours; drugs: %s", ncol(tr), nrow(tr), ncol(te), paste(colnames(A), collapse = ", "))

## ---------------------------------------------------------------------
## 3. predict each drug with the original pipeline
## ---------------------------------------------------------------------
pred <- matrix(NA_real_, ncol(te), ncol(A), dimnames = list(colnames(te), colnames(A)))
for (k in colnames(A)) {
  t0 <- Sys.time(); y <- A[[k]]; okk <- !is.na(y)
  pr <- tryCatch(calcPhenotype_FIXED(trainingExprData = tr[, okk, drop = FALSE], trainingPtype = y[okk],
                   testExprData = te, batchCorrect = "eb", powerTransformPhenotype = FALSE,
                   removeLowVaryingGenes = 0.2, minNumSamples = 10, selection = 1,
                   printOutput = FALSE, removeLowVaringGenesFrom = "homogenizeData"),
                 error = function(e) { msg("  %s FAILED: %s", k, conditionMessage(e)); NULL })
  if (!is.null(pr) && length(pr) == ncol(te)) pred[, k] <- as.numeric(pr)
  msg("  %-12s done in %.1f min (%d cell lines)", k, as.numeric(difftime(Sys.time(), t0, units = "mins")), sum(okk))
}
got <- colnames(pred)[colSums(!is.na(pred)) > 0]
if (!length(got)) stop("No drug could be predicted with the original pipeline - see messages above.")
cell_sd <- apply(A, 2, sd, na.rm = TRUE)

## ---------------------------------------------------------------------
## 4. High vs Low for each grouping
## ---------------------------------------------------------------------
one <- function(isHigh, tag) {
  okg <- !is.na(isHigh)
  do.call(rbind, lapply(got, function(k) {
    x <- pred[okg, k]; g <- isHigh[okg]
    data.frame(drug = k, grouping = tag, delta_SD_op = round((median(x[g]) - median(x[!g])) / cell_sd[[k]], 3),
               more_sensitive_op = ifelse(median(x[g]) < median(x[!g]), "High", "Low"),
               wilcox_p_op = signif(wilcox.test(x ~ g)$p.value, 3))
  }))
}
op_tab <- do.call(rbind, lapply(names(gr), function(n) one(gr[[n]], n)))
op_tab$FDR_op <- signif(ave(op_tab$wilcox_p_op, op_tab$grouping, FUN = function(p) p.adjust(p, "BH")), 3)

## ---------------------------------------------------------------------
## 5. agreement with the Stage 11 (kernel-ridge) implementation
## ---------------------------------------------------------------------
f11 <- file.path(results_dir, "Supplementary_Table_11_CMPIS_drugs.csv")
cmp <- op_tab
if (file.exists(f11)) {
  k11 <- read.csv(f11, check.names = FALSE)
  k11 <- k11[, c("drug", "grouping", "delta_SD", "more_sensitive", "wilcox_p", "FDR", "LOO_r", "reliable")]
  names(k11)[3:6] <- c("delta_SD_kr", "more_sensitive_kr", "wilcox_p_kr", "FDR_kr")
  cmp <- merge(op_tab, k11, by = c("drug", "grouping"), all.x = TRUE)
  cmp$same_direction <- cmp$more_sensitive_op == cmp$more_sensitive_kr
  rho <- suppressWarnings(cor(cmp$delta_SD_op, cmp$delta_SD_kr, method = "spearman", use = "complete.obs"))
  msg("Agreement with Stage 11: Spearman rho of effect sizes = %.2f; same direction in %d/%d drug-grouping pairs",
      rho, sum(cmp$same_direction, na.rm = TRUE), sum(!is.na(cmp$same_direction)))
}
cmp <- cmp[order(cmp$grouping, cmp$wilcox_p_op), ]; rownames(cmp) <- NULL
print(cmp[, intersect(c("grouping","drug","delta_SD_op","more_sensitive_op","wilcox_p_op","FDR_op","delta_SD_kr","more_sensitive_kr","same_direction"), names(cmp))])
utils::write.csv(cmp, file.path(results_dir, "Supplementary_Table_12_drug_oncoPredict.csv"), row.names = FALSE)

## the manuscript's specific claims, in one place
msg("--- Selected drugs (High / CS2 more sensitive?) under the ridge pipeline ---")
for (d in c("BEZ235","PICTILISIB","GSK-2141795","VORINOSTAT","(+)-JQ-1","PACLITAXEL","DOCETAXEL","DABRAFENIB")) {
  r <- cmp[cmp$drug == d & cmp$grouping == "CMPIS_corrected", ]; r2 <- cmp[cmp$drug == d & cmp$grouping == "CMPIS_original", ]
  if (nrow(r)) msg("  %-12s corrected: %s sens. p=%s | original: %s sens. p=%s", d, r$more_sensitive_op, fmtp(r$wilcox_p_op),
                   if (nrow(r2)) r2$more_sensitive_op else "-", if (nrow(r2)) fmtp(r2$wilcox_p_op) else "-")
}

## ---------------------------------------------------------------------
## 6. figure
## ---------------------------------------------------------------------
if (all(c("delta_SD_kr", "delta_SD_op") %in% names(cmp))) {
  pdf(file.path(results_dir, "Figure_S7_oncoPredict_check.pdf"), width = 7.5, height = 7)
  cols <- c(CMPIS_corrected = "#B4232A", CMPIS_original = "#999999", CS2_vs_CS1 = "#1F4E79")
  v <- cmp[!is.na(cmp$delta_SD_kr), ]; lim <- range(c(v$delta_SD_kr, v$delta_SD_op, 0)) * 1.1
  plot(v$delta_SD_kr, v$delta_SD_op, col = cols[v$grouping], pch = ifelse(isTRUE(v$reliable) | v$reliable %in% TRUE, 19, 1), cex = 1.2,
       xlim = lim, ylim = lim, xlab = "Stage 11 kernel ridge: predicted AUC shift High-Low (SD)",
       ylab = "Original oncoPredict/pRRophetic ridge: predicted AUC shift (SD)", main = "Drug-sensitivity prediction: implementation check")
  abline(0, 1, lty = 2); abline(h = 0, v = 0, col = "grey80")
  big <- abs(v$delta_SD_op) > 0.15 | abs(v$delta_SD_kr) > 0.15
  if (any(big)) text(v$delta_SD_kr[big], v$delta_SD_op[big], v$drug[big], pos = 3, cex = 0.6)
  legend("topleft", c(names(cols), "filled = passed LOO QC"), col = c(cols, "black"), pch = c(19, 19, 19, 19), bty = "n", cex = 0.8)
  dev.off(); msg("Wrote Figure_S7_oncoPredict_check.pdf")
}
saveRDS(list(table = cmp, pred = pred), file.path(results_dir, "drug_oncopredict.rds"))
msg("Stage 12 complete.")
