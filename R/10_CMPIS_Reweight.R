## =====================================================================
## Stage 10 — Final CMPIS weights and re-validation
## ---------------------------------------------------------------------
## Final weights (primary scheme):
##   w_g = stepcox_coef_g * max(imp_g, 0) / max(imp)
##   sign and effect size from the stepwise Cox model, emphasis from RSF
##   importance; higher score = higher risk, no inversion.
## Sensitivity schemes: an equal-average ensemble of the standardised Cox
## coefficient and RSF importance (labelled "original" in the code, applied
## with risk_direction = "invert"); StepCox coefficients alone; and a
## multivariable Cox model on all 8 genes.
##
## Cohorts: TCGA-SKCM (training), GSE19234 + GSE98394 (used in the
## model-selection loop), GSE65904 + GSE53118 (independent of selection).
##
## Outputs:
##   skcm_analysis_results/final_model_info_corrected.rds
##   results/Supplementary_Table_6_CMPIS_weights.csv
##   results/Supplementary_Table_7_CMPIS_revalidation.csv
##   results/Figure_S_CMPIS_reweighting.pdf
##   results/cmpis_reweight.rds         (declared stage output)
## =====================================================================

suppressMessages(library(survival))
set.seed(if (exists("cfg") && !is.null(cfg$seed)) cfg$seed else 20240101L)
results_dir <- "results"; if (!dir.exists(results_dir)) dir.create(results_dir, recursive = TRUE)
msg <- function(...) message(sprintf(...))
genes8 <- c("CFLAR","CYBB","CYLD","GYS1","IRF1","NFKBIA","PIM2","VCAM1")

## ---------------------------------------------------------------------
## 1. Training data (identical to SKCM_ML_Algorithms.R::prepare_tcga_data)
## ---------------------------------------------------------------------
ex <- as.matrix(read.table("tcga_skcm_expression.txt", header = TRUE, row.names = 1,
                           sep = "\t", check.names = FALSE))
colnames(ex) <- substr(colnames(ex), 1, 12)
cl <- read.delim("tcga_skcm_clinical.txt", check.names = FALSE)
names(cl) <- tolower(names(cl))
dead <- tolower(cl$vital_status) == "dead"
cl$time   <- ifelse(dead, as.numeric(cl$days_to_death), as.numeric(cl$days_to_last_follow_up))
cl$status <- as.integer(dead)
cl <- cl[!is.na(cl$time) & cl$time > 0 & !is.na(cl$status), ]
sm <- intersect(colnames(ex), cl$submitter_id)
train <- data.frame(time = cl$time[match(sm, cl$submitter_id)],
                    status = cl$status[match(sm, cl$submitter_id)],
                    t(ex[intersect(genes8, rownames(ex)), sm, drop = FALSE]), check.names = FALSE)
msg("TCGA training: %d samples x %d genes, %d events", nrow(train), ncol(train) - 2, sum(train$status))

## ---------------------------------------------------------------------
## 2. Raw StepCox coefficients + RSF importance
##    (a) recover from the saved ML workspace; (b) else refit identically
## ---------------------------------------------------------------------
find_final <- function(x, depth = 0) {
  if (depth > 5 || !is.list(x)) return(NULL)
  if (all(c("stepcox_coefficients", "rsf_importance_selected") %in% names(x))) return(x)
  for (el in x) { r <- tryCatch(find_final(el, depth + 1), error = function(e) NULL)
                  if (!is.null(r)) return(r) }
  NULL
}
cox_raw <- NULL; imp_raw <- NULL; source_of_raw <- NA_character_
if (file.exists("SKCM_ML_Algorithms.RData")) {
  env <- new.env()
  ok <- tryCatch({ load("SKCM_ML_Algorithms.RData", envir = env); TRUE }, error = function(e) FALSE)
  if (ok) for (nm in ls(env)) {
    ff <- find_final(get(nm, envir = env))
    if (!is.null(ff)) { cox_raw <- ff$stepcox_coefficients; imp_raw <- ff$rsf_importance_selected
                        source_of_raw <- sprintf("SKCM_ML_Algorithms.RData (%s)", nm); break }
  }
}
if (is.null(cox_raw)) {
  msg("Raw model components not found in workspace — refitting StepCox + RSF on TCGA.")
  gcols <- setdiff(names(train), c("time","status"))
  null_fit <- coxph(Surv(time, status) ~ 1, data = train)
  up <- as.formula(paste("~", paste(sprintf("`%s`", gcols), collapse = " + ")))
  sfit <- step(null_fit, scope = list(lower = ~1, upper = up), direction = "forward", trace = 0)
  cox_raw <- coef(sfit); names(cox_raw) <- gsub("`", "", names(cox_raw))
  if (requireNamespace("randomForestSRC", quietly = TRUE)) {
    rf <- randomForestSRC::rfsrc(Surv(time, status) ~ ., data = train, ntree = 1000,
                                 importance = TRUE, seed = -20240101)
    imp_raw <- rf$importance
  } else {
    msg("randomForestSRC not installed — RSF importance set to 1 (weights = StepCox only).")
    imp_raw <- setNames(rep(1, length(gcols)), gcols)
  }
  source_of_raw <- "refit on TCGA training data"
}
msg("Raw components from: %s", source_of_raw)

## ---------------------------------------------------------------------
## 3. Weight schemes
## ---------------------------------------------------------------------
orig <- readRDS(file.path("skcm_analysis_results", "final_model_info.rds"))
w_orig <- orig$coefficients; if (is.null(w_orig)) w_orig <- orig$ensemble_coefficients
fill <- function(v) { out <- setNames(rep(0, length(genes8)), genes8)
                      k <- intersect(names(v), genes8); out[k] <- v[k]; out }
cox_f <- fill(cox_raw); imp_f <- fill(imp_raw)
imp_n <- pmax(imp_f, 0); if (max(imp_n) > 0) imp_n <- imp_n / max(imp_n)
w_corr <- cox_f * imp_n                                   # PRIMARY (pre-specified)
w_step <- cox_f                                           # sensitivity 1
mv <- coxph(as.formula(paste("Surv(time, status) ~",
            paste(sprintf("`%s`", intersect(genes8, names(train))), collapse = " + "))), data = train)
w_mvc  <- fill(setNames(coef(mv), gsub("`", "", names(coef(mv)))))   # sensitivity 2

## per-gene univariate HR in training (context for the table)
uni <- sapply(genes8, function(g) if (g %in% names(train))
  exp(coef(coxph(Surv(time, status) ~ scale(train[[g]]), data = train))) else NA)

wt <- data.frame(gene = genes8,
  original_ensemble = round(fill(w_orig), 4),
  original_effective = round(-fill(w_orig), 4),          # after risk_direction = "invert"
  stepcox_coef = round(cox_f, 4), rsf_importance = round(imp_f, 4),
  corrected_weight = round(w_corr, 4),
  univariate_HR_TCGA = round(uni, 3),
  original_direction  = ifelse(-fill(w_orig) > 0, "risk", ifelse(-fill(w_orig) < 0, "protective", "none")),
  corrected_direction = ifelse(w_corr > 0, "risk", ifelse(w_corr < 0, "protective", "none")),
  univariate_direction = ifelse(uni > 1, "risk", "protective"), stringsAsFactors = FALSE)
utils::write.csv(wt, file.path(results_dir, "Supplementary_Table_6_CMPIS_weights.csv"), row.names = FALSE)
print(wt)

## ---------------------------------------------------------------------
## 4. Re-validation: original vs corrected (+ sensitivity) across cohorts
## ---------------------------------------------------------------------
cohorts <- list(TCGA_train = train)
for (f in c("gse19234_cohort.txt","gse98394_cohort.txt","gse65904_cohort.txt","gse53118_cohort.txt")) {
  if (!file.exists(f)) next
  d <- read.delim(f, check.names = FALSE)
  d <- d[!is.na(d$time) & d$time > 0 & !is.na(d$status), ]
  cohorts[[toupper(sub("_cohort.txt", "", f))]] <- d
}
schemes <- list(original  = list(w = fill(w_orig), sign = -1),
                corrected = list(w = w_corr,       sign = +1),
                stepcox   = list(w = w_step,       sign = +1),
                mv_cox    = list(w = w_mvc,        sign = +1))
score <- function(d, sc) {
  g <- intersect(genes8, names(d)); z <- scale(as.matrix(d[, g])); z[is.na(z)] <- 0
  as.numeric(scale(sc$sign * (z %*% sc$w[g])))
}
rows <- list()
for (cn in names(cohorts)) for (sn in names(schemes)) {
  d <- cohorts[[cn]]; s <- score(d, schemes[[sn]])
  if (all(!is.finite(s))) next
  C  <- tryCatch(concordance(Surv(time, status) ~ s, data = d, reverse = TRUE)$concordance, error = function(e) NA)
  cf <- tryCatch(summary(coxph(Surv(d$time, d$status) ~ s)), error = function(e) NULL)
  g2 <- factor(ifelse(s >= median(s), "high", "low"))
  lp <- tryCatch(1 - pchisq(survdiff(Surv(d$time, d$status) ~ g2)$chisq, 1), error = function(e) NA)
  rows[[paste(cn, sn)]] <- data.frame(cohort = cn, scheme = sn, n = nrow(d), events = sum(d$status),
    C_index = round(C, 3),
    HR_per_SD = if (!is.null(cf)) round(cf$conf.int[1,1], 3) else NA,
    HR_lower  = if (!is.null(cf)) round(cf$conf.int[1,3], 3) else NA,
    HR_upper  = if (!is.null(cf)) round(cf$conf.int[1,4], 3) else NA,
    cox_p     = if (!is.null(cf)) signif(cf$coefficients[1,5], 3) else NA,
    KM_logrank_p = signif(lp, 3), stringsAsFactors = FALSE)
}
rv <- do.call(rbind, rows); rownames(rv) <- NULL
utils::write.csv(rv, file.path(results_dir, "Supplementary_Table_7_CMPIS_revalidation.csv"), row.names = FALSE)
print(rv)

## ---------------------------------------------------------------------
## 5. Save corrected model (consumed by R/cmpis_weights.R)
## ---------------------------------------------------------------------
saveRDS(list(coefficients = w_corr, method = "stepcox_coef x normalised RSF importance (sign-preserving)",
             direction = "higher score = higher risk (no inversion)",
             stepcox_coefficients = cox_f, rsf_importance = imp_f, source = source_of_raw,
             sensitivity = list(stepcox = w_step, mv_cox = w_mvc), original = fill(w_orig)),
        file.path("skcm_analysis_results", "final_model_info_corrected.rds"))

## ---------------------------------------------------------------------
## 6. Figure
## ---------------------------------------------------------------------
pdf(file.path(results_dir, "Figure_S_CMPIS_reweighting.pdf"), width = 12, height = 5.5)
op <- par(mfrow = c(1, 2), mar = c(6, 5, 3, 1))
W <- rbind(original = -fill(w_orig) / max(abs(fill(w_orig))), corrected = w_corr / max(abs(w_corr)))
barplot(W, beside = TRUE, las = 2, col = c("grey60", "#1F4E79"), ylim = c(-1.1, 1.1),
        ylab = "effective weight (scaled; + = risk)", main = "A. Gene weights: original vs corrected")
abline(h = 0); legend("topleft", c("original (inverted)", "corrected"), fill = c("grey60", "#1F4E79"), bty = "n", cex = .85)
M <- xtabs(C_index ~ scheme + cohort, rv)[c("original","corrected","stepcox","mv_cox"), , drop = FALSE]
barplot(M, beside = TRUE, las = 2, ylim = c(0.4, max(0.75, max(M, na.rm = TRUE) + .03)), xpd = FALSE,
        col = c("grey60","#1F4E79","#6FA0CF","#C55A11"), ylab = "C-index", main = "B. Prognostic C-index by cohort")
abline(h = 0.5, lty = 2)
legend("topright", rownames(M), fill = c("grey60","#1F4E79","#6FA0CF","#C55A11"), bty = "n", cex = .8)
par(op); dev.off()

saveRDS(list(weights = wt, revalidation = rv, source = source_of_raw, generated = Sys.time()),
        file.path(results_dir, "cmpis_reweight.rds"))
msg("Stage 10 complete.")
