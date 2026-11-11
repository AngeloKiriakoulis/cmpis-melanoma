## =====================================================================
## Stage 14 — Manuscript versions of the CMPIS weight, drug-effect and
##            implementation-agreement figures (drawn from saved tables)
## ---------------------------------------------------------------------
## Inputs : results/Supplementary_Table_6_CMPIS_weights.csv     (Stage 10)
##          results/Supplementary_Table_7_CMPIS_revalidation.csv (Stage 10)
##          results/Supplementary_Table_11_CMPIS_drugs.csv       (Stage 11)
##          results/Supplementary_Table_12_drug_oncoPredict.csv  (Stage 12)
## Outputs: results/Figure_MS_CMPIS_weights.pdf
##          results/Figure_MS_drug_effects.pdf
##          results/Figure_MS_drug_implementations.pdf
## Only the final CMPIS weighting and the CS2-vs-CS1 split are shown.
## Run: source("config.R"); source("R/14_Manuscript_Figures.R")
## =====================================================================
rd <- "results"
rcsv <- function(f) read.csv(file.path(rd, f), check.names = FALSE, stringsAsFactors = FALSE)
nice <- function(v) ifelse(grepl("[0-9]", v), v, tools::toTitleCase(tolower(v)))

## ---- 1. CMPIS gene weights + discrimination by cohort --------------
w  <- rcsv("Supplementary_Table_6_CMPIS_weights.csv")
rv <- rcsv("Supplementary_Table_7_CMPIS_revalidation.csv")
hr <- rcsv("Supplementary_Table_5_CMPIS_modules.csv")           # univariate HR per SD, full TCGA-SKCM
w$univariate_HR_TCGA <- hr$univariate_HR[match(w$gene, hr$gene)]
w  <- w[order(w$corrected_weight), ]
pdf(file.path(rd, "Figure_MS_CMPIS_weights.pdf"), width = 12, height = 5.5)
par(mfrow = c(1, 2), mar = c(6, 5, 3.5, 1))
col <- ifelse(w$corrected_weight > 0, "#B4232A", ifelse(w$corrected_weight < 0, "#1F4E79", "grey70"))
bp <- barplot(w$corrected_weight, names.arg = sprintf("%s\nHR %.2f", w$gene, w$univariate_HR_TCGA), las = 1, cex.names = 0.75,
              col = col, border = NA, ylim = c(min(w$corrected_weight) * 1.1, 0.8),
              ylab = "CMPIS weight (+ = risk)", main = "A. CMPIS gene weights")
abline(h = 0)
legend("bottomright", c("risk weight", "protective weight", "zero weight"), fill = c("#B4232A", "#1F4E79", "grey70"), border = NA, bty = "n", cex = 0.8)
mtext("HR: univariate hazard ratio in TCGA-SKCM", side = 1, line = 4, cex = 0.7)
ord <- c(TCGA_train = "TCGA-SKCM\n(training)", GSE19234 = "GSE19234\n(selection)", GSE98394 = "GSE98394\n(selection)",
         GSE65904 = "GSE65904\n(independent)", GSE53118 = "GSE53118\n(independent)")
rv <- rv[rv$scheme %in% c("corrected", "stepcox"), ]
M <- sapply(names(ord), function(k) sapply(c("corrected", "stepcox"), function(s) rv$C_index[rv$cohort == k & rv$scheme == s][1]))
par(mar = c(6, 5, 3.5, 1))
barplot(M, beside = TRUE, names.arg = ord, las = 1, cex.names = 0.75, col = c("#1F4E79", "#9DB9D5"), border = NA,
        ylim = c(0.45, 0.75), xpd = FALSE, ylab = "C-index", main = "B. Discrimination by cohort")
abline(h = 0.5, lty = 2); box()
legend("topright", c("CMPIS (StepCox x RSF)", "StepCox weights only"), fill = c("#1F4E79", "#9DB9D5"), border = NA, bty = "n", cex = 0.8)
dev.off()

## ---- 2. drug effect sizes (models passing QC) ----------------------
dr <- rcsv("Supplementary_Table_11_CMPIS_drugs.csv")
dr <- dr[dr$grouping %in% c("CMPIS_corrected", "CS2_vs_CS1") & dr$reliable %in% c(TRUE, "TRUE"), ]
eff <- function(g) setNames(ifelse(dr$more_sensitive[dr$grouping == g] %in% c("High"), -abs(dr$delta_SD[dr$grouping == g]), abs(dr$delta_SD[dr$grouping == g])),
                            dr$drug[dr$grouping == g])
e1 <- eff("CMPIS_corrected"); e2 <- eff("CS2_vs_CS1")
drugs <- names(sort(e1)); M <- rbind(`CMPIS High vs Low` = e1[drugs], `CS2 vs CS1` = e2[drugs])
fd <- function(g) setNames(dr$FDR[dr$grouping == g], dr$drug[dr$grouping == g])
f1 <- fd("CMPIS_corrected")[drugs]; f2 <- fd("CS2_vs_CS1")[drugs]
pdf(file.path(rd, "Figure_MS_drug_effects.pdf"), width = 9, height = 0.32 * length(drugs) + 2.2)
par(mar = c(5, 10, 3, 1))
bp <- barplot(M[, rev(drugs)], beside = TRUE, horiz = TRUE, las = 1, names.arg = nice(rev(drugs)), col = c("#B4232A", "#1F4E79"), border = NA,
              xlim = range(c(M, 0), na.rm = TRUE) * 1.15, cex.names = 0.8,
              xlab = "Shift in predicted AUC (SD of cell-line AUC)\nnegative = CMPIS-High / CS2 predicted more sensitive",
              main = "Predicted drug sensitivity (models passing QC)")
abline(v = 0)
sig <- rbind(f1, f2)[, rev(drugs)] < 0.05; sig[is.na(sig)] <- FALSE
xs <- M[, rev(drugs)]; text(ifelse(xs < 0, xs - 0.015, xs + 0.015)[sig], bp[sig], "*", cex = 1.1)
legend("topright", c("CMPIS High vs Low", "CS2 vs CS1", "* FDR < 0.05"), fill = c("#B4232A", "#1F4E79", NA), border = NA, bty = "n", cex = 0.8)
dev.off()

## ---- 3. agreement between the two implementations -------------------
ip <- rcsv("Supplementary_Table_12_drug_oncoPredict.csv")
ip <- ip[ip$grouping %in% c("CMPIS_corrected", "CS2_vs_CS1"), ]
sgn <- function(d, s) ifelse(s == "High", -abs(d), abs(d))
x <- sgn(ip$delta_SD_kr, ip$more_sensitive_kr); y <- sgn(ip$delta_SD_op, ip$more_sensitive_op)
rho <- suppressWarnings(cor(x, y, method = "spearman", use = "complete.obs"))
pdf(file.path(rd, "Figure_MS_drug_implementations.pdf"), width = 7, height = 6.5)
par(mar = c(5, 5, 3, 1))
cl <- ifelse(ip$grouping == "CMPIS_corrected", "#B4232A", "#1F4E79"); rel <- ip$reliable %in% c(TRUE, "TRUE")
lim <- range(c(x, y), na.rm = TRUE) * 1.1
plot(x, y, pch = ifelse(rel, 19, 1), col = cl, xlim = lim, ylim = lim, cex = 1.2,
     xlab = "Kernel ridge model: shift in predicted AUC (SD)", ylab = "oncoPredict ridge workflow: shift in predicted AUC (SD)",
     main = sprintf("Agreement between implementations (Spearman rho = %.2f, %d pairs)", rho, sum(is.finite(x) & is.finite(y))), cex.main = 0.95)
abline(0, 1, lty = 2, col = "grey50"); abline(h = 0, v = 0, col = "grey85")
lab <- abs(x) > 0.1 | abs(y) > 0.1
text(x[lab], y[lab], nice(ip$drug[lab]), pos = ifelse(ip$grouping[lab] == "CS2_vs_CS1", 2, 4), cex = 0.65)
legend("topleft", c("CMPIS High vs Low", "CS2 vs CS1", "open = model failed QC"), col = c("#B4232A", "#1F4E79", "grey40"),
       pch = c(19, 19, 1), bty = "n", cex = 0.8)
dev.off()
message("Stage 14 complete: Figure_MS_CMPIS_weights.pdf, Figure_MS_drug_effects.pdf, Figure_MS_drug_implementations.pdf")
