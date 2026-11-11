## =====================================================================
## Stage 9 — Mechanistic network figure
## ---------------------------------------------------------------------
## Builds a figure linking
##   PCD dysregulation -> transcription factors (STAT1/2 vs ZEB2/E2F4)
##   -> immune state (hot vs cold) -> therapeutic response (ICI vs HDAC/PI3K)
## using the inferred regulatory networks / results already in the study.
##
## Data-grounded where possible (no heavy recompute):
##   * 8-gene PCD modules split by the sign of the CMPIS coefficients
##     (skcm_analysis_results/final_model_info.rds): risk(+) vs protective(-)
##   * key TF directions confirmed from the limma CS-vs-Others DE tables
##     (results/tables/consensusMOIC_..._CS{1,2}_vs_Others.txt)
## The four layers and their CS1/CS2 assignment follow Figs 3 (regulons),
## 4 (immune deconvolution) and 6/S6 (oncogenic signalling, drug sensitivity).
##
## Output: results/Figure_Mechanistic_Network.pdf  (+ mechanistic_network.rds)
## Run: source("config.R"); source("R/09_Mechanistic_Network.R")
## =====================================================================

results_dir <- "results"; if (!dir.exists(results_dir)) dir.create(results_dir, recursive = TRUE)
msg <- function(...) message(sprintf(...))

## ---- 1. data-grounded bits -----------------------------------------
## 8-gene modules — derived from data
## module = sign of each gene's contribution to the FINAL MPCDS score
##          (MPCDS = -sum(coef * expr), i.e. risk_direction = "invert")
## cross-check = the gene's own univariate Cox HR for OS in TCGA-SKCM.
pcd8 <- c("CFLAR","CYBB","CYLD","GYS1","IRF1","NFKBIA","PIM2","VCAM1")
risk_genes <- c("GYS1","CYBB","IRF1","NFKBIA"); prot_genes <- c("VCAM1","CYLD","PIM2","CFLAR") # fallback: manuscript text
mod_tbl <- NULL; co <- NULL
source("R/cmpis_weights.R"); cmpis <- load_cmpis()
if (!is.null(cmpis)) { co <- cmpis$coef; msg("CMPIS weights in use: %s", cmpis$version) }
if (!is.null(co)) msg("CMPIS coefficients: %s", paste(sprintf("%s=%.4f", names(co), co), collapse=", "))
if (!is.null(co) && file.exists("skcm.tcga.RData")) {
  tryCatch({
    suppressMessages(library(survival))
    ev <- new.env(); load("skcm.tcga.RData", envir = ev); tc <- get("skcm.tcga", envir = ev)
    ex <- as.matrix(tc$mRNA.expr); if (max(ex, na.rm = TRUE) > 50) ex <- log2(ex + 1)
    cl <- tc$clin.info; cl <- cl[!is.na(cl$futime) & cl$futime > 0 & !is.na(cl$fustat), , drop = FALSE]
    sm <- intersect(colnames(ex), rownames(cl)); tt <- cl[sm, "futime"]; ee <- cl[sm, "fustat"]
    g8 <- intersect(pcd8, rownames(ex))
    mod_tbl <- do.call(rbind, lapply(g8, function(g) {
      cf <- if (g %in% names(co)) as.numeric(co[g]) else 0
      z  <- as.numeric(scale(ex[g, sm]))
      sfit <- summary(coxph(Surv(tt, ee) ~ z))
      ctb <- cmpis$sign * cf
      data.frame(gene = g, coefficient = signif(cf, 4), contribution_to_MPCDS = signif(ctb, 4),
        module = ifelse(abs(ctb) < 1e-8, "neutral", ifelse(ctb > 0, "risk", "protective")),
        univariate_HR = round(sfit$coef[1,2], 3), HR_lower = round(sfit$conf.int[1,3], 3),
        HR_upper = round(sfit$conf.int[1,4], 3), univariate_p = signif(sfit$coef[1,5], 3),
        HR_direction = ifelse(sfit$coef[1,2] > 1, "risk", "protective"), stringsAsFactors = FALSE)
    }))
    mod_tbl$agrees <- mod_tbl$module == "neutral" | mod_tbl$module == mod_tbl$HR_direction
    utils::write.csv(mod_tbl, file.path(results_dir, "Supplementary_Table_5_CMPIS_modules.csv"), row.names = FALSE)
    print(mod_tbl)
    lab <- function(d) paste0(d$gene, ifelse(d$agrees, "", " *"))
    risk_genes <- lab(mod_tbl[mod_tbl$module == "risk", ])
    prot_genes <- lab(mod_tbl[mod_tbl$module == "protective", ])
    neu <- mod_tbl$gene[mod_tbl$module == "neutral"]
    if (length(neu)) msg("Neutral (coefficient ~0, excluded from modules): %s", paste(neu, collapse=","))
    msg("Data-derived modules: risk=%s | protective=%s  (* = univariate HR direction differs)",
        paste(risk_genes, collapse=","), paste(prot_genes, collapse=","))
  }, error = function(e) msg("Module derivation failed (%s) - using manuscript modules.", conditionMessage(e)))
}
## confirm a few TF directions from DE tables (optional annotation)
tf_dir <- list()
de_dir <- function(file, genes) {
  if (!file.exists(file)) return(NULL)
  d <- tryCatch(read.delim(file, row.names = 1, check.names = FALSE), error = function(e) NULL)
  if (is.null(d)) return(NULL)
  lc <- grep("log2fc", names(d), ignore.case = TRUE, value = TRUE)[1]
  if (is.na(lc)) return(NULL)
  g <- intersect(genes, rownames(d)); setNames(d[g, lc], g)
}
cs1de <- de_dir("results/tables/consensusMOIC_TCGA-SKCM_limma_test_result.CS1_vs_Others.txt",
                c("STAT1","STAT2","IRF1","MEF2A","ZBTB7A"))
cs2de <- de_dir("results/tables/consensusMOIC_TCGA-SKCM_limma_test_result.CS2_vs_Others.txt",
                c("ZEB2","E2F4","MITF"))
if (!is.null(cs1de)) msg("CS1 TF expression log2FC (not regulon): %s", paste(sprintf("%s=%.2f", names(cs1de), cs1de), collapse=", "))
if (!is.null(cs2de)) msg("CS2 TF expression log2FC (not regulon): %s", paste(sprintf("%s=%.2f", names(cs2de), cs2de), collapse=", "))

## PCD genes are placed on a track by SUBTYPE EXPRESSION (CS1 vs Others, limma),
## and annotated with their role in the CMPIS: (R) risk, (P) protective, (0) zero weight.
## A gene's CMPIS weight is multivariable, so it need not match its subtype.
pcd_de <- NULL
f1 <- "results/tables/consensusMOIC_TCGA-SKCM_limma_test_result.CS1_vs_Others.txt"
if (file.exists(f1)) {
  d <- tryCatch(read.delim(f1, row.names = 1, check.names = FALSE), error = function(e) NULL)
  if (!is.null(d)) {
    lc <- grep("log2fc", names(d), ignore.case = TRUE, value = TRUE)[1]
    pc <- grep("padj|fdr", names(d), ignore.case = TRUE, value = TRUE)[1]
    g  <- intersect(pcd8, rownames(d))
    pcd_de <- data.frame(gene = g, CS1_vs_Others_log2FC = round(d[g, lc], 3),
                         CS1_vs_Others_padj = signif(d[g, pc], 3), stringsAsFactors = FALSE)
    pcd_de$enriched_in <- ifelse(pcd_de$CS1_vs_Others_padj >= 0.05, "none",
                          ifelse(pcd_de$CS1_vs_Others_log2FC > 0, "CS1", "CS2"))
  }
}
cs1_genes <- prot_genes; cs2_genes <- risk_genes; ns_genes <- character(0)
cs1_title <- "CS1 PCD module (protective)"; cs2_title <- "CS2 PCD module (risk)"
if (!is.null(mod_tbl) && !is.null(pcd_de)) {
  mod_tbl <- merge(mod_tbl, pcd_de, by = "gene", all.x = TRUE)
  utils::write.csv(mod_tbl, file.path(results_dir, "Supplementary_Table_5_CMPIS_modules.csv"), row.names = FALSE)
  role <- c(risk = "R", protective = "P", neutral = "0")[mod_tbl$module]
  mod_tbl$label <- paste0(mod_tbl$gene, " (", role, ")", ifelse(mod_tbl$agrees, "", "*"))
  pack <- function(x) if (!length(x)) "-" else
    vapply(split(x, ceiling(seq_along(x) / 2)), paste, "", collapse = ", ")
  cs1_genes <- pack(mod_tbl$label[mod_tbl$enriched_in %in% "CS1"])
  cs2_genes <- pack(mod_tbl$label[mod_tbl$enriched_in %in% "CS2"])
  ns_genes  <- mod_tbl$label[!mod_tbl$enriched_in %in% c("CS1", "CS2")]
  cs1_title <- "PCD genes up in CS1"; cs2_title <- "PCD genes up in CS2"
  msg("Subtype placement: CS1=%s | CS2=%s | n.s.=%s",
      paste(mod_tbl$gene[mod_tbl$enriched_in %in% "CS1"], collapse=","),
      paste(mod_tbl$gene[mod_tbl$enriched_in %in% "CS2"], collapse=","), paste(ns_genes, collapse=","))
}

## ---- therapeutic layer: data-driven from Stage 11 ---------------------
## Predicted GDSC2 sensitivity, CS2 vs CS1, drugs whose model predicts held-out cell lines, FDR < 0.05.
chunk2 <- function(v) if (!length(v)) "none (FDR < 0.05)" else
  vapply(split(v, ceiling(seq_along(v) / 2)), paste, "", collapse = ", ")
rx1 <- c("Predicted drug sensitivity: run Stage 11", "ICI: inflamed phenotype; response not", "confirmed in anti-PD-1 cohorts (Fig S9)")
rx2 <- c("Predicted drug sensitivity: run Stage 11")
cf <- file.path(results_dir, "cmpis_clinical.rds")
if (file.exists(cf)) tryCatch({
  dr <- readRDS(cf)$drugs
  if (!is.null(dr)) {
    z <- dr[dr$grouping == "CS2_vs_CS1" & dr$reliable & !is.na(dr$FDR) & dr$FDR < 0.05, , drop = FALSE]
    z <- z[order(-abs(z$delta_SD)), , drop = FALSE]
    nm <- function(v) head(ifelse(grepl("[0-9]", v), v, tools::toTitleCase(tolower(v))), 4)   # keep GSK-2141795, BEZ235 as-is
    c2 <- nm(z$drug[z$more_sensitive == "High"]); c1 <- nm(z$drug[z$more_sensitive == "Low"])
    rx1 <- c("More sensitive in CS1 (in silico):", chunk2(c1), "ICI: inflamed phenotype; response not", "confirmed in anti-PD-1 cohorts (Fig S9)")
    rx2 <- c("More sensitive in CS2 (in silico):", chunk2(c2))
    msg("Therapeutic layer from Stage 11: CS2=%s | CS1=%s", paste(c2, collapse = ","), paste(c1, collapse = ","))
  }
}, error = function(e) msg("Stage 11 drug table unreadable (%s).", conditionMessage(e)))

## ---- 2. figure ------------------------------------------------------
CS1 <- "#1F4E79"; CS2 <- "#B4232A"; CS1f <- "#EAF1F8"; CS2f <- "#FBEBEC"; ink <- "#222222"

node <- function(x, y, w, h, title, lines, edge, fill) {
  rect(x-w/2, y-h/2, x+w/2, y+h/2, col = fill, border = edge, lwd = 2)
  rect(x-w/2, y+h/2-2.4, x+w/2, y+h/2, col = edge, border = edge)          # header bar
  text(x, y+h/2-1.2, title, col = "white", font = 2, cex = 0.82)
  if (length(lines))
    text(x, y+h/2-3.0, paste(lines, collapse = "\n"), col = ink, cex = 0.72, adj = c(.5,1))
}
arrow_between <- function(x1, y, x2, col) arrows(x1, y, x2, y, length = .10, lwd = 2.5, col = col)

pdf(file.path(results_dir, "Figure_Mechanistic_Network.pdf"), width = 13, height = 7.2)
op <- par(mar = c(1,1,3,1)); plot.new(); plot.window(xlim = c(0,100), ylim = c(0,100))

## column centres and layer titles
xc <- c(13, 38, 63, 87); hdr <- c("PCD program", "Transcriptional regulators",
                                  "Immune microenvironment", "Therapeutic vulnerability")
for (i in seq_along(xc)) text(xc[i], 87, hdr[i], font = 2, cex = 0.95, col = ink)
text(50, 96, "PCD dysregulation -> transcription factors -> immune state -> therapeutic response",
     font = 2, cex = 1.05)
y1 <- 67; y2 <- 34; w <- 21; h <- 17                   # CS1 track (top), CS2 track (bottom)

## CS1 (immune-hot) track
node(xc[1], y1, w, h, cs1_title, cs1_genes, CS1, CS1f)
node(xc[2], y1, w, h, "TFs: interferon / homeostasis",
     c("STAT1, STAT2, IRF1","MEF2A, ZBTB7A"), CS1, CS1f)
node(xc[3], y1, w, h, "Immune-HOT",
     c("CD8+ T, NK, memory B","IFN-gamma, antigen presentation","high Immune/Stromal score"), CS1, CS1f)
node(xc[4], y1, w, h, "Predicted therapeutic context", rx1, CS1, CS1f)

## CS2 (immune-cold) track
node(xc[1], y2, w, h, cs2_title, cs2_genes, CS2, CS2f)
node(xc[2], y2, w, h, "TFs: lineage / cell cycle",
     c("MITF, E2F4, ZEB2","(melanocytic, MITF-high)"), CS2, CS2f)
node(xc[3], y2, w, h, "Immune-COLD",
     c("immune-desert","WNT/beta-catenin, VEGFR","T-cell exclusion"), CS2, CS2f)
node(xc[4], y2, w, h, "Predicted therapeutic context", rx2, CS2, CS2f)

## arrows within tracks
for (i in 1:3) { arrow_between(xc[i]+w/2, y1, xc[i+1]-w/2, CS1)
                 arrow_between(xc[i]+w/2, y2, xc[i+1]-w/2, CS2) }

## legend / caveat
legend(38, 17, c("CS1 (immune-hot)", "CS2 (immune-cold)"), fill = c(CS1, CS2),
       border = NA, bty = "n", cex = 0.85, horiz = TRUE)
text(50, 7, paste("Arrows denote directional associations inferred from regulon activity (Fig 3),",
     "immune deconvolution (Fig 4) and drug-sensitivity analyses (Fig S10); associative, not causal."),
     cex = 0.72, col = "#555555")
if (!is.null(mod_tbl) && !is.null(mod_tbl$label))
  text(50, 3, paste0("PCD genes placed by subtype expression (limma FDR < 0.05). (R)/(P)/(0) = risk / protective / zero weight in the CMPIS; ",
       "* = multivariable weight opposes the gene's univariate HR",
       if (length(ns_genes)) paste0("; not subtype-specific: ", paste(ns_genes, collapse = ", ")) else "", "."),
       cex = 0.62, col = "#555555")
par(op); dev.off()
msg("Wrote %s", file.path(results_dir, "Figure_Mechanistic_Network.pdf"))

saveRDS(list(module_table = mod_tbl, cs1_pcd = cs1_genes, cs2_pcd = cs2_genes, risk_module = risk_genes, protective_module = prot_genes,
             cs1_tf_log2fc = cs1de, cs2_tf_log2fc = cs2de, generated = Sys.time()),
        file.path(results_dir, "mechanistic_network.rds"))
msg("Stage 9 complete.")
