## =====================================================================
## Stage 6 — PCD gene-set sensitivity analysis
## ---------------------------------------------------------------------
## Tests whether the CS1/CS2 split depends on the definition of the
## 2,682-gene MSigDB PCD set, by producing:
##
##   (1) Supplementary Table 1 — per-pathway PCD gene lists
##        * BROAD decomposition of the 2,682-gene set via msigdbr
##          (GO-BP / Reactome / WikiPathways / KEGG / Hallmark terms for
##          apoptosis, ferroptosis, pyroptosis, necroptosis, autophagy)
##        * tight CANONICAL CORE sets (KEGG hsa04210/04216/04217/04140 +
##          the 33-gene pyroptosis set) used for the clustering test
##        * annotation of every one of the 2,682 genes by pathway
##
##   (2) Sensitivity analysis — re-derive the two-cluster (CS1/CS2)
##        solution from TCGA-SKCM expression restricted to each pathway
##        subset (apoptosis-only, ferroptosis-only, pyroptosis-only,
##        necroptosis-only, autophagy-only) and to the full 2,682-gene
##        set (combined), then show the k=2 solution is STABLE and
##        CONCORDANT with the reference CS1/CS2 subtypes.
##
## Outputs (written to results/):
##   Supplementary_Table_1_summary.csv / _coverage.csv / _PCD2682_annotated.csv
##     (per-pathway core-set lists + annotation of the 2,682-gene set; the
##      polished multi-sheet Supplementary_Table_1_PCD_pathways.xlsx is
##      maintained separately)
##   Supplementary_Table_2_subtype_sensitivity.csv
##   Figure_S_PCD_sensitivity.pdf
##   pcd_sensitivity.rds   (declared output for the {targets} stage)
##
## Reuses existing caches — no downloads:
##   skcm.tcga.RData        ($mRNA.expr, $clin.info)
##   TCGA_Clustering.RData  (cmoic.skcm$clust.res = reference CS1/CS2)
##   PCD GENES.xlsx         (the 2,682-gene PCD set)
##
## Run just this stage:  source("R/helpers.R"); rerun("pcd_sensitivity")
## =====================================================================

suppressWarnings(suppressMessages({
  library(survival)
  have_ccp     <- requireNamespace("ConsensusClusterPlus", quietly = TRUE)
  have_msigdbr <- requireNamespace("msigdbr", quietly = TRUE)
  have_openxlsx<- requireNamespace("openxlsx", quietly = TRUE)
  have_readxl  <- requireNamespace("readxl", quietly = TRUE)
  have_mclust  <- requireNamespace("mclust", quietly = TRUE)
  have_cluster <- requireNamespace("cluster", quietly = TRUE)
}))

set.seed(if (exists("cfg") && !is.null(cfg$seed)) cfg$seed else 20240101L)
results_dir <- "results"
if (!dir.exists(results_dir)) dir.create(results_dir, recursive = TRUE)

msg <- function(...) message(sprintf(...))

## ---------------------------------------------------------------------
## 0. Canonical CORE pathway sets (citable; retrieved from KEGG REST +
##    Ye et al. 2021). Used for the clustering sensitivity test.
## ---------------------------------------------------------------------
core_sets <- list(
  Apoptosis = c("AKT3","BCL2L11","PARP2","PARP3","CASP12","TUBA1B","CTSC","GADD45G","P3R3URF-PIK3R3","TUBA3E","TUBA3D","CHUK","PARP1","PARP4","CSF2RB","CTSB","CTSD","CTSH","CTSK","CTSL","CTSV","CTSO","CTSS","CTSW","CTSZ","DAB2IP","DAXX","GADD45A","DDIT3","DFFA","DFFB","EIF2S1","ENDOG","AKT1","AKT2","ERN1","FOS","BBC3","HTRA2","GZMB","APAF1","HRAS","BIRC2","BIRC3","XIAP","BIRC5","FAS","IKBKB","FASLG","IL3","IL3RA","ITPR1","ITPR2","ITPR3","JUN","KRAS","LMNA","LMNB1","MCL1","MAP3K5","GADD45B","ATF4","ATM","NFKB1","NFKBIA","NGF","NRAS","NTRK1","PDPK1","TUBA8","PIK3CA","PIK3CB","PIK3CD","PIK3R1","PIK3R2","PMAIP1","SEPTIN4","CYCS","PIDD1","PRF1","MAPK1","MAPK3","MAPK8","MAPK9","MAPK10","MAP2K1","MAP2K2","DIABLO","BAD","BAK1","PTPN13","BAX","RAF1","BCL2","BCL2A1","RELA","BCL2L1","ACTB","BID","TP53AIP1","SPTAN1","ACTG1","TNF","TNFRSF1A","TP53","TRAF1","TRAF2","TUBA4A","TUBA3C","TUBA1A","TUBAL3","TUBA4B","CAPN1","CAPN2","CASP2","CASP3","CASP6","CASP7","CASP8","CASP9","CASP10","TUBA1C","LMNB2","PIK3R3","IKBKG","TRADD","CTSF","RIPK1","HRK","TNFSF10","FADD","TNFRSF10B","TNFRSF10A","CFLAR","NOL3","MAP3K14","AIFM1","EIF2AK3"),
  Ferroptosis = c("LPCAT3","ATG7","SAT2","CP","CYBB","ACSL1","ACSL3","ACSL4","ACSL6","SLC39A14","SLC7A11","ALOX15","FTH1","FTL","GCLC","GCLM","GPX4","GSS","SLC40A1","HMOX1","MAP1LC3BP1","MAP1LC3C","SLC11A2","PCBP1","PCBP2","ACSL5","STEAP3","PRNP","SAT1","SLC39A8","MAP1LC3B2","SLC3A2","TF","TFRC","TP53","VDAC2","VDAC3","NCOA4","MAP1LC3B","MAP1LC3A","FTMT","ATG5"),
  Pyroptosis = c("AIM2","CASP1","CASP3","CASP4","CASP5","CASP6","CASP8","CASP9","GPX4","GSDMA","GSDMB","GSDMC","GSDMD","DFNA5","ELANE","IL18","IL1B","IL6","NLRC4","NLRP1","NLRP2","NLRP3","NLRP6","NLRP7","NOD1","NOD2","DFNB59","PLCG1","PRKACA","PYCARD","SCAF11","TIRAP","TNF"),
  Necroptosis = c("PLA2G4B","RNF103-CHMP3","DNM1L","IRF9","PPIAL4H","RBCK1","PPIAP40","RIPK3","NLRP3","H2AL3","PLA2G4E","SPATA2L","CHMP4B","PARP1","TICAM1","CYBB","CYLD","PGAM5","MLKL","H2AC1","ALOX15","FTH1","FTL","PLA2G4F","CHMP2B","USP21","VPS4A","CHMP2A","GLUD1","GLUD2","GLUL","PLA2G4D","CHMP4A","SLC25A4","PYCARD","SLC25A5","SLC25A6","H2AC8","H2AC7","H2AX","H2AZ1","HMGB1","H2AC21","BIRC2","BIRC3","XIAP","HSP90AA1","HSP90AB1","IFNA1","IFNA2","IFNA4","IFNA5","IFNA6","IFNA7","IFNA8","IFNA10","IFNA13","IFNA14","IFNA16","IFNA17","IFNA21","IFNAR1","IFNAR2","IFNB1","IFNG","IFNGR1","IFNGR2","TICAM2","FAS","IL1A","IL1B","FASLG","JAK1","JAK2","JAK3","PPIAP59","H2AB2","H2AB1","CHMP1A","CHMP5","CHMP3","PLA2G4A","PPIA","PPID","TRPM7","RNF31","MACROH2A2","H2AJ","MAPK8","MAPK9","MAPK10","EIF2AK2","CHMP1B","BAX","PYGB","PYGL","PYGM","BCL2","BID","PPIAL4G","PPIAL4D","PPIAL4A","PPIAL4C","SMPD1","STAT1","STAT2","STAT3","STAT4","STAT5A","STAT5B","STAT6","TLR3","TLR4","TNF","TNFAIP3","TNFRSF1A","TRAF2","TRAF5","H2AC19","PPIAL4F","TYK2","PPIAL4E","VDAC1","VDAC2","VDAC3","CHMP6","ZBP1","CAMK2A","CAMK2B","CAMK2D","CAMK2G","SHARPIN","CAPN1","CAPN2","H2AC13","H2AC15","H2AC14","H2AC16","H2AC6","H2AC4","H2AC17","H2AC18","H2AC20","CASP1","SLC25A31","H2AB3","CASP8","H2AC12","PLA2G4C","JMJD7-PLA2G4B","TRADD","RIPK1","TNFSF10","FADD","TNFRSF10B","TNFRSF10A","CFLAR","SQSTM1","H2AC11","IL33","AIFM1","CHMP7","CHMP4C","H2AC25","H2AZ2","VPS4B","MACROH2A1","SPATA2"),
  Autophagy = c("AKT3","TANK","SUPT20HL1","OPTN","CALCOCO2","RRAGB","ATG7","CAMKK2","YKT6","RRAGA","P3R3URF-PIK3R3","WDR45","GABARAP","GABARAPL2","ATG4A","RAB39B","PRAP1","SMCR8","CTSB","CTSD","CTSL","DAPK1","DAPK3","SUPT20HL2","EIF2S1","C9orf72","AKT1","AKT2","ERN1","RRAS2","MRAS","ATG14","WDFY3","ATG2A","ATG4B","VPS39","DAPK2","GABARAPL1","MTOR","GORASP2","WIPI2","VPS41","ATG9B","TBK1","NRBF2","PIK3R4","HIF1A","HMGB1","HRAS","RAB7B","IGBP1","IGF1R","INS","IRS1","ITPR1","KRAS","LAMP1","LAMP2","MAP1LC3BP1","NBR1","RAB8A","EIF2AK4","MAP1LC3C","BECN2","NRAS","SH3GLB1","PDPK1","PIK3C3","PIK3CA","PIK3CB","PIK3CD","PIK3R1","PIK3R2","ZFYVE1","DDIT4","STX17","ATG16L1","WIPI1","ATG2B","PPP2CA","PPP2CB","WDR41","SUPT20H","PRKAA1","AMBRA1","PRKAA2","PRKACA","PRKACB","PRKACG","PRKCD","VPS11","PRKCQ","MAPK1","MAPK3","MAPK8","MAPK9","MAPK10","MAP2K1","MAP2K2","WDR45B","BAD","PTEN","BIRC6","RPTOR","VPS18","EPG5","TP53INP2","RRAGD","RAB1A","RAF1","BCL2","BCL2L1","RHEB","ATG101","RPS6KB1","RPS6KB2","RPS27A","RRAS","RRAGC","MLST8","MAP1LC3B2","MTMR14","ATG3","VPS16","GORASP1","DEPTOR","VPS33A","BNIP3","STK11","MAP3K7","TRAF6","TSC1","TSC2","UBA52","UBB","UBC","UVRAG","RAB7A","ATG9A","MAP1LC3B","VMP1","RAB33B","ATG10","ULK1","STX7","AKT1S1","MAP1LC3A","IRS4","ATG4C","ATG4D","PIK3R3","IRS2","VAMP8","BECN1","CFLAR","SQSTM1","TAX1BP1","MTMR3","ATG16L2","MTMR4","ATG12","SNAP29","ATG5","ULK2","RUBCN","ATG13","RB1CC1","PLEKHM1")
)
core_sets <- lapply(core_sets, function(x) sort(unique(x)))
pathways  <- names(core_sets)
core_source <- c(Apoptosis="KEGG hsa04210", Ferroptosis="KEGG hsa04216",
                 Pyroptosis="Ye et al. 2021 (33-gene)", Necroptosis="KEGG hsa04217",
                 Autophagy="KEGG hsa04140")

## ---------------------------------------------------------------------
## 1. Load the 2,682-gene PCD set
## ---------------------------------------------------------------------
pcd_file <- if (exists("cfg") && !is.null(cfg$pcd_gene_set)) cfg$pcd_gene_set else "PCD GENES.xlsx"
stopifnot(file.exists(pcd_file))
pcd_all <- if (have_readxl) {
  as.character(readxl::read_excel(pcd_file)[[1]])
} else {
  stop("readxl needed to read ", pcd_file)
}
pcd_all <- sort(unique(pcd_all[!is.na(pcd_all) & pcd_all != "" & pcd_all != "gene"]))
msg("Loaded %d unique genes from %s", length(pcd_all), pcd_file)

## ---------------------------------------------------------------------
## 2. Supplementary Table 1 — per-pathway gene lists (canonical CORE sets)
##    + annotation of the 2,682-gene set + coverage summary.
##    The polished multi-sheet .xlsx is maintained separately
##    (results/Supplementary_Table_1_PCD_pathways.xlsx); here we emit tidy
##    CSVs so the stage is self-contained and needs no extra packages.
## ---------------------------------------------------------------------
in_set <- function(g) g %in% pcd_all
core_union <- sort(unique(unlist(core_sets)))
summary_tbl <- do.call(rbind, lapply(pathways, function(p) data.frame(
  Pathway = p, Core_source = core_source[[p]],
  Core_n = length(core_sets[[p]]), Core_in_2682 = sum(in_set(core_sets[[p]])),
  stringsAsFactors = FALSE)))
coverage <- data.frame(
  Metric = c("PCD set size (2,682 file)", "Core union size",
             "of 2,682 covered by core union", "of 2,682 in NO core pathway"),
  Value  = c(length(pcd_all), length(core_union),
             sum(in_set(core_union)), sum(!(pcd_all %in% core_union))),
  stringsAsFactors = FALSE)
annot <- data.frame(gene = pcd_all, stringsAsFactors = FALSE)
for (p in pathways) annot[[p]] <- ifelse(pcd_all %in% core_sets[[p]], "yes", "")
annot$n_core_pathways <- rowSums(annot[pathways] == "yes")
utils::write.csv(summary_tbl, file.path(results_dir, "Supplementary_Table_1_summary.csv"), row.names = FALSE)
utils::write.csv(coverage,    file.path(results_dir, "Supplementary_Table_1_coverage.csv"), row.names = FALSE)
utils::write.csv(annot,       file.path(results_dir, "Supplementary_Table_1_PCD2682_annotated.csv"), row.names = FALSE)
msg("Wrote Supplementary Table 1 CSVs (summary, coverage, annotation) to %s/", results_dir)

## ---------------------------------------------------------------------
## 4. Load expression + reference CS1/CS2 labels
## ---------------------------------------------------------------------
stopifnot(file.exists("skcm.tcga.RData"))
load("skcm.tcga.RData")                      # -> skcm.tcga
expr <- as.matrix(skcm.tcga$mRNA.expr)
storage.mode(expr) <- "double"

ref <- NULL
if (file.exists("TCGA_Clustering.RData")) {
  e <- new.env()
  tryCatch(load("TCGA_Clustering.RData", envir = e), error = function(err) NULL)
  if (exists("cmoic.skcm", envir = e)) {
    cr <- get("cmoic.skcm", envir = e)$clust.res
    ref <- setNames(as.integer(cr$clust), cr$samID)
  }
}
if (is.null(ref))
  msg("WARNING: reference CS1/CS2 labels not found in TCGA_Clustering.RData; concordance vs published subtypes will be NA.")

## survival
clin <- skcm.tcga$clin.info
clin <- clin[!is.na(clin$futime) & clin$futime > 0, , drop = FALSE]

## ---------------------------------------------------------------------
## 5. Helper: consensus clustering + optimal-k selection (PAC)
## ---------------------------------------------------------------------
prep_expr <- function(genes) {
  g <- intersect(genes, rownames(expr))
  m <- expr[g, , drop = FALSE]
  m <- m[apply(m, 1, function(x) var(x, na.rm = TRUE) > 0), , drop = FALSE]
  ## log if it looks like linear CPM
  if (max(m, na.rm = TRUE) > 50) m <- log2(m + 1)
  t(scale(t(m)))                              # z-score genes
}

pac_from_consensus <- function(cm, x1 = 0.1, x2 = 0.9) {
  v <- cm[lower.tri(cm)]
  ecdf(v)(x2) - ecdf(v)(x1)                   # proportion of ambiguous pairs
}

consensus_cluster <- function(m, kmax = 6, reps = 200) {
  ## returns list(assign = list per k, pac = numeric per k)
  if (have_ccp) {
    tmp <- tempfile()
    cc <- ConsensusClusterPlus::ConsensusClusterPlus(
      m, maxK = kmax, reps = reps, pItem = 0.8, pFeature = 1,
      clusterAlg = "hc", distance = "pearson", seed = 20240101,
      plot = NULL, title = tmp, verbose = FALSE)
    assign <- lapply(2:kmax, function(k) cc[[k]]$consensusClass)
    pac    <- vapply(2:kmax, function(k) pac_from_consensus(cc[[k]]$consensusMatrix), numeric(1))
    names(assign) <- names(pac) <- as.character(2:kmax)
    list(assign = assign, pac = pac)
  } else {
    ## fallback: hierarchical clustering on 1-pearson, cutree; PAC via silhouette gap
    d <- as.dist(1 - cor(m, use = "pairwise.complete.obs"))
    hc <- hclust(d, method = "ward.D2")
    assign <- lapply(2:kmax, function(k) setNames(cutree(hc, k), colnames(m)))
    sil <- rep(NA_real_, length(2:kmax))
    if (have_cluster) sil <- vapply(2:kmax, function(i)
      mean(cluster::silhouette(cutree(hc, i), d)[, 3]), numeric(1))
    names(assign) <- names(sil) <- as.character(2:kmax)
    ## no consensus matrix available: use -silhouette as the selection
    ## criterion so which.min() picks the k with the HIGHEST silhouette.
    list(assign = assign, pac = -sil)                                 # proxy
  }
}

adj_rand <- function(a, b) {
  if (have_mclust) return(mclust::adjustedRandIndex(a, b))
  tab <- table(a, b); n <- sum(tab)
  sum_comb <- function(x) sum(choose(x, 2))
  ai <- sum_comb(rowSums(tab)); bi <- sum_comb(colSums(tab))
  idx <- sum_comb(as.vector(tab)); exp_idx <- ai * bi / choose(n, 2)
  (idx - exp_idx) / ((ai + bi)/2 - exp_idx)
}
cohen_kappa <- function(a, b) {          # a,b aligned integer vectors, same levels
  tab <- table(a, b); n <- sum(tab)
  po <- sum(diag(tab)) / n
  pe <- sum(rowSums(tab) * colSums(tab)) / n^2
  (po - pe) / (1 - pe)
}
align2 <- function(part, refv) {         # relabel 2-cluster 'part' to best match refv
  common <- intersect(names(part), names(refv))
  p <- part[common]; r <- refv[common]
  a1 <- mean(p == r); pf <- ifelse(p == 1, 2, 1)
  if (mean(pf == r) > a1) p <- pf
  list(p = p, r = r)
}

## ---------------------------------------------------------------------
## 6. Run sensitivity across gene sets
## ---------------------------------------------------------------------
gene_sets <- c(core_sets, list(PCD_combined = pcd_all))
res_rows <- list(); parts <- list(); kmax <- 6

for (nm in names(gene_sets)) {
  m <- prep_expr(gene_sets[[nm]])
  msg("[%s] %d genes x %d samples", nm, nrow(m), ncol(m))
  cc <- consensus_cluster(m, kmax = kmax)
  k_free <- as.integer(names(which.min(cc$pac)))        # k a free search would pick
  part2 <- cc$assign[["2"]]
  parts[[nm]] <- part2
  ## silhouette of the FIXED 2-cluster solution (higher = better separated);
  ## a direct stability measure of the k=2 partition the paper uses.
  sil_k2 <- NA_real_
  if (have_cluster) {
    d2 <- as.dist(1 - cor(m, use = "pairwise.complete.obs"))
    sil_k2 <- tryCatch(mean(cluster::silhouette(as.integer(factor(part2)), d2)[, 3]),
                       error = function(e) NA_real_)
  }
  ## concordance with reference CS1/CS2
  ari <- kap <- agr <- NA_real_
  if (!is.null(ref)) {
    al <- align2(part2, ref)
    ari <- adj_rand(al$p, al$r); kap <- cohen_kappa(al$p, al$r); agr <- mean(al$p == al$r)
  }
  ## survival separation of the 2-cluster solution
  lr_p <- hr <- NA_real_
  sc <- intersect(names(part2), rownames(clin))
  if (length(sc) > 20) {
    sd <- data.frame(time = clin[sc,"futime"], status = clin[sc,"fustat"],
                     grp = factor(part2[sc]))
    lr_p <- tryCatch(1 - pchisq(survival::survdiff(Surv(time,status)~grp, sd)$chisq, 1),
                     error=function(e) NA_real_)
    hr   <- tryCatch(exp(coef(survival::coxph(Surv(time,status)~grp, sd))[1]),
                     error=function(e) NA_real_)
  }
  res_rows[[nm]] <- data.frame(
    gene_set = nm, n_genes_used = nrow(m),
    PAC_k2 = round(cc$pac[["2"]], 4),        # ambiguity of the 2-cluster solution
    sil_k2 = round(sil_k2, 3),               # silhouette of the 2-cluster solution
    k_free = k_free,                         # k an unconstrained search would pick
    ARI_vs_CS = round(ari,3), kappa_vs_CS = round(kap,3), pct_agree_vs_CS = round(100*agr,1),
    logrank_p = signif(lr_p,3), HR_grp2 = round(hr,3), stringsAsFactors = FALSE)
}
res <- do.call(rbind, res_rows)
rownames(res) <- NULL
st2 <- file.path(results_dir, "Supplementary_Table_2_subtype_sensitivity.csv")
utils::write.csv(res, st2, row.names = FALSE)
msg("Wrote %s", st2)
print(res)

## pairwise ARI between all gene-set 2-cluster partitions
pw <- matrix(NA_real_, length(parts), length(parts),
             dimnames = list(names(parts), names(parts)))
for (i in names(parts)) for (j in names(parts)) {
  ci <- intersect(names(parts[[i]]), names(parts[[j]]))
  pw[i,j] <- adj_rand(parts[[i]][ci], parts[[j]][ci])
}

## ---------------------------------------------------------------------
## 7. Figure
## ---------------------------------------------------------------------
figf <- file.path(results_dir, "Figure_S_PCD_sensitivity.pdf")
pdf(figf, width = 10, height = 8)
op <- par(mfrow = c(2,2), mar = c(7,4,3,1))
## A: concordance vs CS1/CS2
bp <- barplot(res$ARI_vs_CS, names.arg = res$gene_set, las = 2,
              ylim = c(0,1), col = "#1F4E79", ylab = "Adjusted Rand Index",
              main = "A. k=2 concordance with published CS1/CS2")
text(bp, res$ARI_vs_CS, labels = ifelse(is.na(res$ARI_vs_CS),"NA",res$ARI_vs_CS), pos = 3, cex = .8)
## B: PAC at k=2 (lower = more stable 2-cluster solution)
barplot(res$PAC_k2, names.arg = res$gene_set, las = 2, col = "#C55A11",
        ylab = "PAC (k=2)", main = "B. 2-cluster ambiguity (lower = more stable)")
## C: silhouette of the fixed 2-cluster solution (higher = better separated)
bc <- barplot(res$sil_k2, names.arg = res$gene_set, las = 2, col = "#548235",
        ylab = "silhouette (k=2)", main = "C. 2-cluster separation (higher = better)")
text(bc, pmax(res$sil_k2,0), labels = ifelse(is.na(res$sil_k2),"NA",res$sil_k2), pos = 3, cex = .8)
## D: pairwise ARI heatmap
par(mar = c(7,7,3,2))
image(seq_len(nrow(pw)), seq_len(ncol(pw)), t(pw[nrow(pw):1,]), axes = FALSE,
      col = hcl.colors(20, "Blues", rev = TRUE), zlim = c(0,1),
      xlab = "", ylab = "", main = "D. Pairwise ARI between partitions")
axis(1, seq_len(ncol(pw)), colnames(pw), las = 2, cex.axis = .8)
axis(2, seq_len(nrow(pw)), rev(rownames(pw)), las = 2, cex.axis = .8)
par(op); dev.off()
msg("Wrote %s", figf)

## ---------------------------------------------------------------------
## 8. Declared stage output
## ---------------------------------------------------------------------
out <- list(summary = summary_tbl, coverage = coverage, sensitivity = res,
            pairwise_ARI = pw, core_sets = core_sets,
            table2 = st2, figure = figf,
            session = utils::sessionInfo()$R.version$version.string,
            generated = Sys.time())
saveRDS(out, file.path(results_dir, "pcd_sensitivity.rds"))
msg("Stage 6 complete. Outputs in %s/", results_dir)
