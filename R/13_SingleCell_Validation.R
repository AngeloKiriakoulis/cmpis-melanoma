## =====================================================================
## Stage 13 — Single-cell validation of the CMPIS / PCD programme (Fig 7)
## ---------------------------------------------------------------------
## Datasets (GEO, downloaded into data/sc/ on first run):
##   Tirosh 2016        GSE72056   4.6k cells, malignant + non-malignant
##   Jerby-Arnon 2018   GSE115978  7.2k cells, malignant + non-malignant
##   Sade-Feldman 2018  GSE120575  16k CD45+ cells (immune only), anti-PD-1/CTLA-4
##
## Questions -> where answered
##   1. CMPIS genes across malignant/T/NK/macrophage/fibroblast cells   Fig 7e + Table 13
##   2. PCD score by cell type                                          Fig 7f + Table 13b
##   3. Is the immune-cold subtype (CS2) enriched for dedifferentiated
##      melanoma cells?     tumour-level CS2-vs-CS1 score vs programmes  Table 13c
##   4. Do CS2-like tumours show less CD8+ infiltration?                 Table 13c
##   5. Do WNT/beta-catenin signatures track low PCD activity?           Table 13c
##   +  Cell-cell communication (curated ligand-receptor list)           Table 13d, Fig 7 p.4
##
## Design notes (read before interpreting):
##  * CS1/CS2 signatures = top 100 genes up in each subtype (TCGA limma, CS-vs-Others).
##    A signature that includes immune genes is circular with immune infiltration, so
##    the tumour-level tests use a TUMOUR-INTRINSIC version: only signature genes
##    detected in >=10% of that dataset's malignant cells, scored on malignant cells only.
##  * Tumour-level statistics treat each tumour as one observation (no pseudo-replication);
##    with ~19-30 tumours per dataset they have limited power.
##  * Module score = Seurat-style (signature mean minus expression-matched control genes).
##  * Melanocytic / dedifferentiated lists are curated from the literature (see
##    Supplementary_Table_13e); they are approximations, not a validated classifier.
##
## Run: source("config.R"); source("R/13_SingleCell_Validation.R")
## =====================================================================

suppressMessages(suppressWarnings({ library(Matrix) }))
set.seed(if (exists("cfg") && !is.null(cfg$seed)) cfg$seed else 20240101L)
options(timeout = 7200)
results_dir <- "results"; if (!dir.exists(results_dir)) dir.create(results_dir, recursive = TRUE)
sc_dir <- file.path("data", "sc");   if (!dir.exists(sc_dir)) dir.create(sc_dir, recursive = TRUE)
msg  <- function(...) message(sprintf(...))
has  <- function(p) requireNamespace(p, quietly = TRUE)
fmtp <- function(p) ifelse(is.na(p), "NA", ifelse(p < 1e-4, formatC(p, format = "e", digits = 1), formatC(p, format = "g", digits = 2)))
for (p in c("uwot", "irlba", "readxl")) if (!has(p)) stop("package '", p, "' is required")
source("R/cmpis_weights.R")
cm <- load_cmpis(); if (is.null(cm)) stop("CMPIS weights not found - run Stage 10 first.")
cmpis_genes <- names(cm$coef)
out <- list()

## ---------------------------------------------------------------------
## 0. Gene sets
## ---------------------------------------------------------------------
pcd_file <- if (exists("cfg") && !is.null(cfg$pcd_gene_set)) cfg$pcd_gene_set else "PCD GENES.xlsx"
pcd_all <- as.character(readxl::read_excel(pcd_file)[[1]])
pcd_all <- sort(unique(pcd_all[!is.na(pcd_all) & pcd_all != "" & pcd_all != "gene"]))

cs_top <- function(f, n = 100) {
  x <- read.delim(f, check.names = FALSE, stringsAsFactors = FALSE)
  x <- x[!is.na(x$padj) & x$padj < 0.01 & x$log2fc > 0, , drop = FALSE]
  head(x$id[order(-x$log2fc)], n)
}
tab_dir <- file.path("results", "tables")
f1 <- file.path(tab_dir, "consensusMOIC_TCGA-SKCM_limma_test_result.CS1_vs_Others.txt")
f2 <- file.path(tab_dir, "consensusMOIC_TCGA-SKCM_limma_test_result.CS2_vs_Others.txt")
if (!file.exists(f1) || !file.exists(f2)) stop("MOVICS limma tables not found in results/tables/")
cs_sig <- list(CS1 = cs_top(f1), CS2 = cs_top(f2))
cs_sig$CS1 <- setdiff(cs_sig$CS1, cs_sig$CS2); cs_sig$CS2 <- setdiff(cs_sig$CS2, cs_sig$CS1)
msg("CS signatures: CS1 %d genes, CS2 %d genes (top log2FC, FDR<0.01)", length(cs_sig$CS1), length(cs_sig$CS2))

melanocytic <- c("MITF","SOX10","MLANA","PMEL","TYR","DCT","TYRP1","MLPH","GPR143","SLC45A2","OCA2","IRF4","TRPM1","QPCT")
dedifferentiated <- c("AXL","NGFR","WNT5A","TGFBI","SERPINE1","LOXL2","EGFR","FN1","CCN1","CYR61","CCN2","CTGF",
                      "THBS1","FOSL1","JUN","SPARC","TGFB2","ITGA3","VIM")
wnt_fallback <- c("AXIN1","AXIN2","CCND2","CSNK1E","CTNNB1","CUL1","DKK1","DKK4","DLL1","DVL2","FZD1","FZD8","FRAT1",
                  "GNAI1","HDAC11","HDAC2","HDAC5","HEY1","HEY2","JAG1","JAG2","KAT2A","LEF1","MAML1","MYC","NCOR2",
                  "NCSTN","NKD1","NOTCH1","NOTCH4","NUMB","PPARD","PSEN2","PTCH1","RBPJ","SKP2","TCF7","TP53","TRIB1",
                  "WNT1","WNT5B","WNT6")
wnt <- tryCatch({
  m <- tryCatch(msigdbr::msigdbr(species = "Homo sapiens", collection = "H"),
                error = function(e) msigdbr::msigdbr(species = "Homo sapiens", category = "H"))
  m$gene_symbol[m$gs_name == "HALLMARK_WNT_BETA_CATENIN_SIGNALING"]
}, error = function(e) character(0))
if (length(wnt) < 10) wnt <- wnt_fallback
msg("WNT/beta-catenin set: %d genes", length(wnt))

## curated ligand-receptor pairs (one subunit per receptor complex)
lr <- do.call(rbind, strsplit(c(
 "CXCL9-CXCR3","CXCL10-CXCR3","CXCL11-CXCR3","CCL5-CCR5","CCL4-CCR5","CCL3-CCR5","CCL2-CCR2","CCL19-CCR7","CCL21-CCR7",
 "CXCL12-CXCR4","CXCL13-CXCR5","CXCL8-CXCR2","CXCL16-CXCR6","CX3CL1-CX3CR1",
 "CD274-PDCD1","PDCD1LG2-PDCD1","CD80-CTLA4","CD86-CTLA4","CD80-CD28","CD86-CD28","PVR-TIGIT","NECTIN2-TIGIT","LGALS9-HAVCR2",
 "CD70-CD27","TNFSF9-TNFRSF9","TNFSF4-TNFRSF4","HLA-E-KLRC1","HLA-E-KLRD1","CD40LG-CD40",
 "IFNG-IFNGR1","IFNG-IFNGR2","TNF-TNFRSF1A","TGFB1-TGFBR1","TGFB1-TGFBR2","TGFB2-TGFBR2","IL6-IL6R","IL10-IL10RA",
 "VEGFA-KDR","VEGFA-FLT1","CSF1-CSF1R","IL34-CSF1R","MIF-CD74","MIF-CXCR4","MDK-LRP1","SPP1-CD44","SPP1-ITGAV",
 "WNT5A-FZD5","WNT5A-ROR2","FASLG-FAS","TNFSF10-TNFRSF10B","TNFSF10-TNFRSF10A","ICAM1-ITGB2","ICAM1-ITGAL",
 "VCAM1-ITGA4","VCAM1-ITGB1","GAS6-AXL","PROS1-AXL","HGF-MET","PDGFB-PDGFRB","ANXA1-FPR1","LGALS3-LAG3"), "-(?=[^-]+$)", perl = TRUE))
lr <- data.frame(ligand = lr[, 1], receptor = lr[, 2], stringsAsFactors = FALSE)
lr$pair <- paste0(lr$ligand, "-", lr$receptor)

marker_sets <- list(
  T = c("CD3D","CD3E","CD3G","CD2","TRAC"), NK = c("KLRF1","NCAM1","KLRC1","FCGR3A","KLRB1","TYROBP"),
  B = c("MS4A1","CD79A","CD79B","CD19","CD37"), Plasma = c("MZB1","SDC1","JCHAIN","TNFRSF17","CD38"),
  Myeloid = c("CD14","LYZ","CD68","CSF1R","AIF1","LST1","CD163","FCGR3A"))

gene_union <- unique(c(pcd_all, cmpis_genes, cs_sig$CS1, cs_sig$CS2, melanocytic, dedifferentiated, wnt,
                       lr$ligand, lr$receptor, unlist(marker_sets), "CD8A","CD8B","CD4","IL7R","FOXP3","PTPRC"))

## ---------------------------------------------------------------------
## 1. Download
## ---------------------------------------------------------------------
fetch_geo <- function(acc, patterns) {
  dir_url <- sprintf("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE%snnn/%s/suppl/", substr(acc, 4, nchar(acc) - 3), acc)
  idx <- tryCatch(readLines(dir_url, warn = FALSE), error = function(e) character(0))
  files <- unique(unlist(regmatches(idx, gregexpr('(?<=href=")[^"?/]+(?=")', idx, perl = TRUE))))
  hit <- NA_character_
  for (p in patterns) { h <- grep(p, files, value = TRUE, ignore.case = TRUE); if (length(h)) { hit <- h[1]; break } }
  if (is.na(hit)) { msg("  [%s] no file matching %s in %s (found: %s)", acc, paste(patterns, collapse = " | "), dir_url,
                        paste(head(files, 10), collapse = ", ")); return(NA_character_) }
  dest <- file.path(sc_dir, hit)
  if (file.exists(dest) && file.size(dest) > 1e4) { msg("  have %s", hit); return(dest) }
  msg("  downloading %s ...", hit)
  ok <- tryCatch({ download.file(paste0(dir_url, hit), dest, mode = "wb", quiet = FALSE); TRUE }, error = function(e) { msg("  FAILED: %s", conditionMessage(e)); FALSE })
  if (ok) dest else NA_character_
}
msg("Downloading / locating single-cell files in %s", sc_dir)
files <- list(
  tirosh  = fetch_geo("GSE72056",  c("revised_v2.*\\.txt\\.gz$", "single_cell.*\\.txt\\.gz$")),
  jerby   = fetch_geo("GSE115978", c("tpm.*\\.csv\\.gz$")),
  jerby_a = fetch_geo("GSE115978", c("annotations.*\\.csv\\.gz$")),
  sade    = fetch_geo("GSE120575", c("TPM.*\\.txt\\.gz$")),
  sade_p  = fetch_geo("GSE120575", c("patient.*\\.txt\\.gz$")))

## ---------------------------------------------------------------------
## 2. Streaming reader: gene-by-cell text -> sparse matrix (optionally only some genes)
## ---------------------------------------------------------------------
read_sc <- function(file, sep = "\t", keep = NULL, meta_pat = NULL) {
  con <- gzfile(file, "rt"); on.exit(close(con))
  unq <- function(s) gsub('^"|"$', "", s)
  hdr <- unq(strsplit(readLines(con, n = 1L), sep, fixed = TRUE)[[1]])
  I <- list(); J <- list(); V <- list(); genes <- character(0); meta <- list(); ncell <- NA_integer_; nline <- 0L
  repeat {
    ln <- readLines(con, n = 500L); if (!length(ln)) break
    for (l in ln) {
      nline <- nline + 1L
      p <- regexpr(sep, l, fixed = TRUE); if (p < 1) next
      g <- unq(substr(l, 1, p - 1))
      if (!is.null(keep) && nline > 8L && !(g %in% keep)) next
      if (g %in% genes) next
      f <- strsplit(l, sep, fixed = TRUE)[[1]]
      if (!is.null(meta_pat) && grepl(meta_pat, g, ignore.case = TRUE)) { meta[[g]] <- unq(f[-1]); next }   # numeric annotation rows (Tirosh)
      v <- suppressWarnings(as.numeric(f[-1]))
      if (mean(is.na(v)) > 0.5) { meta[[g]] <- unq(f[-1]); next }              # annotation row
      if (is.na(ncell)) ncell <- length(v)
      if (length(v) != ncell) next
      v[is.na(v)] <- 0
      nz <- which(v != 0); genes <- c(genes, g); k <- length(genes)
      I[[k]] <- rep.int(k, length(nz)); J[[k]] <- nz; V[[k]] <- v[nz]
    }
  }
  cells <- utils::tail(hdr, ncell)
  X <- sparseMatrix(i = unlist(I), j = unlist(J), x = unlist(V), dims = c(length(genes), ncell), dimnames = list(genes, cells))
  if (max(X@x) > 50) X@x <- log2(X@x + 1)
  list(X = X, meta = meta)
}
find_col <- function(df, pats) { for (p in pats) { h <- grep(p, names(df), ignore.case = TRUE); if (length(h)) return(h[1]) }; NA_integer_ }

## ---------------------------------------------------------------------
## 3. Module score (expression-matched control genes) and CD8/CD4 split
## ---------------------------------------------------------------------
make_scorer <- function(X, nbin = 24, ctrl = 5) {
  avg <- Matrix::rowMeans(X)
  bins <- cut(rank(avg, ties.method = "first"), nbin, labels = FALSE); names(bins) <- rownames(X)
  pool <- split(rownames(X), bins)
  function(genes) {
    genes <- intersect(genes, rownames(X)); if (length(genes) < 3) return(rep(NA_real_, ncol(X)))
    cg <- unique(unlist(lapply(genes, function(g) {
      p <- setdiff(pool[[as.character(bins[g])]], genes); if (!length(p)) return(character(0)); sample(p, min(ctrl, length(p))) })))
    s <- Matrix::colMeans(X[genes, , drop = FALSE])
    if (length(cg)) s <- s - Matrix::colMeans(X[cg, , drop = FALSE])
    as.numeric(s)
  }
}
split_t <- function(X, is_t) {                 # CD8 / CD4 / other among T cells, by marker expression
  g <- function(x) if (x %in% rownames(X)) as.numeric(X[x, ]) else rep(0, ncol(X))
  cd8 <- pmax(g("CD8A"), g("CD8B")); cd4 <- pmax(g("CD4"), 0.5 * g("IL7R"))
  ifelse(!is_t, NA, ifelse(cd8 > 0 & cd8 > cd4, "CD8 T cell", ifelse(cd4 > 0, "CD4 T cell", "T cell (other)")))
}

## ---------------------------------------------------------------------
## 4. Dataset loaders -> list(X = sparse genes x cells, meta = data.frame)
## ---------------------------------------------------------------------
load_tirosh <- function(f) {
  r <- read_sc(f, "\t", keep = NULL, meta_pat = "^(tumor|malignant|non-malignant)"); m <- r$meta; X <- r$X; n <- ncol(X)
  gv <- function(pat) { h <- grep(pat, names(m), ignore.case = TRUE)[1]; if (is.na(h)) stop("Tirosh metadata row '", pat, "' not found; rows: ", paste(names(m), collapse = " | ")); as.numeric(m[[h]]) }
  tum <- gv("^tumor"); mal <- gv("^malignant"); nmt <- gv("^non-malignant")
  lab <- c("T cell","B cell","Macrophage","Endothelial","CAF","NK")
  type <- ifelse(mal == 2, "Malignant", lab[ifelse(nmt %in% 1:6, nmt, NA)])
  type[mal == 0 | (mal == 1 & !(nmt %in% 1:6))] <- NA
  meta <- data.frame(cell = colnames(X), tumour = paste0("Mel", tum), type = type, stringsAsFactors = FALSE)
  keep <- !is.na(meta$type); X <- X[, keep]; meta <- meta[keep, ]
  tt <- split_t(X, meta$type == "T cell"); meta$type[meta$type == "T cell"] <- tt[meta$type == "T cell"]
  list(X = X, meta = meta)
}
load_jerby <- function(f, fa) {
  ann <- read.csv(gzfile(fa), check.names = FALSE, stringsAsFactors = FALSE)
  r <- read_sc(f, ",", keep = NULL); X <- r$X
  ic <- find_col(ann, c("^cells?$", "^cell")); if (is.na(ic)) ic <- 1
  tc <- find_col(ann, c("cell.?type")); sc <- find_col(ann, c("^samples?$", "sample"))
  if (is.na(tc) || is.na(sc)) stop("Jerby-Arnon annotation columns not found: ", paste(names(ann), collapse = ", "))
  ann <- ann[match(colnames(X), ann[[ic]]), ]
  if (mean(is.na(ann[[ic]])) > 0.1) stop("Jerby-Arnon: cell ids in the matrix do not match the annotation file")
  tp <- c(`B.cell` = "B cell", CAF = "CAF", `Endo.` = "Endothelial", Macrophage = "Macrophage", Malignant = "Malignant",
          NK = "NK", `T.CD4` = "CD4 T cell", `T.CD8` = "CD8 T cell", `T.cell` = "T cell (other)", Mal = "Malignant")
  ty <- unname(tp[as.character(ann[[tc]])]); ty[is.na(ty)] <- as.character(ann[[tc]])[is.na(ty)]
  ty[ty %in% c("?", "")] <- NA
  meta <- data.frame(cell = colnames(X), tumour = as.character(ann[[sc]]), type = ty, stringsAsFactors = FALSE)
  keep <- !is.na(meta$type) & !is.na(meta$tumour)
  cc <- find_col(ann, c("^cohort$"))                       # the GSE115978 file also contains the Tirosh-2016 cells
  if (!is.na(cc)) { own <- !grepl("Tirosh", as.character(ann[[cc]]), ignore.case = TRUE)
    msg("  Jerby-Arnon file: dropping %d cells that belong to the Tirosh cohort (analysed separately)", sum(keep & !own)); keep <- keep & own }
  list(X = X[, keep], meta = meta[keep, ])
}
load_sade <- function(f, fp) {
  r <- read_sc(f, "\t", keep = gene_union); X <- r$X
  pf <- iconv(readLines(gzfile(fp), warn = FALSE), "latin1", "UTF-8", sub = ""); k <- grep("^\"?Sample name", pf, useBytes = TRUE)[1]
  if (is.na(k)) stop("Sade-Feldman patient file: header not found")
  pf2 <- pf[k:length(pf)]; first <- sub("\t.*", "", pf2); stp <- which(first == "" | grepl("^(PROTOCOLS|#)", first))[1]
  if (!is.na(stp)) pf2 <- pf2[seq_len(stp - 1)]
  pd <- read.delim(text = pf2, check.names = FALSE, stringsAsFactors = FALSE, quote = "")
  idc <- NA
  for (j in seq_len(min(3, ncol(pd)))) if (mean(colnames(X) %in% pd[[j]]) > 0.5) { idc <- j; break }
  if (is.na(idc)) stop("Sade-Feldman: cell ids not matched in patient file")
  pd <- pd[match(colnames(X), pd[[idc]]), ]
  pc <- find_col(pd, c("pati")); rc <- find_col(pd, c("response")); tc <- find_col(pd, c("therapy"))
  meta <- data.frame(cell = colnames(X), sample = as.character(pd[[pc]]),
                     response = as.character(pd[[rc]]), therapy = if (is.na(tc)) NA else as.character(pd[[tc]]), stringsAsFactors = FALSE)
  meta$timepoint <- ifelse(grepl("^pre", meta$sample, ignore.case = TRUE), "Pre", "Post")
  meta$tumour <- meta$sample
  sc <- make_scorer(X)
  S <- sapply(marker_sets, sc)
  tt <- split_t(X, rep(TRUE, ncol(X)))
  ty <- rep("Other", ncol(X))
  S[is.na(S)] <- 0
  Z <- scale(S); Z[is.na(Z)] <- 0                              # lineage scores on a common scale
  best <- colnames(Z)[max.col(Z, ties.method = "first")]; bz <- apply(Z, 1, max)
  lin <- c(T = NA, NK = "NK", B = "B cell", Plasma = "Plasma cell", Myeloid = "Macrophage")
  ty[bz > 0.75] <- ifelse(best[bz > 0.75] == "T", tt[bz > 0.75], lin[best[bz > 0.75]])
  meta$type <- ty
  list(X = X, meta = meta)
}

## ---------------------------------------------------------------------
## 5. Per-dataset scoring + UMAP
## ---------------------------------------------------------------------
embed_umap <- function(X, ngenes = 2000, npc = 30) {
  keep <- rownames(X)[Matrix::rowMeans(X > 0) >= 0.03]
  Xk <- X[keep, , drop = FALSE]
  m <- Matrix::rowMeans(Xk); v <- Matrix::rowMeans(Xk^2) - m^2; disp <- log(v / pmax(m, 1e-6) + 1e-6)
  bin <- cut(rank(m, ties.method = "first"), 20, labels = FALSE)
  z <- ave(disp, bin, FUN = function(x) if (length(x) > 1 && sd(x) > 0) (x - mean(x)) / sd(x) else rep(0, length(x)))
  hv <- keep[order(z, decreasing = TRUE)[seq_len(min(ngenes, length(keep)))]]
  Z <- t(as.matrix(X[hv, , drop = FALSE])); Z <- scale(Z); Z[is.na(Z)] <- 0; Z[Z > 10] <- 10; Z[Z < -10] <- -10
  npc <- min(npc, ncol(Z) - 1, nrow(Z) - 1)
  s <- irlba::irlba(Z, nv = npc); pcs <- s$u %*% diag(s$d, npc)
  um <- uwot::umap(pcs, n_neighbors = 30, min_dist = 0.3, verbose = FALSE)
  colnames(um) <- c("UMAP1", "UMAP2"); um
}
prep_dataset <- function(name, d) {
  X <- d$X; meta <- d$meta
  msg("[%s] %d genes x %d cells; types: %s", name, nrow(X), ncol(X),
      paste(names(table(meta$type)), table(meta$type), sep = "=", collapse = ", "))
  sc <- make_scorer(X)
  meta$PCD <- sc(pcd_all); meta$CS1 <- sc(cs_sig$CS1); meta$CS2 <- sc(cs_sig$CS2); meta$CSdiff <- meta$CS2 - meta$CS1
  meta$Melanocytic <- sc(melanocytic); meta$Dedifferentiated <- sc(dedifferentiated)
  meta$WNT <- sc(wnt)
  g <- intersect(cmpis_genes, rownames(X))
  meta$CMPIS <- as.numeric(score_cmpis(as.matrix(X[g, , drop = FALSE]), cm))
  mal <- meta$type == "Malignant"
  meta$CSdiff_int <- NA_real_
  if (sum(mal) >= 50) {
    det <- Matrix::rowMeans(X[, mal, drop = FALSE] > 0); mn <- Matrix::rowMeans(X[, mal, drop = FALSE])
    ok <- names(det)[det >= 0.10 & mn > 0.1]
    i1 <- intersect(cs_sig$CS1, ok); i2 <- intersect(cs_sig$CS2, ok)
    msg("[%s] tumour-intrinsic signatures: CS1 %d genes, CS2 %d genes", name, length(i1), length(i2))
    if (length(i1) >= 10 && length(i2) >= 10) meta$CSdiff_int <- sc(i2) - sc(i1)
    d$intrinsic <- list(CS1 = i1, CS2 = i2)
  }
  d$umap <- tryCatch(embed_umap(X), error = function(e) { msg("[%s] UMAP failed: %s", name, conditionMessage(e)); NULL })
  d$meta <- meta; d$name <- name; d$cmpis_expr <- as.matrix(X[g, , drop = FALSE]); d
}

## ---------------------------------------------------------------------
## 6. Tables: CMPIS genes and PCD score by cell type
## ---------------------------------------------------------------------
expr_by_type <- function(d) {
  meta <- d$meta; types <- sort(unique(meta$type)); E <- d$cmpis_expr
  do.call(rbind, lapply(types, function(tp) {
    k <- meta$type == tp
    data.frame(dataset = d$name, cell_type = tp, n_cells = sum(k), gene = rownames(E),
               mean_expr = rowMeans(E[, k, drop = FALSE]), pct_detected = 100 * rowMeans(E[, k, drop = FALSE] > 0),
               row.names = NULL)
  }))
}
pcd_by_type <- function(d) {
  meta <- d$meta
  do.call(rbind, lapply(sort(unique(meta$type)), function(tp) { k <- meta$type == tp
    data.frame(dataset = d$name, cell_type = tp, n_cells = sum(k), PCD_mean = mean(meta$PCD[k]), PCD_median = median(meta$PCD[k]),
               CMPIS_mean = mean(meta$CMPIS[k]), CSdiff_mean = mean(meta$CSdiff[k])) })) }

## ---------------------------------------------------------------------
## 7. Tumour-level tests (datasets with malignant cells)
## ---------------------------------------------------------------------
min_mal <- 15; min_nonmal <- 20
tumour_table <- function(d) {
  meta <- d$meta; mal_t <- table(meta$tumour[meta$type == "Malignant"])
  tums <- names(mal_t)[mal_t >= min_mal]
  do.call(rbind, lapply(tums, function(t) {
    mal <- meta$tumour == t & meta$type == "Malignant"; nm <- meta$tumour == t & meta$type != "Malignant"
    nn <- sum(nm)
    data.frame(dataset = d$name, tumour = t, n_malignant = sum(mal), n_nonmalignant = nn,
               CSdiff_int = mean(meta$CSdiff_int[mal]), CMPIS = mean(meta$CMPIS[mal]), PCD = mean(meta$PCD[mal]),
               WNT = mean(meta$WNT[mal]), Melanocytic = mean(meta$Melanocytic[mal]), Dedifferentiated = mean(meta$Dedifferentiated[mal]),
               DedifMinusMelan = mean(meta$Dedifferentiated[mal] - meta$Melanocytic[mal]),
               CD8_frac = if (nn >= min_nonmal) sum(nm & meta$type == "CD8 T cell") / nn else NA,
               Tcell_frac = if (nn >= min_nonmal) sum(nm & grepl("T cell", meta$type)) / nn else NA)
  }))
}
spear <- function(x, y) {
  k <- is.finite(x) & is.finite(y); if (sum(k) < 5) return(c(n = sum(k), rho = NA, p = NA))
  ct <- suppressWarnings(cor.test(x[k], y[k], method = "spearman", exact = FALSE)); c(n = sum(k), rho = unname(ct$estimate), p = ct$p.value)
}
tumour_tests <- function(tt) {
  pairs <- list(c("CSdiff_int","Dedifferentiated"), c("CSdiff_int","Melanocytic"), c("CSdiff_int","DedifMinusMelan"),
                c("CSdiff_int","CD8_frac"), c("CSdiff_int","Tcell_frac"), c("CSdiff_int","WNT"), c("CSdiff_int","PCD"),
                c("CMPIS","CD8_frac"), c("CMPIS","Tcell_frac"), c("CMPIS","Dedifferentiated"), c("CMPIS","Melanocytic"),
                c("WNT","PCD"), c("WNT","CD8_frac"), c("PCD","CD8_frac"), c("CSdiff_int","CMPIS"))
  do.call(rbind, lapply(pairs, function(pp) { r <- spear(tt[[pp[1]]], tt[[pp[2]]])
    data.frame(dataset = tt$dataset[1], level = "tumour", x = pp[1], y = pp[2], n = r["n"], rho = r["rho"], p = r["p"], row.names = NULL) }))
}
cell_tests <- function(d) {                    # within-tumour correlations among malignant cells
  meta <- d$meta; mal <- meta[meta$type == "Malignant", ]
  tums <- names(which(table(mal$tumour) >= 30)); if (!length(tums)) return(NULL)
  pr <- list(c("CSdiff_int","Dedifferentiated"), c("CSdiff_int","Melanocytic"), c("WNT","PCD"), c("CMPIS","PCD"), c("CMPIS","Dedifferentiated"))
  do.call(rbind, lapply(pr, function(pp) {
    rs <- sapply(tums, function(t) { z <- mal[mal$tumour == t, ]; spear(z[[pp[1]]], z[[pp[2]]])["rho"] }); rs <- rs[is.finite(rs)]
    wp <- if (length(rs) >= 5) suppressWarnings(wilcox.test(rs, mu = 0)$p.value) else NA
    data.frame(dataset = d$name, level = "within-tumour (malignant cells)", x = pp[1], y = pp[2], n = length(rs),
               rho = if (length(rs)) median(rs) else NA, p = wp, row.names = NULL) }))
}

## ---------------------------------------------------------------------
## 8. Cell-cell communication (curated L-R, tumour pseudo-bulk)
## ---------------------------------------------------------------------
broad_class <- function(ty) ifelse(grepl("T cell", ty), "T cell", ty)
comm_tests <- function(d, tt) {
  meta <- d$meta; X <- d$X; lg <- intersect(unique(c(lr$ligand, lr$receptor)), rownames(X))
  L <- as.matrix(X[lg, , drop = FALSE]); cls <- broad_class(meta$type)
  partners <- intersect(c("T cell","NK","Macrophage","CAF","Endothelial","B cell"), unique(cls))
  pbm <- function(t, cl, minc) { k <- meta$tumour == t & cls == cl; if (sum(k) < minc) return(NULL); rowMeans(L[, k, drop = FALSE]) }
  rows <- list(); strength <- list()
  for (P in partners) for (dir in c("out", "in")) {
    S <- matrix(NA_real_, nrow(lr), nrow(tt)); colnames(S) <- tt$tumour
    for (j in seq_len(nrow(tt))) {
      m <- pbm(tt$tumour[j], "Malignant", min_mal); p <- pbm(tt$tumour[j], P, 10); if (is.null(m) || is.null(p)) next
      s <- if (dir == "out") m[lr$ligand] * p[lr$receptor] else p[lr$ligand] * m[lr$receptor]; S[, j] <- s
    }
    S[!(lr$ligand %in% lg & lr$receptor %in% lg), ] <- NA
    tot <- colSums(S, na.rm = TRUE); tot[colSums(!is.na(S)) == 0] <- NA
    r <- spear(tt$CSdiff_int[match(colnames(S), tt$tumour)], tot)
    strength[[length(strength) + 1]] <- data.frame(dataset = d$name, direction = ifelse(dir == "out", "malignant -> partner", "partner -> malignant"),
                                                   partner = P, pair = "ALL (summed)", n = r["n"], rho = r["rho"], p = r["p"], row.names = NULL)
    for (i in seq_len(nrow(lr))) { if (all(is.na(S[i, ]))) next
      r <- spear(tt$CSdiff_int[match(colnames(S), tt$tumour)], S[i, ])
      rows[[length(rows) + 1]] <- data.frame(dataset = d$name, direction = ifelse(dir == "out", "malignant -> partner", "partner -> malignant"),
                                            partner = P, pair = lr$pair[i], n = r["n"], rho = r["rho"], p = r["p"], row.names = NULL) }
  }
  list(strength = do.call(rbind, strength), pairs = if (length(rows)) do.call(rbind, rows) else NULL)
}

## ---------------------------------------------------------------------
## 9. Sade-Feldman: immune-cell CMPIS / PCD vs clinical response (pre-treatment)
## ---------------------------------------------------------------------
sade_response <- function(d) {
  meta <- d$meta; pre <- meta[meta$timepoint == "Pre" & meta$response %in% c("Responder", "Non-responder"), ]
  do.call(rbind, lapply(c("CD8 T cell", "CD4 T cell", "Macrophage", "NK", "B cell"), function(tp) {
    z <- pre[pre$type == tp, ]; ns <- table(z$sample); ok <- names(ns)[ns >= 10]; if (length(ok) < 6) return(NULL)
    agg <- do.call(rbind, lapply(ok, function(s) { k <- z[z$sample == s, ]
      data.frame(sample = s, response = k$response[1], CMPIS = mean(k$CMPIS), PCD = mean(k$PCD), CSdiff = mean(k$CSdiff)) }))
    do.call(rbind, lapply(c("CMPIS", "PCD", "CSdiff"), function(v) {
      a <- agg[[v]][agg$response == "Responder"]; b <- agg[[v]][agg$response == "Non-responder"]
      p <- if (length(a) >= 3 && length(b) >= 3) suppressWarnings(wilcox.test(a, b)$p.value) else NA
      data.frame(cell_type = tp, score = v, n_resp = length(a), n_nonresp = length(b),
                 median_resp = median(a), median_nonresp = median(b), wilcox_p = p) }))
  }))
}

## ---------------------------------------------------------------------
## 10. TCGA bulk sanity check: is CS2 melanocytic or dedifferentiated?
## ---------------------------------------------------------------------
tcga_check <- function() {
  if (!file.exists("skcm.tcga.RData") || !file.exists("TCGA_Clustering.RData")) return(NULL)
  load("skcm.tcga.RData"); e <- new.env(); load("TCGA_Clustering.RData", envir = e)
  expr <- as.matrix(skcm.tcga$mRNA.expr); storage.mode(expr) <- "double"; if (max(expr, na.rm = TRUE) > 50) expr <- log2(expr + 1)
  expr <- expr[, !duplicated(colnames(expr)), drop = FALSE]
  cr <- get("cmoic.skcm", envir = e)$clust.res; CS <- setNames(ifelse(cr$clust == 1, "CS1", "CS2"), cr$samID)
  sm <- intersect(colnames(expr), names(CS)); cs <- CS[sm]
  zm <- function(genes) { g <- intersect(genes, rownames(expr)); if (length(g) < 3) return(rep(NA, length(sm)))
    z <- t(scale(t(expr[g, sm, drop = FALSE]))); z[is.na(z)] <- 0; colMeans(z) }
  sc <- list(Melanocytic = zm(melanocytic), Dedifferentiated = zm(dedifferentiated), WNT_beta_catenin = zm(wnt),
             PCD_score = zm(pcd_all), CMPIS = as.numeric(score_cmpis(expr[, sm], cm)))
  tab <- do.call(rbind, lapply(names(sc), function(n) { x <- sc[[n]]
    data.frame(score = n, n_CS1 = sum(cs == "CS1"), n_CS2 = sum(cs == "CS2"),
               median_CS1 = median(x[cs == "CS1"]), median_CS2 = median(x[cs == "CS2"]),
               higher_in = ifelse(median(x[cs == "CS2"]) > median(x[cs == "CS1"]), "CS2", "CS1"),
               wilcox_p = wilcox.test(x[cs == "CS2"], x[cs == "CS1"])$p.value) }))
  list(table = tab, scores = sc, cs = cs)
}

## ---------------------------------------------------------------------
## 11. Run
## ---------------------------------------------------------------------
loaders <- list(
  Tirosh = function() load_tirosh(files$tirosh),
  `Jerby-Arnon` = function() load_jerby(files$jerby, files$jerby_a),
  `Sade-Feldman` = function() load_sade(files$sade, files$sade_p))
req <- list(Tirosh = "tirosh", `Jerby-Arnon` = c("jerby", "jerby_a"), `Sade-Feldman` = c("sade", "sade_p"))
ds <- list()
for (nm in names(loaders)) {
  if (any(is.na(unlist(files[req[[nm]]])))) { msg("== %s: file(s) missing, skipped", nm); next }
  msg("== %s ==", nm)
  r <- tryCatch(prep_dataset(nm, loaders[[nm]]()), error = function(e) { msg("   %s FAILED: %s", nm, conditionMessage(e)); NULL })
  if (!is.null(r)) ds[[nm]] <- r
}
if (!length(ds)) stop("no single-cell dataset could be processed")

t13  <- do.call(rbind, lapply(ds, expr_by_type)); t13b <- do.call(rbind, lapply(ds, pcd_by_type))
write.csv(t13,  file.path(results_dir, "Supplementary_Table_13_sc_CMPIS_expression.csv"), row.names = FALSE)
write.csv(t13b, file.path(results_dir, "Supplementary_Table_13b_sc_PCD_by_celltype.csv"), row.names = FALSE)

tum <- list(); tests <- list(); comm_s <- list(); comm_p <- list()
for (nm in names(ds)) {
  d <- ds[[nm]]
  if (!any(d$meta$type == "Malignant") || all(is.na(d$meta$CSdiff_int))) next
  tt <- tumour_table(d); tt <- tt[is.finite(tt$CSdiff_int), ]
  msg("[%s] tumour-level: %d tumours with >= %d malignant cells", nm, nrow(tt), min_mal)
  if (nrow(tt) < 5) next
  tum[[nm]] <- tt; tests[[nm]] <- rbind(tumour_tests(tt), cell_tests(d))
  cc <- tryCatch(comm_tests(d, tt), error = function(e) { msg("[%s] communication failed: %s", nm, conditionMessage(e)); NULL })
  if (!is.null(cc)) { comm_s[[nm]] <- cc$strength; comm_p[[nm]] <- cc$pairs }
}
t13c <- if (length(tests)) do.call(rbind, tests) else NULL
if (!is.null(t13c)) { t13c$fdr <- NA_real_
  for (g in unique(t13c$level)) { k <- t13c$level == g; t13c$fdr[k] <- p.adjust(t13c$p[k], "BH") }
  write.csv(t13c, file.path(results_dir, "Supplementary_Table_13c_sc_tumour_tests.csv"), row.names = FALSE)
  write.csv(do.call(rbind, tum), file.path(results_dir, "Supplementary_Table_13c2_sc_tumour_level_data.csv"), row.names = FALSE) }
cs_all <- if (length(comm_s)) do.call(rbind, comm_s) else NULL
cp_all <- if (length(comm_p)) do.call(rbind, comm_p) else NULL
if (!is.null(cp_all)) { cp_all$fdr <- ave(cp_all$p, cp_all$dataset, cp_all$direction, FUN = function(x) p.adjust(x, "BH"))
  write.csv(rbind(cbind(cs_all, fdr = NA), cp_all[, c(names(cs_all), "fdr")]),
            file.path(results_dir, "Supplementary_Table_13d_sc_communication.csv"), row.names = FALSE) }
sets_tab <- data.frame(set = rep(c("CS1 (TCGA top100)", "CS2 (TCGA top100)", "Melanocytic", "Dedifferentiated", "WNT/beta-catenin"),
                                 c(length(cs_sig$CS1), length(cs_sig$CS2), length(melanocytic), length(dedifferentiated), length(wnt))),
                       gene = c(cs_sig$CS1, cs_sig$CS2, melanocytic, dedifferentiated, wnt))
write.csv(sets_tab, file.path(results_dir, "Supplementary_Table_13e_signature_genes.csv"), row.names = FALSE)
resp <- if ("Sade-Feldman" %in% names(ds)) tryCatch(sade_response(ds[["Sade-Feldman"]]), error = function(e) { msg("response analysis failed: %s", conditionMessage(e)); NULL }) else NULL
if (!is.null(resp)) write.csv(resp, file.path(results_dir, "Supplementary_Table_13f_sc_response.csv"), row.names = FALSE)
tc <- tryCatch(tcga_check(), error = function(e) { msg("TCGA check failed: %s", conditionMessage(e)); NULL })
if (!is.null(tc)) write.csv(tc$table, file.path(results_dir, "Supplementary_Table_13g_TCGA_programme_check.csv"), row.names = FALSE)

## ---------------------------------------------------------------------
## 12. Figures
## ---------------------------------------------------------------------
type_cols <- c(Malignant = "#B2182B", `CD8 T cell` = "#2166AC", `CD4 T cell` = "#67A9CF", `T cell (other)` = "#92C5DE", NK = "#1B7837",
               `B cell` = "#F4A582", `Plasma cell` = "#FDB863", Macrophage = "#7B3294", CAF = "#8C6D31", Endothelial = "#636363", Other = "grey80")
tcol <- function(ty) { z <- type_cols[ty]; z[is.na(z)] <- "grey70"; unname(z) }
pal <- if (has("viridis")) viridis::viridis(100) else heat.colors(100)
umap_cont <- function(um, v, main) {
  lo <- quantile(v, 0.02, na.rm = TRUE); hi <- quantile(v, 0.98, na.rm = TRUE); vv <- pmin(pmax(v, lo), hi)
  o <- order(vv); col <- pal[cut(vv, seq(lo, hi, length.out = 101), include.lowest = TRUE, labels = FALSE)]
  col[is.na(col)] <- "grey80"
  plot(um[o, 1], um[o, 2], col = col[o], pch = 16, cex = 0.35, axes = FALSE, xlab = "UMAP1", ylab = "UMAP2", main = main, cex.main = 0.9); box()
}
dot_panel <- function(d, types_show) {
  e <- t13[t13$dataset == d$name & t13$cell_type %in% types_show, ]
  tps <- intersect(types_show, unique(e$cell_type)); gs <- cmpis_genes[cmpis_genes %in% e$gene]
  mm <- sapply(tps, function(tp) sapply(gs, function(g) e$mean_expr[e$cell_type == tp & e$gene == g]))
  pp <- sapply(tps, function(tp) sapply(gs, function(g) e$pct_detected[e$cell_type == tp & e$gene == g]))
  z <- t(scale(t(mm))); z[is.na(z)] <- 0
  par(mar = c(7, 6, 3, 4))
  plot(NA, xlim = c(0.5, length(tps) + 0.5), ylim = c(0.5, length(gs) + 0.5), axes = FALSE, xlab = "", ylab = "",
       main = "CMPIS genes (colour: z-scored mean; size: % detected)", cex.main = 0.8)
  rb <- colorRampPalette(c("#2166AC", "white", "#B2182B"))(101)
  for (i in seq_along(gs)) for (j in seq_along(tps))
    points(j, length(gs) + 1 - i, pch = 21, bg = rb[round((max(min(z[i, j], 2), -2) + 2) / 4 * 100) + 1], col = "grey30",
           cex = 0.4 + 2.2 * sqrt(pp[i, j] / 100))
  axis(1, at = seq_along(tps), labels = tps, las = 2, cex.axis = 0.75); axis(2, at = rev(seq_along(gs)), labels = gs, las = 1, cex.axis = 0.8)
  w <- cm$coef[gs]; mtext(sprintf("%+.2f", w), side = 4, at = rev(seq_along(gs)), las = 1, cex = 0.6, line = 0.3)
}
page_dataset <- function(d) {
  meta <- d$meta; um <- d$umap; if (is.null(um)) return(invisible())
  par(mfrow = c(2, 3), mar = c(4, 4, 3, 1))
  ty <- meta$type; plot(um[, 1], um[, 2], col = tcol(ty), pch = 16, cex = 0.35, axes = FALSE, xlab = "UMAP1", ylab = "UMAP2",
                        main = sprintf("%s: cell type (n = %d)", d$name, nrow(meta)), cex.main = 0.9); box()
  tl <- names(sort(table(ty), decreasing = TRUE)); legend("topright", legend = tl, pch = 16, col = tcol(tl), cex = 0.6, bg = "#FFFFFFCC", box.col = NA)
  umap_cont(um, meta$CMPIS, "CMPIS (corrected weights)")
  umap_cont(um, meta$PCD, "PCD score (2,682-gene set)")
  umap_cont(um, meta$CSdiff, "CS2-like minus CS1-like (full signatures)")
  dot_panel(d, intersect(c("Malignant","CD8 T cell","CD4 T cell","T cell (other)","NK","Macrophage","CAF","Endothelial","B cell"), unique(ty)))
  par(mar = c(8, 5, 3, 1)); ord <- names(sort(tapply(meta$PCD, ty, median), decreasing = TRUE))
  boxplot(meta$PCD ~ factor(ty, ord), las = 2, col = tcol(ord), outline = FALSE, ylab = "PCD score", xlab = "", main = "PCD score by cell type", cex.axis = 0.7)
  abline(h = 0, lty = 3)
}
scatter <- function(tt, x, y, xl, yl) {
  r <- spear(tt[[x]], tt[[y]])
  if (sum(is.finite(tt[[x]]) & is.finite(tt[[y]])) < 3) { plot.new(); text(0.5, 0.5, paste("no data:", y)); return(invisible()) }
  plot(tt[[x]], tt[[y]], pch = 19, col = "#B2182B", xlab = xl, ylab = yl,
       main = sprintf("rho = %.2f, p = %s (n = %d)", r["rho"], fmtp(r["p"]), r["n"]), cex.main = 0.85)
  if (is.finite(r["rho"])) abline(lm(tt[[y]] ~ tt[[x]]), lty = 2, col = "grey40")
}
page_tumour <- function(nm) {
  tt <- tum[[nm]]; par(mfrow = c(2, 3), mar = c(4.5, 4.5, 3, 1), oma = c(0, 0, 2, 0))
  scatter(tt, "CSdiff_int", "Dedifferentiated", "Tumour-intrinsic CS2 - CS1", "Dedifferentiated programme")
  scatter(tt, "CSdiff_int", "Melanocytic", "Tumour-intrinsic CS2 - CS1", "Melanocytic programme")
  scatter(tt, "CSdiff_int", "CD8_frac", "Tumour-intrinsic CS2 - CS1", "CD8 T-cell fraction of non-malignant")
  scatter(tt, "CSdiff_int", "WNT", "Tumour-intrinsic CS2 - CS1", "WNT/beta-catenin score")
  scatter(tt, "WNT", "PCD", "WNT/beta-catenin score", "PCD score (malignant)")
  scatter(tt, "CMPIS", "CD8_frac", "CMPIS (malignant mean)", "CD8 T-cell fraction")
  mtext(sprintf("%s: one point = one tumour (>= %d malignant cells)", nm, min_mal), outer = TRUE, line = 0, cex = 0.8); par(oma = c(0, 0, 0, 0))
}
page_comm <- function() {
  if (is.null(cs_all)) return(invisible())
  par(mfrow = c(2, 2), mar = c(7, 4.5, 3, 1)); cols <- c(`malignant -> partner` = "#B2182B", `partner -> malignant` = "#2166AC")
  for (nm in unique(cs_all$dataset)) {
    z <- cs_all[cs_all$dataset == nm, ]; z$rho[!(z$n >= 6)] <- NA; parts <- unique(z$partner); M <- sapply(parts, function(p) sapply(names(cols), function(dn) { r <- z$rho[z$partner == p & z$direction == dn]; if (length(r)) r else NA }))
    barplot(M, beside = TRUE, col = cols, names.arg = parts, las = 2, ylim = c(-1, 1), ylab = "Spearman rho (summed L-R strength vs CS2-CS1)",
            main = nm, cex.main = 0.9, cex.names = 0.8); abline(h = 0); legend("topright", names(cols), fill = cols, cex = 0.6, bty = "n")
  }
  if (!is.null(cp_all)) for (nm in unique(cp_all$dataset)) {
    z <- cp_all[cp_all$dataset == nm & is.finite(cp_all$p) & cp_all$n >= 6, ]; z <- head(z[order(z$p), ], 12); if (!nrow(z)) next
    lab <- sprintf("%s [%s %s]", z$pair, ifelse(grepl("^malig", z$direction), "M>", "<M"), z$partner)
    par(mar = c(4.5, 14, 3, 1)); barplot(rev(z$rho), horiz = TRUE, names.arg = rev(lab), las = 1, cex.names = 0.6, col = ifelse(rev(z$rho) > 0, "#B2182B", "#2166AC"),
            xlab = "Spearman rho with tumour CS2 - CS1", main = paste0(nm, ": top 12 L-R pairs (nominal p)"), cex.main = 0.85); abline(v = 0)
  }
}
page_tcga <- function() {
  if (is.null(tc)) return(invisible()); par(mfrow = c(1, 5), mar = c(4, 4, 3, 1))
  for (n in names(tc$scores)) boxplot(tc$scores[[n]] ~ tc$cs, col = c("#2166AC", "#B2182B"), outline = FALSE, ylab = "z-mean score", xlab = "",
                                      main = sprintf("%s\np = %s", n, fmtp(tc$table$wilcox_p[tc$table$score == n])), cex.main = 0.85)
}
pdf(file.path(results_dir, "Figure_7_singlecell.pdf"), width = 13, height = 8.5)
for (nm in names(ds)) tryCatch(page_dataset(ds[[nm]]), error = function(e) msg("figure page '%s' failed: %s", nm, conditionMessage(e)))
for (nm in names(tum)) tryCatch(page_tumour(nm), error = function(e) msg("tumour page '%s' failed: %s", nm, conditionMessage(e)))
tryCatch(page_comm(), error = function(e) msg("communication page failed: %s", conditionMessage(e)))
tryCatch(page_tcga(), error = function(e) msg("TCGA page failed: %s", conditionMessage(e)))
dev.off()
msg("Wrote results/Figure_7_singlecell.pdf")

## ---------------------------------------------------------------------
## 13. Console summary
## ---------------------------------------------------------------------
cat("\n--- CMPIS-gene expression: share of cells with signal, by cell type (mean % detected over the 6 weighted genes) ---\n")
print(aggregate(pct_detected ~ dataset + cell_type, t13[t13$gene %in% names(cm$coef)[cm$coef != 0], ], mean), digits = 3)
cat("\n--- PCD score by cell type (median) ---\n"); print(t13b[, c("dataset", "cell_type", "n_cells", "PCD_median", "CMPIS_mean")], digits = 3)
if (!is.null(t13c)) { cat("\n--- Tumour-level / within-tumour tests (Spearman) ---\n"); print(t13c, digits = 3, row.names = FALSE) }
if (!is.null(cs_all)) { cat("\n--- Communication: summed L-R strength vs tumour CS2-CS1 ---\n"); print(cs_all, digits = 3, row.names = FALSE) }
if (!is.null(resp)) { cat("\n--- Sade-Feldman: pre-treatment responders vs non-responders ---\n"); print(resp, digits = 3, row.names = FALSE) }
if (!is.null(tc)) { cat("\n--- TCGA bulk: programme scores CS2 vs CS1 ---\n"); print(tc$table, digits = 3, row.names = FALSE) }

out$tables <- list(expr_by_type = t13, pcd_by_type = t13b, tests = t13c, tumour = tum, comm_strength = cs_all, comm_pairs = cp_all,
                   response = resp, tcga = if (!is.null(tc)) tc$table else NULL, signature_genes = sets_tab)
out$meta <- lapply(ds, function(d) d$meta); out$umap <- lapply(ds, function(d) d$umap)
saveRDS(out, if (exists("cfg") && !is.null(cfg$sc_output)) cfg$sc_output else file.path(results_dir, "singlecell.rds"))
msg("Stage 13 complete.")
