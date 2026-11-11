library("MOVICS")
# Load required libraries
library(TCGAbiolinks)
library(DESeq2)
library(ggplot2)
library(Rtsne)
library(umap)

library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)

library(dplyr)
library(viper)
library(dorothea)
library(pheatmap)

library(IOBR)
library(reshape2)
library(limma)

library(GEOquery)

library(devtools)
library(GSVA)
library(msigdbr)

library(lubridate)

# Core
library(survival)
library(survminer)
library(survcomp)
library(illuminaHumanv3.db)
library(org.Hs.eg.db)
library(hgu133plus2.db)

devtools::install_github("IOBR/IOBR")


load("skcm.tcga.RData")

# extract multi-omics data
mo.data   <- skcm.tcga[1:5]

# extract raw count data for downstream analyses
count     <- skcm.tcga$count

# extract fpkm data for downstream analyses
fpkm      <- skcm.tcga$fpkm

# extract maf for downstream analysis
maf       <- skcm.tcga$maf
maf$Tumor_Sample_Barcode <- substr((maf$Tumor_Sample_Barcode), 1, 12)

# Tumor Type categorization
rna_query <- GDCquery(project = 'TCGA-SKCM',
                      data.category = 'Transcriptome Profiling',
                      access = 'open',
                      experimental.strategy = 'RNA-Seq')

rna <- GDCprepare(rna_query)

rna_matrix <- assay(rna, 'unstranded')

dds_input <- DESeqDataSetFromMatrix(countData = rna_matrix,
                                    colData = colData(rna),
                                    design = ~sample_type)

tumor_type <- colData(dds_input)$definition
short_names <- substr(rownames(colData(dds_input)), 1, 12)
names(tumor_type) <- short_names


# extract survival information
surv.info <- skcm.tcga$clin.info
surv.info$tumor_type <- tumor_type[rownames(surv.info)]

elite.mRNA <- getElites(dat       = skcm.tcga$mRNA.expr,
                        method    = "cox",
                        surv.info = surv.info, # must provide survival information with 'futime' and 'fustat'
                        na.action = "impute",
                        p.cutoff  = 0.05)
elite.mRNA$elite.dat <- na.omit(elite.mRNA$elite.dat)

elite.miRNA <- getElites(dat       = skcm.tcga$miRNA.expr,
                        method    = "cox",
                        surv.info = surv.info, # must provide survival information with 'futime' and 'fustat'
                        na.action = "impute",
                        p.cutoff  = 0.05)
elite.miRNA$elite.dat <- na.omit(elite.miRNA$elite.dat)

elite.lncRNA <- getElites(dat       = skcm.tcga$lncRNA.expr,
                        method    = "cox",
                        surv.info = surv.info, # must provide survival information with 'futime' and 'fustat'
                        na.action = "impute",
                        p.cutoff  = 0.05)

elite.meth <- getElites(dat       = skcm.tcga$meth.beta,
                        method    = "cox",
                        surv.info = surv.info, # must provide survival information with 'futime' and 'fustat'
                        na.action = "impute",
                        p.cutoff  = 0.05)

elite.mut <- getElites(dat       = skcm.tcga$mut.status,
                       method    = "freq",
                       na.action = "impute",
                       elite.pct = 0.02)

mo.data <- list(mRNA.expr = elite.mRNA$elite.dat,
                miRNA.expr = elite.miRNA$elite.dat,
                lncRNA.expr = elite.lncRNA$elite.dat,
                meth.beta = elite.meth$elite.dat,
                mut.status = elite.mut$elite.dat)

common.samples <- Reduce(intersect, lapply(mo.data, colnames))

mo.data <- lapply(mo.data, function(x) x[, common.samples, drop = FALSE])

optk.tcga <- getClustNum(data        = mo.data,
                         is.binary   = c(F,F,F,F,T), # note: the 5th data is somatic mutation which is a binary matrix
                         try.N.clust = 2:8, # try cluster number from 2 to 8
                         fig.name    = "CLUSTER NUMBER OF TCGA-SKCM")

moic.res.list.leightweight <- getMOIC(data        = mo.data,
                                       methodslist = list("PINSPlus", "LRAcluster", "ConsensusClustering"), 
                                       N.clust     = 2,
                                       type        = c("gaussian", "gaussian", "gaussian", "gaussian", "binomial"))

moic.res.list.moderate <- getMOIC(data        = mo.data,
                                      methodslist = list("MoCluster", "NEMO", "COCA", "IntNMF"),
                                      N.clust     = 2,
                                      type        = c("gaussian", "gaussian", "gaussian", "gaussian", "binomial"))

moic.res.list.heavy <- getMOIC(data        = mo.data,
                                      methodslist = list("SNF", "CIMLR"),
                                      N.clust     = 2,
                                      type        = c("gaussian", "gaussian", "gaussian", "gaussian", "binomial"))

mo.data_mRNA_reduced <- mo.data$mRNA.expr[order(apply(mo.data$mRNA.expr, 1, var), decreasing = TRUE)[1:500], ]
mo.data_meth_reduced <- mo.data$meth.beta[order(apply(mo.data$meth.beta, 1, var), decreasing = TRUE)[1:1000], ]
mo.data_mut_reduced <- mo.data$mut.status[order(apply(mo.data$mut.status, 1, var), decreasing = TRUE)[1:50], ]

moic.res.list.iCluster <- getMOIC(data        = list(mo.data_mRNA_reduced,
                                                     mo.data$miRNA.expr,
                                                     mo.data$lncRNA.expr,
                                                     mo.data_meth_reduced,
                                                     mo.data_mut_reduced),
                                  methodslist = list("iClusterBayes"),
                                  N.clust     = 2,
                                  type        = c("gaussian", "gaussian", "gaussian", "gaussian", "binomial"))

moic.res.list.iCluster <- list(iClusterBayes = moic.res.list.iCluster)
moic.res.list <- c(moic.res.list.leightweight, 
                   moic.res.list.moderate, 
                   moic.res.list.heavy, 
                   moic.res.list.iCluster)

cmoic.skcm <- getConsensusMOIC(moic.res.list = moic.res.list,
                               fig.name      = "CONSENSUS HEATMAP",
                               distance      = "euclidean",
                               linkage       = "average")

getSilhouette(sil      = cmoic.skcm$sil, # a sil object returned by getConsensusMOIC()
              fig.path = getwd(),
              fig.name = "SILHOUETTE",
              height   = 5.5,
              width    = 5)

indata <- mo.data
indata$meth.beta <- log2(indata$meth.beta / (1 - indata$meth.beta))

# data normalization for heatmap
plotdata <- getStdiz(data       = indata,
                     halfwidth  = c(2,2,2,2,NA), # no truncation for mutation
                     centerFlag = c(T,T,T,T,F), # no center for mutation
                     scaleFlag  = c(T,T,T,T,F)) # no scale for mutation

feat   <- moic.res.list$MoCluster$feat.res
feat1  <- feat[which(feat$dataset == "mRNA.expr"),][1:10,"feature"] 
feat2  <- feat[which(feat$dataset == "lncRNA.expr"),][1:10,"feature"]
feat3  <- feat[which(feat$dataset == "miRNA.expr"),][1:10,"feature"]
feat4  <- feat[which(feat$dataset == "meth.beta"),][1:10,"feature"]
feat5  <- feat[which(feat$dataset == "mut.status"),][1:10,"feature"]
annRow <- list(feat1, feat2, feat3, feat4, feat5)


mRNA.col   <- c("#00FF00", "#008000", "#000000", "#800000", "#FF0000")
lncRNA.col <- c("#6699CC", "white"  , "#FF3C38")
miRNA.col <-  c("yellow", "black"  , "blue")
meth.col   <- c("#0074FE", "#96EBF9", "#FEE900", "#F00003")
mut.col    <- c("grey90" , "black")
col.list   <- list(mRNA.col, lncRNA.col, miRNA.col, meth.col, mut.col)

getMoHeatmap(data          = plotdata,
             row.title     = c("mRNA","lncRNA", "miRNA","Methylation","Mutation"),
             is.binary     = c(F,F,F,F,T), # the 4th data is mutation which is binary
             legend.name   = c("mRNA","lncRNA","miRNA","M value","Mutated"),
             clust.res     = moic.res.list$MoCluster$clust.res, # cluster results
             clust.dend    = NULL, # no dendrogram
             show.rownames = c(F,F,F,F,F), # specify for each omics data
             show.colnames = FALSE, # show no sample names
             annRow        = annRow, # mark selected features
             color         = col.list,
             annCol        = NULL, # no annotation for samples
             annColors     = NULL, # no annotation color
             width         = 10, # width of each subheatmap
             height        = 3, # height of each subheatmap
             fig.name      = "COMPREHENSIVE HEATMAP OF Mo Cluster")

getMoHeatmap(data          = plotdata,
             row.title     = c("mRNA","lncRNA","miRNA","Methylation","Mutation"),
             is.binary     = c(F,F,F,F,T), # the 4th data is mutation which is binary
             legend.name   = c("mRNA","lncRNA","miRNA","M value","Mutated"),
             clust.res     = moic.res.list$COCA$clust.res, # cluster results
             clust.dend    = moic.res.list$COCA$clust.dend, # show dendrogram for samples
             color         = col.list,
             width         = 10, # width of each subheatmap
             height        = 3, # height of each subheatmap
             fig.name      = "COMPREHENSIVE HEATMAP OF COCA")

getMoHeatmap(data          = plotdata,
             row.title     = c("mRNA","lncRNA","miRNA","Methylation","Mutation"),
             is.binary     = c(F,F,F,F,T), # the 4th data is mutation which is binary
             legend.name   = c("mRNA","lncRNA","miRNA","M value","Mutated"),
             clust.res     = moic.res.list$NEMO$clust.res, # cluster results
             clust.dend    = NULL, # show dendrogram for samples
             annRow        = annRow,
             color         = col.list,
             width         = 10, # width of each subheatmap
             height        = 3, # height of each subheatmap
             fig.name      = "COMPREHENSIVE HEATMAP OF NEMO")

getMoHeatmap(data          = plotdata,
             row.title     = c("mRNA","lncRNA","miRNA","Methylation","Mutation"),
             is.binary     = c(F,F,F,F,T), # the 4th data is mutation which is binary
             legend.name   = c("mRNA","lncRNA","miRNA","M value","Mutated"),
             clust.res     = moic.res.list$iClusterBayes$clust.res, # cluster results
             clust.dend    = NULL, # show dendrogram for samples
             annRow        = annRow,
             color         = col.list,
             width         = 10, # width of each subheatmap
             height        = 3, # height of each subheatmap
             fig.name      = "COMPREHENSIVE HEATMAP OF iClusterBayes")

# survival comparison
surv.skcm <- compSurv(moic.res         = cmoic.skcm,
                      surv.info        = surv.info,
                      convt.time       = "m", # convert day unit to month
                      surv.median.line = "h", # draw horizontal line at median survival
                      xyrs.est         = c(5,10), # estimate 5 and 10-year survival
                      fig.name         = "KAPLAN-MEIER CURVE OF CONSENSUSMOIC")

clin.skcm <- compClinvar(moic.res      = cmoic.skcm,
                         var2comp      = surv.info, # data.frame needs to summarize (must has row names of samples)
                         strata        = "Subtype", # stratifying variable (e.g., Subtype in this example)
                         factorVars    = c("pstage","fustat","tumor_type"), # features that are considered categorical variables
                         nonnormalVars = "futime", # feature(s) that are considered using nonparametric test
                         exactVars     = "pstage", # feature(s) that are considered using exact test
                         doWord        = TRUE, # generate .docx file in local path
                         tab.name      = "SUMMARIZATION OF CLINICAL FEATURES")
#> --all samples matched.

# mutational frequency comparison
annCol    <- surv.info[,c("pstage", "age", "tumor_type"), drop = TRUE]

# generate corresponding colors for sample annotation
annColors <- list(age    = circlize::colorRamp2(breaks = c(min(annCol$age, na.rm = TRUE),
                                                           median(annCol$age, na.rm = TRUE),
                                                           max(annCol$age, na.rm = TRUE)), 
                                                colors = c("#0000AA", "#555555", "#AAAA00")),
                  pstage = c("Ti"    = "violet",
                             "T0"    = "white",
                             "T1"    = "green",
                             "T2"    = "blue",
                             "T3"    = "red",
                             "T4"    = "yellow", 
                             "TX"    = "black"),
                  tumor_type=c("Primary solid Tumor" = "#4DAF4A","Additional Metastatic" = "#AAAA00", "Metastatic" = "#E41A1C")  # or pick other contrasting colors
)


# mutational frequency comparison
mut.skcm <- compMut(moic.res     = cmoic.skcm,
                    mut.matrix   = skcm.tcga$mut.status, # binary somatic mutation matrix
                    doWord       = TRUE, # generate table in .docx format
                    doPlot       = TRUE, # draw OncoPrint
                    freq.cutoff  = 0.05, # keep those genes that mutated in at least 5% of samples
                    p.adj.cutoff = 1.1, # keep those genes with adjusted p value < 0.05 to draw OncoPrint
                    innerclust   = TRUE, # perform clustering within each subtype
                    annCol       = annCol, # same annotation for heatmap
                    annColors    = annColors, # same annotation color for heatmap
                    width        = 10, 
                    height       = 10,
                    fig.name     = "ONCOPRINT FOR SIGNIFICANT MUTATIONS",
                    tab.name     = "INDEPENDENT TEST BETWEEN SUBTYPE AND MUTATION")
# compare TMB
tmb.skcm <- compTMB(moic.res     = cmoic.skcm,
                    maf          = maf,
                    rmDup        = TRUE, # remove duplicated variants per sample
                    rmFLAGS      = FALSE, # keep FLAGS mutations
                    exome.size   = 38, # estimated exome size
                    test.method  = "nonparametric", # statistical testing method
                    fig.name     = "DISTRIBUTION OF TMB AND TITV")

# drug sensitivity comparison
drug.brca <- compDrugsen(moic.res    = cmoic.skcm,
                         norm.expr   = fpkm[,cmoic.skcm$clust.res$samID], # double guarantee sample order
                         drugs       = c("Cisplatin", "Paclitaxel"), # a vector of names of drug in GDSC
                         tissueType  = "skin", # choose specific tissue type to construct ridge regression model
                         test.method = "nonparametric", # statistical testing method
                         prefix      = "BOXVIOLIN OF ESTIMATED IC50")

surv.info$pstage <- factor(surv.info$pstage, levels = c("TX","Ti","T0","T1","T2","T3","T4"))

agree.skcm <- compAgree(moic.res  = cmoic.skcm,
                        subt2comp = surv.info[,c("tumor_type","pstage")],
                        doPlot    = TRUE,
                        box.width = 0.2,
                        fig.name  = "AGREEMENT OF CONSENSUSMOIC WITH PSTAGE")


runDEA(dea.method = "limma",
       expr       = fpkm, # raw-scale FPKM
       moic.res   = cmoic.skcm,
       prefix     = "TCGA-SKCM")

bad.samples <- rownames(annCol)[rowSums(is.na(annCol)) > 0]
annCol <- annCol[!rownames(annCol) %in% bad.samples, ]
fpkm   <- fpkm[, colnames(fpkm) %in% rownames(annCol)]

# choose edgeR result to identify subtype-specific up-regulated biomarkers
marker.up <- runMarker(moic.res      = cmoic.skcm,
                       dea.method    = "limma",
                       prefix        = "TCGA-SKCM",
                       dirct         = "up",
                       n.marker      = 1000, # switch to 1000
                       doplot        = TRUE,
                       norm.expr     = fpkm,
                       annCol        = annCol,
                       annColors     = annColors,
                       fig.name      = "UPREGULATED BIOMARKER HEATMAP")

skcm.meta <- read.csv("scaled_expression_data.csv", row.names = 1)

# run NTP in Yau cohort by using up-regulated biomarkers
meta.ntp.pred <- runNTP(expr       = skcm.meta,
                       templates  = marker.up$templates, # the template has been already prepared in runMarker()
                       scaleFlag  = TRUE, # scale input data (by default)
                       centerFlag = TRUE, # center input data (by default)
                       doPlot     = TRUE, # to generate heatmap
                       fig.name   = "NTP HEATMAP FOR META")

meta.pam.pred <- runPAM(train.expr  = fpkm,
                       moic.res    = cmoic.skcm,
                       test.expr   = skcm.meta)

tcga.ntp.pred <- runNTP(expr      = fpkm,
                        templates = marker.up$templates,
                        doPlot    = FALSE)

fpkm <- fpkm[apply(fpkm, 1, sd) != 0, ]
sd_vals <- apply(fpkm_log2, 1, sd)
fpkm <- fpkm[!is.na(sd_vals) & !is.nan(sd_vals) & sd_vals != 0, ]
tcga.pam.pred <- runPAM(train.expr  = fpkm_log2,
                        moic.res    = cmoic.skcm,
                        test.expr   = fpkm_log2)

common_samples <- intersect(cmoic.skcm$clust.res$samID, 
                            tcga.pam.pred$clust.res$samID)

pam_clust <- tcga.pam.pred$clust.res$clust[common_samples]
ntp_clust <- tcga.ntp.pred$clust.res$clust[common_samples]
cmoic_clust <- cmoic.skcm$clust.res$clust[common_samples]

pam_clust_res <- tcga.pam.pred$clust.res[tcga.pam.pred$clust.res$samID %in% common_samples, ]
ntp_clust_res <- tcga.ntp.pred$clust.res[tcga.ntp.pred$clust.res$samID %in% common_samples, ]
cmoic_clust_res <- cmoic.skcm$clust.res[cmoic.skcm$clust.res$samID %in% common_samples, ]

pam_clust_res <- pam_clust_res[order(pam_clust_res$samID), ]
ntp_clust_res <- ntp_clust_res[order(ntp_clust_res$samID), ]
cmoic_clust_res <- cmoic_clust_res[order(cmoic_clust_res$samID), ]

pam_clust <- pam_clust_res$clust
ntp_clust <- ntp_clust_res$clust
cmoic_clust <- cmoic_clust_res$clust

runKappa(subt1     = cmoic_clust,
         subt2     = ntp_clust,
         subt1.lab = "CMOIC",
         subt2.lab = "NTP",
         fig.name  = "CONSISTENCY HEATMAP FOR TCGA between CMOIC and NTP")

# check consistency between current and PAM-predicted subtype in discovery TCGA-SKCM
runKappa(subt1     = cmoic_clust,
         subt2     = pam_clust,
         subt1.lab = "CMOIC",
         subt2.lab = "PAM",
         fig.name  = "CONSISTENCY HEATMAP FOR TCGA between CMOIC and PAM")

# check consistency between NTP and PAM-predicted subtype in validation META-SKCM
runKappa(subt1     = meta.ntp.pred$clust.res$clust,
         subt2     = meta.pam.pred$clust.res$clust,
         subt1.lab = "NTP",
         subt2.lab = "PAM",
         fig.name  = "CONSISTENCY HEATMAP FOR NTP AND PAM")

selected_genes <- unique(unlist(marker.up$unqlist))
expr_subset <- t(fpkm[selected_genes, ])
expr_scaled <- scale(expr_subset)
clusters <- as.factor(cmoic_clust_res$clust)
names(clusters) <- cmoic_clust_res$samID

commonSamples <- intersect(rownames(expr_scaled), names(clusters))
expr_scaled <- expr_scaled[common_samples, ]
clusters <- clusters[common_samples]

pca_res <- prcomp(expr_scaled, center = TRUE, scale. = TRUE)
pca_df <- data.frame(PC1 = pca_res$x[,1],
                     PC2 = pca_res$x[,2],
                     Cluster = clusters)

set.seed(42)
tsne_out <- Rtsne(expr_scaled, dims = 2, perplexity = 30)
tsne_df <- data.frame(tSNE1 = tsne_out$Y[,1],
                      tSNE2 = tsne_out$Y[,2],
                      Cluster = clusters)

umap_out <- umap(expr_scaled)
umap_df <- data.frame(UMAP1 = umap_out$layout[,1],
                      UMAP2 = umap_out$layout[,2],
                      Cluster = clusters)

ggplot(pca_df, aes(PC1, PC2, color = Cluster)) +
  geom_point(size = 2, alpha = 0.8) +
  theme_minimal() +
  labs(title = "PCA of Marker Genes")

ggplot(tsne_df, aes(tSNE1, tSNE2, color = Cluster)) +
  geom_point(size = 2, alpha = 0.8) +
  theme_minimal() +
  labs(title = "t-SNE of Marker Genes")

ggplot(umap_df, aes(UMAP1, UMAP2, color = Cluster)) +
  geom_point(size = 2, alpha = 0.8) +
  theme_minimal() +
  labs(title = "UMAP of Marker Genes")
#-------------------------------------------------------------------------------
# Clinical Annotation

clin_data_full <- GDCquery_clinic("TCGA-SKCM", type = "clinical")
rownames(clin_data_full) <- clin_data_full$submitter_id

cs_vector <- setNames(cmoic.skcm$clust.res$clust, cmoic.skcm$clust.res$samID)
clin_ids <- rownames(clin_data_full)
cs_aligned <- cs_vector[clin_ids]
clin_data_full$CS <- cs_aligned

clin_data_full <- clin_data_full[colnames(elite.mRNA$elite.dat), ]
clin_data_full$ajcc_pathologic_stage <- gsub("Stage ([IV]+)[A-C]?", "Stage \\1", clin_data_full$ajcc_pathologic_stage)
clin_data_full$ajcc_pathologic_t <- gsub("T([0-4])[a-d]?", "T\\1", clin_data_full$ajcc_pathologic_t)
clin_data_full$ajcc_pathologic_t[!clin_data_full$ajcc_pathologic_t %in% c("T1", "T2", "T3", "T4")] <- "unknow"
clin_data_full$ajcc_pathologic_n <- gsub("N([0-3])[a-c]?", "N\\1", clin_data_full$ajcc_pathologic_n)
clin_data_full$ajcc_pathologic_m <- gsub("M([0-1])[a-c]?", "M\\1", clin_data_full$ajcc_pathologic_m)
clin_data_full$age_at_index <- ifelse(clin_data_full$age_at_index <= 65, "<=65", ">65")
clin_data_full$CS <- ifelse(clin_data_full$CS == "1", "CS1", "CS2")

expr_raw <- as.matrix(elite.mRNA$elite.dat)
# 2. Scale values globally to range [-10, 10]

expr_scaled <- t(scale(t(expr_raw)))
# Optional: trim extreme values to avoid color distortion
expr_scaled[expr_scaled > 3] <- 3
expr_scaled[expr_scaled < -3] <- -3

clin_data_filtered <- clin_data_full[!is.na(clin_data_full$CS), ]
clin_data_filtered$CS <- factor(clin_data_filtered$CS, levels = c("CS1", "CS2"))
ordered_samples <- rownames(clin_data_filtered)[order(clin_data_filtered$CS)]

expr_scaled_filtered <- expr_scaled[, ordered_samples]
clin_data_filtered <- clin_data_filtered[order(clin_data_filtered$CS),]

ha <- HeatmapAnnotation(
  Stage = clin_data_filtered$ajcc_pathologic_stage,
  M = clin_data_filtered$ajcc_pathologic_m,
  N = clin_data_filtered$ajcc_pathologic_n,
  T = clin_data_filtered$ajcc_pathologic_t,
  Gender = clin_data_filtered$gender,
  Age = clin_data_filtered$age_at_index,
  CS = clin_data_filtered$CS,
  col = list(
    Stage = c("Stage 0"= "#ffffc8","Stage I"="#fee8c8", "Stage II"="#fdbb84", "Stage III"="#fc8d59", "Stage IV"="#e34a33"),
    M = c("M0"="#66c2a5", "M1"="#fc8d62", "unknow"="gray"),
    N = c("NX" = "#ffafc8", "N0"="#8dd3c7", "N1"="#ffffb3", "N2"="#bebada", "N3"="#fb8072", "unknow"="gray"),
    T = c("T1"="#1f78b4", "T2"="#33a02c", "T3"="#ff7f00", "T4"="#6a3d9a", "unknow"="gray"),
    Gender = c("female"="gold", "male"="green4"),
    Age = c("<=65"="lightgreen", ">65"="darkgreen"),
    CS = c("CS1"="cyan", "CS2"="red")
  ),
  annotation_height = unit(5, "mm"),
  show_legend = TRUE
)

existing_samples <- rownames(skcm.tcga$clin.info)

# Match clin_data_filtered to existing samples
matched_indices <- match(existing_samples, rownames(clin_data_filtered))

# Create matched clinical data (will have NAs for missing samples)
skcm.tcga$clin.info <- data.frame(
  fustat = skcm.tcga$clin.info$fustat,  # Keep existing survival data
  futime = skcm.tcga$clin.info$futime,  # Keep existing survival data
  Stage = clin_data_filtered$ajcc_pathologic_stage[matched_indices],
  M = clin_data_filtered$ajcc_pathologic_m[matched_indices],
  N = clin_data_filtered$ajcc_pathologic_n[matched_indices],
  T = clin_data_filtered$ajcc_pathologic_t[matched_indices],
  Gender = clin_data_filtered$gender[matched_indices],
  Age = as.numeric(clin_data$age_at_index[matched_indices]),
  CS = clin_data_filtered$CS[matched_indices],
  stringsAsFactors = FALSE
)
rownames(skcm.tcga$clin.info) <- existing_samples

save(skcm.tcga, file = "skcm.tcga.RData")

Heatmap(
  expr_scaled_filtered,
  name = "Expression",
  show_row_names = FALSE,
  show_column_names = FALSE,
  top_annotation = ha,
  col = colorRamp2(c(-1, 0, 1), c("cyan", "white", "red")),
  cluster_columns = FALSE,
  cluster_rows = FALSE,
  column_title = "SKCM Cancer Subtypes - Clinical Annotation"
)

# Transcription Factors & Chromatin Remodeling Analysis
# Install if needed
BiocManager::install("viper")
BiocManager::install("dorothea")

# 1. Prepare expression matrix (genes x samples)
# 2. Load dorothea regulon (e.g. human confidence level A-C)
data(dorothea_hs, package = "dorothea")
regulon <- dplyr::filter(dorothea_hs, confidence %in% c("A", "B", "C"))

# 3. Infer TF activity (assumes normalized expression matrix)
tf_activity <- run_viper(expr_scaled_filtered, regulon, options = list(method = "scale"))

# 1. Compute row variances (across samples)
regulon_var <- apply(tf_activity, 1, var)

# 2. Order by variance and get top 20
top20_regulons <- names(sort(regulon_var, decreasing = TRUE))[1:20]

# 3. Subset the tf_activity matrix
tf_top20 <- tf_activity[top20_regulons, ]


chromatin_remodelers <- c('CDKN2B-AS1', 'FAM47E', 'HDAC6', 'HDAC5', 'PARTICL', 'PTENP1-AS', 'TRIM28', 'CDK2', 'PRMT3', 'ZMPSTE24', 'CDKN1C', 'SMYD5', 'PCGF3', 'KLF2', 'PRMT5', 'RBM14', 'ZNHIT1', 'CTCF', 'EHMT2', 'SPIN1', 'APOBEC2', 'CBX1', 'WDR5', 'CHEK1', 'KAT7', 'GLMN', 'BAZ2A', 'BAZ1A', 'PADI2', 'SGF29', 'CBX3', 'UHRF2', 'H2AL3', 'TDRD9', 'METTL23', 'MBD3L2', 'ARB2BP', 'H1-8', 'ASZ1', 'ATF2', 'ZFP92', 'CTCFL', 'PIWIL4', 'SPTY2D1', 'C19orf84', 'TDRD5', 'DDB1', 'APOBEC3H', 'RNF168', 'DNMT1', 'DNMT3A', 'DYRK1A', 'EGR1', 'SDR16C5', 'APOBEC3A', 'APOBEC3F', 'TET3', 'CENPV', 'EP300', 'EZH1', 'EZH2', 'H2AC1', 'KDM1B', 'MTF2', 'MORC2', 'BAHD1', 'SIRT2', 'SCMH1', 'SPEN', 'KDM1A', 'PHF8', 'TASOR', 'UBR2', 'TUT4', 'SMCHD1', 'RRP8', 'SIRT1', 'CBX5', 'SUZ12', 'WBP2', 'H2AP', 'CIZ1', 'DCAF13', 'SIN3A', 'L3MBTL1', 'PHF19', 'NOC2L', 'USP21', 'MORC1', 'APOBEC3C', 'MBD3L5', 'STPG4', 'GPX1', 'HIPK2', 'ATAD2', 'N6AMT1', 'BRD7', 'UHRF1', 'GSK3A', 'RBM15B', 'DNMT3L', 'H1-0', 'H1-2', 'H2AC8', 'H2AC7', 'H2AX', 'H2AZ1', 'H3-3A', 'H3-3B', 'HDAC1', 'HDAC2', 'HELLS', 'HMGB1', 'H2AC21', 'HNRNPK', 'HNRNPU', 'APEX1', 'APOBEC1', 'EZHIP', 'IFI16', 'SAMD7', 'ZFP57', 'IGF2', 'APP', 'PADI6', 'ZNF445', 'DPPA3', 'ING2', 'JARID2', 'H1-9P', 'LBR', 'LMNA', 'LMNB1', 'TEX19', 'KPNA7', 'MIR182', 'MIR29A', 'MIR29B1', 'MIR29C', 'MBD1', 'MECP2', 'MEN1', 'KMT2A', 'ASIP', 'MOV10', 'TRIM37', 'MYC', 'MYCN', 'NDN', 'H2AB2', 'H2AB1', 'OGG1', 'RLIM', 'UBR5', 'TRMT112', 'PPHLN1', 'SIRT7', 'SIRT6', 'HDAC7', 'ARID4B', 'PHB1', 'PHF1', 'PHF2', 'PIK3CA', 'MBD3', 'MIS18A', 'POLE3', 'ATAD2B', 'MOV10L1', 'PRMT7', 'DDX4', 'ATRX', 'MPHOSPH8', 'ALKBH4', 'TET2', 'TASOR2', 'NRDE2', 'PIWIL2', 'MSL2', 'RIF1', 'RESF1', 'MACROH2A2', 'ATF7IP', 'H2AJ', 'LRIF1', 'PCID2', 'HDAC8', 'KMT2E', 'TEX15', 'TDRD1', 'METTL3', 'SPHK2', 'SMARCAD1', 'ZNFX1', 'ZNF304', 'AICDA', 'ARID1B', 'BEND3', 'ZDBF2', 'CREBZF', 'RB1', 'ARID4A', 'KDM5A', 'RBBP5', 'TRIM27', 'RLF', 'BCL6', 'RNF2', 'APOBEC3G', 'ZNF335', 'PRDM14', 'BTBD18', 'RBM15', 'BMI1', 'METTL4', 'MBD3L4', 'MBD3L3', 'SMARCA1', 'SMARCA2', 'SMARCA4', 'SMARCB1', 'SMARCD1', 'SPI1', 'BRCA1', 'SUV39H1', 'VPS72', 'TDG', 'TFAP2C', 'TNP1', 'TP53', 'TPR', 'H2AC19', 'MBD3L2B', 'WT1', 'XIST', 'ZNF91', 'USP7', 'TUT7', 'SUV39H2', 'EHMT1', 'TLE6', 'HDAC11', 'TET1', 'EPC1', 'KMT2D', 'HMGA2', 'ZNF93', 'ARID1A', 'AXIN1', 'BAP1', 'H2AC13', 'H2AC15', 'H2AC16', 'H2AC6', 'H2AC4', 'H2AC17', 'H2AC18', 'H2AC20', 'H3C1', 'H3C4', 'H3C3', 'H3C6', 'H3C11', 'H3C8', 'H3C12', 'H3C10', 'H3C2', 'H2AB3', 'HDAC10', 'ARB2A', 'KAT8', 'PCGF5', 'DOT1L', 'L3MBTL3', 'GLYR1', 'DPY30', 'SMARCA5', 'FKBP6', 'LMNB2', 'PPM1D', 'MAEL', 'HAT1', 'H2AC12', 'CGGBP1', 'MBD3L1', 'EED', 'HDAC3', 'ALKBH1', 'KAT2B', 'MBD2', 'H3C7', 'H2AC11', 'RNF8', 'SAMD1', 'ASH2L', 'DIRAS3', 'SPOCD1', 'MTA1', 'TDRD12', 'YTHDC1', 'MTA2', 'PIWIL1', 'H2AC25', 'TRIP12', 'LHX2', 'MYOCD', 'H2AZ2', 'CDYL', 'MACROH2A1', 'APOBEC3B', 'CTR9', 'HDAC9', 'HDAC4', 'SETDB1')
# Subset and check which remodelers are present
present_remodelers <- intersect(rownames(expr_scaled_filtered), chromatin_remodelers)
# Get the top 30 expressed
remodeler_expr <- expr_scaled_filtered[present_remodelers, ]
remodeler_top30 <- remodeler_expr[order(rowMeans(remodeler_expr), decreasing = TRUE)[1:30], ]
cs_vector <- cs_vector[colnames(tf_activity)]
cs_vector <- factor(cs_vector, levels = c(1, 2), labels = c("CS1", "CS2"))

ha <- HeatmapAnnotation(
  CS = cs_vector,
  col = list(
    CS = c("CS1"="cyan", "CS2"="red")
  ),
  annotation_height = unit(5, "mm"),
  show_legend = TRUE
)


ht1 <- Heatmap(
  tf_top20,
  name = "Regulon",
  show_row_names = TRUE,
  show_column_names = FALSE,
  top_annotation = ha,
  col = colorRamp2(c(-3, 0, 3), c("deepskyblue", "black", "yellow")),
  cluster_columns = FALSE,
  cluster_rows = FALSE
)

ht2 <- Heatmap(
  remodeler_top30,
  name = "Chromatin Remodeling",
  show_row_names = TRUE,
  show_column_names = FALSE,
  col = colorRamp2(c(-1, 0, 1), c("turquoise", "white", "red")),
  cluster_columns = FALSE,
  cluster_rows = FALSE,
  column_title = "Top 30 Chromatin Remodelers"
)

ht1 %v% ht2

immune_checkpoints <- c('ABI1', 'ABL2', 'ABRACL', 'ACAA1', 'ACOT9', 'ACTB', 'ACTR2', 'ACTR3', 'ADA', 'ADRM1', 'AHR', 'AIF1', 'AKAP13', 'ALDH2', 'ALOX5AP', 'ANXA1', 'ANXA11', 'ANXA2', 'ANXA5', 'AP1S2', 'AP2S1', 'ARF6', 'ARFGAP3', 'ARHGDIA', 'ARHGDIB', 'ARL4C', 'ARL6IP5', 'ARL8B', 'ARPC1B', 'ARPC2', 
                        'ARPC3', 'ARPC4', 'ARPC5', 'ARRB2', 'ASAH1', 'ATF5', 'ATOX1', 'ATP13A3', 'ATP1B3', 'ATP2B1', 'ATP5F1E', 'ATP5MF', 'ATP6AP2', 'ATP6V0B', 'ATP6V0D1', 'ATP6V0E1', 'ATP6V1B2', 'ATP6V1F', 'ATP6V1H', 'AURKAIP1', 'B2M', 'B4GALT1', 'BASP1', 'BAX', 'BAZ1A', 'BCL2A1', 'BID', 'BLOC1S1', 'BLOC1S2', 'BLVRB', 'BNIP3L', 'BRI3', 'BTG1', 'BZW1', 'TMEM273', 'ATP5MJ', 'C15orf48', 'C1QA', 'C1orf162', 'STMP1', 'CACYBP', 'CALR', 'CAP1', 'CAPG', 'CAPZA1', 'CAPZB', 'CARD16', 'CASP1', 'CASP4', 'CCR7', 'CD14', 'CD1C', 'CD37', 'CD44', 'CD48', 'CD52', 'CD53', 'CD58', 'CD63', 'CD68', 'CD74', 'CD83', 'CD86', 'CDC37', 'CDC42', 'CDC42EP3', 'CDC42SE1', 'CDKN1A', 'CEBPB', 'CFD', 'CFL1', 'CFLAR', 'CFP', 'CHCHD10', 'CHMP1B', 'CHMP4B', 'CIB1', 'CKLF', 'CLEC10A', 'CLEC7A', 'CLIC1', 'CLIC2', 'CMTM6', 'CNPY3', 'COPE', 'CORO1A', 'CORO1C', 'COTL1', 'COX17', 'COX5A', 'COX6B1', 'COX7B', 'COX8A', 'CPVL', 'CREM', 'CRIP1', 'CSF2RA', 'CST3', 'CSTA', 'CSTB', 'CTSB', 'CTSC', 'CTSD', 'CTSH', 'CTSL', 'CTSS', 'CTSZ', 'CXCL16', 'CXCL2', 'CXCL8', 'CXCR4', 'CYBA', 'CYTIP', 'DAZAP2', 'DBI', 'DDX39A', 'DNTTIP2', 'DOK2', 'DPP7', 'DSE', 'DUSP1', 'EFHD2', 'EHD1', 'ELF1', 'ELOF1', 'EMP3', 'ENO1', 'ETS2', 'EVI2A', 'EVI2B', 'EZR', 'FABP5', 'CALHM6', 'CYRIA', 'CYRIB', 'FCER1A', 'FCER1G', 'FCGR2A', 'FCGR2B', 'FCGRT', 'FCN1', 'FERMT3', 'FGL2', 'FGR', 'FKBP4', 'FLNA', 'FNBP1', 'FNIP2', 'FPR3', 'FTH1', 'FTL', 'FXYD5', 'G0S2', 'GABARAP', 'GAPDH', 'GLA', 'GLIPR1', 'GLIPR2', 'GMFG', 'GNA15', 'GNAI2', 'GNG5', 'GPAT3', 'GPR137B', 'GPR183', 'GPSM3', 'GPX1', 'GPX4', 'GRB2', 'GRINA', 'GRN', 'GSTK1', 'GSTO1', 'MACROH2A1', 'H3-3A', 'HAVCR2', 'HCLS1', 'HCST', 'HES4', 'HEXB', 'HLA-A', 'HLA-B', 'HLA-C', 'HLA-DMA', 'HLA-DMB', 'HLA-DPA1', 'HLA-DPB1', 'HLA-DQA1', 'HLA-DQB1', 'HLA-DRA', 'HLA-DRB1', 'HLA-DRB5', 'HM13', 'HMGA1', 'HMOX1', 'HSBP1', 'HSP90AA1', 'HSPA1A', 'HSPA6', 'HSPB1', 'HSPE1', 'IER3', 'IFI30', 'IFNGR1', 'IFNGR2', 'IGSF6', 'IL10RA', 'IL2RG', 'IL4I1', 'INSIG1', 'IQGAP1', 'ISG15', 'ISG20', 'ITGB2', 'JAML', 'JARID2', 'KYNU', 'LACTB', 'LAMTOR2', 'LAMTOR4', 'LAP3', 'LAPTM5', 'LCP1', 'LDHA', 'LGALS1', 'LGALS2', 'LGALS3', 'LGALS9', 'LGMN', 'LIMD2', 'LIMS1', 'CYTOR', 'LITAF', 'LMNA', 'LPXN', 'LRRFIP1', 'LSP1', 'LST1', 'LY86', 'LY96', 'LYZ', 'M6PR', 'MAFB', 'MAP2K1', 'MAP2K3', 'MAP3K8', 'MARCKS', 'MBP', 'MFSD1', 'MGAT1', 'MIR22HG', 'MIR4435-2HG', 'MMP9', 'MNDA', 'MOB1A', 'MPP1', 'MRPL18', 'MS4A6A', 'MS4A7', 'MSN', 'MT-CO1', 'MT-CO2', 'MTPN', 'MYL6', 'MYO9B', 'NAA50', 'NABP1', 'NAMPT', 'NAPA', 'NCOA4', 'NDUFA13', 'NDUFB1', 'NDUFV2', 'NEAT1', 'NFKB1', 'NFKBIA', 'NINJ1', 'NOP10', 'NPC2', 'NUMB', 'OAZ1', 'OLR1', 'OSBPL8', 'OSTF1', 'PABPC1', 'PCBP1', 'PDE4DIP', 'PEA15', 'PET100', 'PFN1', 'PGK1', 'PGLS', 'PHACTR1', 'PIM3', 'PKIB', 'PKM', 'PLAUR', 'PLEK', 'PLEKHB2', 'PLIN2', 'POLR2E', 'POLR2L', 'POMP', 'PPIF', 'PPP1CA', 'PPP1R18', 'PPP4C', 'PPT1', 'PRDX1', 'PRR13', 'PSAP', 'PSMB3', 'PSMB9', 'PSME2', 'PTP4A2', 'PTPN1', 'PTPRC', 'PTPRE', 'PYCARD', 'RAB10', 'RAB5C', 'RAB7A', 'RAB8A', 'RAB8B', 'RALA', 'RAP1A', 'RAP1B', 'RBM47', 'RBX1', 'REL', 'RGCC', 'RGS1', 'RGS10', 'RGS2', 'RHOA', 'RHOG', 'RILPL2', 'RIT1', 'RNASE6', 'RNASET2', 'RNF13', 'RNF130', 'RNF145', 'RNF181', 'RNH1', 'S100A10', 'S100A11', 'S100A4', 'S100A6', 'S100A9', 'SAMSN1', 'SAT1', 'SDCBP', 'SDS', 'SEC11A', 'SEC61B', 'SERF2', 'SERPINA1', 'SERPINB1', 'SERPINB9', 'SFT2D1', 'SH3BGRL3', 'SLC16A3', 'SNAP23', 'SNN', 'SNX3', 'SOD2', 'SPAG9', 'SPCS3', 'SPHK1', 'SPI1', 'SPNS1', 'SQSTM1', 'SRGN', 'SSR4', 'STK4', 'STX11', 'STX4', 'SYNGR2', 'TAGLN2', 'TALDO1', 'TANK', 'TAOK3', 'TAPBP', 'TES', 'THEMIS2', 'TIMP1', 'TMEM123', 'TMEM167A', 'TMEM9B', 'TMSB4X', 'TNFAIP8', 'TNFRSF1B', 'TNFSF13B', 'TNIP1', 'TPD52L2', 'TPM3', 'TPP1', 'TRAPPC1', 'TSPO', 'TUBA1C', 'TUBB4B', 'TXN', 'TXNRD1', 'TYMP', 'TYROBP', 'UBA52', 'UBE2L3', 'UBE2R2', 'UPP1', 'UQCR11', 'ATP5MK', 'VAMP8', 'VASP', 'VMO1', 'VPS29', 'VSIG4', 'WARS1', 'WAS', 'WIPF1', 'WTAP', 'YBX1', 'YWHAH', 'YWHAZ', 'ZEB2', 'ZFAND2A', 'ZFYVE16', 'ZNF706')
present_checkpoints <- intersect(rownames(expr_scaled_filtered), immune_checkpoints)
# Get the top 30 expressed
checkpoints_expr <- expr_scaled_filtered[present_checkpoints, ]
checkpoints_top30 <- checkpoints_expr[order(rowMeans(checkpoints_expr), decreasing = TRUE)[1:30], ]


# Use CIBERSORT to estimate immune infiltration
res_cibersort <- deconvo_tme(eset = expr_scaled_filtered, method = "cibersort", arrays = FALSE, perm = 200)
res_cibersort$ID <- NULL
res_cibersort <- apply(res_cibersort, 2, as.numeric)
rownames(res_cibersort) <- colnames(expr_scaled_filtered)  # or sample IDs
res_cibersort <- t(res_cibersort)
res_cibersort <- t(scale(t(as.matrix(res_cibersort))))
rownames(res_cibersort) <- substr(rownames(res_cibersort), 1, nchar(rownames(res_cibersort)) - 10)

estimate<-deconvo_tme(eset = expr_scaled_filtered, method = "estimate")

ha <- HeatmapAnnotation(
  ImmuneScore = estimate$ImmuneScore_estimate,
  StromalScore = estimate$StromalScore_estimate,
  CS = cs_vector,
  col = list(
    ImmuneScore = colorRamp2(c(-1000, 1000, 3000), c("blue", "white", "black")),
    StromalScore = colorRamp2(c(-1000, 0, 1000), c("blue", "white", "black")),
    CS = c("CS1"="cyan", "CS2"="red")
  )
)

ht3 <- Heatmap(
  checkpoints_top30,
  name = "Immune Checkpoint",
  show_row_names = TRUE,
  show_column_names = FALSE,
  top_annotation = ha,
  col = colorRamp2(c(-2, 0, 2), c("deepskyblue", "black", "yellow")),
  cluster_columns = FALSE,
  cluster_rows = FALSE
)

ht4 <- Heatmap(
  res_cibersort,
  name = "CIBERSORT",
  show_row_names = TRUE,
  show_column_names = FALSE,
  col = colorRamp2(c(-1, 0, 1), c("turquoise", "white", "red")),
  cluster_columns = FALSE,
  cluster_rows = FALSE
)

ht3 %v% ht4


# This works with msigdbr v10.0.0+
kegg_msigdb <- msigdbr(species = "Homo sapiens", collection = "C2", subcollection = "CP:KEGG_LEGACY")

# Format as named list for GSVA
kegg_gene_sets <- split(kegg_msigdb$gene_symbol, kegg_msigdb$gs_name)

common_samples <- intersect(colnames(expr_scaled_filtered), names(cs_vector))
cs1_samples <- names(cs_vector[cs_vector == "CS1"])
cs2_samples <- names(cs_vector[cs_vector == "CS2"])

expr_cs1 <- expr_scaled_filtered[, colnames(expr_scaled_filtered) %in% cs1_samples]
expr_cs2 <- expr_scaled_filtered[, colnames(expr_scaled_filtered) %in% cs2_samples]

gsvapar_1 <- gsvaParam(expr_cs1, kegg_gene_sets, assay="logCPM", minSize=5,
                     maxSize=300)

gsvapar_2 <- gsvaParam(expr_cs2, kegg_gene_sets, assay="logCPM", minSize=5,
                       maxSize=300)

gsva_1 <- gsva(gsvapar_1)
gsva_2 <- gsva(gsvapar_2)

rownames(gsva_1) <- substr(rownames(gsva_1), 6, nchar(rownames(gsva_1)))
rownames(gsva_2) <- substr(rownames(gsva_2), 6, nchar(rownames(gsva_2)))
gsva_1_top30 <- gsva_1[order(rowMeans(gsva_1), decreasing = TRUE)[1:70], ]
gsva_2_top30 <- gsva_2[order(rowMeans(gsva_2), decreasing = TRUE)[1:70], ]

common_pathways <- intersect(rownames(gsva_1_top30), rownames(gsva_2_top30))

gsva_1_top30 <- gsva_1_top30[common_pathways, , drop = FALSE]
gsva_2_top30 <- gsva_2_top30[common_pathways, , drop = FALSE]

ha <- HeatmapAnnotation(
  CS = cs_vector[cs_vector == "CS1"],
  col = list(
    CS = c("CS1"="cyan", "CS2"="red")
  ),
  annotation_height = unit(5, "mm"),
  show_annotation_name = FALSE,
  show_legend = TRUE
)


ht5 <- Heatmap(
  gsva_1_top30,
  name = "Pathways",
  show_row_names = FALSE,
  show_column_names = FALSE,
  top_annotation = ha,
  col = colorRamp2(c(-1, 0, 1), c("turquoise", "white", "red")),
  cluster_columns = FALSE,
  cluster_rows = FALSE,
  row_names_gp = gpar(fontsize = 6)
)

ha <- HeatmapAnnotation(
  CS = cs_vector[cs_vector == "CS2"],
  col = list(
    CS = c("CS1"="cyan", "CS2"="red")
  ),
  annotation_height = unit(5, "mm"),
  show_legend = TRUE
)

ht6 <- Heatmap(
  gsva_2_top30,
  show_row_names = TRUE,
  show_column_names = FALSE,
  show_heatmap_legend = FALSE,
  top_annotation = ha,
  col = colorRamp2(c(-1, 0, 1), c("turquoise", "white", "red")),
  cluster_columns = FALSE,
  cluster_rows = FALSE,
  row_names_gp = gpar(fontsize = 6)
)

draw(ht5 + ht6, merge_legend = TRUE)

melanoma_msigdb <- msigdbr(species = "Homo sapiens", collection = "C2", subcollection = "CP:KEGG_LEGACY") %>% filter(gs_name == "KEGG_MELANOMA")
melanoma_gene_sets <- split(melanoma_msigdb$gene_symbol, melanoma_msigdb$gs_name)

melanoma_gsvapar <- gsvaParam(expr_scaled_filtered, melanoma_gene_sets, assay="logCPM", minSize=5,
                          maxSize=300)

melanoma_gsva <- gsva(melanoma_gsvapar)

#------------------------------------------
melanogenesis_msigdb <- msigdbr(species = "Homo sapiens", collection = "C2", subcollection = "CP:KEGG_LEGACY") %>% filter(gs_name == "KEGG_MELANOGENESIS")
melanogenesis_gene_sets <- split(melanogenesis_msigdb$gene_symbol, melanogenesis_msigdb$gs_name)

melanogenesis_gsvapar <- gsvaParam(expr_scaled_filtered, melanogenesis_gene_sets, assay="logCPM", minSize=5,
                              maxSize=300)

melanogenesis_gsva <- gsva(melanogenesis_gsvapar)
#------------------------------------------
mapk_msigdb <- msigdbr(species = "Homo sapiens", collection = "C2", subcollection = "CP:KEGG_LEGACY") %>% filter(gs_name == "KEGG_MAPK_SIGNALING_PATHWAY")
mapk_gene_sets <- split(mapk_msigdb$gene_symbol, mapk_msigdb$gs_name)

mapk_gsvapar <- gsvaParam(expr_scaled_filtered, mapk_gene_sets, assay="logCPM", minSize=5,
                           maxSize=300)

mapk_gsva <- gsva(mapk_gsvapar)

#------------------------------------------
vegfr_msigdb <- msigdbr(species = "Homo sapiens", collection = "C2", subcollection = "CP:WIKIPATHWAYS") %>% filter(gs_name == "WP_VEGFAVEGFR2_SIGNALING")
vegfr_gene_sets <- split(vegfr_msigdb$gene_symbol, vegfr_msigdb$gs_name)

vegfr_common_samples <- intersect(colnames(expr_scaled_filtered), names(cs_vector))

vegfr_gsvapar <- gsvaParam(expr_scaled_filtered, vegfr_gene_sets, assay="logCPM", minSize=5,
                       maxSize=300)

vegfr_gsva <- gsva(vegfr_gsvapar)

#------------------------------------------
fgfr3_msigdb <- msigdbr(species = "Homo sapiens", collection = "C2", subcollection = "CP:REACTOME") %>% filter(gs_name == "REACTOME_SIGNALING_BY_FGFR3")
fgfr3_gene_sets <- split(fgfr3_msigdb$gene_symbol, fgfr3_msigdb$gs_name)

fgfr3_common_samples <- intersect(colnames(expr_scaled_filtered), names(cs_vector))

fgfr3_gsvapar <- gsvaParam(expr_scaled_filtered, fgfr3_gene_sets, assay="logCPM", minSize=5,
                           maxSize=300)

fgfr3_gsva <- gsva(fgfr3_gsvapar)
#------------------------------------------

#-----------------------------------------
pparg_msigdb <- msigdbr(species = "Homo sapiens", collection = "C3", subcollection = "TFT:GTRD") %>% filter(gs_name == "PPARGC1A_TARGET_GENES")
pparg_gene_sets <- split(pparg_msigdb$gene_symbol, pparg_msigdb$gs_name)

pparg_common_samples <- intersect(colnames(expr_scaled_filtered), names(cs_vector))

pparg_gsvapar <- gsvaParam(expr_scaled_filtered, pparg_gene_sets, assay="logCPM", minSize=5,
                           maxSize=300)

pparg_gsva <- gsva(pparg_gsvapar)
#-----------------------------------------
wnt_msigdb <- msigdbr(species = "Homo sapiens", collection = "H") %>% filter(gs_name == "HALLMARK_WNT_BETA_CATENIN_SIGNALING")
wnt_gene_sets <- split(wnt_msigdb$gene_symbol, wnt_msigdb$gs_name)

wnt_common_samples <- intersect(colnames(expr_scaled_filtered), names(cs_vector))

wnt_gsvapar <- gsvaParam(expr_scaled_filtered, wnt_gene_sets, assay="logCPM", minSize=5,
                           maxSize=300)

wnt_gsva <- gsva(wnt_gsvapar)
#-----------------------------------------
hypoxia_msigdb <- msigdbr(species = "Homo sapiens", collection = "H") %>% filter(gs_name == "HALLMARK_HYPOXIA")
hypoxia_gene_sets <- split(hypoxia_msigdb$gene_symbol, hypoxia_msigdb$gs_name)

hypoxia_common_samples <- intersect(colnames(expr_scaled_filtered), names(cs_vector))

hypoxia_gsvapar <- gsvaParam(expr_scaled_filtered, hypoxia_gene_sets, assay="logCPM", minSize=5,
                         maxSize=300)

hypoxia_gsva <- gsva(hypoxia_gsvapar)

#-----------------------------------------
cell_msigdb <- msigdbr(species = "Homo sapiens", collection = "C2", subcollection = "CP:KEGG_LEGACY") %>% filter(gs_name == "KEGG_CELL_CYCLE")
cell_gene_sets <- split(cell_msigdb$gene_symbol, cell_msigdb$gs_name)

cell_common_samples <- intersect(colnames(expr_scaled_filtered), names(cs_vector))

cell_gsvapar <- gsvaParam(expr_scaled_filtered, cell_gene_sets, assay="logCPM", minSize=5,
                             maxSize=300)

cell_gsva <- gsva(cell_gsvapar)
#-----------------------------------------
dna_msigdb <- msigdbr(species = "Homo sapiens", collection = "C2", subcollection = "CP:BIOCARTA") %>% filter(gs_name == "BIOCARTA_CDC25_PATHWAY")
dna_gene_sets <- split(dna_msigdb$gene_symbol, dna_msigdb$gs_name)

dna_gsvapar <- gsvaParam(expr_scaled_filtered, dna_gene_sets, assay="logCPM", minSize=5,
                          maxSize=300)

dna_gsva <- gsva(dna_gsvapar)

rownames(melanoma_gsva) <- "Melanoma"
rownames(melanogenesis_gsva) <- "Melanogenesis"
rownames(mapk_gsva) <- "MAPK"
rownames(vegfr_gsva) <- "VEGFR"
rownames(fgfr3_gsva) <- "FGFR3"
rownames(pparg_gsva) <- "PPARG"
rownames(wnt_gsva) <- "WNT-beta catenin"
rownames(hypoxia_gsva) <- "Hypoxia"
rownames(cell_gsva) <- "Cell_cycle"
rownames(dna_gsva) <- "DNA_Replication"

combined_gsva <- rbind(melanoma_gsva, melanogenesis_gsva, mapk_gsva, vegfr_gsva, fgfr3_gsva, pparg_gsva, wnt_gsva, hypoxia_gsva, cell_gsva)

ha <- HeatmapAnnotation(
  CS = cs_vector,
  col = list(
    CS = c("CS1"="cyan", "CS2"="red")
  ),
  annotation_height = unit(5, "mm"),
  show_legend = TRUE
)

row_group <- c(
  rep("Oncogenic Signaling", 3),
  rep("EGFR Network", 4),
  rep("Radiotherapy Predicted Pathways", 2)
)

left_ann <- rowAnnotation(
  Group = row_group,
  name = "Path",
  col = list(Group = 
               c(
                 "Oncogenic Signaling" = "darkorange",
                 "EGFR Network" = "violet",
                 "HNSCC Signatures" = "magenta",
                 "Radiotherapy Predicted Pathways" = "black"
               )),
  annotation_label = "Pathway Enrichment",
  show_annotation_name = FALSE
)

ht7 <- Heatmap(
  combined_gsva,
  name = "Pathway Type",
  show_row_names = TRUE,
  show_column_names = FALSE,
  top_annotation = ha,
  col = colorRamp2(c(-0.5, 0, 0.5), c("deepskyblue", "black", "yellow")),
  cluster_columns = FALSE,
  cluster_rows = FALSE,
  row_split = row_group,
  row_gap = unit(2, "mm"),
  row_title = NULL
  )

left_ann + ht7

#--------------------------------------------------------------------

expr <- skcm.tcga$mRNA.expr[, names(cs_vector)]
expr <- scale(expr)
group <- factor(cs_vector)
design <- model.matrix(~0 + group)
colnames(design) <- levels(group)

fit <- lmFit(expr, design)
fit <- eBayes(fit)

upregulated_list <- list()

for (cs in levels(group)) {
  others <- setdiff(levels(group), cs)
  contrast_formula <- paste(cs, "-", paste(others, collapse = "-/length(others)"), sep = "")
  
  contrast_matrix <- makeContrasts(contrasts = contrast_formula, levels = design)
  fit2 <- contrasts.fit(fit, contrast_matrix)
  fit2 <- eBayes(fit2)
  
  top_genes <- topTable(fit2, number = Inf, adjust.method = "fdr")
  top_upregulated <- top_genes %>%
    filter(logFC > 0, adj.P.Val < 0.05) %>%
    arrange(-logFC) %>%
    head(2000)
  
  upregulated_list[[cs]] <- rownames(top_upregulated)
}

pcdrgs <- readxl::read_excel(path = "PCD GENES.xlsx")$gene
candidate_pcdrgs_list <- lapply(upregulated_list, function(cs_genes) {
  intersect(cs_genes, pcdrgs)
})

all_candidates <- unique(unlist(candidate_pcdrgs_list))

# Genes shared across subtypes (optional)
shared_genes <- Reduce(intersect, candidate_pcdrgs_list)

cat("Total unique candidate PCDRGs:", length(all_candidates), "\n")
cat("Shared across all CSs:", length(shared_genes), "\n")
tcga.expr <- expr[all_candidates,]
tcga.expr <- scale(tcga.expr)

elite.tcga_pcdrgs <- getElites(dat       = tcga.expr,
                        method    = "cox",
                        surv.info = surv.info, # must provide survival information with 'futime' and 'fustat'
                        na.action = "impute",
                        p.cutoff  = 0.05)

elite.tcga_pcdrgs$elite.dat <- na.omit(elite.tcga_pcdrgs$elite.dat)

meta.cleaned <- skcm.meta[!is.na(skcm.meta[,1]) & skcm.meta[,1] != "", ]


common.META <- Reduce(intersect, list(
  all_candidates,
  rownames(meta.cleaned)
))

meta.pcdrg.expr <- meta.cleaned[common.META,]

gse15605 <- getGEO("GSE15605", GSEMatrix = TRUE)[[1]]
pdata_15605 <- pData(gse15605)
gse22153 <- getGEO("GSE22153", GSEMatrix = TRUE)[[1]]
pdata_22153 <- pData(gse22153)

os_22153 <- data.frame(
  sample_id = rownames(pdata_22153),
  os_time = as.numeric(gsub(".*: ", "", pdata_22153$`characteristics_ch1.1`)),
  os_event = as.numeric(gsub(".*: ", "", pdata_22153$`characteristics_ch1.2`))
)
gse46517 <- getGEO("GSE46517", GSEMatrix = TRUE)[[1]]
pdata_46517 <- pData(gse46517)
resection_date <- dmy(paste0("01-", pdata_46517$`date of resection:ch1`), format = "%d-%b-%y")

dates_raw <- pdata_46517$`last clinical status:ch1`

# Replace German → English
months_map <- c("Dez" = "Dec", "Mai" = "May", "Mar" = "Mar", "Okt" = "Oct", 
                "Jul" = "Jul", "Jun" = "Jun", "Feb" = "Feb", "Jan" = "Jan", 
                "Sep" = "Sep", "Aug" = "Aug", "Nov" = "Nov", "Apr" = "Apr")

# Replace each manually
for (ger in names(months_map)) {
  dates_raw <- gsub(ger, months_map[ger], dates_raw)
}

# Handle entries with just the year
dates_raw[grepl("^([0-9]{4}).*$", dates_raw)] <- paste0("Jan ", dates_raw[grepl("^([0-9]{4}).*$", dates_raw)])

# Now prepend "01-" to make day explicit
dates_standardized <- paste0("01-", dates_raw)

last_update_date <- ifelse(grepl(",", dates_standardized),
       sub(",.*$", "", dates_standardized),
       NA)

last_update_date <- dmy(last_update_date)

os_time_days <- as.numeric(difftime(last_update_date, resection_date, units = "days"))

os_46517 <- data.frame(
  sample_id = rownames(pdata_46517),
  os_time = os_time_days[2:122],
  os_event = ifelse(
    grepl("alive", pdata_46517$`last clinical status`),
    0,
    1
  )
)

gse54467 <- getGEO("GSE54467", GSEMatrix = TRUE)[[1]]
pdata_54467 <- pData(gse54467)

os_54467 <- data.frame(
  sample_id = rownames(pdata_54467),
  os_time = as.integer(30.4375*as.numeric(gsub(".*: ", "", pdata_54467$`characteristics_ch1.4`))),
  os_event = ifelse(
    grepl("Alive", pdata_54467$`characteristics_ch1.5`),
    0,
    1
  )
)

surv.meta <- bind_rows(
  os_22153,
  os_46517,
  os_54467
)


colnames(surv.meta)[colnames(surv.meta) == "os_time"] <- "futime"
colnames(surv.meta)[colnames(surv.meta) == "os_event"] <- "fustat"

surv.meta <- surv.meta %>%
  filter(!is.na(futime), !is.na(fustat))
rownames(surv.meta) <- surv.meta$sample_id
surv.meta$sample_id <- NULL  # optional cleanup

common_samples <- intersect(colnames(meta.pcdrg.expr), rownames(surv.meta))

meta.pcdrg.expr <- meta.pcdrg.expr[, common_samples]
surv.meta <- surv.meta[common_samples, ]

elite.meta <- getElites(dat = meta.pcdrg.expr,
                        method    = "cox",
                        surv.info = surv.meta,
                        na.action = "impute",
                        p.cutoff  = 0.05)

elite.meta$elite.dat <- na.omit(elite.meta$elite.dat)

tcga.sig <- subset(elite.tcga_pcdrgs$unicox.res, pvalue < 0.05)
meta.sig <- subset(elite.meta$unicox.res, pvalue < 0.05)

common.genes <- intersect(tcga.sig$gene, meta.sig$gene)

# Δημιουργία data.frame για εύκολη σύγκριση
tcga.df <- tcga.sig[tcga.sig$gene %in% common.genes, ]
meta.df <- meta.sig[meta.sig$gene %in% common.genes, ]

# Ταξινόμηση για να ταιριάξουν οι γραμμές
tcga.df <- tcga.df[order(tcga.df$gene), ]
meta.df <- meta.df[order(meta.df$gene), ]

same.dir <- (tcga.df$HR > 1 & meta.df$HR > 1) |
  (tcga.df$HR < 1 & meta.df$HR < 1)

final.pcdrgs <- tcga.df$gene[same.dir]

final.expr.tcga <- tcga.expr[final.pcdrgs, , drop = FALSE]
final.expr.meta <- skcm.meta[final.pcdrgs, , drop = FALSE]

#-------------------output format for ML models---------------------------
convert_expression_matrix <- function(expr_matrix, output_file) {
  
  # Ensure genes are row names
  if(is.null(rownames(expr_matrix))) {
    stop("Expression matrix must have gene names as row names")
  }
  
  # Write to file with gene names as first column
  write.table(expr_matrix, file = output_file, 
              sep = "\t", quote = FALSE, 
              col.names = NA, row.names = TRUE)
  
  cat("Expression matrix saved to:", output_file, "\n")
  cat("Dimensions:", nrow(expr_matrix), "genes x", ncol(expr_matrix), "samples\n")
}

convert_meta_cohort <- function(combined_matrix, time_col, status_col, 
                                gene_cols = NULL, output_file) {
  
  # If gene_cols not specified, assume all non-survival columns are genes
  if(is.null(gene_cols)) {
    gene_cols <- setdiff(colnames(combined_matrix), c(time_col, status_col))
  }

  # Reorder columns: time, status, then genes
  ordered_cols <- c(time_col, status_col, gene_cols)
  formatted_data <- combined_matrix[, ordered_cols]
  
  # Standardize column names
  colnames(formatted_data)[1:2] <- c("time", "status")
  
  # Write to file
  write.table(formatted_data, file = output_file, 
              sep = "\t", quote = FALSE, 
              col.names = TRUE, row.names = FALSE)
  
  cat("META cohort saved to:", output_file, "\n")
  cat("Samples:", nrow(formatted_data), ", Genes:", length(gene_cols), "\n")
}

convert_clinical_data <- function(clinical_df, sample_id_col, time_col = NULL, 
                                  status_col = NULL, output_file) {
  
  # Standardize column names
  colnames(clinical_df)[colnames(clinical_df) == sample_id_col] <- "submitter_id"
  
  # Handle different time column formats
  if(!is.null(time_col)) {
    colnames(clinical_df)[colnames(clinical_df) == time_col] <- "days_to_death"
  }
  
  # Handle status column
  if(!is.null(status_col)) {
    colnames(clinical_df)[colnames(clinical_df) == status_col] <- "vital_status"
    # Convert to standard format if needed
    if(is.numeric(clinical_df$vital_status)) {
      clinical_df$vital_status <- ifelse(clinical_df$vital_status == 1, "Dead", "Alive")
    }
  }
  
  # Write to file
  write.table(clinical_df, file = output_file, 
              sep = "\t", quote = FALSE, 
              col.names = TRUE, row.names = FALSE)
  
  cat("Clinical data saved to:", output_file, "\n")
  cat("Samples:", nrow(clinical_df), "\n")
}

convert_gene_list <- function(gene_vector, output_file) {
  
  gene_df <- data.frame(gene_name = gene_vector)
  
  write.table(gene_df, file = output_file, 
              sep = "\t", quote = FALSE, 
              col.names = TRUE, row.names = FALSE)
  
  cat("Gene list saved to:", output_file, "\n")
  cat("Total genes:", length(gene_vector), "\n")
}

tcga_expression <- final.expr.tcga[,intersect(colnames(final.expr.tcga), rownames(clin_data_full))]
convert_expression_matrix(tcga_expression, "tcga_skcm_expression.txt")

meta_final <- merge(surv.meta, t(final.expr.meta), by=0, all=TRUE)
rownames(meta_final)=meta_final$Row.names
meta_final$Row.names <- NULL

convert_meta_cohort(
  combined_matrix = meta_final, 
  time_col = "futime",
  status_col = "fustat",
  output_file = "meta_cohort.txt")

clinical_df <- clin_data_full[, c("submitter_id", "days_to_death","days_to_last_follow_up", "vital_status")]

convert_clinical_data(clinical_df = clinical_df, sample_id_col = "submitter_id", time_col = "days_to_death", status_col = "vital_status", output_file = "tcga_skcm_clinical.txt")

convert_gene_list(final.pcdrgs, "pcdrg_genes.txt")


gse65904 <- getGEO("GSE65904", GSEMatrix = TRUE)[[1]]
pdata_65904 <- pData(gse65904)
expr_65904 <- exprs(gse65904)
probe_ids <- rownames(expr_65904)
gene_symbols <- mapIds(illuminaHumanv3.db,
                       keys = probe_ids,
                       column = "SYMBOL",
                       keytype = "PROBEID",
                       multiVals = "first")
rownames(expr_65904) <- gene_symbols
expr_65904 <- expr_65904[final.pcdrgs,]
expr_65904 <- t(scale(t(expr_65904)))

os_65904 <- data.frame(
  sample_id = rownames(pdata_65904),
  os_time = as.integer(sub(".*: ", "", pdata_65904$`characteristics_ch1.6`)),
  os_event = ifelse(
    grepl(0, sub(".*: ", "", pdata_65904$`characteristics_ch1.7`)),
    0,
    1
  )
)
rownames(os_65904) <- os_65904$sample_id
os_65904$sample_id <- NULL

final_65904 <- merge(os_65904, t(expr_65904), by=0, all=TRUE)
rownames(final_65904)=final_65904$Row.names
final_65904$Row.names <- NULL

convert_meta_cohort(
  combined_matrix = final_65904, 
  time_col = "os_time",
  status_col = "os_event",
  output_file = "gse65904_cohort.txt")

gse98394 <- getGEO("GSE98394", GSEMatrix = TRUE)[[1]]
pdata_98394 <- pData(gse98394)
expr_98394 <- read.table("GSE98394_expression.txt", header = TRUE, row.names = 1, sep = "\t")
colnames(expr_98394) <- rownames(pdata_98394)

# Get the Ensembl gene IDs from your data
ensembl_ids <- rownames(expr_98394)

# Use mapIds to convert the IDs
gene_symbols <- mapIds(illuminaHumanv3.db,
                       keys = ensembl_ids,
                       column = "SYMBOL",
                       keytype = "ENSEMBL",
                       multiVals = "first")
valid_genes <- !is.na(gene_symbols)
expr_98394 <- expr_98394[valid_genes, ]
gene_symbols_filtered <- gene_symbols[valid_genes]
expr_98394 <- as.data.frame(expr_98394)
expr_98394 <- limma::avereps(expr_98394, ID = gene_symbols_filtered)
expr_98394 <- expr_98394[final.pcdrgs,]
expr_98394 <- t(scale(t(expr_98394)))

os_98394 <- data.frame(
  sample_id = rownames(pdata_98394),
  os_time = as.integer(sub(".*: ", "", pdata_98394$`characteristics_ch1.8`)),
  os_event = ifelse(
    grepl(0, sub(".*: ", "", pdata_98394$`characteristics_ch1.8`)),
    0,
    1
  )
)

rownames(os_98394) <- os_98394$sample_id
os_98394$sample_id <- NULL

final_98394 <- merge(os_98394, t(expr_98394), by=0, all=TRUE)
rownames(final_98394)=final_98394$Row.names
final_98394$Row.names <- NULL

convert_meta_cohort(
  combined_matrix = final_98394, 
  time_col = "os_time",
  status_col = "os_event",
  output_file = "gse98394_cohort.txt")


gse19234 <- getGEO("GSE19234", GSEMatrix = TRUE)[[1]]
pdata_19234 <- pData(gse19234)
expr_19234 <- exprs(gse19234)
ensembl_ids <- rownames(expr_19234)

gene_symbols <- select(hgu133plus2.db,
                       keys = ensembl_ids,
                       columns = c("SYMBOL", "GENENAME"),
                       keytype = "PROBEID")

probe2symbol <- gene_symbols$SYMBOL
names(probe2symbol) <- gene_symbols$PROBEID

expr_19234 <- expr_19234[rownames(expr_19234) %in% names(probe2symbol), ]
rownames(expr_19234) <- probe2symbol[rownames(expr_19234)]
expr_19234 <- expr_19234[final.pcdrgs,]
expr_19234 <- t(scale(t(expr_19234)))

os_19234 <- data.frame(
  sample_id = rownames(pdata_19234),
  os_time = as.integer(sub(".*: ", "", pdata_19234$`characteristics_ch1.1`)),
  os_event = ifelse(
    grepl(0, sub(".*: ", "", pdata_19234$`characteristics_ch1.5`)),
    0,
    1
  )
)

rownames(os_19234) <- os_19234$sample_id
os_19234$sample_id <- NULL

final_19234 <- merge(os_19234, t(expr_19234), by=0, all=TRUE)
rownames(final_19234)=final_19234$Row.names
final_19234$Row.names <- NULL

convert_meta_cohort(
  combined_matrix = final_19234, 
  time_col = "os_time",
  status_col = "os_event",
  output_file = "gse19234_cohort.txt")

gse53118 <- getGEO("GSE53118", GSEMatrix = TRUE)[[1]]
pdata_53118 <- pData(gse53118)
expr_53118 <- exprs(gse53118)

ensembl_ids <- rownames(expr_53118)

# Use mapIds to convert the IDs
gene_symbols <- mapIds(illuminaHumanv3.db,
                       keys = ensembl_ids,
                       column = "SYMBOL",
                       keytype = "PROBEID",
                       multiVals = "first")
rownames(expr_53118) <- gene_symbols
expr_53118 <- expr_53118[final.pcdrgs,]
expr_53118 <- t(scale(t(expr_53118)))

os_53118 <- data.frame(
  sample_id = rownames(pdata_53118),
  os_time = as.integer(as.integer(sub(".*: ", "", pdata_53118$`characteristics_ch1.4`))*30.4375),
  os_event = ifelse(
    grepl("Alive", pdata_53118$`characteristics_ch1.5`, ignore.case = TRUE),
    0,
    1
  )
)

rownames(os_53118) <- os_53118$sample_id
os_53118$sample_id <- NULL

final_53118 <- merge(os_53118, t(expr_53118), by=0, all=TRUE)
rownames(final_53118)=final_53118$Row.names
final_53118$Row.names <- NULL

convert_meta_cohort(
  combined_matrix = final_53118, 
  time_col = "os_time",
  status_col = "os_event",
  output_file = "gse53118_cohort.txt")

gse_83583 <- getGEO("GSE83583", GSEMatrix = TRUE)
gse_data <- gse_83583[[1]] 
clinical_data <- pData(gse_data)

# Print the column names to look for survival information
colnames(clinical_data)
