### Load our libraries
library(TCGAbiolinks)
library(SummarizedExperiment)
library(tidyverse)
library(DESeq2)
library(umap)
library(survival)
library(survminer)

### Setup working directory
# run from the project root (PCD.Rproj)

### Load the transcription datasets.
rna_query <- GDCquery(project = 'TCGA-SKCM',
                      data.category = 'Transcriptome Profiling',
                      access = 'open',
                      experimental.strategy = 'RNA-Seq')
                      #sample.type = c('Primary Tumor', 'Solid Tissue Normal'))

output_rna_query <- getResults(rna_query)

rna <- GDCprepare(rna_query)

rna_matrix <- assay(rna, 'unstranded')

dds_input <- DESeqDataSetFromMatrix(countData = rna_matrix,
                                    colData = colData(rna),
                                    design = ~sample_type)

keep <- rowSums(counts(dds_input)) >= 10

dds_input <- dds_input[keep,]

dds <- DESeq(dds_input)

### Variance stable transformation
vsd <- vst(dds, blind = FALSE)

### Dimensionality reduction by PCA
plotPCA(vsd, intgroup = 'sample_type')

### Plot UMAP
vsd_count <- assay(vsd)

set.seed(24)

umap <- umap(t(vsd_count), 
             n_neighbors = 15)

umap_res <- as.data.frame(umap$layout)

ggplot(data = umap_res,
       aes(x = V1, y = V2)) +
  geom_point() +
  xlab('umap1') +
  ylab('umap2')

### Load umap clusters by different criteria
umap_merged <- cbind(umap_res, colData(vsd))

ggplot(data = umap_merged, aes(x = V1, y= V2, 
                               color = sample_type)) +
  geom_point() +
  xlab('umap1') +
  ylab('umap2')

umap_merged <- cbind(umap_res, colData(vsd))

ggplot(data = umap_merged, aes(x = V1, y= V2, 
                               color = gender)) +
  geom_point() +
  xlab('umap1') +
  ylab('umap2')

ggplot(data = umap_merged, aes(x = V1, y= V2, 
                               color = primary_diagnosis)) +
  geom_point() +
  xlab('umap1') +
  ylab('umap2')

### Retrieve Clinical Data
clinical_query <- GDCquery_clinic('TCGA-SKCM')

### Neat our clinical datasets
all(c('vital_status', 'days_to_last_follow_up', 'days_to_death') %in% colnames(clinical_query))

which(colnames(clinical_query) %in% c('vital_status', 'days_to_last_follow_up', 'days_to_death'))

View(clinical_query[,c(24, 40, 46)])

table(clinical_query$vital_status)

clinical_query$deceased <- ifelse(clinical_query$vital_status == 'Alive', FALSE, TRUE)

table(clinical_query$deceased)

clinical_query$overall_survival <- ifelse(clinical_query$deceased==FALSE, clinical_query$days_to_last_follow_up, clinical_query$days_to_death)

### Integrate with our Transcription (TXN) data
gene_metadata <- as.data.frame(rowData(rna))

norm_count <- counts(dds, normalized = T) %>%
  as.data.frame() %>%
  rownames_to_column(var = 'gene_id') %>%
  gather(key = 'case_id', value = 'counts', -'gene_id') %>%
  left_join(., gene_metadata, by = 'gene_id')

### Pick our gene of interest
goi <- norm_count %>%
  filter(gene_name=='CDKN1A')

### Filter out Transcription Data from solid tissue normals
tumor_indeces <- grep('-0[1-9][A-Z]-.*', goi$case_id)

tumor_goi <- goi[tumor_indeces,]

### Seperate our tumor group into 2 based on their median Transcription level

median_goi <- median(tumor_goi$counts)
median_goi

tumor_goi$strata <- ifelse(tumor_goi$counts >= median_goi, 'High', 'Low')

table(tumor_goi$strata)

### Merge clinical data
tumor_goi$submitter_id <- gsub('-0[1-9][A-Z]-.*', '', tumor_goi$case_id)

tumor_goi <- merge(tumor_goi, clinical_query, by = 'submitter_id')

### Fitting our survival curve
fit <- survfit(Surv(overall_survival, deceased) ~strata, data = tumor_goi)

ggsurvplot(fit = fit,
           data = tumor_goi,
           pval = T,
           risk.table = T)