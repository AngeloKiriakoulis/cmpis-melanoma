### Install necessary packages

install.packages('BiocManager')
library(BiocManager)

BiocManager::install('remotes')
BiocManager::install('BioinformaticsFMRP/TCGAbiolinksGUI.data')
BiocManager::install('BioinformaticsFMRP/TCGAbiolinks')

library(TCGAbiolinks)

BiocManager::install('tidyverse')
BiocManager::install('SummarizedExperiment')

BiocManager::install('maftools')
BiocManager::install('DESeq2')
BiocManager::install('EnhancedVolcano')
BiocManager::install('umap')
BiocManager::install('survival')
BiocManager::install('survminer')

### Load necessary packages
library(TCGAbiolinks)
library(SummarizedExperiment)
library(tidyverse)
library(maftools)
library(DESeq2)
library(EnhancedVolcano)

### Set up working directory
# run from the project root (PCD.Rproj)

### Get table of contents from TCGAbiolinks

gdc_projects <- getGDCprojects()
getProjectSummary('TCGA-SKCM')

### Build a query, and filter only for Masked Somatic Mutations.
snv_query <- GDCquery(project = 'TCGA-SKCM',
         data.category = 'Simple Nucleotide Variation', access = 'open', data.type = 'Masked Somatic Mutation')
output_snv_query <- getResults(snv_query)

### Retrieve datasets
GDCdownload(snv_query)

### Load our downloaded data into R
maf <-GDCprepare(snv_query)

### Visualize mutational datasets by maftools package
maftools_maf <- read.maf(maf = maf)
getSampleSummary(maftools_maf)
plotmafSummary(maf = maftools_maf,
               rmOutlier = TRUE,
               dashboard = TRUE)

oncoplot(maf = maftools_maf,
         top = 20)

titv(maf = maftools_maf,
     plot = TRUE,
     useSyn = TRUE)

#### Second example
rna_query <- GDCquery(project = 'TCGA-SKCM',
                     data.category = 'Transcriptome Profiling',
                     access = 'open',
                     experimental.strategy = 'RNA-Seq')

output_rna_query <- getResults(rna_query)

GDCdownload(rna_query)

rna <- GDCprepare(rna_query)

rna_matrix <- assay(rna, 'unstranded')

### Differential expression analysis by DESeq2
dds_input <- DESeqDataSetFromMatrix(countData = rna_matrix,
                                    colData = colData(rna),
                                    design = ~sample_type)

### Filter low counts
keep <- rowSums(counts(dds_input)) >= 10

dds_input <- dds_input[keep,]

dds <- DESeq(dds_input)

res <- results(dds)

summary(res)

res

write.csv(res, file = 'de.csv')

### Transfer gene ID
gene_metadata <- as.data.frame(rowData(rna))

res_df <- as.data.frame(res)

res_df$gene_id <- rownames(res_df)

res_df <- left_join(x = res_df, y = gene_metadata, by = 'gene_id')

### Data visualization
EnhancedVolcano(toptable = res_df,
                lab = res_df$gene_name,
                x = 'log2FoldChange',
                y = 'padj',
                pCutoff = 0.05,
                FCcutoff = 1)

