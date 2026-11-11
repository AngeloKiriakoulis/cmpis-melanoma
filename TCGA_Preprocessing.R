# Load required libraries
library(TCGAbiolinks)    # query/download GDC data:contentReference[oaicite:5]{index=5}
library(SummarizedExperiment)
library(biomaRt)        # gene annotation
library(maftools)       # mutation file processing
library(dplyr)          # data manipulation
library(EDASeq)         # for normalization (used by TCGAbiolinks)
library(DESeq2)         # for gene expression
library("org.Hs.eg.db") # gene annotation
library("AnnotationDbi") # gene annotation

# 1. Query and download RNA-seq (gene expression) data
rna_query <- GDCquery(
  project = 'TCGA-SKCM',
  data.category = 'Transcriptome Profiling',
  data.type = 'Gene Expression Quantification',
  workflow.type = 'STAR - Counts',
  access = 'open',
  experimental.strategy = 'RNA-Seq'
)

GDCdownload(rna_query)

rna_data <- GDCprepare(rna_query)

# 2. Extract raw count matrix and FPKM data
counts <- assay(rna_data, 'unstranded')

fpkm <- assay(rna_data, 'fpkm_unstrand')

# 3. Filter out genes with low expression (remove genes with ~0 counts in all samples)
counts_keep_genes <- rowSums(counts > 1) > 0   # at least one sample with counts>1
counts <- counts[counts_keep_genes, ]

fpkm_keep_genes <- rowSums(fpkm > 1) > 0   # at least one sample with FPKM>1
fpkm <- fpkm[fpkm_keep_genes, ]

# 4. Annotate genes with biotype and symbol
counts <- as.data.frame(counts)

# Remove version numbers from Ensembl IDs
rownames(counts) <- sub("\\..*$", "", rownames(counts))

# Map Ensembl IDs to gene symbols
gene_symbols <- mapIds(
  org.Hs.eg.db,
  keys = rownames(counts),
  keytype = "ENSEMBL",
  column = "SYMBOL"
)

# Replace rownames with gene symbols
# Keep only entries with non-NA gene symbols
valid_rows <- !is.na(gene_symbols)
counts <- counts[valid_rows, ]
rownames(counts) <- make.unique(gene_symbols[valid_rows])

fpkm <- as.data.frame(fpkm)

# Remove version numbers from Ensembl IDs
rownames(fpkm) <- sub("\\..*$", "", rownames(fpkm))

# Map Ensembl IDs to gene symbols
gene_symbols <- mapIds(
  org.Hs.eg.db,
  keys = rownames(fpkm),
  keytype = "ENSEMBL",
  column = "SYMBOL"
)

# Replace rownames with gene symbols
# Keep only entries with non-NA gene symbols
valid_rows <- !is.na(gene_symbols)
fpkm <- fpkm[valid_rows, ]
rownames(fpkm) <- make.unique(gene_symbols[valid_rows])

colnames(counts) <- substr(colnames(counts), 1, 12)
colnames(fpkm) <- substr(colnames(fpkm), 1, 12)

mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# Assuming you have the gene names stored in a vector or a column of a dataframe
gene_names <- rownames(counts)  # Example assuming the gene names are in the rownames

# Retrieve the annotation for each gene name
annot <- getBM(attributes = c("external_gene_name", "gene_biotype"),
               filters = "external_gene_name", values = gene_names, mart = mart)

# To differentiate mRNA and lncRNA, we can filter by gene biotype
mRNA_genes <- annot[annot$gene_biotype == "protein_coding", ]
lncRNA_genes <- annot[annot$gene_biotype == "lncRNA", ]

mrna_counts <- counts[mRNA_genes$external_gene_name,]
lncrna_counts <- counts[lncRNA_genes$external_gene_name,]

# Normalize using TMM (Trimmed Mean of M-values) from edgeR package
library(edgeR)

# Create a DGEList object
dge <- DGEList(counts = mrna_counts)

# Normalize using TMM
dge <- calcNormFactors(dge)
mRNA_counts_tmm <- cpm(dge, normalized.lib.sizes = TRUE)
mRNA.expr <- mRNA_counts_tmm

# Repeat for lncRNA counts if needed
dge_lncRNA <- DGEList(counts = lncrna_counts)
dge_lncRNA <- calcNormFactors(dge_lncRNA)
lncRNA_counts_tmm <- cpm(dge_lncRNA, normalized.lib.sizes = TRUE)
lncRNA.expr <- lncRNA_counts_tmm

# 5. Query and download DNA methylation data (HM450K β-values)
query.met <- GDCquery(
  project = "TCGA-SKCM",
  data.category = "DNA Methylation",
  data.type = "Methylation Beta Value",
  platform = "Illumina Human Methylation 450"
)
GDCdownload(query.met)
se.met <- GDCprepare(query.met, summarizedExperiment = TRUE)
# Extract beta-value matrix (probes x samples)
beta_mat <- assay(se.met)

# (Optional) remove probes with all-NA across samples
beta_mat <- beta_mat[rowSums(is.na(beta_mat)) < ncol(beta_mat), ]
# Load required packages
library(minfi)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
# Load annotation data
annotation <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
# Match probes to genes
gene_annotation <- annotation[, c("Name", "UCSC_RefGene_Name")]
gene_annotation <- na.omit(gene_annotation)
gene_annotation$Gene <- sapply(strsplit(gene_annotation$UCSC_RefGene_Name, ";"), `[`, 1)  # Handle multiple genes per probe

# Subset to probes present in your beta matrix
probes_in_matrix <- intersect(rownames(beta_mat), rownames(gene_annotation))
# beta_subset <- beta_mat[probes_in_matrix, ]
gene_map <- gene_annotation[probes_in_matrix, ]

beta_df <- as.data.frame(beta_mat)
beta_df$Gene <- gene_map$Gene

beta_gene_matrix <- beta_df %>%
  group_by(Gene) %>%
  summarise(across(where(is.numeric), mean, na.rm = TRUE)) %>%
  as.data.frame()

beta_gene_matrix <- beta_gene_matrix[!is.na(beta_gene_matrix$Gene) & beta_gene_matrix$Gene != "", ]

# Set rownames to gene names
rownames(beta_gene_matrix) <- beta_gene_matrix$Gene
beta_gene_matrix$Gene <- NULL
colnames(beta_gene_matrix) <- substr(colnames(beta_gene_matrix), 1, 12)
meth.beta <- beta_gene_matrix

### 6. Mutation data
library(maftools)
snv_query <- GDCquery(
  project = "TCGA-SKCM",
  data.category = "Simple Nucleotide Variation",
  data.type = "Masked Somatic Mutation",
  access = "open"
)

GDCdownload(snv_query)

maf <- GDCprepare(snv_query)
maf_object <- read.maf(maf = maf)

library(MAFDash)
mutation_matrix <- createOncoMatrix(maf_object,g=(gene_map@listData$Gene))$numericMatrix

# Binary Table
mutation_matrix <- ifelse(!is.na(mutation_matrix) & mutation_matrix > 0, 1, 0)
mutation_matrix <- as.data.frame(mutation_matrix)
colnames(mutation_matrix) <- substr(colnames(mutation_matrix), 1, 12)

# miRNA Data
library(SummarizedExperiment)
query_mirna <- GDCquery(project = "TCGA-SKCM",
                        data.category = "Transcriptome Profiling",
                        data.type = "miRNA Expression Quantification")
GDCdownload(query_mirna)
miRNA_data <- GDCprepare(query_mirna, summarizedExperiment = FALSE)

# Step 1: Keep miRNA_ID + every 3rd column starting from the 3rd
rpm_columns <- c(1, seq(3, ncol(miRNA_data), by = 3))  # 1 = miRNA_ID, rest are RPMs
# Step 2: Subset the dataframe
miRNA_rpm_data <- miRNA_data[, rpm_columns]

# Step 3: Clean the column names
# Remove the prefix 'reads_per_million_miRNA_mapped_' from each column name
colnames(miRNA_rpm_data) <- c(
  "miRNA_ID",
  gsub("^reads_per_million_miRNA_mapped_", "", colnames(miRNA_rpm_data)[-1])
)

# Step 4: Set miRNA_ID as rownames and convert to matrix
rownames(miRNA_rpm_data) <- miRNA_rpm_data$miRNA_ID
miRNA_rpm_matrix <- as.matrix(miRNA_rpm_data[ , -1])

# Optional: log2 transform
miRNA.expr <- log2(miRNA_rpm_matrix + 1)
colnames(miRNA.expr)<- substr(colnames(miRNA.expr), 1, 12)

clin_data <- GDCquery_clinic("TCGA-SKCM", type = "clinical")
clin.info <- data.frame(
  # Extract patient IDs (typically in format "TCGA-XX-XXXX")
  fustat = ifelse(!is.na(clin_data$vital_status) & clin_data$vital_status == "Dead", 1, 0),
  futime = ifelse(!is.na(clin_data$vital_status) & clin_data$vital_status == "Dead", 
                   clin_data$days_to_death, 
                   clin_data$days_to_last_follow_up),
  pstage = substr(clin_data$ajcc_pathologic_t, start = 1, stop = 2),
  age = as.numeric(round(clin_data$age_at_diagnosis/365))
)
rownames(clin.info) <- clin_data$submitter_id

skcm.tcga <- list(
  mRNA.expr  = as.data.frame(mRNA.expr),
  lncRNA.expr= as.data.frame(lncRNA.expr),
  miRNA.expr = as.data.frame(miRNA.expr),
  meth.beta  = meth.beta,
  mut.status = mutation_matrix,
  count   = counts,
  fpkm = fpkm,
  maf = maf,
  clin.info = clin.info
)

save(skcm.tcga, file = "skcm.tcga.RData")