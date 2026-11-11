# Install if not already
install.packages("BiocManager")
BiocManager::install(c("GEOquery", "sva", "limma", "Biobase"))

# Load libraries
library(GEOquery)
library(sva)
library(limma)
library(Biobase)
library(ggplot2)
library(dplyr)
library(sva)
library(ggplot2)

# Function to map probe IDs to gene symbols and collapse probes
get_gene_mapped_exprs <- function(eset, gpl_id) {
  # Get expression matrix
  expr <- exprs(eset)
  
  # Download platform annotation
  gpl <- getGEO(gpl_id, AnnotGPL = TRUE)
  gpl_table <- Table(gpl)
  
  # Find possible gene symbol column
  possible_cols <- grep("gene.*symbol|symbol", colnames(gpl_table), ignore.case = TRUE, value = TRUE)
  gene_col <- possible_cols[1]
  if (is.na(gene_col)) stop(paste("No gene symbol column found in", gpl_id))
  
  message("Using gene symbol column: ", gene_col)
  
  # Ensure 'ID' is a column
  if (!"ID" %in% colnames(gpl_table)) {
    gpl_table$ID <- rownames(gpl_table)
  }
  
  # Create mapping
  mapping <- gpl_table[, c("ID", gene_col)]
  colnames(mapping) <- c("ID", "GeneSymbol")
  mapping <- mapping[!is.na(mapping$GeneSymbol) & mapping$GeneSymbol != "", ]
  mapping <- mapping[!duplicated(mapping$ID), ]
  
  # Subset expression matrix
  expr <- expr[rownames(expr) %in% mapping$ID, ]
  
  # Merge and collapse
  expr_df <- as.data.frame(expr)
  expr_df$ProbeID <- rownames(expr_df)
  expr_annot <- merge(expr_df, mapping, by.x = "ProbeID", by.y = "ID")
  
  expr_collapsed <- expr_annot %>%
    select(-ProbeID) %>%
    group_by(GeneSymbol) %>%
    summarise_all(mean, na.rm = TRUE)
  
  expr_mat <- as.matrix(expr_collapsed[,-1])
  rownames(expr_mat) <- expr_collapsed$GeneSymbol
  return(expr_mat)
}

# Download datasets
gse15605 <- getGEO("GSE15605", GSEMatrix = TRUE)[[1]]
gse22153 <- getGEO("GSE22153", GSEMatrix = TRUE)[[1]]
gse46517 <- getGEO("GSE46517", GSEMatrix = TRUE)[[1]]
gse54467 <- getGEO("GSE54467", GSEMatrix = TRUE)[[1]]

# Map each to gene symbols and collapse
exprs_list <- list(
  GSE15605 = get_gene_mapped_exprs(gse15605, "GPL570"),
  GSE22153 = get_gene_mapped_exprs(gse22153, "GPL6102"),
  GSE46517 = get_gene_mapped_exprs(gse46517, "GPL96"),
  GSE54467 = get_gene_mapped_exprs(gse54467, "GPL6884")
)

exprs_list[["GSE22153"]] <- log2(exprs_list[["GSE22153"]] + abs(min(exprs_list[["GSE22153"]])) + 1)

# Apply log2 if needed
# exprs_list <- lapply(exprs_list, function(mat) {
#   if (max(mat, na.rm = TRUE) > 100 && all(mat >= 0, na.rm = TRUE)) {
#     log2(mat + 1)
#   } else {
#     mat
#   }
# })

# Find common genes

common_genes <- Reduce(intersect, lapply(exprs_list, rownames))
exprs_list <- lapply(exprs_list, function(mat) mat[common_genes, ])
combined_expr <- do.call(cbind, exprs_list)  # genes x all samples
batch_vector <- rep(names(exprs_list), sapply(exprs_list, ncol))

pca_res <- prcomp(t(combined_expr), scale. = TRUE)

pca_df <- data.frame(PC1 = pca_res$x[,1],
                     PC2 = pca_res$x[,2],
                     Batch = batch_vector)

hulls <- pca_df %>%
  group_by(Batch) %>%
  slice(chull(PC1, PC2))

ggplot(pca_df, aes(x = PC1, y = PC2, color = Batch)) +
  geom_point(size = 2) +
  geom_polygon(data = hulls, aes(fill = Batch), alpha = 0.2, color = NA) +
  theme_minimal() +
  labs(title = "PCA before ComBat batch correction")


combat_expr <- ComBat(dat = as.matrix(combined_expr),
                      batch = batch_vector,
                      par.prior = TRUE,
                      prior.plots = FALSE)

pca_res <- prcomp(t(combat_expr), scale. = TRUE)

pca_df <- data.frame(PC1 = pca_res$x[,1],
                     PC2 = pca_res$x[,2],
                     Batch = batch_vector)

hulls <- pca_df %>%
  group_by(Batch) %>%
  slice(chull(PC1, PC2))

ggplot(pca_df, aes(x = PC1, y = PC2, color = Batch)) +
  geom_point(size = 2) +
  geom_polygon(data = hulls, aes(fill = Batch), alpha = 0.2, color = NA) +
  theme_minimal() +
  labs(title = "PCA after ComBat batch correction")

expr.z <- t(scale(t(combat_expr)))
write.csv(as.data.frame(expr.z), file = "scaled_expression_data.csv", row.names = TRUE)
# ------------------------------------------------------------------