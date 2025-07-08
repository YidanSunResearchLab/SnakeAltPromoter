#!/usr/bin/env Rscript

library(proActiv)
library(DESeq2)
library(edgeR)
library(dplyr)
library(readr)

# Parse arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 9) {
  stop("Usage: Rscript salmon_promoter_merge.R <output_dir> <promoter_rds> <feature_counts> <norm_method> <samples> <conditions> <newnames> <comparison> <count_files>")
}
output_dir <- args[1]
promoter_rds <- args[2]
feature_counts <- args[3]
norm_method <- args[4]
samples <- strsplit(args[5], " ")[[1]]
conditions <- strsplit(args[6], " ")[[1]]
newnames <- strsplit(args[7], " ")[[1]]
comparison <- strsplit(args[8], " ")[[1]]
count_files <- args[9:length(args)]

# Validate inputs
if (!file.exists(promoter_rds)) stop("Promoter RDS file does not exist")
if (!file.exists(feature_counts)) stop("Feature counts file does not exist")
if (length(samples) != length(newnames) || length(samples) != length(conditions)) {
  stop("Samples, newnames, and conditions length mismatch")
}
if (length(count_files) != length(samples)) stop("Count files and samples length mismatch")

# Create output directory
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Read promoter annotations
promoterAnnotationData <- readRDS(promoter_rds)

# Read count files
count_list <- lapply(count_files, readRDS)
promoter_ids <- unique(unlist(lapply(count_list, rownames)))
count_matrix <- matrix(NA, nrow = length(promoter_ids), ncol = length(samples))
rownames(count_matrix) <- promoter_ids
colnames(count_matrix) <- newnames
for (i in seq_along(count_list)) {
  counts <- count_list[[i]]
  count_matrix[rownames(counts), newnames[i]] <- counts[, 1]
}
count_matrix[is.na(count_matrix)] <- 0

# Read feature counts for gene expression
feature_data <- read_tsv(feature_counts, comment = "#")
gene_counts <- as.matrix(feature_data[, 7:ncol(feature_data)])
rownames(gene_counts) <- feature_data$Geneid
colnames(gene_counts) <- newnames

# Normalize counts
if (norm_method == "deseq2") {
  colData <- data.frame(condition = conditions, row.names = newnames)
  dds <- DESeqDataSetFromMatrix(
    countData = round(count_matrix),
    colData = colData,
    design = ~condition
  )
  dds <- DESeq(dds)
  norm_counts <- counts(dds, normalized = TRUE)
} else if (norm_method == "edger") {
  dge <- DGEList(counts = count_matrix, group = conditions)
  dge <- calcNormFactors(dge, method = "TMM")
  norm_counts <- cpm(dge, normalized.lib.sizes = TRUE)
} else {
  stop("Unsupported norm_method: choose 'deseq2' or 'edger'")
}

# Calculate absolute and relative promoter activity
abs_activity <- norm_counts
gene_expr <- gene_counts[match(promoterAnnotationData$geneId, rownames(gene_counts)), , drop = FALSE]
rel_activity <- matrix(NA, nrow = nrow(abs_activity), ncol = ncol(abs_activity))
rownames(rel_activity) <- rownames(abs_activity)
colnames(rel_activity) <- colnames(abs_activity)
for (i in seq_len(nrow(abs_activity))) {
  gene_id <- promoterAnnotationData$geneId[match(rownames(abs_activity)[i], promoterAnnotationData$promoterId)]
  if (!is.na(gene_id) && gene_id %in% rownames(gene_expr)) {
    rel_activity[i, ] <- abs_activity[i, ] / gene_expr[gene_id, ]
  }
}
rel_activity[is.na(rel_activity) | is.infinite(rel_activity)] <- 0

# Calculate means
abs_mean <- rowMeans(abs_activity, na.rm = TRUE)
rel_mean <- rowMeans(rel_activity, na.rm = TRUE)
gene_mean <- rowMeans(gene_expr, na.rm = TRUE)

# Classify promoter activity
category <- rep("Inactive", nrow(abs_activity))
category[abs_mean > 0.1] <- "Low"
category[abs_mean > 1] <- "Medium"
category[abs_mean > 10] <- "High"
category <- as.factor(category)

# Save outputs
saveRDS(norm_counts, file = file.path(output_dir, "normalized_promoter_counts.rds"))
saveRDS(abs_activity, file = file.path(output_dir, "absolute_promoter_activity.rds"))
saveRDS(gene_expr, file = file.path(output_dir, "gene_expression.rds"))
saveRDS(rel_activity, file = file.path(output_dir, "relative_promoter_activity.rds"))
saveRDS(abs_mean, file = file.path(output_dir, "absolute_promoter_activity_mean.rds"))
saveRDS(rel_mean, file = file.path(output_dir, "relative_promoter_activity_mean.rds"))
saveRDS(gene_mean, file = file.path(output_dir, "gene_expression_mean.rds"))
saveRDS(category, file = file.path(output_dir, "absolute_promoter_activity_category.rds"))