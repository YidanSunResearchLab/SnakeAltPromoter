#!/usr/bin/env Rscript

library(proActiv)
library(DESeq2)
library(edgeR)
library(reshape2)

# Custom normalization function
normalizePromoterReadCounts <- function(promoterReadCounts, sizefactor = NULL) {
  if (ncol(promoterReadCounts) == 1) return(promoterReadCounts)
  activePromoters <- which(!is.na(promoterReadCounts[, 1]))
  colData <- data.frame(sampleLabels = colnames(promoterReadCounts))
  rownames(colData) <- colnames(promoterReadCounts)
  message("Calculating normalized read counts...")
  dds <- DESeq2::DESeqDataSetFromMatrix(countData = promoterReadCounts[activePromoters, ], colData = colData, design = ~1)
  dds <- DESeq2::estimateSizeFactors(dds)
  if (!is.null(sizefactor)) {
    DESeq2::sizeFactors(dds) <- 1 / sizefactor
    message("Using external size factors")
  }
  promoterReadCounts[activePromoters, ] <- DESeq2::counts(dds, normalized = TRUE)
  return(promoterReadCounts)
}

# Parse arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 8) {
  stop("Usage: Rscript proactiv_merge.R <output_dir> <promoter_rds> <sizefactor> <norm_method> <samples> <conditions> <newnames> <count_files>")
}
output_dir <- args[1]
promoter_rds <- args[2]
sizefactor_file <- args[3]
norm_method <- args[4]
samples <- strsplit(args[5], " ")[[1]]
conditions <- strsplit(args[6], " ")[[1]]
newnames <- strsplit(args[7], " ")[[1]]
count_files <- args[8:length(args)]

# Validate inputs
if (!file.exists(promoter_rds)) stop("Promoter RDS file does not exist")
if (!file.exists(sizefactor_file)) stop("Size factor file does not exist")
if (length(samples) != length(conditions) || length(samples) != length(newnames)) stop("Samples, conditions, and newnames length mismatch")
if (length(count_files) != length(samples)) stop("Count files and samples length mismatch")
if (!(norm_method %in% c("deseq2", "edger"))) stop("norm_method must be 'deseq2' or 'edger'")

# Create output directory
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Read promoter annotation
promoterAnnotationData <- readRDS(promoter_rds)

# Sample info: may need to comment all
sample_info <- data.frame(samplename = samples, condition = conditions, newnamegender = newnames, stringsAsFactors = FALSE)
rownames(sample_info) <- samples

# Read and merge counts
count_list <- lapply(count_files, readRDS)
promoterCounts <- do.call(cbind, count_list)
colnames(promoterCounts) <- sample_info$newnamegender
saveRDS(promoterCounts, file = file.path(output_dir, "promoter_counts.rds"))

# Read size factors: may need to comment all
sizefactor <- read.delim(sizefactor_file, sep = ":", stringsAsFactors = FALSE, row.names = 1)
colnames(sizefactor) <- "sizeFactor"
sizefactor <- sizefactor[sample_info$newnamegender, , drop = FALSE]

# Normalize counts
if (norm_method == "deseq2") {
  normalizedPromoterCounts <- normalizePromoterReadCounts(promoterCounts, sizefactor = sizefactor$sizeFactor)
} else {
  cds <- DGEList(counts = as.matrix(promoterCounts[which(!is.na(promoterCounts[, 1])), ]), group = factor(1:ncol(promoterCounts)))
  cds <- calcNormFactors(cds, method = "TMM")
  cds$samples$norm.factors <- 1 / sizefactor[rownames(cds$samples), 1]
  normalizedPromoterCounts <- cpm(cds, normalized.lib.sizes = TRUE)
}
saveRDS(normalizedPromoterCounts, file = file.path(output_dir, "normalized_promoter_counts.rds"))

# Calculate absolute promoter activity
absolutePromoterActivity <- getAbsolutePromoterActivity(normalizedPromoterCounts, promoterAnnotationData)
saveRDS(absolutePromoterActivity, file = file.path(output_dir, "absolute_promoter_activity.rds"))

# Calculate gene expression
geneExpression <- getGeneExpression(absolutePromoterActivity)
saveRDS(geneExpression, file = file.path(output_dir, "gene_expression.rds"))

# Calculate relative promoter activity
relativePromoterActivity <- getRelativePromoterActivity(absolutePromoterActivity, geneExpression)
saveRDS(relativePromoterActivity, file = file.path(output_dir, "relative_promoter_activity.rds"))

# Merge replicates with mean
abs_mean <- matrix(NA, nrow = nrow(absolutePromoterActivity), ncol = length(unique(conditions)))
colnames(abs_mean) <- unique(conditions)
for (cond in unique(conditions)) {
  cond_samples <- sample_info$newnamegender[sample_info$condition == cond]
  abs_mean[, cond] <- rowMeans(absolutePromoterActivity[, cond_samples, drop = FALSE], na.rm = TRUE)
}
saveRDS(abs_mean, file = file.path(output_dir, "absolute_promoter_activity_mean.rds"))

rel_mean <- matrix(NA, nrow = nrow(relativePromoterActivity), ncol = length(unique(conditions)))
colnames(rel_mean) <- unique(conditions)
for (cond in unique(conditions)) {
  cond_samples <- sample_info$newnamegender[sample_info$condition == cond]
  rel_mean[, cond] <- rowMeans(relativePromoterActivity[, cond_samples, drop = FALSE], na.rm = TRUE)
}
saveRDS(rel_mean, file = file.path(output_dir, "relative_promoter_activity_mean.rds"))

gene_mean <- matrix(NA, nrow = nrow(geneExpression), ncol = length(unique(conditions)))
colnames(gene_mean) <- unique(conditions)
for (cond in unique(conditions)) {
  cond_samples <- sample_info$newnamegender[sample_info$condition == cond]
  gene_mean[, cond] <- rowMeans(geneExpression[, cond_samples, drop = FALSE], na.rm = TRUE)
}
saveRDS(gene_mean, file = file.path(output_dir, "gene_expression_mean.rds"))

# Promoter classification
absolutePromoterActivity[is.na(absolutePromoterActivity)] <- 0
absolutePromoterActivity$AverageActivity <- rowMeans(absolutePromoterActivity[, sample_info$newnamegender], na.rm = TRUE)

major_minor_promoter_classification <- function(x) {
  gene_id <- absolutePromoterActivity$geneId
  gene_split <- split(x, gene_id)
  
  major_minor_value_change_per_gene <- function(z) {
    z <- as.numeric(z)
    if (max(z) < 0.25) {
      z[z < 0.25] <- 0
    } else {
      z[z < 0.25] <- 0
      z[z == max(z)] <- 100000
      z[z >= 0.25 & z < max(z)] <- 1000
    }
    z[z == 100000] <- "Major"
    z[z == 1000] <- "Minor"
    z[z == 0] <- "Inactive"
    return(z)
  }
  
  category <- lapply(gene_split, major_minor_value_change_per_gene)
  return(as.character(reshape2::melt(category)[, 1]))
}

absolutePromoterActivity$PromoterCategory <- major_minor_promoter_classification(absolutePromoterActivity$AverageActivity)
absolutePromoterActivityCategory <- absolutePromoterActivity[, c("promoterId", "geneId", "AverageActivity", "PromoterCategory")]
saveRDS(absolutePromoterActivityCategory, file = file.path(output_dir, "absolute_promoter_activity_category.rds"))