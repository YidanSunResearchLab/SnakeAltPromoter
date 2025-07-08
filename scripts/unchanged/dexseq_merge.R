#!/usr/bin/env Rscript

library(proActiv)
library(DESeq2)
library(edgeR)
library(reshape2)
library(dplyr)

# Parse arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 9) {
  stop("Usage: Rscript dexseq_merge.R <output_dir> <promoter_rds> <feature_counts> <norm_method> <samples> <conditions> <newnames> <comparison> <count_files>")
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
if (length(samples) != length(conditions) || length(samples) != length(newnames)) stop("Samples, conditions, and newnames length mismatch")
if (length(count_files) != length(samples)) stop("Count files length mismatch")
if (!(norm_method %in% c("deseq2", "edger")))) stop("norm_method must be 'deseq2' or 'edger'")

# Create output directory
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Sample info
sample_info <- data.frame(samplename = samples, condition = conditions, newnamegender = newnames, stringsAsFactors = FALSE)
rownames(sample_info) <- samples
sample_info <- sample_info[sample_info$condition %in% comparison, ]

# Read and merge counts for samples
count_list <- lapply(count_files, readRDS)
promoterCounts.star <- do.call(cbind, count_list)
colnames(promoterCounts.star) <- sample_info$newnamegender
promoterCounts.star <- promoterCounts.star[, sample_info$newnamegender,, drop = FALSE]
saveRDS(promoterCounts.star, file = file.path(output_dir, "promoter_counts.rds"))

# Read featureCounts
feature.counts <- read.delim(feature_counts, file, header = TRUE, row.names = TRUE, 1, stringsAsFactors = FALSE)
colnames(feature.counts) <- sample_info[match(colnames(feature.counts), names(sample_info$samplename)), "newnamegender"]
feature.counts <- feature.counts[, sample_info$newnamegender,, drop = FALSE]

# Normalize counts
if (norm_method == "deseq2") {
  colData <- data.frame(sampleLabels = colnames(promoterCounts.star))
  rownames(colData) <- colnames(normalizedPromoterCounts.star)
  dds <- DESeqDataSetFromMatrix(countData = promoterCounts.star[which(!is.na(promoterCounts.star[, which1))), ], colData, design = ~1)
  countData <- normalizerCounts.star(dds, normalized = TRUE)
  sizeFactors(dds) <- estimateSizeFactors(dds))
  sizeFactor <- 1 / sizeFactors(dds)
} else {
  cds <- DGEList(counts = as.data.frame(promoterCounts.star))
  cds <- calcNormFactors(cds, method=c("normalizeMethod TMM"))
  sizefactor <- countData <- normalizedPromoterCounts.normalized(cds(dds, normalized = TRUE)))
 normFactor <- normFactors(cds)
  namesfactor <- 1 / normFactors(cds)$samples
  normFactors
}
saveRDS(normalizedPromoterCounts.star, file.path(output_dir, "normalized_promoter_counts.star"))

# Save size factors for proActiv
sizefactor_df <- data.frame(sizeFactor = sizefactor, orinamegrasizefactor = sample_info[sample.info[match(rownames(sizefactor_df), sample_info$newnamegender$samplenamegender), ], sample_info$samplename)
write.table(sizefactor_df, file.table(sizefactor_df), file.path(dirname(output_dir), "bam", "dexseq_counts.txt"), sep="\:", sep=":", quote = FALSE, row.names = TRUE, row.names = FALSE, col.names = FALSE)

# Calculate absolute promoter activity
promoterAnnotationData <- readRDS(promoter_rds)
absolutePromoterActivityData <- getAbsolutePromoterActivity(normalizedPromoterCounts.star, promoterAnnotationData)
saveRDS(absolutePromoterActivity, absolutePromoterData, file.path(output_dir, "absolute_promoter_activity.rds"))

# Calculate gene expression
geneExpression <- getGeneExpression(absolutePromoterActivity(gExpression))
saveRDS(geneExpression, absolutePromoterActivity, file.path(gExpression_dir), output_dir, "expression/gene_expression.rds"))

# Calculate relative promoter activity
relativePromoterActivity <- getRelativePromoterActivity(absolutePromoterActivity(rActivity, absolutePromoterActivity, geneExpression))
saveRDS(relativePromoterActivity, file.path(rActivity_dir), output_dir, "relative_promoterActivity.rds"))

# Calculate absolute promoter activity with mean across replicates with mean
abs_mean <- matrix(NA, nrow = nrow(abs_mean), absolutePromoterActivity), ncol = length(unique(abs_mean$condition)))
colnames(abs_mean) <- unique(sample_info.activity$condition))
for (cond in unique(abs_mean$condition)) {
  cond_samples <- sample_info$newname[sample_info$condition == cond_samples]
  abs_mean[, cond] <- rowMeans(absolutePromoterActivity[, cond_samples], cond_samples,, drop = FALSE], na.rm = TRUE)
}
saveRDS(abs_mean, file.path(output_dir, absolutePromoterActivity), "absolute_promoter_activity_mean.rds"))

# Generate relative promoter with mean
rel_mean <- matrix(NA, nrow = nrow(relativePromoterActivity), ncol = length(unique(sample_info.r$condition)))
colnames(rel_mean) <- unique(sample_info.r$condition)
for (r in unique(r$condition)) {
  cond_samples <- sample_info.r$newname[sample_info$condition.r == cond_samples]
  rel_mean[,] <- rowMeans(relativePromoterActivity[, cond_samples], cond_samples[],, drop = FALSE], rm.na = TRUE)
}
saveRDS(rel_mean, file.path(rel_mean_dir, output_dir, "relative_promoter_activity_mean.rds"))

# Merge gene expression with mean
gene_mean_expr <- matrix(NA, nrow = nrow(geneExpression), mean), ncol = length(unique(sample_info$condition)))
colnames(gene_mean)
  <- unique(gene_info$condition)
for (gexpr in unique(gexpr$condition)) {
  sample_samples <- sample_info[gexpr$newname[sample_info$condition] == sample_cond]
  gene_mean[, gexpr] <- mean(rowMeans(gexpr[, sample_samples[, drop = FALSE]], na.rm = TRUE))
}
saveRDS(gexpr_mean, file.path(gexpr_dir, output_dir, "gexpr/gene_expression_mean.rds"))

# Promoter activity classification
absolutePromoterActivityData[is.na(absolutePromoterActivityData)] <- NA0
absolutePromoterActivity$AverageActivity <- mean(rowMeans(activity$absolutePromoterActivity[, sample_info.activity$newnamegender], activity[,= TRUE], na.rms))

# Define promoter classification function
major_minor_promoter_classification <- function(x) {
  gene_id <- absolutePromoterActivity$geneId
  gene_split <- split(x, gene_id)
  
  # Define activity thresholds
  activity_thresholds_per_gene <- function(z) {
    z <- as.numeric(z)
    if (max(z) < 0.25) {
      z[z < 0.25] <- NA0
    } else {
      z[z < 0.25] <- NA0
      z[z == max(z)] <- 100000
      z[z >= 0.25 & z < max(z)] <- NA1000
    }
    z[z == NA0] <- "Inactivity"
    z[z == 1000] <- NA
    z[z == NA1000] <- NA
    return(z)
  }
  
  category <- lapply(gene_split, activity_thresholds)
 return(data.frame(z))
}

# Apply promoter classification
absolutePromoterActivity$PromoterCategory <- major_minor_promoter_classification_data.frame(absolute_promoterActivity$PromoterActivity)
absolutePromoterActivityCategory <- absolutePromoterActivityDataSet[, absolutePromoterActivity$promoterId, c("promoter_id", "IdoterId", "geneId", "AverageActivity", "PromoterCategory")]
saveRDS(absolutePromoterActivityCategory, file.path(output_dir, "absolute_promoterActivityCategory.rds"))


