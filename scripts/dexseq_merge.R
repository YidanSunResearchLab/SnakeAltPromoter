#!/usr/bin/env Rscript

# -------------------------
# DEXSeq Promoter Merge Script
# -------------------------
# Usage:
# Rscript dexseq_promoter_merge.R <out_dir> <promoter_rds> <norm_method> <samples> <conditions> <newnames> <comparison> <count_files...>
# -------------------------

library(proActiv)  # v1.16
library(DESeq2)
library(edgeR)
library(dplyr)

# -------------------------
# Parse arguments
# -------------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 10) {
  stop("Usage: Rscript dexseq_promoter_merge.R <out_dir> <promoter_rds> <norm_method> <samples> <conditions> <newnames> <comparison> <batch> <fit_script> <count_files...>")
}

out_dir     <- args[1]
prom_rds    <- args[2]
norm_method <- args[3]
samples     <- strsplit(args[4], " ")[[1]]
conditions_str  <- args[5]
newnames    <- strsplit(args[6], " ")[[1]]
comparison  <- strsplit(args[7], ",")[[1]]
batch_str       <- args[8]
fit_script <- args[9]  # Path to the dexseq_fit script
count_files <- args[10:length(args)]


conditions  <- strsplit(conditions_str, ",")[[1]]
batch       <- strsplit(batch_str, ",")[[1]]
# -------------------------
# Validate
# -------------------------
stopifnot(length(samples) == length(conditions),
          length(samples) == length(newnames),
          length(samples) == length(count_files),
          file.exists(prom_rds),
          all(file.exists(count_files)))

# -------------------------
# Sample information
# -------------------------
sample_info <- data.frame(
  sample = samples,
  condition = conditions,
  newname = newnames,
  file = count_files,
  stringsAsFactors = FALSE
)

# Keep conditions for comparison
sample_info <- sample_info[sample_info$condition %in% comparison, ]
stopifnot(nrow(sample_info) > 0)
parent_dir <- dirname(dirname(out_dir))
bam_dir <- file.path(parent_dir, "bam")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# -------------------------
# Load promoter annotation and count data
# -------------------------
promAnno <- readRDS(prom_rds)
# Read each RDS file in count_files into a list of matrices
cnt_list <- lapply(sample_info$file, readRDS)

# Extract all unique promoter IDs across all count matrices
prom_ids <- unique(unlist(lapply(cnt_list, rownames)))

# Construct an empty matrix to hold promoter counts
# Rows: all unique promoter IDs
# Columns: one for each sample, using names from newnames
cnt_mat <- matrix(0, nrow = length(prom_ids), ncol = nrow(sample_info),
                  dimnames = list(prom_ids, sample_info$newname))

# Populate the matrix with values from each sample's count object
# rownames(cnt_list[[i]]) gives the promoter IDs present in that sample.
# cnt_list[[i]][, 1] extracts the count values (first column).
# These are assigned into the corresponding rows and column in cnt_mat.
for (i in seq_along(cnt_list)) {
  cnt_mat[rownames(cnt_list[[i]]), sample_info$newname[i]] <- cnt_list[[i]][, 1]
}

saveRDS(cnt_mat, file.path(out_dir, "promoter_counts.rds"))

# -------------------------
# Normalization
# -------------------------
if (norm_method == "deseq2") {
  dds <- DESeqDataSetFromMatrix(round(cnt_mat),
                                colData = data.frame(row.names = colnames(cnt_mat)),
                                design = ~ 1)
  dds <- DESeq(dds, quiet = TRUE)
  # Export normalized counts
  norm_counts <- counts(dds, normalized = TRUE)
  saveRDS(norm_counts, file.path(out_dir, "normalized_promoter_counts.rds")) 
  # Export size factors
  sizeFactors <- sizeFactors(dds) 
  saveRDS(sizeFactors, file.path(out_dir, "size_factors.rds"))  
  # Export estimated mu (fitted mean values for each promoter in each sample)
  mu_mat <- assays(dds)[["mu"]]
  saveRDS(mu_mat, file.path(out_dir, "promoter_fitted_mu.rds"))
  # Export dispersion
  disp_vec <- mcols(dds)$dispersion
  saveRDS(disp_vec, file.path(out_dir, "dispersion_per_gene.rds"))

  # Export size for NB
  size_vec <- 1 / disp_vec
  size_vec[is.na(size_vec) | !is.finite(size_vec) | size_vec <= 0] <- NA
  saveRDS(size_vec, file.path(out_dir, "NB_size_per_gene.rds"))
} else if (norm_method == "edger") {
 library(edgeR)
  dge <- DGEList(counts = cnt_mat)
  dge <- calcNormFactors(dge)
  # Get dispersion and mu
  design <- model.matrix(~1, data = dge$samples)
  dge <- estimateDisp(dge, design)
  fit <- glmFit(dge, design)
  # mu: Row is gene, column is sample
  mu_mat <- fit$fitted.values
  saveRDS(mu_mat, file.path(out_dir, "promoter_fitted_mu.rds"))
  # Normalized counts
  norm_counts <- cpm(dge, normalized.lib.sizes = TRUE)
  saveRDS(norm_counts, file.path(out_dir, "normalized_promoter_counts.rds"))
  # size factors
  sizeFactors <- dge$samples$norm.factors
  names(sizeFactors) <- colnames(cnt_mat)
  saveRDS(sizeFactors, file.path(out_dir, "size_factors.rds"))
  # Dispersion
  disp_vec <- dge$tagwise.dispersion
  saveRDS(disp_vec, file.path(out_dir, "dispersion_per_gene.rds"))
  # NB size = 1 / dispersion
  size_vec <- 1 / disp_vec
  size_vec[is.na(size_vec) | !is.finite(size_vec) | size_vec <= 0] <- NA
  saveRDS(size_vec, file.path(out_dir, "NB_size_per_gene.rds"))

} else stop("norm_method must be deseq2 or edger")

# -------------------------
# Save sizefactor
# -------------------------

# Conver to data.frame before merging
sizefactor_df <- data.frame(
  newname = names(sizeFactors),
  sizeFactor = as.numeric(sizeFactors),
  stringsAsFactors = FALSE
)
# Each sample has both size factor and sample
sizefactor_df <- merge(sizefactor_df, sample_info[, c("sample", "newname")], by = "newname")
write.table(sizefactor_df[, c("sample", "sizeFactor")],
            file = file.path(bam_dir, "dexseq_counts.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)


# -------------------------
# Promoter → Gene Expression
# -------------------------
# proActiv 1.16 helper returns data.frame: columns = promoterId / geneId / samples
# getAbsolutePromoterActivity is from proActiv 1.18 but version used is 1.16 so should do proActiv:: instead of calling directly
# cannot use proActiv 1.18 because conda cannot install that version (time attempted is 2025/6/18)
# getting gene expression using this method may result in shortened file (dim 7*7)
# Dataframe with promoterId, geneId, and sample columns
absDF  <- proActiv:::getAbsolutePromoterActivity(norm_counts, promAnno)
# Returns all sample columns
numCol <- setdiff(names(absDF), c("promoterId","geneId"))
# Matrix has rows as promoterId and columns as samples, values are activity
absMat <- as.matrix(absDF[ , numCol])
rownames(absMat) <- absDF$promoterId
saveRDS(absMat, file.path(out_dir, "absolute_promoter_activity.rds"))

absMat <- readRDS(file.path(out_dir, "absolute_promoter_activity.rds"))

# promoterId, geneId, transcriptName
map <- promAnno@promoterIdMapping
# Remove duplicate promoterId entries if any as duplicated may cause issues in mapping
n_dup <- sum(duplicated(map$promoterId))
if (n_dup > 0) {
  warning(sprintf("Found %d duplicated promoterId entries in promoterIdMapping. Only the first geneId was kept for each.", n_dup))
}
# Keeping the first occurrence, because exon bins are sorted by their position strand-aware
map_unique <- map[!duplicated(map$promoterId), c("promoterId", "geneId")]
# Drop version (ENSG000001.1 → ENSG000001)
map_unique$geneClean <- sub("\\..*$", "", map_unique$geneId)
# Keep only promoters that are present in absMat and add average activity for promoters
map_unique <- map_unique[map_unique$promoterId %in% rownames(absMat), ]
avgAct <- rowMeans(absMat, na.rm = TRUE)
map_unique$AverageActivity <- avgAct[match(map_unique$promoterId, rownames(absMat))]

# Match promoterId to geneId
gene_vec <- map_unique$geneClean[match(rownames(absMat), map_unique$promoterId)]
if (anyNA(gene_vec)) {
  stop("Some promoterId in absMat not found in map, check mapping.")
}
# Sum promoters of the same gene and create gene expression matrix of gene by sample
# The order of rows is that of absMat
gene_sum <- rowsum(absMat, group = gene_vec)
geneExpr <- gene_sum[match(gene_vec, rownames(gene_sum)), , drop = FALSE]
rownames(geneExpr) <- rownames(absMat)
saveRDS(geneExpr, file.path(out_dir, "gene_expression.rds"))

# -------------------------
# Relative Activity
# -------------------------
relMat <- absMat / geneExpr
relMat[!is.finite(relMat)] <- NA
saveRDS(relMat, file.path(out_dir, "relative_promoter_activity.rds"))

# -------------------------
# Call dexseq_fit.R and pass arguments
# -------------------------

system2("Rscript", args = c(
    fit_script, 
    out_dir, 
    file.path(out_dir, "promoter_counts.rds"),
    file.path(out_dir, "NB_size_per_gene.rds"),
    file.path(out_dir, "promoter_fitted_mu.rds"),
    conditions_str, 
    batch_str))

# -------------------------
# Condition Means
# -------------------------
mean_by_cond <- function(mat)
  sapply(unique(sample_info$condition), function(cd)
    rowMeans(mat[, sample_info$newname[sample_info$condition == cd], drop = FALSE], na.rm = TRUE))

saveRDS(mean_by_cond(absMat),  file.path(out_dir, "absolute_promoter_activity_mean.rds"))
saveRDS(mean_by_cond(relMat),  file.path(out_dir, "relative_promoter_activity_mean.rds"))
saveRDS(mean_by_cond(geneExpr), file.path(out_dir, "gene_expression_mean.rds"))

# -------------------------
# Classification
# -------------------------

# Cell-line specific classification (same logic as merge, per geneClean)
threshold_minor <- 0.25

classify_per_gene <- function(df) {
  # Retrieve the average promoter activity values for all promoters of this gene
  v <- df$AverageActivity
  # If no activity is above the threshold, maintain default as "Inactive" and exit function
  if (all(is.na(v)) || max(v, na.rm = TRUE) < threshold_minor) {
    df$PromoterCategory <- "Inactive"
  } else {
    # Default all promoters to "Inactive"
    df$PromoterCategory <- "Inactive"
    # Only one active promoter per gene with the highest activity is classified as "Major"
    major_idx <- which.max(v)
    df$PromoterCategory[major_idx] <- "Major"
    # The rest of the promoters with activity above the threshold are classified as "Minor"
    minor_idx <- which(v >= threshold_minor & seq_along(v) != major_idx)
    df$PromoterCategory[minor_idx] <- "Minor"
  }
  return(df)
}

# Apply classification to each group of promoters sharing the same geneClean identifier 
# Classify all promoters from activity of its corresponding gene
map_classified <- do.call(rbind, by(map_unique, map_unique$geneClean, classify_per_gene))

saveRDS(map_classified[, c("promoterId", "geneId", "AverageActivity", "PromoterCategory")],
        file.path(out_dir, "absolute_promoter_activity_category.rds"))
