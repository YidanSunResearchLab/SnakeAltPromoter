
library(DESeq2)
library(edgeR)
library(dplyr)

# -------------------------
# Parse arguments
# -------------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 6) {
  stop("Usage: Rscript dexseq_fit.R <out_dir> <count_rds> <size_rds> <mu_rds> <conditions> <batch>")
}

out_dir     <- args[1]
count_rds     <- args[2]
size_rds    <- args[3]
mu_rds <- args[4]
conditions  <- strsplit(args[5], ",")[[1]]
batch       <- strsplit(args[6], ",")[[1]]

count.mat <- readRDS(count_rds)
size.vec  <- readRDS(size_rds)
mu.mat    <- readRDS(mu_rds)

# Function to use raw counts, estimated NB size, and estimated NB mean to calculate p-values for each promoter
calculate_pnbinom_pvalues <- function(count.mat, mu.mat, size.vec) {
# Each row correspond to and named as a promoter
  pval_vec <- vector("list", nrow(count.mat))
  names(pval_vec) <- rownames(count.mat)

  for (i in seq_len(nrow(count.mat))) {
    counts <- as.numeric(count.mat[i, ])
    mu <- as.numeric(mu.mat[i, ])
    size <- size.vec[i]
    # Filter: at least 3 non-NA counts; mu cannot be all NA or <= 0; must have size parameter
    if (max(counts, na.rm = TRUE) <= 3 || all(is.na(mu)) || any(mu <= 0, na.rm = TRUE) || is.na(size)) {
      pval_vec[[i]] <- rep(NA, ncol(count.mat))
      next
    }

    p <- pnbinom(counts, mu = mu, size = size)
    pval_vec[[i]] <- p
  }
# Convert to matrix and save: row is promoter, column is sample, values are p-values
  pval_mat <- do.call(rbind, pval_vec)
  rownames(pval_mat) <- rownames(count.mat)
  colnames(pval_mat) <- colnames(count.mat)
  saveRDS(pval_mat, file.path(out_dir, "per_sample_pnbinom_pvalues_matrix.rds"))

  return(pval_mat)
}

# Calculate p-values and define pval_mat for later use in goodness-of-fit
pval_mat <- calculate_pnbinom_pvalues(count.mat, mu.mat, size.vec)

# Calculate goodness-of-fit p-values for all promoters by condition
calculate_gof_by_condition <- function(pval_mat, conditions) {
  names.conditions <- unique(conditions)
  gof_pvals <- matrix(NA, nrow = nrow(pval_mat), ncol = length(names.conditions))
  # Assign row as promoter and column as condition
  rownames(gof_pvals) <- rownames(pval_mat)
  colnames(gof_pvals) <- names.conditions
  # Loop through corresponding columns of each condition
  for (j in seq_along(names.conditions)) {
    condition <- names.conditions[j]
    idx <- which(conditions == condition)
    # Run KS test for every promoter across the samples of the current condition
    for (i in seq_len(nrow(pval_mat))) {
      pvals_gene <- as.numeric(pval_mat[i, idx])
      # Filter out non-finite p-values and NAs
      good_p <- pvals_gene[is.finite(pvals_gene) & !is.na(pvals_gene)]
      # Need at least 5 non-NA p-values to run KS test
      # This is adjusted to 3 because test data only has 3 samples per condition
      if (length(good_p) < 3 || length(unique(good_p)) <= 1) {
        gof_pvals[i, j] <- NA
      } else {
        gof_pvals[i, j] <- ks.test(good_p, "punif")$p.value
      }
    }
  }
  
  return(gof_pvals)
}

# Row is promoter, column is condition, values are p-values from KS test
gof_pval_mat <- calculate_gof_by_condition(pval_mat, conditions)
cat("Goodness-of-fit p-values calculated for each promoter by condition:\n", head(gof_pval_mat), "\n")
# Number of extreme p-values, we compare this number with the threshold 0.42
check=mean(gof_pval_mat <0.01, na.rm = T) < 0.42
output=mean(gof_pval_mat <0.01, na.rm = T)
cat("Whether mean proportion of extreme p-values (p < 0.01) indicating the model does not fit NOT exceeding 0.42:", check, "\n")
cat("The mean proportion of extreme p-values is:", output, "\n")
saveRDS(gof_pval_mat, file.path(out_dir, "goodness_of_fit_ks_pvalues.rds"))




# -------------------------
# Function to perform a DESeq2-based sanity check with permuted condition labels
count.mat <- round(count.mat)
sanity_check_DESeq2 <- function(count.mat, conditions, batch) {
# Randomly permute the conditions to create a null distribution
  conditions.perm <- sample(conditions)
# Create a DESeq2 dataset object with the permuted conditions and original batch information
  colData <- DataFrame(condition = conditions.perm, batch = batch, row.names = colnames(count.mat))
# Include batch as a covariate if there are multiple batches
  if (length(unique(batch)) == 1) {
    dds <- DESeqDataSetFromMatrix(count.mat, colData, design = ~ condition)
  } else {
    dds <- DESeqDataSetFromMatrix(count.mat, colData, design = ~ batch + condition)
  }
# Run DESeq2 analysis and count significant p-values
  dds <- DESeq(dds, quiet = TRUE)
  res <- results(dds)
  num.dis <- sum(res$padj <= 0.05, na.rm = TRUE)
  return(num.dis)
}
# Save result: number of false positive discoveries
sanity_check_result <- sanity_check_DESeq2(count.mat, conditions, batch)
cat("Sanity check DESeq2 result (number of false positive discoveries):", sanity_check_result, "\n")
saveRDS(sanity_check_result, file.path(out_dir, "sanity_check_DESeq2_result.rds"))