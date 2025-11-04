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
count_rds   <- args[2]
size_rds    <- args[3]
mu_rds      <- args[4]
conditions  <- strsplit(args[5], ",")[[1]]
batch       <- strsplit(args[6], ",")[[1]]

count.mat <- readRDS(count_rds)
size.vec  <- readRDS(size_rds)
mu.mat    <- readRDS(mu_rds)

# Function to use raw counts, estimated NB size, and estimated NB mean to calculate p-values for each promoter
calculate_pnbinom_pvalues <- function(count.mat, mu.mat, size.vec) {
  pval_vec <- vector("list", nrow(count.mat))
  names(pval_vec) <- rownames(count.mat)

  for (i in seq_len(nrow(count.mat))) {
    counts <- as.numeric(count.mat[i, ])
    mu <- as.numeric(mu.mat[i, ])
    size <- size.vec[i]
    if (max(counts, na.rm = TRUE) <= 3 || all(is.na(mu)) || any(mu <= 0, na.rm = TRUE) || is.na(size)) {
      pval_vec[[i]] <- rep(NA, ncol(count.mat))
      next
    }
    p <- pnbinom(counts, mu = mu, size = size)
    pval_vec[[i]] <- p
  }
  pval_mat <- do.call(rbind, pval_vec)
  rownames(pval_mat) <- rownames(count.mat)
  colnames(pval_mat) <- colnames(count.mat)
  saveRDS(pval_mat, file.path(out_dir, "per_sample_pnbinom_pvalues_matrix.rds"))
  return(pval_mat)
}

pval_mat <- calculate_pnbinom_pvalues(count.mat, mu.mat, size.vec)

calculate_gof_by_condition <- function(pval_mat, conditions) {
  names.conditions <- unique(conditions)
  gof_pvals <- matrix(NA, nrow = nrow(pval_mat), ncol = length(names.conditions))
  rownames(gof_pvals) <- rownames(pval_mat)
  colnames(gof_pvals) <- names.conditions
  for (j in seq_along(names.conditions)) {
    condition <- names.conditions[j]
    idx <- which(conditions == condition)
    for (i in seq_len(nrow(pval_mat))) {
      pvals_gene <- as.numeric(pval_mat[i, idx])
      good_p <- pvals_gene[is.finite(pvals_gene) & !is.na(pvals_gene)]
      if (length(good_p) < 3 || length(unique(good_p)) <= 1) {
        gof_pvals[i, j] <- NA
      } else {
        gof_pvals[i, j] <- ks.test(good_p, "punif")$p.value
      }
    }
  }
  return(gof_pvals)
}

gof_pval_mat <- calculate_gof_by_condition(pval_mat, conditions)
cat("Goodness-of-fit p-values calculated for each promoter by condition:\n", head(gof_pval_mat), "\n")
check = mean(gof_pval_mat < 0.01, na.rm = TRUE) < 0.42
output = mean(gof_pval_mat < 0.01, na.rm = TRUE)
cat("Whether mean proportion of extreme p-values (p < 0.01) indicating the model does not fit NOT exceeding 0.42:", check, "\n")
cat("The mean proportion of extreme p-values is:", output, "\n")
saveRDS(gof_pval_mat, file.path(out_dir, "goodness_of_fit_ks_pvalues.rds"))

# -------------------------
# Function to perform a DESeq2-based sanity check with permuted condition labels
count.mat <- round(count.mat)
count.mat <- count.mat[complete.cases(count.mat), ]
sanity_check_DESeq2 <- function(count.mat, conditions, batch) {
  conditions.perm <- sample(conditions)
  colData <- DataFrame(condition = conditions.perm, batch = batch, row.names = colnames(count.mat))
  if (length(unique(batch)) == 1) {
    dds <- DESeqDataSetFromMatrix(count.mat, colData, design = ~ condition)
  } else {
    dds <- DESeqDataSetFromMatrix(count.mat, colData, design = ~ batch + condition)
  }
  dds <- DESeq(dds, quiet = TRUE)
  res <- results(dds)
  num.dis <- sum(res$padj <= 0.05, na.rm = TRUE)
  return(num.dis)
}

# -------------------------
# Perform pairwise comparisons for DESeq2 sanity check and padj extraction
unique_conditions <- unique(conditions)
condition_pairs <- combn(unique_conditions, 2, simplify = FALSE)

for (pair in condition_pairs) {
  cond1 <- pair[1]
  cond2 <- pair[2]
  label <- paste0(cond1, "_vs_", cond2)
  cat("\nPerforming pairwise comparison:", label, "\n")

  idx <- which(conditions %in% c(cond1, cond2))
  cnt_sub <- count.mat[, idx]
  cond_sub <- conditions[idx]
  batch_sub <- batch[idx]

  sanity_res <- sanity_check_DESeq2(cnt_sub, cond_sub, batch_sub)
  cat("Sanity check false positive (", label, "):", sanity_res, "\n")
  saveRDS(sanity_res, file.path(out_dir, paste0("sanity_check_DESeq2_result_", label, ".rds")))

  if (!all(c(cond1, cond2) %in% colnames(gof_pval_mat))) {
    warning("One or both conditions not found in gof_pval_mat. Skipping: ", label)
    next
  }

  padj_vec <- apply(gof_pval_mat[, c(cond1, cond2), drop = FALSE], 1, min, na.rm = TRUE)
  res_df <- data.frame(promoterId = rownames(gof_pval_mat), padj = padj_vec)

  write.table(res_df, file = file.path(out_dir, paste0("deseq2_results_", label, ".tsv")),
              sep = "\t", quote = FALSE, row.names = FALSE)
}
