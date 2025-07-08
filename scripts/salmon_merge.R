#!/usr/bin/env Rscript

# ---------------------
# Salmon merge script
# Usage:
# Rscript salmon_promoter_merge.R <out_dir> <promoter_rds> <norm_method> <samples> <conditions> <newnames> <count_files...>
# ---------------------

library(proActiv)   # 1.16
library(DESeq2)
library(edgeR)
library(dplyr)

# --------------------
# Parse arguments
# --------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 9)
  stop("Usage: Rscript salmon_promoter_merge.R <out_dir> <promoter_rds> <norm_method> <samples> <conditions> <newnames> <batch_str> <fit_script> <count_files...>")

out_dir      <- args[1]
prom_rds     <- args[2]              # promoter annotation
norm_method  <- args[3]
samples      <- strsplit(args[4], " ")[[1]]
conditions_str   <- args[5]
newnames     <- strsplit(args[6], " ")[[1]]
batch_str       <- args[7]
fit_script <- args[8]  # Path to the dexseq_fit script
count_files  <- args[9:length(args)]

conditions   <- strsplit(conditions_str, ",")[[1]]
batch       <- strsplit(batch_str, ",")[[1]]

stopifnot(length(samples) == length(conditions),
          length(samples) == length(newnames),
          length(samples) == length(count_files),
          file.exists(prom_rds),
          all(file.exists(count_files)))

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# --------------------
# Annotation & raw counts
# --------------------
promAnno <- readRDS(prom_rds)
# Read each RDS file in count_files into a list of matrices
cnt_list <- lapply(count_files, readRDS)
# Extract all unique promoter IDs across all count matrices
prom_ids <- unique(unlist(lapply(cnt_list, rownames)))
# Construct an empty matrix to hold promoter counts
# Rows: all unique promoter IDs
# Columns: one for each sample, using names from newnames
cnt_mat  <- matrix(0, nrow = length(prom_ids), ncol = length(samples),
                   dimnames = list(prom_ids, newnames))
# Populate the matrix with values from each sample's count object
# rownames(cnt_list[[i]]) gives the promoter IDs present in that sample.
# cnt_list[[i]][, 1] extracts the count values (first column).
# These are assigned into the corresponding rows and column in cnt_mat.
for (i in seq_along(cnt_list))
  cnt_mat[rownames(cnt_list[[i]]), newnames[i]] <- cnt_list[[i]][, 1]
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
  # Export size for NB fit test
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

# --------------------
# Absolute promoter activity
# --------------------
# proActiv 1.16 helper returns data.frame: columns = promoterId / geneId / samples
# getAbsolutePromoterActivity is from proActiv 1.18 but version used is 1.16 so should do proActiv:: instead of calling directly
# Cannot use proActiv 1.18 because conda cannot install that version (time attempted is 2025/6/18) - 
# getting gene expression using this method may result in shortened file (dim 7*7)
absDF  <- proActiv:::getAbsolutePromoterActivity(norm_counts, promAnno)
numCol <- setdiff(names(absDF), c("promoterId","geneId"))
absMat <- as.matrix(absDF[ , numCol])
rownames(absMat) <- absDF$promoterId
saveRDS(absMat, file.path(out_dir, "absolute_promoter_activity.rds"))

# --------------------
# promoter → gene (sum)
# --------------------
# promoterId, geneId, transcriptName
map <- promAnno@promoterIdMapping
# drop version (ENSG000001.1 → ENSG000001)
map$geneClean <- sub("\\..*$", "", map$geneId)
# length = nPromoter
gene_vec <- map$geneClean[match(rownames(absMat), map$promoterId)]
# sum promoters → gene matrix
gene_sum <- rowsum(absMat, group = gene_vec)
geneExpr <- gene_sum[match(gene_vec, rownames(gene_sum)), , drop = FALSE]
rownames(geneExpr) <- rownames(absMat)
saveRDS(geneExpr, file.path(out_dir, "gene_expression.rds"))

# relative activity
relMat <- absMat / geneExpr
relMat[!is.finite(relMat)] <- NA
saveRDS(relMat, file.path(out_dir, "relative_promoter_activity.rds"))

# -------------------------
# Call dexseq_fit2.R and pass arguments
# -------------------------


system2("Rscript", args = c(
    fit_script, 
    out_dir, 
    file.path(out_dir, "promoter_counts.rds"),
    file.path(out_dir, "NB_size_per_gene.rds"),
    file.path(out_dir, "promoter_fitted_mu.rds"),
    conditions_str, 
    batch_str))


# --------------------
# Condition means
# --------------------
mean_by_cond <- function(mat)
  sapply(unique(conditions), function(cd)
    rowMeans(mat[, newnames[conditions == cd], drop = FALSE], na.rm = TRUE))

saveRDS(mean_by_cond(absMat),  file.path(out_dir, "absolute_promoter_activity_mean.rds"))
saveRDS(mean_by_cond(relMat),  file.path(out_dir, "relative_promoter_activity_mean.rds"))
saveRDS(mean_by_cond(geneExpr), file.path(out_dir, "gene_expression_mean.rds"))

# --------------------
# Classification Major / Minor / Inactive
# --------------------
#proActiv result only includes detected promoters, so must keep map the same length
#This ensures we do not classify promoters not detected, thus NA, as inactive
map <- map[map$promoterId %in% rownames(absMat), ]


threshold_minor <- 0.25   # proActiv's threshold for active gene

# compute average promoter activity across all samples to determine which promoters are active
avgAct <- rowMeans(absMat, na.rm = TRUE)
# add average activity to the map, matching the corresponding promoter by promoterId
map$AverageActivity <- avgAct[match(map$promoterId, rownames(absMat))]

# for promoters belonging to the same gene:
# classify all as inactive if all values are na or max is inactive
# the highest activity promoter is major, ensured there is only one per gene using which.max
# the rest of active promoters are minor promoters
# assign not classified as lab's default, inactive
classify <- function(v) {
  if (all(is.na(v)) || max(v, na.rm = TRUE) < threshold_minor)
    return(rep("Inactive", length(v)))
  lab <- rep("Inactive", length(v))
  major_idx <- which.max(v)          
  lab[major_idx] <- "Major"
  lab[-major_idx][v[-major_idx] >= threshold_minor] <- "Minor"
  lab
}

# group promoters by gene and classify them
map$PromoterCategory <- unlist(tapply(map$AverageActivity,
                                      map$geneClean, classify))

saveRDS(map[, c("promoterId","geneId","AverageActivity","PromoterCategory")],
        file.path(out_dir, "absolute_promoter_activity_category.rds"))

message("Done. Output in ", out_dir)


