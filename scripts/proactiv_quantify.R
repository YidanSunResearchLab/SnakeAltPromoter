#!/usr/bin/env Rscript

# ---------------------
# proActiv quantification script
# Usage:
# Rscript proactiv_counts.R <out_dir> <promoter_rds> <sj_files_str> <conditions_str> <reference_condition>
# ---------------------

library(proActiv)

# --------------------- #
# Parse and validate inputs
# --------------------- #

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 7) {
  stop("Usage: Rscript proactiv_counts.R <out_dir> <promoter_rds> <sj_files_str> <conditions_str> <reference_condition (can be empty string)> <batch_str> <fit_script>")
}
out_dir <- args[1]
promoter_rds <- args[2]
sj_files_str <- args[3]
conditions_str <- args[4]
reference_condition <- args[5]
batch_str       <- args[6]
fit_script <- args[7]  # Path to the dexseq_fit script

# Ensure input types are character strings
if (!is.character(sj_files_str)) stop("`sj_files_str` must be a character string.")
if (!is.character(conditions_str)) stop("`conditions_str` must be a character string.")

# Parse file paths and conditions
sj_files <- strsplit(sj_files_str, " ")[[1]]
condition <- strsplit(conditions_str, ",")[[1]]
batch <- strsplit(batch_str, ",")[[1]]

# Validate file existence
if (!file.exists(promoter_rds)) {
  stop("Promoter RDS file does not exist: ", promoter_rds)
}


missing_sj <- sj_files[!file.exists(sj_files)]
if (length(missing_sj) > 0) {
  message("The following SJ files are missing:")
  print(missing_sj)
  stop("Some SJ files do not exist.")
}

# Check consistency between SJ files and conditions
if (length(sj_files) != length(condition)) {
  stop("Number of SJ files does not match number of conditions.")
}

# Create out directory if needed
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# --------------------- #
# Run proActiv
# --------------------- #

# Load promoter annotation
promoterAnnotationData <- readRDS(promoter_rds)

# Estimate promoter activity and summarize results across conditions
result <- proActiv(
  files = sj_files,
  promoterAnnotation = promoterAnnotationData,
  condition = condition
)

# Filter out promoters with NA counts (e.g., from single-exon transcripts)
result <- result[complete.cases(assays(result)$promoterCounts), ]


# -------------------------
# Size factor estimation and dispersion fitting
# -------------------------

library(DESeq2)

dds <- DESeqDataSetFromMatrix(countData = assays(result)$promoterCounts,
                              colData = DataFrame(condition = condition),
                              design = ~ condition)
dds <- estimateSizeFactors(dds)
dds <- estimateDispersions(dds)
saveRDS(assays(result)$promoterCounts, file = file.path(out_dir, "raw_promoter_counts.rds"))
mu_mat <- assays(dds)[["mu"]]
saveRDS(mu_mat, file.path(out_dir, "promoter_fitted_mu.rds"))
disp_vec <- mcols(dds)$dispersion
size_vec <- 1 / disp_vec
size_vec[is.na(size_vec) | !is.finite(size_vec) | size_vec <= 0] <- NA
saveRDS(size_vec, file.path(out_dir, "NB_size_per_gene.rds"))

# -------------------------
# Call dexseq_fit2.R and pass arguments
# -------------------------


system2("Rscript", args = c(
    fit_script, 
    out_dir, 
    file.path(out_dir, "raw_promoter_counts.rds"),
    file.path(out_dir, "NB_size_per_gene.rds"),
    file.path(out_dir, "promoter_fitted_mu.rds"),
    conditions_str, 
    batch_str))



# --------------------- #
# Save outs
# --------------------- #

saveRDS(result, file = file.path(out_dir, "proactiv_result.rds"))
saveRDS(assays(result), file = file.path(out_dir, "assays_list.rds"))
saveRDS(rowData(result), file = file.path(out_dir, "rowData.rds"))
saveRDS(colData(result), file = file.path(out_dir, "colData.rds"))
saveRDS(assays(result)$normalizedPromoterCounts, file = file.path(out_dir, "normalized_promoter_counts.rds"))
saveRDS(assays(result)$absolutePromoterActivity, file = file.path(out_dir, "absolute_promoter_activity.rds"))
saveRDS(assays(result)$relativePromoterActivity, file = file.path(out_dir, "relative_promoter_activity.rds"))
saveRDS(assays(result)$geneExpression, file = file.path(out_dir, "gene_expression.rds"))

# --------------------- #
# Identify alternative promoters if reference is specified
# --------------------- #

if (reference_condition %in% condition && reference_condition != "") {
  alt_promoters <- getAlternativePromoters(result = result, referenceCondition = reference_condition)
  saveRDS(alt_promoters$upReg, file = file.path(out_dir, "alternative_promoters_upReg.rds"))
  saveRDS(alt_promoters$downReg, file = file.path(out_dir, "alternative_promoters_downReg.rds"))
  saveRDS(alt_promoters, file = file.path(out_dir, "full_alternative_promoters.rds"))
} else {
  message("No valid reference condition provided. Skipping alternative promoter analysis.")
}
