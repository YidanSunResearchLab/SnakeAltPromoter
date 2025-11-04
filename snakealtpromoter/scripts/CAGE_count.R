#!/usr/bin/env Rscript

# -------------------------
# CAGE Promoter Merge Script
# -------------------------
# Usage:
# Rscript CAGE_merge.R <out_dir> <promoter_rds> <norm_method> <samples> <conditions> <newnames> <comparison> <count_files...>
# -------------------------

library(proActiv)  # v1.16
library(DESeq2)
library(edgeR)
library(dplyr)

# -------------------------
# Parse arguments
# -------------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5) {
  stop("Usage: Rscript dexseq_promoter_merge.R <out_dir> <promoter_rds> <samples> <conditions> <fc_dir> <cell_line> <fit_script> <batch>")
}

out_dir     <- args[1]
prom_rds    <- args[2]
samples     <- strsplit(args[3], " ")[[1]]
conditions_str <- args[4]
fc_dir <- args[5]


conditions  <- strsplit(conditions_str, ",")[[1]]

# -------------------------
# Validate
# -------------------------

cat("fc_dir passed in:", fc_dir, "\n")
cat("Samples:", paste(samples, collapse = ", "), "\n")
cat("Conditions:", paste(conditions, collapse = ", "), "\n")


stopifnot(length(samples) == length(conditions),
#        dir.exists(fc_dir),
        file.exists(prom_rds))

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# -------------------------
# Load promoter annotation and count data
# -------------------------
prom_anno <- readRDS(prom_rds)
pc <- proActiv:::promoterCoordinates(prom_anno)
all_fc_files <- list.files(fc_dir, pattern = "_promoter_counts\\.txt$", full.names = TRUE)

# Get all promoter-level output files
message("Getting promoter-level output files...")
fc_list <- setNames(lapply(all_fc_files, function(f) {
  df <- read.table(f, header = TRUE, sep = "\t",
                   comment.char = "#", stringsAsFactors = FALSE,
                   check.names = FALSE)
  sample_name <- sub("_promoter_counts\\.txt$", "", basename(f))
  prom_id <- sub("_$", "", df[[1]])  # Extract promoter ID
  counts <- df[[ncol(df)]]  # Extract count values
  names(counts) <- prom_id
  return(counts)
}), nm = sub("_promoter_counts\\.txt$", "", basename(all_fc_files)))

message("Number of promoter‑count files: ", length(all_fc_files))
message("Files:\n", paste(all_fc_files, collapse = "\n"))
# Build promoter-by-sample count matrix
message("Building promoter-by-sample count matrix...")
all_ids <- unique(unlist(lapply(fc_list, names)))
mat <- sapply(fc_list, function(x) x[all_ids])
rownames(mat) <- all_ids
colnames(mat) <- names(fc_list)

# Convert to numeric matrix, fill NAs with 0
mat <- as.matrix(mat)
mode(mat) <- "numeric"
mat[is.na(mat)] <- 0
cage_counts <- mat
saveRDS(cage_counts, file.path(out_dir, "cage_counts.rds"))