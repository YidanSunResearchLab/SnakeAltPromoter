
#!/usr/bin/env Rscript

# ---------------------
# proActiv quantification script
# Usage:
# Rscript proactiv_counts.R <out_dir> <promoter_rds> <sj_files_str> <conditions_str> <reference_condition>
# ---------------------

library(proActiv)
library(S4Vectors)
library(GenomicRanges)
library(SummarizedExperiment)
library(dplyr)
library(edgeR)
library(DESeq2)

# --------------------- #
# Parse and validate inputs
# --------------------- #

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) {
  stop("Usage: Rscript proactiv_counts.R <out_dir> <promoter_rds> <sj_files_str> <conditions_str>")
}
out_dir <- args[1]
promoter_rds <- args[2]
sj_files_str <- args[3]
conditions_str <- args[4]

# Ensure input types are character strings
if (!is.character(sj_files_str)) stop("`sj_files_str` must be a character string.")
if (!is.character(conditions_str)) stop("`conditions_str` must be a character string.")

# Parse file paths and conditions
sj_files <- strsplit(sj_files_str, " ")[[1]]
condition <- strsplit(conditions_str, ",")[[1]]


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
  stop("Number of SJ files does not match number of conditions. Consider re-running after all SJ files are completed.")
}

# Create out directory if needed
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# --------------------- #
# Run proActiv
# --------------------- #

# Load promoter annotation
promoterAnnotationData <- readRDS(promoter_rds)
prom_anno <- promoterAnnotationData

# Estimate promoter activity and summarize results across conditions
result <- proActiv(
  files = sj_files,
  promoterAnnotation = promoterAnnotationData,
  condition = condition
)

junction_counts <- assays(result)$promoterCounts
junctionReadCounts <- junction_counts
#junction_counts <- junction_counts[,order(colnames(junction_counts))]
saveRDS(junction_counts, file = file.path(out_dir, "proactiv_raw_counts.rds"))