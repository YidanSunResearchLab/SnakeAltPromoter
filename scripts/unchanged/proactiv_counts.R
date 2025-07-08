#!/usr/bin/env Rscript

library(proActiv)

# Parse arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) {
  stop("Usage: Rscript proactiv_counts.R <output_dir> <promoter_rds> <sj_file> <sample>")
}
output_dir <- args[1]
promoter_rds <- args[2]
sj_file <- args[3]
sample <- args[4]

# Validate inputs
if (!file.exists(promoter_rds)) stop("Promoter RDS file does not exist")
if (!file.exists(sj_file)) stop("SJ file does not exist")

# Create output directory
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Read promoter annotation
promoterAnnotationData <- readRDS(promoter_rds)
promoterCoordinates(promoterAnnotationData)[is.na(promoterCoordinates(promoterAnnotationData)$internalPromoter)]$internalPromoter <- TRUE
promoterCoordinates(promoterAnnotationData) <- promoterCoordinates(promoterAnnotationData)[!promoterCoordinates(promoterAnnotationData)$internalPromoter]

# Count promoter reads
promoterCounts <- calculatePromoterReadCounts(
  promoterAnnotationData,
  junctionFilePaths = sj_file,
  junctionFileLabels = sample,
  junctionType = "star",
  numberOfCores = parallel::detectCores() - 1
)
promoterCounts <- promoterCounts[as.character(promoterCoordinates(promoterAnnotationData)$promoterId), , drop = FALSE]

# Save counts
saveRDS(promoterCounts, file = file.path(output_dir, paste0(sample, "_promoter_counts.rds")))