#!/usr/bin/env Rscript

library(tximport)
library(proActiv)
library(dplyr)

# Parse arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5) {
  stop("Usage: Rscript salmon_promoter_counts.R <output_dir> <promoter_rds> <quant_file> <sample> <newname>")
}
output_dir <- args[1]
promoter_rds <- args[2]
quant_file <- args[3]
sample <- args[4]
newname <- args[5]

# Validate inputs
if (!file.exists(promoter_rds)) stop("Promoter RDS file does not exist")
if (!file.exists(quant_file)) stop("Quant file does not exist")

# Create output directory
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Read promoter annotations
promoterAnnotationData <- readRDS(promoter_rds)
tx2promoter <- data.frame(
  TXNAME = promoterAnnotationData$transcriptId,
  PROMOTERID = promoterAnnotationData$promoterId,
  stringsAsFactors = FALSE
)

# Import Salmon quantification
txi <- tximport(
  files = quant_file,
  type = "salmon",
  txOut = TRUE,
  ignoreTxVersion = TRUE
)

# Aggregate to promoter level
promoter_counts <- txi$counts
rownames(promoter_counts) <- sub("\\..*", "", rownames(promoter_counts))  # Remove version numbers
promoter_counts_df <- data.frame(
  TXNAME = rownames(promoter_counts),
  counts = promoter_counts[, 1],
  stringsAsFactors = FALSE
)
colnames(promoter_counts_df)[2] <- newname

# Merge with promoter annotations
promoter_counts_merged <- merge(tx2promoter, promoter_counts_df, by = "TXNAME", all.x = FALSE)
promoter_counts_agg <- promoter_counts_merged %>%
  group_by(PROMOTERID) %>%
  summarise(across(all_of(newname), sum, na.rm = TRUE)) %>%
  as.data.frame()

# Convert to matrix
promoter_counts_mat <- as.matrix(promoter_counts_agg[, -1, drop = FALSE])
rownames(promoter_counts_mat) <- promoter_counts_agg$PROMOTERID

# Filter to match proActiv promoters
promoter_ids <- promoterAnnotationData$promoterId
promoter_counts_mat <- promoter_counts_mat[rownames(promoter_counts_mat) %in% promoter_ids, , drop = FALSE]

# Save promoter counts
saveRDS(promoter_counts_mat, file = file.path(output_dir, paste0(sample, "_promoter_counts.rds")))