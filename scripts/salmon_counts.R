
#!/usr/bin/env Rscript

# ---------------------
# Salmon counts script
# Usage:
# Rscript salmon_promoter_counts.R <output_dir> <promoter_rds> <quant_file> <sample>
# ---------------------


library(tximport)
library(proActiv)
library(dplyr)


# Parse arguments
args <- commandArgs(trailingOnly = TRUE)

cat("Received arguments:\n")
print(args)


if (length(args) != 4) {
  stop("Usage: Rscript salmon_promoter_counts.R <output_dir> <promoter_rds> <quant_file> <sample>")
}
output_dir <- args[1]
promoter_rds <- args[2]
quant_file <- args[3]
sample <- args[4]

# Validate inputs
if (!file.exists(promoter_rds)) stop("Promoter RDS file does not exist")
if (!file.exists(quant_file)) stop("Quant file does not exist")

# Create output directory
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Read promoter annotations
promoterAnnotationData <- readRDS(promoter_rds)
# Use @ to access slot in S4 class. Convert to data frame
tx2promoter <- as.data.frame(promoterAnnotationData@promoterIdMapping)
# Extract transcript name to promoter id mapping
tx2promoter <- tx2promoter[, c("transcriptName", "promoterId")]


# ---- Detect transcripts with multiple versions in annotation ----
tx2promoter_versions_checked <- as.data.frame(promoterAnnotationData@promoterIdMapping)

# Extract transcript base and version
tx2promoter_versions_checked <- tx2promoter_versions_checked %>%
  mutate(
    transcript_base = sub("\\..*", "", transcriptName),
    version = suppressWarnings(as.integer(sub(".*\\.", "", transcriptName)))
  )

# Only keep the latest version of each transcript base
tx2promoter_versions_checked <- tx2promoter_versions_checked %>%
  group_by(transcript_base) %>%
  filter(version == max(version, na.rm = TRUE)) %>%
  ungroup()
tx2promoter <- tx2promoter_versions_checked[, c("transcriptName", "promoterId")]

# Number of versions per base
version_count_per_base <- tx2promoter_versions_checked %>%
  group_by(transcript_base) %>%
  summarise(n_versions = n_distinct(version), .groups = "drop") %>%
  filter(n_versions > 1)

# Warning if multiple versions of the same transcript is present
if (nrow(version_count_per_base) > 0) {
  warning("The annotation contains multiple versions of the same transcript base.\n",
          "Deduplicated by using the newest version.\n",
          "Affected transcript bases (up to 5 shown): ",
          paste(head(version_count_per_base$transcript_base, 5), collapse = ", "), " ...")
}

#------------------------------------------


# Rename column names and remove version numbers for transcript
colnames(tx2promoter) <- c("TXNAME", "PROMOTERID")
tx2promoter$TXNAME <- sub("\\..*", "", tx2promoter$TXNAME)

# Import Salmon quantification
# Return transcript-level counts and ignores version number in transcripts
txi <- tximport(
  files = quant_file,
  type = "salmon",
  txOut = TRUE,
  ignoreTxVersion = TRUE
)

# Aggregate to promoter level
# Extract transcript-level count from txi
# Remove version numbers in rownames to match annotation
# Convert to dataframe
promoter_counts <- txi$counts
rownames(promoter_counts) <- sub("\\..*", "", rownames(promoter_counts))  # Remove version numbers
promoter_counts_df <- data.frame(
  TXNAME = rownames(promoter_counts),
  counts = promoter_counts[, 1],
  stringsAsFactors = FALSE
)
# Change column name from counts to newname
colnames(promoter_counts_df)[2] <- sample

cat("Number of transcripts in tx2promoter:", nrow(tx2promoter), "\n")
cat("Number of transcripts in promoter_counts_df:", nrow(promoter_counts_df), "\n")
# Keep only transcripts existing in both quant and annotation
coord <- promoterCoordinates(promoterAnnotationData)
# Filter out internal promoters
# If internalPromoter is NA, set it to TRUE
coord$internalPromoter[is.na(coord$internalPromoter)] <- TRUE
# Extract promoter IDs that are not internal promoters
promoter_ids <- coord$promoterId[!coord$internalPromoter]
# Ensure promoter_ids does not contain duplicates
promoter_ids <- unique(promoter_ids)
# Only keep annotated promoter_ids
tx2promoter <- tx2promoter[tx2promoter$PROMOTERID %in% promoter_ids, , drop = FALSE]
common_tx <- intersect(tx2promoter$TXNAME, promoter_counts_df$TXNAME)
tx2promoter_filtered <- tx2promoter[tx2promoter$TXNAME %in% common_tx, , drop = FALSE]
promoter_counts_df_filtered <- promoter_counts_df[promoter_counts_df$TXNAME %in% common_tx, , drop = FALSE]
cat("Number of overlapping transcripts:", length(common_tx), "\n")
cat("Number of transcripts in promoter_counts_df_filtered:", nrow(promoter_counts_df_filtered), "\n")
# Merge
promoter_counts_merged <- merge(
  tx2promoter_filtered,
  promoter_counts_df_filtered,
  by = "TXNAME",
  all.x = FALSE,
  all.y = FALSE
)
cat("Number of rows after merge: ", nrow(promoter_counts_merged), "\n")
# Print to check for row numbers, promoter IDs, and overlaps
cat("After merge: ", nrow(promoter_counts_merged), "rows\n")
message("Unique promoter IDs    : ", n_distinct(promoter_counts_merged$PROMOTERID))
tx_only <- setdiff(tx2promoter$TXNAME, promoter_counts_df$TXNAME)
count_only <- setdiff(promoter_counts_df$TXNAME, tx2promoter$TXNAME)
cat("Transcripts in annotation but not in quant:", length(tx_only), "\n")
cat("Transcripts in quant but not in annotation:", length(count_only), "\n")
cat("Transcripts in quant but not in annotation (some examples):\n")
print(head(count_only, 20))
if(length(tx_only) > 0){
  cat("Sample transcripts in annotation missing in quant:", head(tx_only, 10), "\n")
}
if(length(count_only) > 0){
  cat("Sample transcripts in quant missing in annotation:", head(count_only, 10), "\n")
}

# Sum the expression levels across all transcripts for every promoter (per row), dropping NA values
promoter_counts_agg <- promoter_counts_merged %>%
  group_by(PROMOTERID) %>%
  summarise(across(all_of(sample), sum, na.rm = TRUE)) %>%
  as.data.frame()

# Convert to matrix, dropping PROMOTERID column and use it as rownames
promoter_counts_mat <- as.matrix(promoter_counts_agg[, -1, drop = FALSE])
rownames(promoter_counts_mat) <- promoter_counts_agg$PROMOTERID

# Filter promoters so the ones kept are annotated by proActiv
promoterIdMapping <- promoterAnnotationData@promoterIdMapping
promoter_ids <- unique(promoterIdMapping$promoterId)
promoter_counts_mat <- promoter_counts_mat[rownames(promoter_counts_mat) %in% promoter_ids, , drop = FALSE]
promoter_counts_mat <- promoter_counts_mat[!is.na(rownames(promoter_counts_mat)), , drop = FALSE]
if (any(duplicated(rownames(promoter_counts_mat)))) {
  warning("Duplicated promoterId found. Aggregating counts.")
  promoter_counts_mat <- rowsum(promoter_counts_mat, group = rownames(promoter_counts_mat))
}
cat("Final promoterId rows: ", nrow(promoter_counts_mat), "\n")
cat("Expected from annotation: ", length(unique(promoterAnnotationData@promoterIdMapping$promoterId)), "\n")
cat("Number of transcripts in quant.sf: ", nrow(promoter_counts_df), "\n")
cat("Number of transcripts in tx2promoter: ", nrow(tx2promoter), "\n")
cat("Number of overlapping TXNAMEs: ", length(intersect(promoter_counts_df$TXNAME, tx2promoter$TXNAME)), "\n")
cat("Unique promoters with counts: ",
    dplyr::n_distinct(promoter_counts_merged$PROMOTERID), "\n")
cat("Rows after aggregation      : ",
    nrow(promoter_counts_agg), "\n")
# Save promoter counts
saveRDS(promoter_counts_mat, file = file.path(output_dir, paste0(sample, "_promoter_counts.rds")))

# Print warnings -----------------------------------
if (length(warnings()) > 0) {
  cat("\n=== WARNINGS ===\n")
  print(warnings())
}