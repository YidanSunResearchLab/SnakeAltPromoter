#!/usr/bin/env Rscript

# ---------------------------------------------
# DEXSeq promoter count extraction script
# Usage:
# Rscript dexseq_counts.R <output_dir> <feature_gtf> <counts_file> <promoter_rds> <samples> <newnames>
# ---------------------------------------------

# Load required libraries
library(GenomicRanges)
library(rtracklayer)
library(proActiv)
library(dplyr)

# Parse command-line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 6) {
  stop("Usage: Rscript dexseq_counts.R <output_dir> <feature_gtf> <counts_file> <promoter_rds> <samples> <newnames>")
}

# Assign arguments
output_dir   <- args[1]
feature_gtf  <- args[2]
counts_file  <- args[3]
promoter_rds <- args[4]
samples      <- strsplit(args[5], " ")[[1]]
newnames     <- strsplit(args[6], " ")[[1]]

# Validate inputs
if (!file.exists(feature_gtf)) stop("Feature GTF file does not exist")
if (!file.exists(counts_file)) stop("Counts file does not exist")
if (!file.exists(promoter_rds)) stop("Promoter RDS file does not exist")
if (length(samples) != length(newnames)) stop("Samples and newnames length mismatch")

# Create output directory if not exists
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Create sample info table
sample_info <- data.frame(samplename = samples, newnamegender = newnames, stringsAsFactors = FALSE)
rownames(sample_info) <- samples

# ---------------------------------------------
# Read GTF feature file (flattened exons)
# ---------------------------------------------
# Read the feature GTF file as a data frame
aggregates <- read.delim(feature_gtf, stringsAsFactors = FALSE, header = FALSE)
# Assign column names
colnames(aggregates) <- c("chr", "source", "class", "start", "end", "ex", "strand", "ex2", "attr")
# GFF categorizes '.' as unstranded, while GenomicRanges uses '*'.
# Replace '.' strand values with '*' for compatibility with GenomicRanges.
aggregates$strand <- gsub("\\.", "*", aggregates$strand)
# Keep only exon or, exonic_part features (from flattened GFF) for mapping with promoters
aggregates <- aggregates[aggregates$class %in% c("exon", "exonic_part"), ]
# Clean up attributes for easier parsing
aggregates$attr <- gsub("\"|=|;", "", aggregates$attr)
# Extract gene_id from attributes by matching the pattern gene_id space ID
aggregates$gene_id <- sub(".*gene_id\\s(\\S+).*", "\\1", aggregates$attr)
# Limit gene_id to 255 characters to ensure compatibility with downstream tools.
aggregates$gene_id <- substr(aggregates$gene_id, 1, 255)

# Extract exon_number from attributes
exonids <- gsub(".*exon_number\\s(\\S+).*", "\\1", aggregates$attr)

# Create GRanges object for exons with chromosome number, start and end position, and strand direction
exoninfo <- GRanges(aggregates$chr, IRanges(start = aggregates$start, end = aggregates$end), strand = aggregates$strand)
if (length(exoninfo) == 0) stop("exoninfo is empty. Check if your feature GTF file contains exon features.")

# Extract transcript IDs
transcripts <- gsub(".*transcripts\\s(\\S+).*", "\\1", aggregates$attr)
# Parse transcripts that are corresponding to the same exon and are connected by + to a list
transcripts <- strsplit(transcripts, "\\+")


# ---------------------------------------------
# Load promoter annotation and overlap
# ---------------------------------------------
promoterAnnotationData <- readRDS(promoter_rds)
promoter.annotation <- promoterCoordinates(promoterAnnotationData)

# Set consistent seqlevels style
seqlevelsStyle(promoter.annotation) <- "UCSC"
seqlevelsStyle(exoninfo) <- "UCSC"

# Identify overlaps between promoters and exons
dexseq.exoninfo.promoter.overlap <- findOverlaps(promoter.annotation, exoninfo)


# ---------------------------------------------
# Load featureCounts file
# ---------------------------------------------
# Read the featureCounts counts file into a data frame. Skip header comments.
counts.table <- read.delim(counts_file, comment.char = "#", stringsAsFactors = FALSE)
# Rename the first six columns to standardized names used for genomic features.
colnames(counts.table)[1:6] <- c("Geneid", "Chr", "Start", "End", "Strand", "Length")
# Remove sufixes
colnames(counts.table) <- gsub(".sorted.bam", "", colnames(counts.table))
colnames(counts.table) <- gsub("^X\\.*", "", colnames(counts.table))
colnames(counts.table) <- gsub(".*bam\\.", "", colnames(counts.table))

gene_ids <- as.character(counts.table$Geneid)
# Number exons per gene using run-length encoding to construct unique exon IDs.
# Exons from the same gene ordered consecutively from the featureCounts output.
exon_numbers <- sequence(rle(gene_ids)$lengths)
# Assign the exon names to exon info.
names(exoninfo) <- sprintf("%s:E%03d", gene_ids , exon_numbers)
# Transcripts and exoninfo are both in the order of aggregate. Construct promoter-exon-transcript mapping.
names(transcripts) <- names(exoninfo)
# Save exon BED and metadata
saveRDS(transcripts, file = file.path(output_dir, "dexseq.transcripts.rds"))
rtracklayer::export(exoninfo, con = file.path(output_dir, "dexseq.exon.bed"))
saveRDS(exoninfo, file = file.path(output_dir, "dexseq.exoninfo.rds"))

# Extract counts matrix
# Read from the 7th column onwards, which contains the counts for each sample.
dcounts <- counts.table[, 7:ncol(counts.table), drop = FALSE]
# Unique exon IDs for every row
rownames(dcounts) <- sprintf("%s:E%03d", gene_ids, exon_numbers)
# Remove rows that are irrelevant for analysis
dcounts <- dcounts[!grepl("^_", rownames(dcounts)), , drop = FALSE]
# Keep only sample names for columns
colnames(dcounts) <- gsub("\\..*$", "", colnames(dcounts))
colnames(dcounts) <- sample_info$newnamegender

# ---------------------------------------------
# Match promoter to exon (DEXSeq)
# ---------------------------------------------
# dexseq.exoninfo.promoter.overlap is a Hits object returned by findOverlaps()
# Each overlap result links a promoter to a DEXSeq exon.

# Indices of the overlapping promoters in promoter.annotation.
qhits <- queryHits(dexseq.exoninfo.promoter.overlap)
# Indices of the overlapping exons in exoninfo.
shits <- subjectHits(dexseq.exoninfo.promoter.overlap)

# Get metadata for overlapping promoters
promoter_mcols <- as.data.frame(mcols(promoter.annotation))[qhits, , drop = FALSE]
# Get corresponding exon IDs
geneId_col <- if ("geneId" %in% colnames(promoter_mcols)) promoter_mcols$geneId else rep(NA, nrow(promoter_mcols))

# Construct mapping table
dexseq.exoninfo.promoter.correspond <- data.frame(
  promoterId       = promoter_mcols$promoterId,
  geneId           = geneId_col,
  internalPromoter = promoter_mcols$internalPromoter,
  dexseq           = names(exoninfo)[shits],
  stringsAsFactors = FALSE
)
# Keep geneId for later merging
dexseq.exoninfo.promoter.correspond$geneId <- sub(":E\\d+$", "", dexseq.exoninfo.promoter.correspond$dexseq)
saveRDS(dexseq.exoninfo.promoter.correspond, file = file.path(output_dir, "dexseq.exoninfo.promoter.correspond.rds"))
write.table(dexseq.exoninfo.promoter.correspond,
            file = file.path(output_dir, "dexseq.exoninfo.promoter.correspond.txt"),
            row.names = FALSE, quote = FALSE, sep = "\t")

# ---------------------------------------------
# Merge counts by promoter
# ---------------------------------------------
# Merge promoter exon mapping with DEXSeq counts by exon Id.
dcounts.promoter <- merge(
  dexseq.exoninfo.promoter.correspond[, c("promoterId", "dexseq")],
  dcounts,
  by.x = "dexseq", by.y = "row.names"
)

cat("Before aggregation: NAs in promoterId = ", sum(is.na(dcounts.promoter$promoterId)), "\n")
# Aggregate exon counts into promoter counts by summing counts exon bins from the same promoter.
dcounts.promoter.merged <- aggregate(. ~ promoterId,
                                     data = dcounts.promoter[, c("promoterId", colnames(dcounts))],
                                     FUN = sum)

# Set row names to promoterId and remove promoterId column so that the matrix has only counts.
rownames(dcounts.promoter.merged) <- dcounts.promoter.merged$promoterId
dcounts.promoter.merged <- dcounts.promoter.merged[, -1, drop = FALSE]
na_like <- grep("^NA", rownames(dcounts.promoter.merged), value = TRUE)
cat("NA-like promoter rownames after aggregation:", length(na_like), "\n")
print(head(na_like))


# ---------------------------------------------
# Filter and save per-sample output
# ---------------------------------------------
# Flag NAs in internalPromoter column as TRUE to ensure only external promoters are kept.
na_mask <- is.na(promoterCoordinates(promoterAnnotationData)$internalPromoter)
if (any(na_mask)) {
#  promoterCoordinates(promoterAnnotationData)$internalPromoter[na_mask] <- TRUE
  # Set internalPromoter NA to FALSE to keep promoters of transcript with only one exon
  promoterCoordinates(promoterAnnotationData)$internalPromoter[na_mask] <- FALSE
}

# Filter to keep only external promoters in promoterAnnotationData
promoterCoordinates(promoterAnnotationData) <- promoterCoordinates(promoterAnnotationData)[
  !promoterCoordinates(promoterAnnotationData)$internalPromoter
]

# Extract promoter counts from dcounts in the order of promoterAnnotationData.
# This ensures that the counts are aligned with the promoter IDs and are only consisting of external promoters.

requested_ids <- as.character(promoterCoordinates(promoterAnnotationData)$promoterId)

available_ids <- rownames(dcounts.promoter.merged)

missing_ids <- setdiff(requested_ids, available_ids)
cat("Number of promoterIds not found in dcounts.promoter.merged:", length(missing_ids), "\n")
print(head(missing_ids))

fill_matrix <- matrix(0,
                      nrow = length(missing_ids),
                      ncol = ncol(dcounts.promoter.merged),
                      dimnames = list(missing_ids, colnames(dcounts.promoter.merged)))

full_counts <- rbind(dcounts.promoter.merged, fill_matrix)
full_counts <- full_counts[requested_ids, , drop = FALSE]

promoterCounts.star <- full_counts

cat("==== Sanity check ====\n")
cat("Total promoters in annotation:", length(requested_ids), "\n")
cat("Rows in final promoterCounts.star:", nrow(promoterCounts.star), "\n")
if (all(rownames(promoterCounts.star) == requested_ids)) {
  cat("promoterCounts.star rownames exactly match promoterAnnotation order\n")
} else {
  cat("promoterCounts.star rownames do NOT match promoterAnnotation order\n")
  mismatch_idx <- which(rownames(promoterCounts.star) != requested_ids)
  cat("First mismatches at:\n")
  print(head(data.frame(row = mismatch_idx,
                        expected = requested_ids[mismatch_idx],
                        actual = rownames(promoterCounts.star)[mismatch_idx])))
}

zero_promoters <- rowSums(promoterCounts.star) == 0
cat("Number of promoters with total count = 0:", sum(zero_promoters), "\n")


# Save one RDS file per sample
for (i in seq_along(newnames)) {
  saveRDS(promoterCounts.star[, i, drop = FALSE],
          file = file.path(output_dir, newnames[i], paste0(newnames[i], "_promoter_counts.rds")))
}
