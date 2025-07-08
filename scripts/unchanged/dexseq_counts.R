#!/usr/bin/env Rscript

library(GenomicRanges)
library(rtracklayer)
library(proActiv)
library(dplyr)

# Parse arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 6) {
  stop("Usage: Rscript dexseq_process_counts.R <output_dir> <feature_gtf> <counts_file> <promoter_rds> <samples> <newnames>")
}
output_dir <- args[1]
feature_gtf <- args[2]
counts_file <- args[3]
promoter_rds <- args[4]
samples <- strsplit(args[5], " ")[[1]]
newnames <- strsplit(args[6], " ")[[1]]

# Validate inputs
if (!file.exists(feature_gtf)) stop("Feature GTF file does not exist")
if (!file.exists(counts_file)) stop("Counts file does not exist")
if (!file.exists(promoter_rds)) stop("Promoter RDS file does not exist")
if (length(samples) != length(newnames)) stop("Samples and newnames length mismatch")

# Create output directory
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Sample info
sample_info <- data.frame(samplename = samples, newnamegender = newnames, stringsAsFactors = FALSE)
rownames(sample_info) <- samples

# Read feature GTF
aggregates <- read.delim(feature_gtf, stringsAsFactors = FALSE, header = FALSE)
colnames(aggregates) <- c("chr", "source", "class", "start", "end", "ex", "strand", "ex2", "attr")
aggregates$strand <- gsub("\\.", "*", aggregates$strand)
aggregates <- aggregates[aggregates$class == "exon", ]
aggregates$attr <- gsub("\"|=|;", "", aggregates$attr)
aggregates$gene_id <- sub(".*gene_id\\s(\\S+).*", "\\1", aggregates$attr)
aggregates$gene_id <- substr(aggregates$gene_id, 1, 255)
exonids <- gsub(".*exon_number\\s(\\S+).*", "\\1", aggregates$attr)
exoninfo <- GRanges(aggregates$chr, IRanges(start = aggregates$start, end = aggregates$end), strand = aggregates$strand)
names(exoninfo) <- paste(aggregates$gene_id, exonids, sep = ":E")
transcripts <- gsub(".*transcripts\\s(\\S+).*", "\\1", aggregates$attr)
transcripts <- strsplit(transcripts, "\\+")
names(transcripts) <- names(exoninfo)

# Export exon BED and save RDS
rtracklayer::export(exoninfo, con = file.path(output_dir, "dexseq.exon.bed"))
saveRDS(exoninfo, file = file.path(output_dir, "dexseq.exoninfo.rds"))
saveRDS(transcripts, file = file.path(output_dir, "dexseq.transcripts.rds"))

# Read promoter annotation
promoterAnnotationData <- readRDS(promoter_rds)
promoter.annotation <- promoterCoordinates(promoterAnnotationData)
seqlevelsStyle(promoter.annotation) <- "Ensembl"

# Find promoter-exon overlaps
dexseq.exoninfo.promoter.overlap <- findOverlaps(promoter.annotation, exoninfo)
dexseq.exoninfo.promoter.correspond <- data.frame(
  promoterId = mcols(promoter.annotation[queryHits(dexseq.exoninfo.promoter.overlap)])$promoterId,
  geneId = mcols(promoter.annotation[queryHits(dexseq.exoninfo.promoter.overlap)])$geneId,
  internalPromoter = mcols(promoter.annotation[queryHits(dexseq.exoninfo.promoter.overlap)])$internalPromoter,
  dexseq = names(exoninfo[subjectHits(dexseq.exoninfo.promoter.overlap)]),
  stringsAsFactors = FALSE
)
saveRDS(dexseq.exoninfo.promoter.correspond, file = file.path(output_dir, "dexseq.exoninfo.promoter.correspond.rds"))
write.table(dexseq.exoninfo.promoter.correspond, file = file.path(output_dir, "dexseq.exoninfo.promoter.correspond.txt"), 
            row.names = FALSE, quote = FALSE, sep = "\t")

# Read featureCounts
counts.table <- read.delim(counts_file, comment.char = "#", stringsAsFactors = FALSE)
colnames(counts.table) <- gsub(".sorted.bam", "", gsub(".*/", "", colnames(counts.table)))
dcounts <- counts.table
id <- as.character(dcounts$Geneid)
n <- id
split(n, id) <- lapply(split(n, id), seq_along)
rownames(dcounts) <- sprintf("%s%s%03.f", id, ":E", as.numeric(n))
dcounts <- dcounts[substr(rownames(dcounts), 1, 1) != "_", ]
dcounts <- dcounts[, 7:ncol(dcounts)]

# Rename columns
dcounts <- dcounts[, sample_info$samplename]
colnames(dcounts) <- sample_info$newnamegender

# Extract promoter counts
promoterCounts.star <- readRDS(file.path(output_dir, "../proactiv/merged/promoter_counts.rds"))
dcounts.promoter <- merge(dexseq.exoninfo.promoter.correspond[, c("promoterId", "dexseq")], 
                         dcounts, by.x = "dexseq", by.y = "row.names")
rownames(dcounts.promoter) <- dcounts.promoter$promoterId
dcounts.promoter <- dcounts.promoter[, -1:-2]
dcounts.promoter <- dcounts.promoter[, colnames(promoterCounts.star)]

# Filter internal promoters
promoterCoordinates(promoterAnnotationData)[is.na(promoterCoordinates(promoterAnnotationData)$internalPromoter)]$internalPromoter <- TRUE
promoterCoordinates(promoterAnnotationData) <- promoterCoordinates(promoterAnnotationData)[!promoterCoordinates(promoterAnnotationData)$internalPromoter]
promoterCounts.star <- dcounts.promoter[as.character(promoterCoordinates(promoterAnnotationData)$promoterId), ]
saveRDS(promoterCounts, file = file.path(output_dir, paste0(sample, "_promoter_counts.rds")))
saveRDS(promoterCounts.star, file = file.path(output_dir, paste0(sample, "_promoter_counts.rds")))
