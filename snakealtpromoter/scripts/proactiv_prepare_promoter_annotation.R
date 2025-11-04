#!/usr/bin/env Rscript

# Load required libraries
library(GenomicFeatures)
library(proActiv)
library(rtracklayer)

# Parse command-line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5) {
  stop("Usage: Rscript prepare_promoter_annotation.R <input_gtf> <output_dir> <species> <promoter_rds> <promoter_bed>")
}
input_gtf <- args[1]
output_dir <- args[2]
species <- args[3]
promoter_rds <- args[4]
promoter_bed <- args[5]

# Create output directories
anno_dir <- file.path(output_dir, "Annotation")
dir.create(anno_dir, showWarnings = FALSE, recursive = TRUE)

# Filter GTF to remove retained_intron and nonsense_mediated_decay
filtered_gtf <- file.path(anno_dir, "proActiv_protein_coding.gtf")
system(paste("grep -v retained_intron", input_gtf, "| grep -v nonsense_mediated_decay >", filtered_gtf))

# Create and save TxDb
txdb <- makeTxDbFromGFF(filtered_gtf)
txdb_file <- file.path(anno_dir, "proActiv_txdb.sqlite")
saveDb(txdb, txdb_file)

# Load TxDb and prepare promoter annotation
txdb <- loadDb(txdb_file)
promoterAnnotationData <- preparePromoterAnnotation(txdb = txdb, species = species)

# Save promoter annotation
saveRDS(promoterAnnotationData, file = promoter_rds)

# Export promoter coordinates as BED
promoter_anno <- promoterAnnotationData@promoterCoordinates
promoter_anno$name <- paste(promoter_anno$promoterId, promoter_anno$geneId, sep = "_")
export(promoter_anno, con = promoter_bed, format = "BED")

# Clean up
gc()