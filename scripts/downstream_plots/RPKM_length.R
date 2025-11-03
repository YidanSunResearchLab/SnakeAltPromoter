library(ggplot2)
library(MASS)
library(viridis)
library(GenomicRanges)
library(S4Vectors)
library(data.table)
library(rtracklayer)


# ---------------------------
# Config
# ---------------------------
cell_lines <- c("Healthy", "Failure")
out_dir    <- "/mnt/citadel2/research/syidan/Projects/SnakeAltPromoterResult/"
prom_anno_rds <- "/mnt/citadel2/research/shared/AltPromoterFlow/genome/organisms/hg38/Annotation/proActiv_promoter_annotation.rds"
gtf_path      <- "/mnt/citadel2/research/shared/SnakeAltPromoter_processed/genome/organisms/hg38/Annotation.gtf"

# ---------------------------
# Helpers
# ---------------------------

normalize_tx_id <- function(x) sub("\\.\\d+$", "", as.character(x))


# ---------------------------
# Annotation: keep NON-internal promoters
# ---------------------------
promoter_anno <- readRDS(prom_anno_rds)
coord <- promoter_anno@promoterCoordinates
meta  <- mcols(coord)

non_internal_ok <- !meta$internalPromoter %in% TRUE #!(isTRUE(meta$internalPromoter))
valid_ids <- as.character(meta$promoterId[non_internal_ok])
coord <- coord[coord$promoterId %in% valid_ids]
meta  <- mcols(coord)

# ---------------------------
# DEXSeq effective length per promoter (firstExon* if available, else TSS/end fallback)
# ---------------------------
pid <- as.character(meta$promoterId)
tss <- ifelse(as.character(strand(coord)) == "+", start(coord), end(coord))
# estimates first exon length per promoter by subtracting the promoter’s TSS coordinate from the end of its first exon
dex_len_bp <- setNames(abs(as.numeric(meta$firstExonEnd) - as.numeric(tss)) + 1, pid)
message("[DEXSeq] Estimated first-exon length from firstExonEnd + strand-aware TSS.")

# ---------------------------
# prom2tx from promoterIdMapping
# ---------------------------
pam_df <- as.data.frame(promoter_anno@promoterIdMapping, stringsAsFactors = FALSE)
stopifnot(all(c("transcriptName","promoterId") %in% colnames(pam_df)))
pam_df$promoterId     <- as.character(pam_df$promoterId)
# Normalize transcript IDs (remove version numbers)
pam_df$transcriptName <- normalize_tx_id(pam_df$transcriptName)
pam_df <- pam_df[pam_df$promoterId %in% valid_ids, , drop = FALSE]
# split into list, where each promoter ID maps to a vector of unique transcript IDs
split_tx <- split(pam_df$transcriptName, pam_df$promoterId)
prom2tx  <- lapply(split_tx, function(v) unique(v[!is.na(v) & nzchar(v)]))
prom2tx <- prom2tx[valid_ids]


# ---------------------------
# transcript exon lengths from GTF → tx_len_bp (Bioc reduce)
# ---------------------------

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
gtf_md5  <- tryCatch(as.character(tools::md5sum(gtf_path)), error = function(e) NA_character_)
gtf_time <- tryCatch(as.numeric(file.info(gtf_path)$mtime), error = function(e) NA_real_)
txlen_cache <- file.path(
  out_dir,
  sprintf("tx_len_bp_cache_%s_%s.rds",
          ifelse(is.na(gtf_md5), "nohash", substr(gtf_md5,1,12)),
          ifelse(is.na(gtf_time), "notime", as.integer(gtf_time)))
)
# Load from cache if available
if (file.exists(txlen_cache)) {
  message("[tx_len_bp] Loading cache: ", txlen_cache)
  tx_len_bp <- readRDS(txlen_cache)

} else {
  message("[tx_len_bp] Importing GTF with rtracklayer + reduce ...")
  gr <- rtracklayer::import(gtf_path)  # GRanges

  mc <- mcols(gr)
  stopifnot("type" %in% colnames(mc))
  tx_col <- ifelse("transcript_id" %in% colnames(mc), "transcript_id",
                 ifelse("transcriptId" %in% colnames(mc), "transcriptId",
                        stop("No transcript_id/transcriptId column in imported GTF.")))
  # Extract exons
  exon_gr <- gr[mcols(gr)$type == "exon"]
  tx_id   <- mcols(exon_gr)[[tx_col]]
  # Define helper function to remove transcript version numbers
  normalize_tx_id <- function(x) sub("\\.\\d+$", "", as.character(x))
  tx_id <- normalize_tx_id(tx_id)
  # Keep only exons with valid transcript IDs
  keep <- !is.na(tx_id) & nzchar(tx_id)
  exon_gr <- exon_gr[keep]
  tx_id   <- tx_id[keep]
  # Attach as a new column to exon_gr
  mcols(exon_gr)$tx <- tx_id
  # Split per transcript → reduce (merge overlapping/touching exons) → width → sum
  ex_list   <- split(exon_gr, mcols(exon_gr)$tx, drop = TRUE)
  ex_red    <- reduce(ex_list, ignore.strand = FALSE)
  tx_len_bp <- sum(width(ex_red))
  tx_len_bp <- setNames(as.numeric(tx_len_bp), names(ex_red))
  # Guards
  med <- median(tx_len_bp[is.finite(tx_len_bp) & tx_len_bp > 0], na.rm = TRUE)
  tx_len_bp[!is.finite(tx_len_bp) | tx_len_bp <= 0] <- med
  saveRDS(tx_len_bp, txlen_cache)
  message("[tx_len_bp] Saved cache: ", txlen_cache)
}


# ---------------------------
# salmon_len_bp from prom2tx + tx_len_bp
# ---------------------------
# repeats each promoter ID for number of transcripts it maps
pid_vec <- rep(names(prom2tx), lengths(prom2tx))
# flattens the list of transcript IDs into a single vector
# match order with pid_vec, one transcript per promoter–transcript pair
tx_vec  <- unlist(prom2tx, use.names = FALSE)
# table with promoter ID, transcript ID, transcript length
dt <- data.table(pid = pid_vec, tx = tx_vec)
dt[, len := tx_len_bp[tx]]
dt <- dt[is.finite(len)]
# group by promoter id, computes mean transcript length per promoter
res <- dt[, .(mean_len = mean(len)), by = pid]
salmon_len_bp <- setNames(rep(NA_real_, length(prom2tx)), names(prom2tx))
salmon_len_bp[res$pid] <- res$mean_len
med <- median(salmon_len_bp[is.finite(salmon_len_bp) & salmon_len_bp > 0], na.rm = TRUE)
salmon_len_bp[!is.finite(salmon_len_bp) | salmon_len_bp <= 0] <- med

save(salmon_len_bp, dex_len_bp, file = file.path(out_dir, "lengths_promoters.RData"))








# # ---------------------------
# # Convert Salmon/DEXSeq to RPKM
# # ---------------------------
# salmon_rpkm <- to_rpkm_by_len(salmon_counts, salmon_len_bp)
# dexseq_rpkm <- to_rpkm_by_len(dexseq_counts, dex_len_bp)

# # ---------------------------
# # Make SIX plots (method × condition)
# # ---------------------------
# for (cl in cell_lines) {
#   # proActiv (counts) vs CAGE (counts)
#   df1 <- make_df_one(cage_counts, proactiv_counts, cl, "proActiv (counts)")
#   if (!is.null(df1)) plot_one(df1, sprintf("scatter_proactiv_%s.pdf", tolower(cl)), y_is_rpkm = FALSE)

#   # Salmon (RPKM) vs CAGE (counts)
#   df2 <- make_df_one(cage_counts, salmon_rpkm, cl, "Salmon (RPKM)")
#   if (!is.null(df2)) plot_one(df2, sprintf("scatter_salmon_rpkm_%s.pdf", tolower(cl)), y_is_rpkm = TRUE)

#   # DEXSeq (RPKM) vs CAGE (counts)
#   df3 <- make_df_one(cage_counts, dexseq_rpkm, cl, "DEXSeq (RPKM)")
#   if (!is.null(df3)) plot_one(df3, sprintf("scatter_dexseq_rpkm_%s.pdf", tolower(cl)), y_is_rpkm = TRUE)
# }

# message("\nDone. 6 plots written to: ", out_dir)
