#!/usr/bin/env Rscript

# Usage:
# Rscript merge_promoter_counts.R <counts_dir> <output_file>

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 6) {
  stop("Usage: Rscript merge_promoter_counts.R <counts_dir> <output_file>")
}

out_dir <- args[1]
#counts_files  <- args[2] 
prom_rds    <- args[2]
samples     <- strsplit(args[3], " ")[[1]]
conditions_str <- args[4]
comparison <- strsplit(args[5], " ")[[1]]
count_files <- args[6:length(args)]

condition  <- strsplit(conditions_str, ",")[[1]]

# -------------------------
# Sample information
# -------------------------
sample_info <- data.frame(
  sample = samples,
  condition = condition,
  newname = samples,
  file = count_files,
  stringsAsFactors = FALSE
)

# Keep conditions for comparison
sample_info <- sample_info[sample_info$condition %in% comparison, ]
print(sample_info$condition)
print(comparison)
stopifnot(nrow(sample_info) > 0)
parent_dir <- dirname(dirname(out_dir))
bam_dir <- file.path(parent_dir, "bam")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# -------------------------
# Load promoter annotation and count data
# -------------------------
promAnno <- readRDS(prom_rds)

# Read each RDS file in count_files into a list of matrices
cnt_list <- lapply(sample_info$file, function(f) {
  m <- readRDS(f)
  m[is.na(m)] <- 0 
  m
})

for (i in seq_along(cnt_list)) {
  rn <- rownames(cnt_list[[i]])
  missing_idx <- which(is.na(rn) | rn == "")
  if (length(missing_idx) > 0) {
    cat("File with invalid rownames (promoterId):", sample_info$file[i], "\n")
    cat("Number of invalid rownames:", length(missing_idx), "\n")
    cat("Example rownames:", head(rn[missing_idx]), "\n")
  }
}


# Extract all unique promoter IDs across all count matrices
prom_ids <- unique(unlist(lapply(cnt_list, rownames)))

# Construct an empty matrix to hold promoter counts
# Rows: all unique promoter IDs
message("Combining promoter counts...")
# Columns: one for each sample, using names from newnames
cnt_mat <- matrix(0, nrow = length(prom_ids), ncol = nrow(sample_info),
                  dimnames = list(prom_ids, sample_info$newname))
# Populate the matrix with values from each sample's count object
# rownames(cnt_list[[i]]) gives the promoter IDs present in that sample.
# cnt_list[[i]][, 1] extracts the count values (first column).
# These are assigned into the corresponding rows and column in cnt_mat.
dim(cnt_mat)
for (i in seq_along(cnt_list)) {
  counts <- cnt_list[[i]][, 1]
  counts[is.na(counts)] <- 0
  cnt_mat[rownames(cnt_list[[i]]), sample_info$newname[i]] <- counts
}
cnt_mat[is.na(cnt_mat)] <- 0
dim(cnt_mat)
na_like_rows <- grep("^NA", rownames(cnt_mat), value = TRUE)
cat("Number of NA-like promoterIds:", length(na_like_rows), "\n")
cat("Example NA-like promoterIds:", head(na_like_rows), "\n")
for (i in seq_along(cnt_list)) {
  bad_rows <- grep("^NA", rownames(cnt_list[[i]]), value = TRUE)
  if (length(bad_rows) > 0) {
    cat("Sample with bad promoterId:", sample_info$file[i], "\n")
    cat("PromoterId(s):", paste(head(bad_rows), collapse = ", "), "\n\n")
  }
}

# Save combined matrix
saveRDS(cnt_mat, file = file.path(out_dir, "merged_promoter_counts.rds"))

