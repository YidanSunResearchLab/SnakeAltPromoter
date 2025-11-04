#!/usr/bin/env Rscript

# -------------------------
# CAGE Promoter Merge Script
# -------------------------
# Usage:
# Rscript CAGE_merge.R <out_dir> <promoter_rds> <norm_method> <samples> <condition> <newnames> <comparison> <count_files...>
# -------------------------

library(proActiv)  # v1.16
library(DESeq2)
library(edgeR)
library(dplyr)

# -------------------------
# Parse arguments
# -------------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 9) {
  stop("Usage: Rscript dexseq_promoter_merge.R <out_dir> <promoter_rds> <samples> <condition> <fc_dir> <cell_line> <fit_script> <batch>")
}

out_dir     <- args[1]
prom_rds    <- args[2]
samples     <- strsplit(args[3], " ")[[1]]
conditions_str <- args[4]
counts <- args[5]
comparison <- strsplit(args[6], " ")[[1]]
print(comparison)
fit_script <- args[7]  # Path to the dexseq_fit script
batch_str <- args[8]
norm_method <- args[9]

condition  <- strsplit(conditions_str, ",")[[1]]

# -------------------------
# Validate
# -------------------------

cat("Samples:", paste(samples, collapse = ", "), "\n")
cat("Conditions:", paste(condition, collapse = ", "), "\n")


stopifnot(length(samples) == length(condition),
#        dir.exists(fc_dir),
        file.exists(prom_rds))

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

readCounts <- readRDS(counts)
prom_anno <- readRDS(prom_rds)


print(sum(rownames(readCounts) %in% mcols(promoterCoordinates(prom_anno))$promoterId))
print(nrow(readCounts))
print(length(mcols(promoterCoordinates(prom_anno))$promoterId))

# Update promoterCoordinates to keep only external promoters
message("Filtering internal promoters...")
pc <- proActiv:::promoterCoordinates(prom_anno)
internal_flag <- mcols(pc)$internalPromoter
internal_flag[is.na(internal_flag)] <- FALSE
mcols(pc)$internalPromoter <- internal_flag
message("class of internal_flag:"); print(class(internal_flag))
message("typeof internal_flag:"); print(typeof(internal_flag))
keep <- !mcols(pc)$internalPromoter
promoterCoordinates(prom_anno) <- pc[keep]


# Subset readCounts to external promoters only
#keep_ids <- mcols(promoterCoordinates(prom_anno))$promoterId
#any(is.na(keep_ids)) 
#sum(keep_ids %in% rownames(readCounts))
#length(keep_ids)
#readCounts <- readCounts[as.character(keep_ids), , drop = FALSE]

keep_ids <- mcols(pc)$promoterId[keep]
valid_ids <- intersect(rownames(readCounts), as.character(keep_ids))
readCounts <- readCounts[valid_ids, , drop = FALSE]

# Re-add geneId to promoterCoordinates after filtering
pid_map <- proActiv:::promoterIdMapping(prom_anno)
pc_ext <- promoterCoordinates(prom_anno)
mcols(pc_ext)$geneId <- pid_map$geneId[match(mcols(pc_ext)$promoterId, pid_map$promoterId)]

# Rank promoter position by strand direction
message("Ranking promoter positions by strand direction...")
df_pos <- data.frame(
  promoterId = mcols(pc_ext)$promoterId,
  geneId     = mcols(pc_ext)$geneId,
  strand     = as.character(strand(pc_ext))
)

df_pos <- df_pos %>%
  group_by(geneId) %>%
  arrange(geneId, if_else(strand == "+", promoterId, -promoterId), .by_group = TRUE) %>%
  mutate(promoterPosition = row_number()) %>%
  ungroup()

# -------------------------
# Normalization
# -------------------------
message("Normalizing promoter counts...")
cnt_mat <- readCounts   
saveRDS(cnt_mat, file.path(out_dir, "raw_promoter_counts.rds"))      
cnt_mat[is.na(cnt_mat)] <- 0
if (norm_method == "deseq2") {
  dds <- DESeqDataSetFromMatrix(round(cnt_mat),
                                colData = data.frame(row.names = colnames(cnt_mat)),
                                design = ~ 1)
  dds <- DESeq(dds, quiet = TRUE)
  # Export normalized counts
  norm_counts <- counts(dds, normalized = TRUE)
  saveRDS(norm_counts, file.path(out_dir, "normalized_promoter_counts.rds")) 
  # Export size factors
  sizeFactors <- sizeFactors(dds) 
  saveRDS(sizeFactors, file.path(out_dir, "size_factors.rds"))  
  # Export estimated mu (fitted mean values for each promoter in each sample)
  mu_mat <- assays(dds)[["mu"]]
  saveRDS(mu_mat, file.path(out_dir, "promoter_fitted_mu.rds"))
  # Export dispersion
  disp_vec <- mcols(dds)$dispersion
  saveRDS(disp_vec, file.path(out_dir, "dispersion_per_gene.rds"))

  # Export size for NB
  size_vec <- 1 / disp_vec
  size_vec[is.na(size_vec) | !is.finite(size_vec) | size_vec <= 0] <- NA
  saveRDS(size_vec, file.path(out_dir, "NB_size_per_gene.rds"))
} else if (norm_method == "edger") {
 library(edgeR)
  dge <- DGEList(counts = cnt_mat)
  dge <- calcNormFactors(dge)
  # Get dispersion and mu
  design <- model.matrix(~1, data = dge$samples)
  dge <- estimateDisp(dge, design)
  fit <- glmFit(dge, design)
  # mu: Row is gene, column is sample
  mu_mat <- fit$fitted.values
  saveRDS(mu_mat, file.path(out_dir, "promoter_fitted_mu.rds"))
  # Normalized counts
  norm_counts <- cpm(dge, normalized.lib.sizes = TRUE)
  saveRDS(norm_counts, file.path(out_dir, "normalized_promoter_counts.rds"))
  # size factors
  sizeFactors <- dge$samples$norm.factors
  names(sizeFactors) <- colnames(cnt_mat)
  saveRDS(sizeFactors, file.path(out_dir, "size_factors.rds"))
  # Dispersion
  disp_vec <- dge$tagwise.dispersion
  saveRDS(disp_vec, file.path(out_dir, "dispersion_per_gene.rds"))
  # NB size = 1 / dispersion
  size_vec <- 1 / disp_vec
  size_vec[is.na(size_vec) | !is.finite(size_vec) | size_vec <= 0] <- NA
  saveRDS(size_vec, file.path(out_dir, "NB_size_per_gene.rds"))

} else stop("norm_method must be deseq2 or edger")

# -------------------------
# Save sizefactor
# -------------------------
message("Saving size factor...")
# Conver to data.frame before merging
sizefactor_df <- data.frame(
  newname = names(sizeFactors),
  sizeFactor = as.numeric(sizeFactors),
  stringsAsFactors = FALSE
)
# Each sample has both size factor and sample
sample_info <- data.frame(
  sample     = samples,
  condition  = condition,
  stringsAsFactors = FALSE
)

sizefactor_df$sample <- sizefactor_df$newname

write.table(sizefactor_df[, c("sample", "sizeFactor")],
            file = file.path(out_dir, "size_factor.tsv"),
            sep  = "\t", quote = FALSE, row.names = FALSE)

abs_activity <- proActiv:::getAbsolutePromoterActivity(norm_counts, prom_anno)
gene_expression <- proActiv:::getGeneExpression(abs_activity)
#gene_expression <- proActiv:::getGeneExpression(abs_activity)
rel_activity <- proActiv:::getRelativePromoterActivity(abs_activity,gene_expression)


abs_mat   <- as.matrix(abs_activity[, -(1:2)])
gene_mat  <- as.matrix(gene_expression[, -(1)])
rel_mat <- as.matrix(rel_activity[, -(1:2)])

gene_id_vec <- abs_activity$geneId
gene_expr_df <- gene_expression
gene_expr_mat <- as.matrix(gene_expr_df[, -1])


idx <- match(gene_id_vec, gene_expr_df$geneId)
if (any(is.na(idx))) {
  stop("Some promoter geneIds not found in gene_expression$geneId")
}

gene_mat <- gene_expr_mat[idx, , drop = FALSE]
rownames(gene_mat) <- abs_activity$promoterId



# Ensure sample order consistency
stopifnot(identical(colnames(abs_mat), colnames(readCounts)))
stopifnot(identical(colnames(gene_mat), colnames(readCounts)))
stopifnot(identical(rownames(cnt_mat), rownames(gene_mat)))

# -------------------------
# Construct SummarizedExperiment
# -------------------------

cat("Dimensions check:\n")
cat("cnt_mat: ", dim(cnt_mat), "\n")
cat("norm_counts: ", dim(norm_counts), "\n")
cat("abs_mat: ", dim(abs_mat), "\n")
head(abs_mat)
cat("gene_mat: ", dim(gene_mat), "\n")
head(gene_mat)
cat("rel_mat: ", dim(rel_mat), "\n")
head(rel_mat)

se <- SummarizedExperiment(
  assays = list(
    promoterCounts            = cnt_mat,
    normalizedPromoterCounts  = norm_counts,
    absolutePromoterActivity  = abs_mat,
    geneExpression            = gene_mat,
    relativePromoterActivity  = rel_mat
  ),
  
  rowData = DataFrame(
    promoterId       = abs_activity$promoterId,
    geneId           = abs_activity$geneId,
    promoterPosition = df_pos$promoterPosition[match(abs_activity$promoterId, df_pos$promoterId)],
    internalPromoter = FALSE
  ),
  
  colData = DataFrame(
    condition = condition,
    row.names = colnames(cnt_mat)
  )
)

se$sampleName <- colnames(se)

saveRDS(se, file.path(out_dir, "Promoter_activity_SE.rds"))


# Set condition labels from sample names
sample_names <- colnames(se)
#condition <- condition

# Collapse samples to per-condition averages and classify promoter status
se_sum <- proActiv:::summarizeAcrossCondition(se, condition)

df <- as.data.frame(rowData(se_sum))
saveRDS(df, file.path(out_dir, "Summary_classified.rds"))

# -------------------------
# Print promoter class distributions
# -------------------------
class_cols <- grep("\\.class$", colnames(df), value = TRUE)
invisible(lapply(class_cols, function(col) {
  cat("\n==", col, "==\n")
  print(table(df[[col]], useNA = "always"))
}))

# -------------------------
# Fix low-expression Major promoters (< 0.25) to Inactive before gene classification
# -------------------------
comparison <- make.names(comparison) 
message("Fixing major promoters...")
cat("Available columns:", colnames(df), "\n")
cat("Comparison values:", comparison, "\n")

message("Fixing major promoters...")
for (cond in comparison) {
  class_col <- paste0(cond, ".class")
  mean_col  <- paste0(cond, ".mean")
  df[[class_col]] <- as.character(df[[class_col]])

  cat("==", cond, "==\n")
  cat("class_col type:", class(df[[class_col]]), "\n")
  cat("mean_col range:\n")
  print(summary(df[[mean_col]]))
  
  cat("table of class:\n")
  print(table(df[[class_col]]))
  
  cat("Number of 'Major' with mean < 0.25:\n")
  #sel <- df[[class_col]] == "Major" & df[[mean_col]] < 0.25
  sel <- df[[class_col]] == "Major" & (is.na(df[[mean_col]]) | df[[mean_col]] < 0.25)
  sel[is.na(sel)] <- TRUE
  print(sum(sel))
  df[[class_col]][sel] <- "Inactive"
}

#rowData(se_sum) <- df
message("Before:", class(rowData(se_sum)))
message("After:", class(df)) 
rowData(se_sum) <- S4Vectors::DataFrame(df)

for (cond in comparison) {
  class_col <- paste0(cond, ".class")
  cat("\n== After correction:", class_col, "==\n")
  print(table(rowData(se_sum)[[class_col]], useNA = "always"))
}
# -------------------------
# Classify genes based on promoter class
# -------------------------
active_promoter_classification <- function(z) {
  if (any(z %in% "Major")) {
    if (any(z %in% "Minor")) {
      category <- "Multipromoter.Multiactive"
    } else if (!any(z %in% "Minor") & any(z %in% "Inactive")) {
      category <- "Multipromoter.Singleactive"
    } else {
      category <- "Singlepromoter.Singleactive"
    }
  } else {
    category <- "Inactive"  
  }
  return(category)
}


# Classify each condition separately
for (cond in comparison) {

  class_col <- paste0(cond, ".class")
  prom_df <- df[, c("geneId", class_col)]
  colnames(prom_df)[2] <- "promoterClass"
  
  gene_cat <- sapply(split(prom_df$promoterClass, prom_df$geneId),
                     active_promoter_classification)
  
  cat("====", cond, "====\n")
  print(table(gene_cat))
  df[[paste0("geneCategory_", cond)]] <- gene_cat[as.character(df$geneId)]
  df_sub <- df[!is.na(df[[paste0("geneCategory_", cond)]]), ]
  promoter_count_per_category <- table(df_sub[[paste0("geneCategory_", cond)]])

  cat("---- Promoter counts per category ----\n")
  print(promoter_count_per_category)
}

# -------------------------
# Save result files
# -------------------------
saveRDS(df, file.path(out_dir, "Summary_classified_rowData.rds"))


# Extract all condition from *.class and geneCategory_* columns
cond_class_cols <- grep("\\.class$", colnames(df), value = TRUE)
cond_category_cols <- grep("^geneCategory_", colnames(df), value = TRUE)

for (i in seq_along(cond_class_cols)) {
  cond <- sub("\\.class$", "", cond_class_cols[i])
  class_col <- cond_class_cols[i]
  cat_col <- cond_category_cols[i]
  
  message("Processing condition: ", cond)

  # Save Major and Minor promoterId lists
  major_ids <- df$promoterId[df[[class_col]] == "Major"]
  minor_ids <- df$promoterId[df[[class_col]] == "Minor"]
  saveRDS(major_ids, file.path(out_dir, paste0("Major_promoterId_", cond, ".rds")))
  saveRDS(minor_ids, file.path(out_dir, paste0("Minor_promoterId_", cond, ".rds")))

  # Save geneId lists per gene category
  gene_cat <- df[[cat_col]]
  if (!(cat_col %in% colnames(df))) {
    message("  [Warning] Missing column: ", cat_col, "; skipping gene/promoter category outputs.")
    next
  }

  gene_cat <- df[[cat_col]]
  if (length(gene_cat) == 0) {
    message("  [Warning] Empty gene_cat for ", cond, "; skipping.")
    next
  }
  message("class_col: ", class_col, "; cat_col: ", cat_col)
  message("gene_cat length: ", length(df[[cat_col]]))
  head(gene_cat)
  head(df)
  names(gene_cat) <- df$geneId
  for (catname in unique(gene_cat)) {
    gene_list <- unique(df$geneId[gene_cat == catname])
    saveRDS(gene_list, file.path(out_dir, paste0("GeneId_", catname, "_", cond, ".rds")))
  }

  # Save promoterId lists per gene category
  for (catname in unique(gene_cat)) {
    prom_list <- df$promoterId[gene_cat == catname]
    saveRDS(prom_list, file.path(out_dir, paste0("PromoterId_", catname, "_", cond, ".rds")))
  }
}


# -------------------------
# Call dexseq_fit2.R and pass arguments
# -------------------------

message("getwd() = ", getwd())
message("fit_script = ", fit_script, " ; exists = ", file.exists(fit_script))

system2("Rscript", args = c(
    fit_script, 
    out_dir, 
    file.path(out_dir, "raw_promoter_counts.rds"),
    file.path(out_dir, "NB_size_per_gene.rds"),
    file.path(out_dir, "promoter_fitted_mu.rds"),
    conditions_str, 
    batch_str))
