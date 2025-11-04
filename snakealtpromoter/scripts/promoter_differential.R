#!/usr/bin/env Rscript

# ---------------------
# Salmon merge script
# Usage:
# Rscript salmon_promoter_merge.R <out_dir> <promoter_rds> <norm_method> <samples> <conditions> <newnames> <count_files...>
# ---------------------

library(proActiv)   # 1.16
library(DESeq2)
library(edgeR)
library(dplyr)
library(matrixStats)

# --------------------
# Parse arguments
# --------------------

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 8)
  stop("Usage: Rscript salmon_promoter_merge.R <out_dir> <promoter_rds> <norm_method> <samples> <conditions> <cell_lines> <newnames> <batch_str> <fit_script> <reference_condition> <count_files...>")

out_dir      <- args[1]
prom_rds     <- args[2] 
samples      <- strsplit(args[3], " ")[[1]]
conditions_str   <- args[4]
cnt_mat_file <- args[5]  # Path to the count matrix file
se_file <- args[6]
baseline_condition <- args[7]  # Baseline condition for comparison
comparison_condition <- args[8]
pFC_raw <- args[9]
gFC_raw <- args[10]
lfcshrink <- as.logical(args[11])
batch_str <- args[12]


cond <- strsplit(conditions_str, ",")[[1]]
head(cond)

batch <- strsplit(batch_str, ",")[[1]]
batch <- factor(batch)
print(batch)

stopifnot(file.exists(se_file),
          file.exists(prom_rds),
          file.exists(cnt_mat_file))

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)


## ---------------------------------------
## Differential promoter activity using DESeq2
## ---------------------------------------

# Round for DESeq2
cnt_mat <- as.matrix(readRDS(cnt_mat_file))
cnt_mat <- round(cnt_mat)
cnt_mat <- cnt_mat[complete.cases(cnt_mat), ]
storage.mode(cnt_mat) <- "integer"
head(cnt_mat)


keep <- cond %in% c(baseline_condition, comparison_condition)
cnt_mat <- cnt_mat[, keep]
head(cnt_mat)
cond <- cond[keep]
condition <- factor(cond, levels = c(baseline_condition, comparison_condition))
colData <- data.frame(condition = condition)
print("==== Check condition factor ====")
print(condition)
print("==== Check colData ====")
print(colData)
if (length(unique(batch)) == 1) {
  print("Differential analysis without batch effect correction.")
  dds <- DESeqDataSetFromMatrix(countData = round(cnt_mat),
                              colData   = data.frame(condition = condition),
                              design    = ~ condition)
} else {
  print("Differential analysis with batch effect correction.")
  dds <- DESeqDataSetFromMatrix(countData = round(cnt_mat),
                              colData   = data.frame(condition = condition, batch = batch),
                              design    = ~ condition + batch)
  }

dds <- DESeq(dds)
# Result is HF vs Healthy, >0 means HF upreg
res <- results(dds, contrast = c("condition",comparison_condition,baseline_condition))
head(res)
# shrink LFC using the normal method
if (lfcshrink){
  print("Using lfcShrink to shrink log2 fold changes.")
  resLFC <- lfcShrink(dds,coef = paste0("condition_", comparison_condition, "_vs_", baseline_condition),type = "normal")
  res <- resLFC
}
# Convert to dataframe to adjust rows and columes
act_tab <- as.data.frame(res)
act_tab$promoterId <- rownames(act_tab)

# Add corresponding geneId of the promoter
promAnno <- readRDS(prom_rds)
map <- promAnno@promoterIdMapping
#map$geneClean <- sub("\\..*$", "", map$geneId)
act_tab$geneId <- map$geneId[ match(act_tab$promoterId, map$promoterId) ]

# Change column names
colnames(act_tab)[colnames(act_tab)=="padj"] <- "FDR"
colnames(act_tab)[colnames(act_tab)=="log2FoldChange"] <- "logFC"
# Filter padj (FDR) and keep only significant results
deseq_deg <- subset(act_tab, FDR < 0.05)

##-------------------
## Save and print output
##-------------------
saveRDS(act_tab, file.path(out_dir, "Promoter_differential_activity.rds"))
saveRDS(deseq_deg, file.path(out_dir, "Promoter_differential_activity_FDR0_05.rds"))
cat("Number of differential promoters unfiltered:", nrow(act_tab), "\n")
cat("Number of differential promoters with FDR<0.05:", nrow(deseq_deg), "\n")

comparison_label <- paste0(comparison_condition, "_vs_", baseline_condition)
archive_dir <- file.path(out_dir, "comparisons_", comparison_label)
dir.create(archive_dir, recursive = TRUE, showWarnings = FALSE)
file.copy(
  from = file.path(out_dir, "Promoter_differential_activity.rds"),
  to   = file.path(archive_dir, "Promoter_differential_activity.rds"),
  overwrite = TRUE
)
file.copy(
  from = file.path(out_dir, "Promoter_differential_activity_FDR0_05.rds"),
  to   = file.path(archive_dir, "Promoter_differential_activity_FDR0_05.rds"),
  overwrite = TRUE
)

## ---------------------------------------
## Differential promoter usage with proActiv (multi-promoter only)
## ---------------------------------------
se <- readRDS(se_file)
# example: keep only A and B
se <- se[, se$condition %in% c(comparison_condition, baseline_condition)]

cat("Running getAlternativePromoters with default promoterFC = 2.0, geneFC = 1.5. \n")

alt <- getAlternativePromoters(
  result             = se,
  referenceCondition = comparison_condition,
  promoterFC         = 2.0,
  geneFC             = 1.5
)
up_count <- if (!is.null(alt$upReg)) length(rownames(alt$upReg)) else 0
down_count <- if (!is.null(alt$downReg)) length(rownames(alt$downReg)) else 0
cat(sprintf("Number of upregulated promoters for promoterFC = 2.0 & geneFC = 1.5: %d \n", up_count))
cat(sprintf("Number of downregulated promoters for promoterFC = 2.0 & geneFC = 1.5: %d \n", down_count))
outfile <- file.path(out_dir, sprintf("Differential_promoter_usage_pFC2_gFC1_5.rds"))
saveRDS(alt, outfile)
file.copy(
  from = file.path(out_dir, sprintf("Differential_promoter_usage_pFC2_gFC1_5.rds")),
  to   = file.path(archive_dir, sprintf("Differential_promoter_usage_pFC2_gFC1_5.rds")),
  overwrite = TRUE
)

# Run getAlternativePromoters
if (pFC_raw != 2.0 & gFC_raw != 1.5) {
  pFC <- as.numeric(pFC_raw)
  gFC <- as.numeric(gFC_raw)
  cat("Running getAlternativePromoters with specified promoterFC =", pFC, "geneFC =", gFC, "\n")
  # Run getAlternativePromoters
  alt <- getAlternativePromoters(
    result             = se,
    referenceCondition = comparison_condition,
    promoterFC         = pFC,
    geneFC             = gFC
  )
  up_count <- if (!is.null(alt$upReg)) length(rownames(alt$upReg)) else 0
  down_count <- if (!is.null(alt$downReg)) length(rownames(alt$downReg)) else 0
  cat(sprintf("Number of upregulated promoters for promoterFC=%.2f & geneFC=%.2f: %d\n", 
              pFC, gFC, up_count))
  cat(sprintf("Number of downregulated promoters for promoterFC=%.2f & geneFC=%.2f: %d\n", 
              pFC, gFC, down_count))
  outfile <- file.path(out_dir, sprintf("Differential_promoter_usage_pFC%.2f_gFC%.2f.rds", pFC, gFC))
  saveRDS(alt, outfile)
  file.copy(
    from = file.path(out_dir, sprintf("Differential_promoter_usage_pFC%.2f_gFC%.2f.rds", pFC, gFC)),
    to   = file.path(archive_dir, sprintf("Differential_promoter_usage_pFC%.2f_gFC%.2f.rds", pFC, gFC)),
    overwrite = TRUE
  )


} else {

  message(
    "No additional minimum promoter fold change and maximum gene fold change thresholds provided. \n",
    "Please specify them using --max_gFC and --min_pFC \n",
    "If you would like to perform differential promoter usage analysis on additional thresholds. \n"
  )

  thresholds <- list(
    c(1.3, 1.1),
    c(1.5, 1.3),
    c(1.8, 1.5),
    c(2.2, 2.0),
    c(2.5, 2.0),
    c(3.0, 2.5),
    c(3.5, 3.0),
    c(4.0, 3.9)
  )

  for (th in thresholds) {
    pFC <- th[[1]]
    gFC <- th[[2]]
    cat("Running getAlternativePromoters with promoterFC =", pFC, "geneFC =", gFC, "\n")
    alt <- getAlternativePromoters(result = se,
                                  referenceCondition = comparison_condition,
                                  promoterFC = pFC,
                                  geneFC = gFC,
                                  maxPval = 1)
  up_count <- if (!is.null(alt$upReg)) length(rownames(alt$upReg)) else 0
  down_count <- if (!is.null(alt$downReg)) length(rownames(alt$downReg)) else 0
  cat(sprintf("Number of upregulated promoters for promoterFC=%.2f & geneFC=%.2f: %d\n", 
              pFC, gFC, up_count))
  cat(sprintf("Number of downregulated promoters for promoterFC=%.2f & geneFC=%.2f: %d\n", 
              pFC, gFC, down_count))
  outfile <- file.path(out_dir, sprintf("Differential_promoter_usage_pFC%.2f_gFC%.2f.rds", pFC, gFC))
  saveRDS(alt, outfile)
  file.copy(
    from = file.path(out_dir, sprintf("Differential_promoter_usage_pFC%.2f_gFC%.2f.rds", pFC, gFC)),
    to   = file.path(archive_dir, sprintf("Differential_promoter_usage_pFC%.2f_gFC%.2f.rds", pFC, gFC)),
    overwrite = TRUE
  )
  cat("\n")
  }

}