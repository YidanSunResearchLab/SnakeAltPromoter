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

cat("---- Received Arguments ----\n", file=stderr())
for (i in seq_along(args)) {
  cat(paste0("[", i, "] ", args[[i]], "\n"), file=stderr())
}
cat("----------------------------\n", file=stderr())

if (length(args) < 10)
  stop("Usage: Rscript salmon_promoter_merge.R <out_dir> <promoter_rds> <norm_method> <samples> <conditions> <cell_lines> <newnames> <batch_str> <fit_script> <reference_condition> <count_files...>")

out_dir      <- args[1]
prom_rds     <- args[2] 
samples      <- strsplit(args[3], " ")[[1]]
conditions_str   <- args[4]
cnt_mat_file <- args[5]  # Path to the count matrix file
baseline_condition <- args[6]  # Baseline condition for comparison
reference_condition <- args[7]
deltaPU_threshold <- as.numeric(args[8])
usage_min         <- as.numeric(args[9])
result_file <- args[10]  # Path to the result file

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

conditions <- strsplit(conditions_str, ",")[[1]]
stopifnot(length(samples) == length(conditions),
          file.exists(prom_rds),
          file.exists(cnt_mat_file))

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)


## ---------------------------------------
## Differentially expressed genes using DESeq2
## ---------------------------------------

##-------------------
## Generate raw gene expression prepare DESeq2 analysis for change between conditions
##-------------------
promAnno <- readRDS(prom_rds)
map <- promAnno@promoterIdMapping
map$geneClean <- sub("\\..*$", "", map$geneId)

cnt_mat <- as.matrix(readRDS(cnt_mat_file))

cond  <- factor(conditions)  # Group: Healthy / Heart_Failure
print(cond)
# Find gene id without version number that corresponds to promoters in cnt_mat using map
gene_vec <- map$geneClean[ match(rownames(cnt_mat), map$promoterId) ]
# Sum promoter counts to gene
keep_promoters <- rowSums(cnt_mat) >= 5
cnt_mat_filt <- cnt_mat[keep_promoters, ]
gene_vec_filt <- map$geneClean[ match(rownames(cnt_mat_filt), map$promoterId) ]
gene_sum_raw <- rowsum(cnt_mat_filt, group = gene_vec_filt)
#gene_sum_raw <- rowsum(cnt_mat, group = gene_vec) 
contrast_levels <- c(baseline_condition, reference_condition)
cond <- factor(conditions, levels = contrast_levels)
stopifnot(ncol(gene_sum_raw) == length(cond))
#Round gene expression raw counts because DESeq2 only takes integers
gene_sum_raw_gene <- round(gene_sum_raw)
# Add 1 to avoid zero counts, DESeq2 will log2 transform the counts
gene_sum_raw_gene <- gene_sum_raw_gene + 1
storage.mode(gene_sum_raw_gene) <- "integer" 


hist(rowSums(gene_sum_raw_gene), breaks = 100, main = "Gene-level sum (counts+1)")
summary(rowSums(gene_sum_raw_gene))
table(rowSums(gene_sum_raw_gene) == 0)
#summary(gene_sum_raw_gene$baseMean)
cat("Samples per condition:\n")
print(table(cond))
cat("Total gene counts:\n")
print(summary(rowSums(gene_sum_raw_gene)))
cat("Promoter count matrix:\n")
print(dim(cnt_mat_filt))
print(summary(rowSums(cnt_mat_filt)))

##-------------------
## Perform DESeq2 analysis for fold change, p values and adjusted p values for gene expression change between conditions
##-------------------
dds <- DESeqDataSetFromMatrix(countData = round(cnt_mat),
                              colData   = data.frame(condition = cond),
                              design    = ~ condition)

dds <- DESeq(dds, 
            quiet = TRUE, 
            useT = FALSE, 
            minReplicatesForReplace = Inf)

# Result is HF vs Healthy, >0 means HF upreg
res <- results(dds, contrast = c("condition","Heart_Failure","Healthy"))
#res <- results(dds, 
#              contrast = c("condition", "Heart_Failure", "Healthy"),
#              lfcThreshold = 0, 
#              altHypothesis = "greaterAbs", 
#              independentFiltering = FALSE,
#              cooksCutoff = FALSE,
#              minmu = 1e-6)
# Convert to dataframe to adjust rows and columes
gene_tab <- as.data.frame(res)
table(gene_tab$pvalue < 0.05)
gene_tab$promoterId <- rownames(gene_tab)
colnames(gene_tab)[colnames(gene_tab)=="padj"] <- "FDR"
colnames(gene_tab)[colnames(gene_tab)=="log2FoldChange"] <- "logFC"
#Add column to assess overall gene expression level
gene_tab$logCPM <- log2(gene_tab$baseMean + 0.1)
# Filter padj (FDR) and keep only significant results
deseq_deg <- subset(gene_tab, FDR < 0.05)

##-------------------
## Save and print output
##-------------------
saveRDS(gene_tab, file.path(out_dir, "DESeq2_gene_allResults_keepall.rds"))
saveRDS(deseq_deg, file.path(out_dir, "DESeq2_gene_DE_FDR<0.05_keepall.rds"))
cat("DESeq2 promoter tested :", nrow(gene_tab), "\n")
cat("DESeq2 promoter FDR<0.05:", nrow(deseq_deg), "\n")
summary(gene_tab$logCPM)
hist(gene_tab$pvalue, breaks = 50, main = "Gene-level raw p-value")




# --------------------- #
# Identify alternative promoters if reference is specified
# --------------------- #
result <- readRDS(result_file)
#pu_mat <- assays(result)$promoterUsage
#summary(as.numeric(pu_mat))
#hist(as.numeric(pu_mat), breaks=100)
summary(result)
#Convert condition to factor with levels
colData(result)$condition <- factor(result$condition, levels = c("Healthy", "Heart_Failure"))
# Fix: convert condition to 2-level factor (reference vs. other)
ref <- reference_condition
cond <- colData(result)$condition
cond <- ifelse(cond == ref, ref, "other")
colData(result)$condition <- factor(cond, levels = c(ref, "other"))

# Use proActiv to identify alternative promoters
if (reference_condition %in% conditions && reference_condition != "") {
  alt_promoters <- getAlternativePromoters(result = result, referenceCondition = reference_condition, minAbs = 0, minRel = 0.01, promoterFC = 1, geneFC = Inf)
}
print(paste0("Number of alternative promoters: ", union(nrow(alt_promoters$upReg), nrow(alt_promoters$downReg))))
saveRDS(alt_promoters, paste0(out_dir, "altPromoter_relativeUsage_delta10pct_allSamples1pct.rds"))



## ---------------------------------------
## Relative promoter usage with DESeq2 offset (multi-promoter only)
## ---------------------------------------
# --------------------------
## Set parameters for usage difference and minimal usage
# --------------------------
deltaPU_threshold <- deltaPU_threshold      # absolute usage difference cutoff
usage_min         <- usage_min      # per-sample minimal retained usage
contrast_levels   <- c("Healthy","Heart_Failure")  # order: baseline, comparison
save_prefix       <- file.path(out_dir, "altPromoter_relativeUsage_delta10pct_allSamples1pct")

# --------------------------
## Set conditions to compare and filter for multi-promoter genes
# --------------------------
cond <- factor(conditions, levels = contrast_levels)
stopifnot(!any(is.na(cond)), ncol(cnt_mat) == length(cond))
# Map promoter in count matrix to corresponding gene
prom2gene <- map$geneClean[ match(rownames(cnt_mat), map$promoterId) ]
stopifnot(length(prom2gene) == nrow(cnt_mat))
# Number of time every gene id appear in prom2gene
tab <- table(prom2gene)
# Keep only genes that appear more than once, thus with more than one promoters
multi_genes <- names(tab[tab >= 2])
keep_multi  <- prom2gene %in% multi_genes
pc_raw      <- cnt_mat[keep_multi, ]
p2g         <- prom2gene[keep_multi]

cat("Promoters (multi-gene subset):", nrow(pc_raw), 
    " | Genes:", length(unique(p2g)), "\n")

# --------------------------
## Calculate gene expression counts and offset
# --------------------------
# Sum promoter counts to gene
gene_sum_raw <- rowsum(pc_raw, group = p2g)        # gene x samples
gene_total_for_prom <- gene_sum_raw[p2g, , drop=FALSE]
# Use promoter id as rowname for gene counts matrix
rownames(gene_total_for_prom) <- rownames(pc_raw)

# Ensure dimnames for offset
if (is.null(rownames(gene_total_for_prom))) rownames(gene_total_for_prom) <- rownames(pc_raw)
if (is.null(colnames(gene_total_for_prom))) colnames(gene_total_for_prom) <- colnames(pc_raw)


# Offset is used to consider promoter usage change
# Set offset as the overall gene expression level of that sample
log_offset <- log(gene_total_for_prom + 1)

# Ensure dimensions align
stopifnot(nrow(log_offset) == nrow(pc_raw),
          ncol(log_offset) == ncol(pc_raw),
          all(rownames(log_offset) == rownames(pc_raw)),
          all(colnames(log_offset) == colnames(pc_raw)))

#--------------------------
## Perform DESeq2 with offset to calculate p values and adjusted p values
#--------------------------
dds_usage <- DESeqDataSetFromMatrix(
  countData = round(pc_raw),                # ensure integer-like
  colData   = data.frame(condition = cond),
  design    = ~ condition
)

# Turn off normalization
sizeFactors(dds_usage) <- rep(1, ncol(dds_usage))
# Promoter usage is the ratio of promoter counts to gene expression 
# Set gene expression offset
# Use promoter count minus log of gene expression for DESeq2 to compute log of promoter usage to find significant change
assays(dds_usage, withDimnames = FALSE)[["offset"]] <- log_offset

# Use T distribution for small sample size, do not process outliers
dds_usage <- DESeq(dds_usage,
                   useT = FALSE,
                   minReplicatesForReplace = Inf,
                   quiet = TRUE)
# Convert to dataframe and contrast conditions
res_usage <- as.data.frame(
  results(dds_usage,
          contrast = c("condition", contrast_levels[2], contrast_levels[1]))
)
res_usage$promoterId <- rownames(res_usage)

# --------------------------
## Compute raw usage matrix, filter low usage and low usage change
# --------------------------
# Compute promoter usage for averaging and using linear filters
usage_mat <- pc_raw / gene_sum_raw[p2g, ]
usage_mat[!is.finite(usage_mat)] <- NA

# Mask promoter with usage ratio smaller than threshold as NA
# A promoter must have all samples satisfying this condition to be kept
usage_masked <- usage_mat
usage_masked[usage_masked < usage_min] <- NA
keep_allSamples <- rowSums(is.na(usage_masked)) == 0
usage_keep <- usage_masked[keep_allSamples, ]
p2g_keep   <- p2g[keep_allSamples]
res_usage_aligned <- res_usage[match(rownames(usage_keep), res_usage$promoterId), ]
stopifnot(all(rownames(usage_keep) == res_usage_aligned$promoterId))
cat("After per-sample usage ≥", usage_min,
    ": promoters =", nrow(usage_keep),
    " genes =", length(unique(p2g_keep)), "\n")

# Average promoter usage for samples within the same condition and compare with promoter usage average with the other
usage_baseline <- rowMeans2(usage_keep[, cond == contrast_levels[1], drop=FALSE], na.rm = FALSE)
usage_case     <- rowMeans2(usage_keep[, cond == contrast_levels[2], drop=FALSE], na.rm = FALSE)
deltaPU        <- usage_case - usage_baseline

# Only promoters with promoter usage change larger than threshold are kept
hit <- which(!is.na(deltaPU) & abs(deltaPU) >= deltaPU_threshold)

# Sanity check for padj >= pvalue
stopifnot(sum(res_usage_aligned$padj < res_usage_aligned$pvalue, na.rm=TRUE) == 0)

# --------------------------
## Organize and save output
# --------------------------
# Put results into dataframe
ap_paper <- data.frame(
  promoterId    = rownames(usage_keep)[hit],
  geneId        = p2g_keep[hit],
  usage_baseline= round(usage_baseline[hit], 4),
  usage_case    = round(usage_case[hit], 4),
  deltaPU       = round(deltaPU[hit], 3),
  log2FC_usage  = round(res_usage_aligned$log2FoldChange[hit], 3),
  pvalue_usage  = res_usage_aligned$pvalue[hit],
  padj_usage    = res_usage_aligned$padj[hit],
  direction     = ifelse(deltaPU[hit] > 0,
                         paste0(contrast_levels[2], ">", contrast_levels[1]),
                         paste0(contrast_levels[2], "<", contrast_levels[1])),
  stringsAsFactors = FALSE
)

# Print number of promoters passing check in change of significant usage / threshold in usage change
xt <- table(ap_paper$pvalue_usage < 0.05,
            ap_paper$padj_usage  < 0.05,
            useNA = "ifany")
print(xt)
cat("ΔPU threshold:", deltaPU_threshold, "\n")
cat("AP promoters Δ≥ and padj < 0.05: ", length(ap_paper$promoterId[ap_paper$deltaPU >= deltaPU_threshold & ap_paper$pvalue_usage < 0.05]), "\n", sep="")
cat("AP promoters Δ≥ and pval < 0.05: ", length(ap_paper$promoterId[ap_paper$deltaPU >= deltaPU_threshold & ap_paper$padj_usage < 0.05]), "\n", sep="")

# Save results
saveRDS(ap_paper, paste0(save_prefix, ".rds"))
write.csv(ap_paper, paste0(save_prefix, ".csv"), row.names = FALSE)

# Save full aligned table for all promoters with promoter usage larger than 1%
full_aligned <- data.frame(
  promoterId    = rownames(usage_keep),
  geneId        = p2g_keep,
  usage_baseline= round(usage_baseline,4),
  usage_case    = round(usage_case,4),
  deltaPU       = round(deltaPU,3),
  log2FC_usage  = round(res_usage_aligned$log2FoldChange,3),
  pvalue_usage  = res_usage_aligned$pvalue,
  padj_usage    = res_usage_aligned$padj,
  stringsAsFactors = FALSE
)
saveRDS(full_aligned, paste0(save_prefix, "_fullAligned.rds"))