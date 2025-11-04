#/mnt/citadel2/research/shared/SnakeAltPromoter_processed/genome/.snakemake_conda/ac686adccb9fe18c6effcb34e4ad27cb_


#!/usr/bin/env Rscript


# ----- Libraries ------------------------------------------------
library(proActiv)
library(ggplot2)
library(dplyr)
library(scales)
library(tidyr)
library(Rtsne)
# future.apply is used for parallel processing, speeding up for large datasets.
library(future.apply)
library(ggpubr)

# --------------------- #
# Parse and validate inputs
# --------------------- #
output_dir <- "/mnt/citadel2/research/shared/SnakeAltPromoter_processed"
promoter_rds <- "/mnt/citadel2/research/shared/AltPromoterFlow/genome/organisms/hg38/Annotation/proActiv_promoter_annotation.rds"
row_data     <- "/mnt/citadel2/research/shared/AltPromoterFlow/RNA_NEW/proactiv/merge_tot/Summary_classified_rowData.rds"
se_file      <- "/mnt/citadel2/research/shared/AltPromoterFlow/RNA_NEW/proactiv/merge_tot/Promoter_activity_SE.rds"
cell_lines_raw <- strsplit("GM12878 H1-hESC K562", " ")[[1]]
cell_lines <- make.names(cell_lines_raw)
names(cell_lines_raw) <- cell_lines 
print(cell_lines)
print(cell_lines_raw)
condition <- strsplit("GM12878,GM12878,GM12878,H1-hESC,H1-hESC,H1-hESC,K562,K562,K562", ",")[[1]]
print(condition)

stopifnot(file.exists(promoter_rds),
          file.exists(row_data))

if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# --------------------- #
# Load data and set theme
# --------------------- #

# ----- Load Global data -------------------
# Data frame to access subset
rowData <- as.data.frame(readRDS(row_data))
cat("Reading se_file:", se_file, "\n")
se <- readRDS(se_file)
print(class(se))
activity <- assays(se)$absolutePromoterActivity

# ----- Theme and Color for all plots -------------------------------------------
plot_theme <- theme_classic() + theme(
  plot.title = element_text(hjust = .5, family = "Helvetica"),
  axis.text  = element_text(colour = "black", family = "Helvetica"),
  axis.title = element_text(colour = "black", family = "Helvetica"))
cols_type <- c("Major" = "#ff1e56", "Minor" = "#ffac41", "Inactive" = "#323232")



# Plot 9: t-SNE -----------------------------------

# A t-SNE plot visualizing promoter activity across all samples (every point is a sample).
# X axis is tSNE1; y axis is tSNE2.
# Ideally, samples from the same cell line cluster together, indicating similar promoter activity patterns.
message("Generating global t-SNE plot …")
# Load the complete result output
# Extract promoter activity matrix from the assay
mat <- assays(se)$absolutePromoterActivity
# Filter out promoters with no activity across all samples.
mat <- mat[rowSums(mat) > 0, ]

# Transpose the matrix: now rows = samples (observation), columns = promoters (feature).
# tsne expects one observation per row.
tsne_df <- data.frame(t(mat))

# Create a sample label vector (cell line name for each row/sample)
#cond <- factor(rep(cell_lines, length.out = nrow(tsne_df)))
# Create a sample label vector (cell line name for each row/sample)
#cond <- factor(condition)
condition <- make.names(condition) 
cond <- factor(condition, levels = cell_lines)
# Add a column for sample labels
tsne_df$Sample <- cond

# Set random seed for reproducibility
set.seed(40)
# Run t-SNE on expression matrix, excluding the Sample column.
#pca_res <- prcomp(tsne_df[ , !names(tsne_df) %in% "Sample"], scale.=TRUE)
#Y <- Rtsne(pca_res$x[, 1:2], perplexity=1)$Y
cat("Input matrix for t-SNE has", nrow(tsne_df), "samples (rows)\n")
print(rownames(tsne_df))
# Y is a 2D matrix with t-SNE coordinates for each sample
Y <- Rtsne(as.matrix(tsne_df[ , !names(tsne_df) %in% "Sample"]), perplexity = 1)$Y
rownames(Y) <- rownames(tsne_df)
print(Y)
# Define colors for each cell line and match by name.
cell_cols <- setNames(c("#ff1e56", "#ffac41", "#323232")[seq_along(cell_lines)],
                      cell_lines)

pdf(file.path(output_dir, "promoter_activity_tsne_plot_RNANEW_test_proactiv.pdf"), width = 12/2.54, height = 10/2.54)
# Allow drawing outside the plot region for legend.
par(xpd = NA)
# Set shape, filled color, border color, size, and aspect ratio for the t-SNE plot.
# ,1 is tSNE1, x axis, the first dimension after dimentionally reduced to 2D.
plot(Y[,1], Y[,2], pch = 24, bg = cell_cols[as.character(cond)], col = "black", cex = 1.4, asp = 1,
     xlab = "tSNE1", ylab = "tSNE2", main = "t-SNE plot (promoters active >=1 sample)")
# legend matches point shape, filled colors, and border color of cell lines.
legend("topright", title = "Cell Lines", legend = levels(cond), pch = 24,
       pt.bg = cell_cols[levels(cond)], col = "black", bty = "n", cex = 0.9)
dev.off()

message("t-SNE plot saved as promoter_activity_tsne_plot.pdf")

message("All plots (1-9) generated successfully.")

# -- End of script ------------------------------------------------