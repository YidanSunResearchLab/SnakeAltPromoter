# ---------------------------
# Load Required Libraries
# ---------------------------
library(ggvenn)
library(ggplot2)
library(gridExtra)
library(grid)
library(dplyr)
# ---------------------------
# Parse Command Line Arguments
# ---------------------------
args <- commandArgs(trailingOnly = TRUE)
out_dir <- args[1]

# Load input files
CAGE_gene     <- readRDS(args[2])
CAGE_promoter <- readRDS(args[3])
CAGE_minor    <- readRDS(args[4])
CAGE_major    <- readRDS(args[5])

salmon_gene     <- readRDS(args[6])
salmon_promoter <- readRDS(args[7])
salmon_minor    <- readRDS(args[8])
salmon_major    <- readRDS(args[9])

dexseq_gene     <- readRDS(args[10])
dexseq_promoter <- readRDS(args[11])
dexseq_minor    <- readRDS(args[12])
dexseq_major    <- readRDS(args[13])

proactiv_gene     <- readRDS(args[14])
proactiv_promoter <- readRDS(args[15])
rowData   <- as.data.frame(readRDS(args[16]))
proactiv_major <- rowData$promoterId[rowData$Healthy.class == "Major"]
proactiv_minor <- rowData$promoterId[rowData$Healthy.class == "Minor"]
proactiv_major_ids <- proactiv_major[!is.na(proactiv_major)]
proactiv_minor_ids <- proactiv_minor[!is.na(proactiv_minor)]








#--------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------
# ---------------------------
# Extract promoter IDs for major minor category plot
# ---------------------------
get_ids <- function(obj, class = "Major") {
  if (is.data.frame(obj)) {
    ids <- obj$promoterId[obj[["Healthy.class"]] == class]
  } else {
    ids <- obj[["Healthy"]]$promoterId
  }
  ids[!is.na(ids)]
}

CAGE_major_ids     <- get_ids(CAGE_major)
CAGE_minor_ids     <- get_ids(CAGE_minor, class = "Minor")
salmon_major_ids   <- get_ids(salmon_major)
salmon_minor_ids   <- get_ids(salmon_minor, class = "Minor")
dexseq_major_ids   <- get_ids(dexseq_major)
dexseq_minor_ids   <- get_ids(dexseq_minor, class = "Minor")
#proactiv_major_ids <- get_ids(proactiv_major)
#proactiv_minor_ids <- get_ids(proactiv_minor, class = "Minor")

# ---------------------------
# Output shared results
# ---------------------------
write_shared_sets <- function(name, x, y) {
  shared <- intersect(x, y)
  write.table(shared, file = file.path(out_dir, paste0(name, ".txt")), 
              quote = FALSE, row.names = FALSE, col.names = FALSE)
}

shared_all_major <- Reduce(intersect, list(salmon_major_ids, proactiv_major_ids, dexseq_major_ids, CAGE_major_ids))
shared_all_minor <- Reduce(intersect, list(salmon_minor_ids, proactiv_minor_ids, dexseq_minor_ids, CAGE_minor_ids))

write_shared_sets("cage_salmon_major", CAGE_major_ids, salmon_major_ids)
write_shared_sets("cage_proactiv_major", CAGE_major_ids, proactiv_major_ids)
write_shared_sets("cage_dexseq_major", CAGE_major_ids, dexseq_major_ids)
write_shared_sets("cage_all3_major", CAGE_major_ids, shared_all_major)

write_shared_sets("cage_salmon_minor", CAGE_minor_ids, salmon_minor_ids)
write_shared_sets("cage_proactiv_minor", CAGE_minor_ids, proactiv_minor_ids)
write_shared_sets("cage_dexseq_minor", CAGE_minor_ids, dexseq_minor_ids)
write_shared_sets("cage_all3_minor", CAGE_minor_ids, shared_all_minor)

# ---------------------------
# Plot Venn diagrams
# ---------------------------
fill_colors <- c("#b3cde3", "#fbb4ae")

# Generate and return a ggvenn plot
make_venn_plot <- function(list_obj, title_text = NULL) {
  p <- ggvenn(
    list_obj,
    fill_color = fill_colors,
    stroke_size = 0.5,
    text_size = 5,
    show_percentage = FALSE,
    set_name_size = 4    
  )
  if (!is.null(title_text)) {
    p <- p + labs(title = title_text) +
      theme(plot.title = element_text(hjust = 0.5))
  }
    p <- p + 
    theme(
      plot.title = element_text(hjust = 0.5, size = 10), 
      text = element_text(size = 8),              
      legend.text = element_text(size = 8)
    )
  return(p)
}

# Prepare Venn plot objects
venn_plots <- list(
  sc_major = make_venn_plot(list(Salmon = salmon_major_ids, CAGE = CAGE_major_ids), "Salmon vs CAGE (Major)"),
  dc_major = make_venn_plot(list(DEXSeq = dexseq_major_ids, CAGE = CAGE_major_ids), "DEXSeq vs CAGE (Major)"),
  pc_major = make_venn_plot(list(proActiv = proactiv_major_ids, CAGE = CAGE_major_ids), "proActiv vs CAGE (Major)"),
  sc_minor = make_venn_plot(list(Salmon = salmon_minor_ids, CAGE = CAGE_minor_ids), "Salmon vs CAGE (Minor)"),
  dc_minor = make_venn_plot(list(DEXSeq = dexseq_minor_ids, CAGE = CAGE_minor_ids), "DEXSeq vs CAGE (Minor)"),
  pc_minor = make_venn_plot(list(proActiv = proactiv_minor_ids, CAGE = CAGE_minor_ids), "proActiv vs CAGE (Minor)")
)

# Save individual plots
for (name in names(venn_plots)) {
  ggsave(file.path(out_dir, paste0(name, "_venn_ggvenn.pdf")), venn_plots[[name]], width = 6, height = 6)
}

# Combine into major and minor panels
pdf(file.path(out_dir, "major_promoter_venn_combined.pdf"), width = 15, height = 4)
grid.arrange(venn_plots$sc_major, venn_plots$dc_major, venn_plots$pc_major, ncol = 3,
             top = textGrob("Major Promoter Overlap Across Methods vs CAGE", gp = gpar(fontsize = 10, fontface = "bold")))
dev.off()

pdf(file.path(out_dir, "minor_promoter_venn_combined.pdf"), width = 15, height = 4)
grid.arrange(venn_plots$sc_minor, venn_plots$dc_minor, venn_plots$pc_minor, ncol = 3,
             top = textGrob("Minor Promoter Overlap Across Methods vs CAGE", gp = gpar(fontsize = 10, fontface = "bold")))
dev.off()


# ---------------------------
# Plot pie charts for promoter category proportions
# ---------------------------

# Prepare summary data
pie_data <- data.frame(
  Method = rep(c("CAGE", "Salmon", "DEXSeq", "proActiv"), each = 2),
  Category = rep(c("Major", "Minor"), times = 4),
  Count = c(length(CAGE_major_ids), length(CAGE_minor_ids),
            length(salmon_major_ids), length(salmon_minor_ids),
            length(dexseq_major_ids), length(dexseq_minor_ids),
            length(proactiv_major_ids), length(proactiv_minor_ids))
)

# Make percentage labels
pie_data <- pie_data |>
  dplyr::group_by(Method) |>
  dplyr::mutate(Fraction = Count / sum(Count),
                Label = paste0(Category, "\n", round(Fraction * 100), "%"))

# Function to plot one pie chart
plot_pie <- function(df) {
  ggplot(df, aes(x = "", y = Count, fill = Category)) +
    geom_bar(stat = "identity", width = 1) +
    coord_polar("y") +
    geom_text(aes(label = Label), position = position_stack(vjust = 0.5), size = 3) +
    scale_fill_manual(values = c(Major = "#b3cde3", Minor = "#fbb4ae")) +
    ggtitle(df$Method[1]) +
    theme_void() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 10),
      legend.position = "none"
    )
}

# Split by method and generate plots
pie_plots <- lapply(split(pie_data, pie_data$Method), plot_pie)

# Arrange into 2x2 grid and save
pdf(file.path(out_dir, "promoter_category_pie_plots.pdf"), width = 8, height = 6)
gridExtra::grid.arrange(grobs = pie_plots, ncol = 2,
                        top = grid::textGrob("Proportion of Major and Minor Promoters per Method",
                                             gp = grid::gpar(fontsize = 12, fontface = "bold")))
dev.off()






#--------------------------------------------------------------------------------------------------------
##log2FC violin plots of differential promoter counts by CAGE
#--------------------------------------------------------------------------------------------------------
plot_violin_logFC_by_CAGE <- function(all_methods, cage_df,
                                      fdr_col = "FDR", logfc_col = "logFC",
                                      prefix = "log2FC_violin_byCAGE") {
  
  # Extract CAGE up/down regulated promoterId
  # (FDR < 0.05 & logFC > 0 for up, < 0 for down)
  cage_up_ids   <- unique(cage_df$promoterId[cage_df[[fdr_col]] < 0.05 & cage_df[[logfc_col]] > 0])
  print(length(cage_up_ids))
  cage_down_ids <- unique(cage_df$promoterId[cage_df[[fdr_col]] < 0.05 & cage_df[[logfc_col]] < 0])
  print(length(cage_down_ids))

 # Loop over both directions: up and down
  for (dir in c("up", "down")) {
    # Use promoter IDs based on regulation direction
    message("Now processing direction: ", dir)
    promoter_ids <- if (dir == "up") cage_up_ids else cage_down_ids
    message("Promoter count: ", length(promoter_ids))

    # get promoter log FC for each method
    violin_data_list <- lapply(names(all_methods), function(method) {
      df <- all_methods[[method]]
      if (is.null(df) || nrow(df) == 0) return(NULL)
      # Filter to keep only rows (promoters) in the method same with promoterId in the CAGE set
      df <- df[df$promoterId %in% promoter_ids, , drop = FALSE]
      # Ensure logFC column exists and is finite
      df <- df[is.finite(df[[logfc_col]]), , drop = FALSE]
      if (nrow(df) == 0) return(NULL)
      # Return a data frame with method and logFC
      data.frame(method = method, logFC = df[[logfc_col]])
    })

    # Remove NULL and combine data
    violin_data <- do.call(rbind, Filter(Negate(is.null), violin_data_list))
    if (is.null(violin_data) || nrow(violin_data) == 0) {
      message("No data for ", dir, " promoters. Skipping.")
      next
    }


    p <- ggplot(violin_data, aes(x = method, y = logFC, fill = method)) +
      geom_violin(trim = TRUE, scale = "width", width = 0.8, alpha = 0.6) +
      geom_boxplot(width = 0.1, outlier.shape = NA) +
      stat_summary(fun = median, geom = "point", shape = 21, fill = "black", size = 2) +
      theme_bw(base_size = 11) +
      ylab("log2 Fold Change") + xlab("") +
      ggtitle(paste0("Promoter log2FC (CAGE ", dir, "-regulated)")) +
      theme(plot.title = element_text(hjust = 0.5),
            axis.text.x  = element_text(angle = 45, hjust = 1),
            legend.position = "none")

    ggsave(file.path(out_dir, paste0(prefix, "_", dir, ".pdf")), p, width = 6, height = 4)
  }
}

# Create a named list of gene-level results from each method
gene_fdr <- list(CAGE = CAGE_gene,
                 Salmon = salmon_gene,
                 DEXSeq = dexseq_gene,
                 proActiv = proactiv_gene)

# Ensure all methods have the same promoterId column
lapply(gene_fdr, function(x) head(colnames(x)))
lapply(gene_fdr, nrow)

# Run the violin plot function
plot_violin_logFC_by_CAGE(gene_fdr, cage_df = gene_fdr$CAGE)
list.files(out_dir, pattern = "violin", full.names = TRUE)





## --------------------------------------------------
## MA Plot on the log2FC comparing methods with CAGE
## --------------------------------------------------
library(DESeq2)
library(ggplot2)
library(apeglm)


cage_df   <- gene_fdr$CAGE
fdr_col   <- "FDR"
logfc_col <- "logFC"
cage_up_ids   <- unique(cage_df$promoterId[cage_df[[fdr_col]] < 0.05 & cage_df[[logfc_col]] > 0])
cage_down_ids <- unique(cage_df$promoterId[cage_df[[fdr_col]] < 0.05 & cage_df[[logfc_col]] < 0])

# Read in promoter-level raw count matrices
cage_counts     <- readRDS("/mnt/citadel2/research/shared/AltPromoterFlow/CAGE_HEART/salmon/merged/promoter_counts.rds")
salmon_counts   <- readRDS("/mnt/citadel2/research/shared/AltPromoterFlow/RNA_HEART/salmon/merged/promoter_counts.rds")
dexseq_counts   <- readRDS("/mnt/citadel2/research/shared/AltPromoterFlow/RNA_HEART/dexseq/merge/promoter_counts.rds")
proactiv_counts <- readRDS("/mnt/citadel2/research/shared/AltPromoterFlow/RNA_HEART/proactiv/quantify/raw_promoter_counts.rds")
#out_dir <- "/mnt/citadel2/research/shared/AltPromoterFlow/RNA_HEART/comparison2"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# Only plot healthy sample counts
get_healthy <- function(mat) mat[, grepl("Healthy", colnames(mat)), drop = FALSE]

# Compare CAGE counts to other methods using DESeq2
compare_to_CAGE <- function(cage_counts, other_counts,
                            method_name, out_dir,
                            subset_ids   = NULL,
                            fc_ylim      = c(-5,5),
                            alpha_cut    = 0.1,
                            lfc_thresh   = 1.2) {


  subset_label <- if (is.null(subset_ids)) "" else if (identical(subset_ids, cage_up_ids)) "_CAGEup" else if (identical(subset_ids, cage_down_ids)) "_CAGEdown" else "_subset"

  # Subset counts to healthy samples
  cage_h  <- get_healthy(cage_counts)
  other_h <- get_healthy(other_counts)

 # Get common promoters across both methods
  common_ids <- intersect(rownames(cage_h), rownames(other_h))
  if (!is.null(subset_ids)) common_ids <- intersect(common_ids, subset_ids)
  if (!length(common_ids)) {
    message("No promoters after intersect for ", method_name, subset_label)
    return(invisible(NULL))
  }

  # Combine counts from CAGE and other method
  mat <- cbind(cage_h[common_ids, , drop = FALSE],
               other_h[common_ids, , drop = FALSE])

  # Filter extremely low expression
  keep <- rowMeans(mat) > 1
  mat  <- mat[keep, , drop = FALSE]
  if (!nrow(mat)) {
    message("Nothing left after mean>1 filter for ", method_name, subset_label)
    return(invisible(NULL))
  }

  # DESeq2
  cond <- factor(c(rep("CAGE", ncol(cage_h)), rep(method_name, ncol(other_h))))
  dds  <- DESeqDataSetFromMatrix(countData = round(mat),
                                 colData   = data.frame(condition = cond),
                                 design    = ~ condition)
  dds  <- DESeq(dds, quiet = TRUE)


  #res <- lfcShrink(dds, contrast = c("condition", method_name, "CAGE"),
  #                 type = "apeglm")
  # Shrink lfc
  coef_name <- paste0("condition_", method_name, "_vs_CAGE")
  res <- lfcShrink(dds, coef = coef_name, type = "apeglm")
  res$promoterId <- rownames(res)

  # Gray points of padj>=alpha_cut or NA and |log2FC| <= lfc_thresh
  gray_mask <- (is.na(res$padj) | res$padj >= alpha_cut) &
               (is.na(res$log2FoldChange) | abs(res$log2FoldChange) <= lfc_thresh)
  gray_ids <- res$promoterId[gray_mask]
  writeLines(gray_ids,
             file.path(out_dir, paste0("gray_promoters_", method_name, "_vs_CAGE", subset_label, ".txt")))

  # Save MA plot
  pdf(file.path(out_dir, paste0("MAplot_", method_name, "_vs_CAGE", subset_label, ".pdf")),
      width = 6, height = 4)
  plotMA(res, ylim = fc_ylim,
         main  = paste0("MA: ", method_name, " vs CAGE", subset_label),
         alpha = alpha_cut)
  abline(h = c(-lfc_thresh, lfc_thresh), col = "red", lty = 2)
  dev.off()


  # res_df <- as.data.frame(res)
  # res_df$status <- ifelse(gray_mask, "gray", "signif")
  # gg <- ggplot(res_df, aes(baseMean, log2FoldChange, color = status)) +
  #   geom_point(alpha = 0.6, size = 0.7) +
  #   scale_x_log10() +
  #   scale_color_manual(values = c(gray = "gray70", signif = "#E64B35")) +
  #   geom_hline(yintercept = c(-lfc_thresh, lfc_thresh), lty = 2, col = "red") +
  #   theme_bw(base_size = 11) +
  #   labs(title = paste0("MA (ggplot): ", method_name, " vs CAGE", subset_label),
  #        x = "baseMean (log10)", y = "log2FC (shrunk)") +
  #   coord_cartesian(ylim = fc_ylim)
  # ggsave(file.path(out_dir, paste0("MAplot_", method_name, "_vs_CAGE", subset_label, "_ggplot.pdf")),
  #        gg, width = 6.5, height = 5)

  message(method_name, subset_label, " gray=", length(gray_ids), " / total=", nrow(res))
  invisible(list(res = res, gray_ids = gray_ids))
}

methods <- list(
  salmon   = salmon_counts,
  proactiv = proactiv_counts,
  dexseq   = dexseq_counts
)


for (m in names(methods)) {
  compare_to_CAGE(cage_counts, methods[[m]], m, out_dir)
  compare_to_CAGE(cage_counts, methods[[m]], m, out_dir, subset_ids = cage_up_ids)
  compare_to_CAGE(cage_counts, methods[[m]], m, out_dir, subset_ids = cage_down_ids)
}














  ## --------------------------------------------------
  ## Use featurecounts to plot scatter plot
  ## --------------------------------------------------
  fc_dir   <- "/mnt/citadel2/research/shared/AltPromoterFlow/CAGE_HEART/featureCounts"
  fc_files <- list.files(fc_dir, pattern = "_promoter_counts\\.txt$", full.names = TRUE)

  # Get all promoter level output files
  fc_list <- setNames(lapply(fc_files, function(f) {
    df <- read.table(f, header = TRUE, sep = "\t",
                    comment.char = "#", stringsAsFactors = FALSE,
                    check.names = FALSE)
    sample_name <- sub("_promoter_counts\\.txt$", "", basename(f))
    # Extract promoter id
    prom_id <- sub("_$", "", df[[1]])
    # Extract count value
    counts <- df[[ncol(df)]]
    # Use promoter id as the name of vector
    names(counts) <- prom_id
    return(counts)
  }), nm = sub("_promoter_counts\\.txt$", "", basename(fc_files)))


  # Build promoter by sample count matrix
  all_ids <- unique(unlist(lapply(fc_list, names)))
  mat <- sapply(fc_list, function(x) x[all_ids])
  rownames(mat) <- all_ids
  colnames(mat) <- names(fc_list)

  # Convert to numeric, fill NAs with 0, remove promoters with zero counts
  mat <- as.matrix(mat)
  mode(mat) <- "numeric"
  mat[is.na(mat)] <- 0
  mat <- mat[rowSums(mat) > 0, , drop = FALSE]
  cage_counts <- mat



anno <- readRDS("/mnt/citadel2/research/shared/AltPromoterFlow/genome/organisms/hg38/Annotation/proActiv_promoter_annotation.rds")
promoterCoordinates <- anno@promoterCoordinates
coord_df <- as.data.frame(mcols(anno@promoterCoordinates))
sum(is.na(coord_df$internalPromoter)) 

#coord_df$internalPromoter[is.na(coord_df$internalPromoter)] <- FALSE
internal_ids <- as.character(coord_df$promoterId[ coord_df$internalPromoter %in% TRUE ])
cage_counts <- cage_counts[ !(rownames(cage_counts) %in% internal_ids) , , drop = FALSE ]


## --------------------------------------------------
## Inputs assumed already in env:
##   cage_counts, salmon_counts, dexseq_counts, proactiv_counts
##   out_dir
## --------------------------------------------------

library(DESeq2)
library(ggplot2)

# Compare counts across different methods using avearage from Healthy samples
# Extract columns corresponding to Healthy samples
get_healthy <- function(mat) mat[, grepl("Healthy", colnames(mat)), drop = FALSE]

# Normalize and return rowMeans of normalized counts
#norm_mean <- function(mat) {
#  cond <- factor(rep("A", ncol(mat)))       # dummy condition
#  dds  <- DESeqDataSetFromMatrix(round(mat), data.frame(cond), ~ cond)
#  dds  <- estimateSizeFactors(dds)
#  rowMeans(counts(dds, normalized = TRUE))
#}

joint_norm <- function(cage_mat, other_mat, tag) {
  # Find overlap between promoter ids in CAGE and other method
  common <- intersect(rownames(cage_mat), rownames(other_mat))
  if (length(common) == 0L)
      stop(sprintf("No common promoters for %s vs CAGE", tag))

  # Subset matrix to only common promoters
  cage_mat  <- round(cage_mat [common, , drop = FALSE])
  other_mat <- round(other_mat[common, , drop = FALSE])
  all_mat   <- cbind(cage_mat, other_mat)

  # sample metadata for DESeq2 input
  colData <- data.frame(
    tech = factor(c(rep("CAGE" , ncol(cage_mat)),
                    rep(tag      , ncol(other_mat)))),
    row.names = colnames(all_mat)
  )

  # no group comparison, only size factor estimation to normalize counts of methods
  dds <- DESeqDataSetFromMatrix(all_mat, colData = colData, design = ~1)
  dds <- estimateSizeFactors(dds)

  cat("\n[Joint normalisation]  tag =", tag, "\n")
  print(sizeFactors(dds))

  #get normalized count for all samples
  norm_all <- counts(dds, normalized = TRUE)
  list(
    cage_norm  = norm_all[ , colnames(cage_mat),  drop = FALSE],
    other_norm = norm_all[ , colnames(other_mat), drop = FALSE]
  )
}


# Build a comparison data frame for one method vs CAGE
# Select healthy samples and common promoters
build_comp_df <- function(cage_mat, other_mat, method_name) {
  # subset to only healthy samples
  cage_h   <- get_healthy(cage_mat)
  other_h  <- get_healthy(other_mat)
  common   <- intersect(rownames(cage_h), rownames(other_h))

# Normalize counts for both matrices that are compared
 # cage_nm  <- norm_mean(cage_h[common, , drop = FALSE], method_name = paste0(method_name, "_CAGE"))
 # other_nm <- norm_mean(other_h[common, , drop = FALSE], method_name = method_name)
 # other_nm_rescaled <- rescale_to_cage(other_nm, cage_nm, method_name)

  # joint DESeq2 normalization for CAGE and the other method that is compared
  nj <- joint_norm(cage_h, other_h, tag = method_name)
  # Calculate rowmean of each promoter within each method
  cage_norm  <- rowMeans(nj$cage_norm)
  other_norm <- rowMeans(nj$other_norm)

# Compute log2 fold change of normalized counts and create a data frame
  df <- data.frame(
    #promoterId = common,
    #cage_norm  = cage_nm,
    #other_norm = other_nm_rescaled,
    #log2FC     = log2((other_nm + 1) / (cage_nm + 1)),
    promoterId = names(cage_norm),
    cage_norm  = cage_norm,
    other_norm = other_norm,
    log2FC     = log2((other_norm + 1) / (cage_norm + 1)),
    method     = method_name,
    stringsAsFactors = FALSE
  )
  df
}

# Input expression matrices from different methods
methods <- list(
  Salmon   = salmon_counts,
  proActiv = proactiv_counts,
  DEXSeq   = dexseq_counts
)

# Combine all methods' data frames
all_df <- do.call(rbind, lapply(names(methods), function(m){
  build_comp_df(cage_counts, methods[[m]], m)
}))


# Filter out promoters with low counts in both methods (<10) and subset based on log2FC
cnt_ok <- all_df$cage_norm > 0 & all_df$other_norm > 0
set1 <- subset(all_df, cnt_ok & abs(log2FC) < 1)
set2 <- subset(all_df, cnt_ok)

# Save IDs
write.table(set1$promoterId, file.path(out_dir, "set1_ids_cnt>=10_absLFC<1.txt"),
            row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(set2$promoterId, file.path(out_dir, "set2_ids_cnt>=10.txt"),
            row.names = FALSE, col.names = FALSE, quote = FALSE)

# Plot function with correlation label
plot_density_scatter <- function(df, title, outfile){
  if (nrow(df) == 0) {
    message("No rows for ", outfile)
    return(invisible(NULL))
  }

  # Correlation (Pearson and Spearman)
  #cor_labels <- do.call(rbind, lapply(split(df, df$method), function(sub){
  #  data.frame(
  #    method = unique(sub$method),
  #    x = min(sub$cage_norm, na.rm = TRUE),
  #    y = max(sub$other_norm, na.rm = TRUE),
  #    label = sprintf("Pearson = %.3f\nSpearman = %.3f",
  #                    cor(sub$cage_norm, sub$other_norm, method = "pearson"),
  #                    cor(sub$cage_norm, sub$other_norm, method = "spearman"))
  #  )
  #}))

   # Correlation (Pearson and Spearman)
  cor_labels <- do.call(rbind, lapply(split(df, df$method), function(sub) {
    cor_pearson  <- cor.test(sub$cage_norm, sub$other_norm, method = "pearson")
    cor_spearman <- cor.test(sub$cage_norm, sub$other_norm, method = "spearman")

      # Format p-values with threshold
    pearson_p_str <- if (cor_pearson$p.value < 1e-300) {
      "p < 1e-300"
    } else {
      sprintf("p = %.2g", cor_pearson$p.value)
    }

    spearman_p_str <- if (cor_spearman$p.value < 1e-300) {
      "p < 1e-300"
    } else {
      sprintf("p = %.2g", cor_spearman$p.value)
    }

    data.frame(
      method = unique(sub$method),
      x = min(sub$cage_norm, na.rm = TRUE),
      y = max(sub$other_norm, na.rm = TRUE),
      label = sprintf("Pearson r = %.3f (p = %.2g)\nSpearman \u03C1 = %.3f (p = %.2g)",
                      cor_pearson$estimate,  cor_pearson$p.value,
                      cor_spearman$estimate, cor_spearman$p.value)
    )
  }))

  # Create scatter plot, use log10 to scale axes
  #  coord_fixed(ratio = 1) +
  p <- ggplot(df, aes(x = cage_norm, y = other_norm)) +
    geom_point(alpha = 0.5, size = 0.6) +
    geom_abline(slope = 1, intercept = 0, color = "gray40", linetype = "dashed") +
    scale_x_log10() + scale_y_log10() +
    facet_wrap(~ method, nrow = 1, scales = "free") +
    geom_text(data = cor_labels, aes(x = x, y = y, label = label),
              hjust = 0, vjust = 1, size = 3.3, inherit.aes = FALSE) +
    labs(x = "CAGE normalized count", y = "Other method normalized count",
         title = title) +
    theme_bw(base_size = 11) +
    theme(plot.title = element_text(hjust = 0.5))

  ggsave(file.path(out_dir, outfile), p, width = 12, height = 4)
}

#plot_density_scatter(set1, "Scatter (counts ≥10 & |log2FC| < 1)", "scatter_set1_cnt10_absLFC1.pdf")
#plot_density_scatter(set2, "Scatter (counts ≥10)",                "scatter_set2_cnt10.pdf")

library(MASS)
library(viridis)



plot_density_scatter <- function(df, title, outfile, use_density = TRUE, n_grid = 200) {
  if (nrow(df) == 0) {
    message("No rows for ", outfile); return(invisible(NULL))
  }

  df$logC <- log10(df$cage_norm  + 1)
  df$logO <- log10(df$other_norm + 1)
  #df$logC <- df$cage_norm
  #df$logO <- df$other_norm
  

  df$density <- NA_real_
  if (use_density) {
    keep <- is.finite(df$logC) & is.finite(df$logO)
    dens <- with(df[keep, ], MASS::kde2d(logC, logO, n = n_grid))
    ix <- findInterval(df$logC[keep], dens$x)
    iy <- findInterval(df$logO[keep], dens$y)
    df$density[keep] <- log10(dens$z[cbind(ix, iy)] + 1e-8)
  }

  cor_labels <- do.call(rbind, lapply(split(df, df$method), function(sub) {
    #x <- sub$logC; y <- sub$logO
    x <- sub$cage_norm
    y <- sub$other_norm
    cp <- cor.test(x, y, method = "pearson")
    cs <- cor.test(x, y, method = "spearman")
    data.frame(
      method = unique(sub$method),
      x = min(sub$logC, na.rm = TRUE),
      y = max(sub$logO, na.rm = TRUE),
      #x = log10(20),
      #y = log10(4e6), 
      label = sprintf("Pearson r = %.3f (p %.1g)\nSpearman ρ = %.3f (p %.1g)",
                      cp$estimate, cp$p.value, cs$estimate, cs$p.value)
    )
  }))


  lims <- range(c(df$logC, df$logO), finite = TRUE)

  p <- ggplot(df, aes(x = cage_norm, y = other_norm)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey50") +
    scale_x_log10(limits = 10^lims, breaks = scales::log_breaks()) +
    scale_y_log10(limits = 10^lims, breaks = scales::log_breaks()) +
    facet_wrap(~method, nrow = 1) +
    geom_text(data = cor_labels, aes(x = 10^x, y = 10^y, label = label),
              hjust = 0, vjust = 1, size = 3.2, inherit.aes = FALSE) +
    labs(x = "CAGE normalized count",
         y = "Other method normalized count",
         title = title) +
    theme_bw(base_size = 11) +
    theme(plot.title = element_text(hjust = 0.5))

  if (use_density) {
    p <- p + geom_point(aes(color = density), size = 1, alpha = 0.85) +
      scale_color_viridis_c(option = "A", name = "log10(density)")
  } else {
    p <- p + geom_point(alpha = 0.5, size = 0.7)
  }

  ggsave(file.path(out_dir, outfile), p, width = 12, height = 4)
}






## --------------------------------------------------
## Use featurecounts to plot scatter plot
## --------------------------------------------------

#Load and standarize promoter id row names across methods
methods <- list(
  Salmon   = salmon_counts,
  proActiv = proactiv_counts,
  DEXSeq   = dexseq_counts
)

#methods <- lapply(methods, function(m) { rownames(m) <- as.character(rownames(m)); m })
methods <- lapply(methods, function(m){
  m <- m[ !(rownames(m) %in% internal_ids) , , drop = FALSE ]
  rownames(m) <- as.character(rownames(m))
  m
})


# Find promoter ids shared by CAGE and other method used
#valid_ids <- Reduce(intersect, c(list(rownames(cage_counts)),
#                                 lapply(methods, rownames)))

#cat("common promoterId:", length(valid_ids), "\n")
#if (length(valid_ids) == 0) stop("No common promoterIds. Check your ID harmonization.")
# Subset matrixes to common promoters
#cage_counts <- cage_counts[valid_ids, , drop = FALSE]
#methods     <- lapply(methods, function(m) m[valid_ids, , drop = FALSE])
#sum(rownames(cage_counts) %in% no_intron_ids) # should be smaller than previous but should not be 0


# Compute normalized expression, log2fold change, and combine to dataframe
dfs <- lapply(names(methods), function(m) {
  message("Processing ", m)
  build_comp_df(cage_counts, methods[[m]], m)
})
sizes <- sapply(dfs, nrow); print(sizes)


# Filter low expression and subset log2fold change
dfs <- Filter(function(x) nrow(x) > 0, dfs)
if (!length(dfs)) stop("All methods returned 0 rows after harmonization.")
all_df <- do.call(rbind, dfs)
cnt_ok <- all_df$cage_norm > 0 & all_df$other_norm > 0
set1   <- subset(all_df, cnt_ok & abs(log2FC) < 1)
set2   <- subset(all_df, cnt_ok)

write.table(set1$promoterId, file.path(out_dir, "featureCounts_set1_ids_cnt>=10_absLFC<1_IP.txt"),
            row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(set2$promoterId, file.path(out_dir, "featureCounts_set2_ids_cnt>=10_IP.txt"),
            row.names = FALSE, col.names = FALSE, quote = FALSE)

plot_density_scatter(set1, "Scatter (counts ≥10 & |log2FC| < 1)", "scatter_set1_featureCounts_density.pdf")
plot_density_scatter(set2, "Scatter (counts ≥10)",                "scatter_set2_featureCounts_density.pdf")




##############################################################################
##  Category-wise scatter plot
##############################################################################

#---------------
# CAGE junction reads raw counts processing
#----------------
#prom_anno_rds <- anno
anno <- readRDS("/mnt/citadel2/research/shared/AltPromoterFlow/genome/organisms/hg38/Annotation/proActiv_promoter_annotation.rds")
pa <- anno

fc_dir   <- "/mnt/citadel2/research/shared/AltPromoterFlow/CAGE_HEART/featureCounts"
all_fc_files <- list.files(fc_dir, pattern = "_promoter_counts\\.txt$", full.names = TRUE)

healthy_fc_files <- all_fc_files[grepl("Healthy", all_fc_files)]

length(healthy_fc_files)   
basename(healthy_fc_files) 


# Get all promoter level output files
fc_list <- setNames(lapply(healthy_fc_files, function(f) {
  df <- read.table(f, header = TRUE, sep = "\t",
                   comment.char = "#", stringsAsFactors = FALSE,
                   check.names = FALSE)
  sample_name <- sub("_promoter_counts\\.txt$", "", basename(f))
  # Extract promoter id
  prom_id <- sub("_$", "", df[[1]])
  # Extract count value
  counts <- df[[ncol(df)]]
  # Use promoter id as the name of vector
  names(counts) <- prom_id
  return(counts)
}), nm = sub("_promoter_counts\\.txt$", "", basename(healthy_fc_files)))


# Build promoter by sample count matrix
all_ids <- unique(unlist(lapply(fc_list, names)))
mat <- sapply(fc_list, function(x) x[all_ids])
rownames(mat) <- all_ids
colnames(mat) <- names(fc_list)

# Convert to numeric, fill NAs with 0, remove promoters with zero counts
mat <- as.matrix(mat)
mode(mat) <- "numeric"
mat[is.na(mat)] <- 0
#mat <- mat[rowSums(mat) > 0, , drop = FALSE]
cage_counts <- mat


norm_counts <- proActiv:::normalizePromoterReadCounts(cage_counts)

abs_activity <- proActiv:::getAbsolutePromoterActivity(norm_counts, pa)

gene_expression <- proActiv:::getGeneExpression(abs_activity)

activity_mean <- rowMeans(abs_activity[, -c(1, 2)], na.rm = TRUE)


abs_activity$meanActivity <- rowMeans(abs_activity[, -c(1, 2)], na.rm = TRUE)

saveRDS(abs_activity,
        file.path(out_dir,"CAGE_abs_activity.rds"))

activity_threshold <- 0.25

major_minor_promoter_classification <- function(df,
                                                activity_col = "meanActivity",
                                                activity_threshold = 0.25) {

  stopifnot("geneId" %in% colnames(df), activity_col %in% colnames(df))

  gene_split <- split(seq_len(nrow(df)), df$geneId)

  classes <- unlist(lapply(gene_split, function(idx) {
    z <- as.numeric(df[[activity_col]][idx]) 

    if (max(z, na.rm = TRUE) < activity_threshold) {
      z[z < activity_threshold] <- -1
    } else {
      z[z < activity_threshold]                    <- -1          # Inactive
      z[ z == max(z, na.rm = TRUE) ]               <- -2          # Major
      z[ z >= activity_threshold &
         z <  max(z, na.rm = TRUE) ]               <- -3          # Minor
    }

    z[z == -1] <- "Inactive"
    z[z == -2] <- "Major"
    z[z == -3] <- "Minor"

    return(z)
  }))

  return(classes)
}


major_minor_promoter_classification <- function(df,
                                                activity_col = "meanActivity",
                                                activity_threshold = 0.25) {
  stopifnot("geneId" %in% colnames(df), activity_col %in% colnames(df))

  gene_split <- split(seq_len(nrow(df)), df$geneId)

  classes <- unlist(lapply(gene_split, function(idx) {
    z <- as.numeric(df[[activity_col]][idx])
    
    if (all(is.na(z))) {
      return(rep(NA, length(z)))
    }

    labels <- rep("Inactive", length(z))

    z_max <- max(z, na.rm = TRUE)
    if (z_max < activity_threshold) {
      return(labels)
    }
    major_idx <- which.max(z)
    if (!is.na(z[major_idx]) && z[major_idx] >= activity_threshold) {
      labels[major_idx] <- "Major"
    }
    other_idx <- setdiff(seq_along(z), major_idx)
    labels[other_idx[z[other_idx] >= activity_threshold]] <- "Minor"

    return(labels)
  }))

  return(classes)
}


#groups <- split(seq_len(nrow(prom_df)), prom_df$geneId)
#abs_activity$promoterClass <- unlist(lapply(groups, major_minor_call))[ order(unlist(groups)) ]
abs_activity$promoterClass <- major_minor_promoter_classification(abs_activity)
#abs_activity$promoterClass <- active_promoter_category(abs_activity)


anno <- readRDS("/mnt/citadel2/research/shared/AltPromoterFlow/genome/organisms/hg38/Annotation/proActiv_promoter_annotation.rds")
pc <- proActiv:::promoterCoordinates(anno)
internal_map <- setNames(mcols(pc)$internalPromoter, mcols(pc)$promoterId)


abs_activity$internalPromoter <- internal_map[as.character(abs_activity$promoterId)]


abs_activity$promoterClass[abs_activity$internalPromoter == TRUE] <- NA

table(abs_activity$promoterClass, useNA = "always")

prom_df <- abs_activity




#--------------------
# Fit genes into promoter category
#--------------------
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

# Use the promoters grouped by gene to categorize gene type
gene_cat <- sapply(split(prom_df$promoterClass, prom_df$geneId),
                   active_promoter_classification)
table(gene_cat)

#--------------------
# Save results
#--------------------

## Major/Minor promoterId
saveRDS(prom_df$promoterId[prom_df$promoterClass=="Major"],
        file.path(out_dir,"CAGE_major_promoterId.rds"))
saveRDS(prom_df$promoterId[prom_df$promoterClass=="Minor"],
        file.path(out_dir,"CAGE_minor_promoterId.rds"))
saveRDS(prom_df$promoterId[prom_df$promoterClass=="Inactive"],
        file.path(out_dir,"CAGE_inactive_promoterId.rds"))

# Group promoters by gene category and write promoterId list for every category
cat_list <- split(prom_df$promoterId, gene_cat[prom_df$geneId])
invisible(lapply(names(cat_list), \(cc){
  saveRDS(cat_list[[cc]],
          file.path(out_dir, paste0("CAGE_promoterId_", cc, ".rds")))
}))

#Extract all gene ids from each category
invisible(lapply(unique(gene_cat), \(cc){
  saveRDS(names(gene_cat)[gene_cat==cc],
          file.path(out_dir, paste0("CAGE_geneId_", cc, ".rds")))
}))

message("Saved results for category-wise scatter plot at:", out_dir)




#-------------
# Scatter plot of promoter in gene category
#---------------
mm <- readRDS("/mnt/citadel2/research/shared/AltPromoterFlow/RNA_HEART/comparison3/CAGE_promoterId_Multipromoter.Multiactive.rds")
ms <- readRDS("/mnt/citadel2/research/shared/AltPromoterFlow/RNA_HEART/comparison3/CAGE_promoterId_Multipromoter.Singleactive.rds")
ss <- readRDS("/mnt/citadel2/research/shared/AltPromoterFlow/RNA_HEART/comparison3/CAGE_promoterId_Singlepromoter.Singleactive.rds")
sapply(list(mm=mm, ms=ms, ss=ss), length)
sapply(list(mm=mm, ms=ms, ss=ss), function(x) length(unique(x)))


methods <- list(
  Salmon   = salmon_counts,
  proActiv = proactiv_counts,
  DEXSeq   = dexseq_counts
)


promoter_sets <- list(
  "Multipromoter.Multiactive"     = unique(mm),
  "Multipromoter.Singleactive"    = unique(ms),
  "Singlepromoter.Singleactive"   = unique(ss)
)


for (method_name in names(methods)) {
  message("Processing method: ", method_name)
  method_counts <- methods[[method_name]]

  df <- build_comp_df(cage_counts, method_counts, method_name)
  # Plot for every subset
  for (cat in names(promoter_sets)) {
    cat_ids <- promoter_sets[[cat]]
    df_sub  <- subset(df, promoterId %in% cat_ids)
    df_sub <- subset(df_sub, cage_norm > 0 & other_norm > 0)

    title <- paste0("Scatter - ", method_name, " [", cat, "]")
    filename <- paste0("scatter_", method_name, "_", cat, ".pdf")

    plot_density_scatter(df_sub, title = title, outfile = filename)
  }
}

#promoter_sets <- list(
#  "Multipromoter.Multiactive"     = mm,
#  "Multipromoter.Singleactive"    = ms,
#  "Singlepromoter.Singleactive"   = ss
#)

for (cat in names(promoter_sets)) {
  message("→ Processing category: ", cat)
  cat_ids <- promoter_sets[[cat]]

  # Extract all data of all methods on promoter in each category
  dfs <- lapply(names(methods), function(method_name) {
    method_counts <- methods[[method_name]]
    df <- build_comp_df(cage_counts, method_counts, method_name)
    df <- subset(df, promoterId %in% cat_ids & cage_norm > 0 & other_norm > 0)
    return(df)
  })

  all_df <- do.call(rbind, dfs)

  if (nrow(all_df) == 0) {
    message("  No valid rows for category: ", cat)
    next
  }

  outfile <- paste0("scatter_", cat, ".pdf")
  title   <- paste0("Scatter - ", cat)

  plot_density_scatter(all_df, title = title, outfile = outfile)
}





#-------------
# Scatter plot of promoter in activity category
#---------------

major <- readRDS("/mnt/citadel2/research/shared/AltPromoterFlow/RNA_HEART/comparison3/CAGE_major_promoterId.rds")
minor <- readRDS("/mnt/citadel2/research/shared/AltPromoterFlow/RNA_HEART/comparison3/CAGE_minor_promoterId.rds")
inactive <- readRDS("/mnt/citadel2/research/shared/AltPromoterFlow/RNA_HEART/comparison3/CAGE_inactive_promoterId.rds")
intersect <- readRDS("/mnt/citadel2/research/shared/AltPromoterFlow/RNA_HEART/comparison/intersect_major_minor.rds")

methods <- list(
  Salmon   = salmon_counts,
  proActiv = proactiv_counts,
  DEXSeq   = dexseq_counts
)


promoter_sets <- list(
  "Major_promoters"     = major,
  "Minor_promoters"    = minor,
#  "Inactive_promoters" = inactive,
  "Intersect" = intersect
)


for (method_name in names(methods)) {
  message("Processing method: ", method_name)
  method_counts <- methods[[method_name]]

  df <- build_comp_df(cage_counts, method_counts, method_name)
  # Plot for every subset
  for (cat in names(promoter_sets)) {
    cat_ids <- promoter_sets[[cat]]
    df_sub  <- subset(df, promoterId %in% cat_ids)
    before_n <- nrow(df_sub)
    df_sub   <- subset(df_sub, cage_norm > 0 & other_norm > 0)
    after_n  <- nrow(df_sub)
    message(sprintf("  %-6s | %-30s : kept %d / %d",
                    method_name, cat, after_n, before_n))
    title <- paste0("Scatter - ", method_name, " [", cat, "]")
    filename <- paste0("scatter_", method_name, "_", cat, ".pdf")

    plot_density_scatter(df_sub, title = title, outfile = filename)
  }
}

for (cat in names(promoter_sets)) {
  message("→ Processing category: ", cat)
  cat_ids <- promoter_sets[[cat]]

  # Extract all data of all methods on promoter in each category
  dfs <- lapply(names(methods), function(method_name) {
    method_counts <- methods[[method_name]]
    df <- build_comp_df(cage_counts, method_counts, method_name)
    df <- subset(df, promoterId %in% cat_ids & cage_norm > 0 & other_norm > 0)
    return(df)
  })

  all_df <- do.call(rbind, dfs)

  if (nrow(all_df) == 0) {
    message("  No valid rows for category: ", cat)
    next
  }

  outfile <- paste0("scatter_", cat, ".pdf")
  title   <- paste0("Scatter - ", cat)

  plot_density_scatter(all_df, title = title, outfile = outfile)
}






# ---------------------------
# Barplot: Single stacked bar per method (Up/Down segments)
# ---------------------------
cage_counts     <- readRDS("/mnt/citadel2/research/shared/AltPromoterFlow/CAGE_HEART/salmon/merged/promoter_counts.rds")

library(dplyr)
library(ggplot2)

count_up_down_by_method <- function(data_list, logfc_col, pval_col, pval_thresh = 0.05) {
  do.call(rbind, lapply(names(data_list), function(method) {
    df <- data_list[[method]]
    df <- df[is.finite(df[[logfc_col]]) & is.finite(df[[pval_col]]), ]
    df <- df[df[[pval_col]] < pval_thresh, ]
    up <- sum(df[[logfc_col]] > 0)
    down <- sum(df[[logfc_col]] < 0)
    tibble(
      method = method,
      direction = c("Up", "Down"),
      count_raw = c(up, down),
      count = c(up, -down),
      vjust = c(-0.3, 1.3)
    )
  }))
}

plot_stacked_updown_bar <- function(count_df, title, outfile) {
  count_df$method <- factor(count_df$method, levels = c("CAGE", "Salmon", "proActiv", "DEXSeq"))
  count_df$direction <- factor(count_df$direction, levels = c("Up", "Down"))

  p <- ggplot(count_df, aes(x = method, y = count, fill = direction)) +
    geom_bar(stat = "identity", width = 0.6, color = "black") +
    geom_hline(yintercept = 0) +
    geom_text(aes(label = abs(count_raw), vjust = vjust), size = 3) +
    scale_fill_manual(values = c("Up" = "#fbb4ae", "Down" = "#b3cde3")) +
    labs(x = NULL, y = "Significant count", title = title) +
    theme_bw() +
    theme(legend.position = "right",
          plot.title = element_text(hjust = 0.5, face = "bold"))

  ggsave(file.path(out_dir, outfile), p, width = 6, height = 4.5)
}

pdfs <- list.files(out_dir, pattern = "^MAplot_.*\\.pdf$", full.names = TRUE)
system(paste("pdfunite", paste(shQuote(pdfs), collapse = " "), 
             shQuote(file.path(out_dir, "ALL_MAplots_combined.pdf"))))

# ---------------------------
# Prepare data lists
# ---------------------------
gene_list <- list(
  CAGE     = CAGE_gene,
  Salmon   = salmon_gene,
  proActiv = proactiv_gene,
  DEXSeq   = dexseq_gene
)

# ---------------------------
# Count Up/Down
# ---------------------------
gene_fdr_df      <- count_up_down_by_method(gene_list,      "logFC",         "FDR",          0.05)

# ---------------------------
# Plot stacked barplots
# ---------------------------
plot_stacked_updown_bar(gene_fdr_df,      "Promoter Up/Down Regulate Count (FDR < 0.05)",          "promoter_count_FDR.pdf")

# ---------------------------
# Pie chart helper
# ---------------------------
get_overlap_pie <- function(all_ids, method_ids, method_name) {
  data <- data.frame(
    category = c("Overlap", "CAGE-only"),
    count = c(length(intersect(all_ids, method_ids)),
              length(setdiff(all_ids, method_ids)))
  )
  data$percent <- paste0(round(100 * data$count / sum(data$count), 1), "%")

  ggplot(data, aes(x = "", y = count, fill = category)) +
    geom_bar(stat = "identity", width = 1, color = "white") +
    coord_polar(theta = "y") +
    scale_fill_manual(values = c("Overlap" = "#fec44f", "CAGE-only" = "#fff7bc")) +
    labs(title = paste0(method_name, " vs CAGE"), x = NULL, y = NULL) +
    theme_void() +
    theme(plot.title = element_text(hjust = 0.5, size = 10),
          legend.title = element_blank(),
          legend.text = element_text(size = 8)) +
    geom_text(aes(label = percent), position = position_stack(vjust = 0.5), size = 3.5)
}
# ---------------------------
# Gene pie charts (FDR)
# ---------------------------
#cage_genes <- CAGE_gene$geneId[CAGE_gene$FDR < 0.05]
#salmon_gene_ids <- salmon_gene$geneId[salmon_gene$FDR < 0.05]
#dexseq_gene_ids <- dexseq_gene$geneId[dexseq_gene$FDR < 0.05]
#proactiv_gene_ids <- proactiv_gene$geneId[proactiv_gene$FDR < 0.05]
cage_promoter_ids <- CAGE_gene$promoterId[CAGE_gene$FDR < 0.05]
salmon_promoter_ids <- salmon_gene$promoterId[salmon_gene$FDR < 0.05]
dexseq_promoter_ids <- dexseq_gene$promoterId[dexseq_gene$FDR < 0.05]
proactiv_promoter_ids <- proactiv_gene$promoterId[proactiv_gene$FDR < 0.05]

salmon_gene_pie   <- get_overlap_pie(cage_promoter_ids, salmon_promoter_ids,   "Salmon")
dexseq_gene_pie   <- get_overlap_pie(cage_promoter_ids, dexseq_promoter_ids,   "DEXSeq")
proactiv_gene_pie <- get_overlap_pie(cage_promoter_ids, proactiv_promoter_ids, "proActiv")

# ---------------------------
# Save pie plots
# ---------------------------

pdf(file.path(out_dir, "combined_promoter_overlap_pie.pdf"), width = 10, height = 4)
grid.arrange(salmon_gene_pie, dexseq_gene_pie, proactiv_gene_pie,
             ncol = 3,
             top = textGrob("CAGE Overlap of DE Genes (FDR < 0.05)", gp = gpar(fontsize = 12, fontface = "bold")))
dev.off()



# ---------------------------
# Venn diagram for gene overlap (FDR)
# ---------------------------
#cage_gene_ids <- cage_genes
venn_fdr_gene <- list(
  CAGE = cage_promoter_ids,
  Salmon = salmon_promoter_ids,
  DEXSeq = dexseq_promoter_ids,
  proActiv = proactiv_promoter_ids
)
#venn_fdr_gene <- list(
#  CAGE = cage_gene_ids,
#  Salmon = salmon_gene$geneId[salmon_gene$FDR < 0.05],
#  DEXSeq = dexseq_gene$geneId[dexseq_gene$FDR < 0.05],
#  proActiv = proactiv_gene$geneId[proactiv_gene$FDR < 0.05]
#)


# ---------------------------
# Individual method vs CAGE (Gene FDR)
# ---------------------------
venn_gene_salmon   <- ggvenn(list(Salmon = salmon_promoter_ids,
                                  CAGE = cage_promoter_ids),
                             fill_color = c("#b3cde3", "#fbb4ae"),
                             stroke_size = 0.5, text_size = 4, show_percentage = FALSE, set_name_size = 4) +
                      labs(title = "Salmon vs CAGE (DEG)") +
                      theme(plot.title = element_text(hjust = 0.5, size = 10))

venn_gene_dexseq   <- ggvenn(list(DEXSeq = dexseq_promoter_ids,
                                  CAGE = cage_promoter_ids),
                             fill_color = c("#b3cde3", "#fbb4ae"),
                             stroke_size = 0.5, text_size = 4, show_percentage = FALSE, set_name_size = 4) +
                      labs(title = "DEXSeq vs CAGE (DEG)") +
                      theme(plot.title = element_text(hjust = 0.5, size = 10))

venn_gene_proactiv <- ggvenn(list(proActiv = proactiv_promoter_ids,
                                  CAGE = cage_promoter_ids),
                             fill_color = c("#b3cde3", "#fbb4ae"),
                             stroke_size = 0.5, text_size = 4, show_percentage = FALSE, set_name_size = 4) +
                      labs(title = "proActiv vs CAGE (DEG)") +
                      theme(plot.title = element_text(hjust = 0.5, size = 10))

# Combine gene Venns
pdf(file.path(out_dir, "combined_promoter_overlap_venn.pdf"), width = 10, height = 4)
grid.arrange(
  venn_gene_salmon, venn_gene_dexseq, venn_gene_proactiv,
  ncol = 3,
  top = textGrob("Venn Diagram: Promoter Overlap (FDR < 0.05)", gp = gpar(fontsize = 12, fontface = "bold"))
)
dev.off()


# ---- UPREG ----
cage_up_promoter     <- CAGE_gene$promoterId[ CAGE_gene$FDR < 0.05 & CAGE_gene$logFC > 0 ]
salmon_up_promoter   <- salmon_gene$promoterId[ salmon_gene$FDR < 0.05 & salmon_gene$logFC > 0 ]
dexseq_up_promoter   <- dexseq_gene$promoterId[ dexseq_gene$FDR < 0.05 & dexseq_gene$logFC > 0 ]
proactiv_up_promoter <- proactiv_gene$promoterId[ proactiv_gene$FDR < 0.05 & proactiv_gene$logFC > 0 ]

# ---- DOWNREG ----
cage_down_promoter     <- CAGE_gene$promoterId[ CAGE_gene$FDR < 0.05 & CAGE_gene$logFC < 0 ]
salmon_down_promoter   <- salmon_gene$promoterId[ salmon_gene$FDR < 0.05 & salmon_gene$logFC < 0 ]
dexseq_down_promoter   <- dexseq_gene$promoterId[ dexseq_gene$FDR < 0.05 & dexseq_gene$logFC < 0 ]
proactiv_down_promoter <- proactiv_gene$promoterId[ proactiv_gene$FDR < 0.05 & proactiv_gene$logFC < 0 ]

library(ggvenn)

# ---- UPREG overlap with CAGE ----
venn_up_salmon <- ggvenn(list(Salmon = salmon_up_promoter, CAGE = cage_up_promoter),
                         fill_color = c("#b3cde3", "#fbb4ae"), stroke_size = 0.5,
                         text_size = 4, show_percentage = FALSE, set_name_size = 4) +
                  labs(title = "Upregulated Promoters (Salmon vs CAGE)") +
                  theme(plot.title = element_text(hjust = 0.5, size = 10))

venn_up_dexseq <- ggvenn(list(DEXSeq = dexseq_up_promoter, CAGE = cage_up_promoter),
                         fill_color = c("#b3cde3", "#fbb4ae"), stroke_size = 0.5,
                         text_size = 4, show_percentage = FALSE, set_name_size = 4) +
                  labs(title = "Upregulated Promoters (DEXSeq vs CAGE)") +
                  theme(plot.title = element_text(hjust = 0.5, size = 10))

venn_up_proactiv <- ggvenn(list(proActiv = proactiv_up_promoter, CAGE = cage_up_promoter),
                           fill_color = c("#b3cde3", "#fbb4ae"), stroke_size = 0.5,
                           text_size = 4, show_percentage = FALSE, set_name_size = 4) +
                    labs(title = "Upregulated Promoters (proActiv vs CAGE)") +
                    theme(plot.title = element_text(hjust = 0.5, size = 10))


# ---- DOWNREG overlap with CAGE ----
venn_down_salmon <- ggvenn(list(Salmon = salmon_down_promoter, CAGE = cage_down_promoter),
                           fill_color = c("#b3cde3", "#fbb4ae"), stroke_size = 0.5,
                           text_size = 4, show_percentage = FALSE, set_name_size = 4) +
                    labs(title = "Downregulated Promoters (Salmon vs CAGE)") +
                    theme(plot.title = element_text(hjust = 0.5, size = 10))

venn_down_dexseq <- ggvenn(list(DEXSeq = dexseq_down_promoter, CAGE = cage_down_promoter),
                           fill_color = c("#b3cde3", "#fbb4ae"), stroke_size = 0.5,
                           text_size = 4, show_percentage = FALSE, set_name_size = 4) +
                    labs(title = "Downregulated Promoters (DEXSeq vs CAGE)") +
                    theme(plot.title = element_text(hjust = 0.5, size = 10))

venn_down_proactiv <- ggvenn(list(proActiv = proactiv_down_promoter, CAGE = cage_down_promoter),
                             fill_color = c("#b3cde3", "#fbb4ae"), stroke_size = 0.5,
                             text_size = 4, show_percentage = FALSE, set_name_size = 4) +
                      labs(title = "Downregulated Promoters (proActiv vs CAGE)") +
                      theme(plot.title = element_text(hjust = 0.5, size = 10))

pdf(file.path(out_dir, "venn_up_down_promoters.pdf"), width = 10, height = 8)
grid.arrange(
  venn_up_salmon, venn_up_dexseq, venn_up_proactiv,
  venn_down_salmon, venn_down_dexseq, venn_down_proactiv,
  ncol = 3,
  top = textGrob("Venn Diagram: Up/Downregulated Promoter Overlap (FDR < 0.05)",
                 gp = gpar(fontsize = 12, fontface = "bold"))
)
dev.off()
