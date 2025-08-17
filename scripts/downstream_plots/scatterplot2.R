
#/genome/.snakemake_conda/671ce29470b1e8f250516402ed04d982_

# ---------------------------
# Plot: feature counts and no intron
# ---------------------------

library(DESeq2)
library(ggplot2)
library(MASS)
library(viridis)


cell_lines <- c("GM12878", "H1hESC", "K562")

cell_lines <- c("Healthy", "Failure")
out_dir <- "/mnt/citadel2/research/syidan/Projects/SnakeAltPromoterResult"


clean_mat <- function(mat) {
  mat[is.na(mat)] <- 0
  mat
}

joint_norm <- function(cage_mat, other_mat, tag) {
  common <- intersect(rownames(cage_mat), rownames(other_mat))
  if (length(common) == 0L)
      stop(sprintf("No common promoters for %s vs CAGE", tag))
  cage_mat  <- round(cage_mat[common, , drop = FALSE])
  other_mat <- round(other_mat[common, , drop = FALSE])
  all_mat   <- cbind(cage_mat, other_mat)
  all_mat[is.na(all_mat)] <- 0
  colData <- data.frame(
    tech = factor(c(rep("CAGE", ncol(cage_mat)), rep(tag, ncol(other_mat)))),
    row.names = colnames(all_mat)
  )
  dds <- DESeqDataSetFromMatrix(all_mat, colData = colData, design = ~1)
  dds <- estimateSizeFactors(dds)
  cat("\n[Joint normalisation]  tag =", tag, "\n")
  print(sizeFactors(dds))
  norm_all <- counts(dds, normalized = TRUE)
  list(
    cage_norm  = norm_all[ , colnames(cage_mat),  drop = FALSE],
    other_norm = norm_all[ , colnames(other_mat), drop = FALSE]
  )
}

make_df_one_method_one_cl <- function(cage_mat, other_mat, cl, method_name) {
  sel_cols <- function(mat) {
    if (any(grepl("^CAGE_", colnames(mat)))) {
      mat[, grepl(paste0("^CAGE_", toupper(cl)), colnames(mat)), drop = FALSE]
    } else {
      mat[, grepl(paste0("^RNA_", toupper(cl)), colnames(mat)), drop = FALSE]
    }
  }

  cage_cl   <- sel_cols(cage_mat)
  other_cl  <- sel_cols(other_mat)
  if (ncol(other_cl) == 0 | ncol(cage_cl) == 0)
      return(NULL)
  nj <- joint_norm(cage_cl, other_cl, tag = paste0(method_name, "_", cl))
  cage_norm  <- rowMeans(nj$cage_norm)
  other_norm <- rowMeans(nj$other_norm)
  data.frame(
    promoterId = names(cage_norm),
    cage_norm  = cage_norm,
    other_norm = other_norm,
    log2FC     = log2((other_norm + 1) / (cage_norm + 1)),
    method     = method_name,
    cell_line  = cl,
    stringsAsFactors = FALSE
  )
}



plot_density_scatter <- function(df, title = "", outfile, use_density = TRUE, n_grid = 200) {
  if (nrow(df) == 0) {
    message("No rows for ", outfile); return(invisible(NULL))
  }

  # Use log10 to calculate coordinates
  df$logC <- log10(df$cage_norm + 1)
  df$logO <- log10(df$other_norm + 1)
  df$density <- NA_real_

  # Use log 10 to calculate density
  if (use_density) {
    keep <- is.finite(df$logC) & is.finite(df$logO)
    dens <- with(df[keep, ], MASS::kde2d(logC, logO, n = n_grid))
    ix <- findInterval(df$logC[keep], dens$x)
    iy <- findInterval(df$logO[keep], dens$y)
    df$density[keep] <- log10(dens$z[cbind(ix, iy)] + 1e-8)
  }

  # Correlation using original value
  cor_labels <- do.call(rbind, lapply(split(df, df$method), function(sub) {
    cp <- cor.test(sub$cage_norm, sub$other_norm, method = "pearson")
    cs <- cor.test(sub$cage_norm, sub$other_norm, method = "spearman")
    pval <- cs$p.value
    pval_str <- ifelse(is.na(pval), "NA", ifelse(pval < 1e-16, "< 1e-16", sprintf("%.1e", pval)))
    data.frame(
      method = unique(sub$method),
      x = min(sub$logC, na.rm = TRUE) + 0.2,
      y = max(sub$logO, na.rm = TRUE) - 0.2,
      #label = sprintf("Pearson r = %.3f\nSpearman ρ = %.3f\n(p = %s)", cp$estimate, cs$estimate, pval_str)
      label = sprintf("Pearson r = %.3f\n(p = %s)", cp$estimate,pval_str)
      #label = sprintf("Spearman = %.3f \n (p = %s)", cs$estimate, pval_str)
    )
  }))

  # Coordinate range
  lims <- range(c(df$logC, df$logO), finite = TRUE)

  # Plot
  p <- ggplot(df, aes(x = logC, y = logO)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey50") +
    scale_x_continuous(limits = lims) +
    scale_y_continuous(limits = lims) +
    facet_wrap(~method, nrow = 1) +
    geom_text(data = cor_labels,
              aes(x = x, y = y, label = label),
              inherit.aes = FALSE,
              hjust = 0, vjust = 1, size = 3.2) +
    labs(x = "CAGE (log10 normalized)", y = "Other method (log10 normalized)") +
    theme_bw(base_size = 11) +
    theme(plot.title = element_blank())

  if (use_density) {
    p <- p + geom_point(aes(color = density), size = 1, alpha = 0.85) +
      scale_color_viridis_c(option = "A", name = "Density")
  } else {
    p <- p + geom_point(alpha = 0.5, size = 0.7)
  }

  ggsave(file.path(out_dir, outfile), p, width = 5, height = 4)
}


#--------------------------------------------------------
#Plot scatterplot of intronless promoter counts of Salmon and DEXseq vs CAGE counts
promoter_anno <- readRDS("/mnt/citadel2/research/shared/AltPromoterFlow/genome/organisms/hg38/Annotation/proActiv_promoter_annotation.rds")
coord <- promoter_anno@promoterCoordinates
valid_ids <- as.character(coord$promoterId[!coord$internalPromoter %in% TRUE &
                                            lengths(coord$intronId) == 0])

cage_counts     <- clean_mat(readRDS("/mnt/citadel2/research/shared/AltPromoterFlow/Heart/cage/merge/raw_promoter_counts.rds"))
salmon_counts   <- clean_mat(readRDS("/mnt/citadel2/research/shared/AltPromoterFlow/Heart/salmon/merge/raw_promoter_counts.rds"))
dexseq_counts   <- clean_mat(readRDS("/mnt/citadel2/research/shared/AltPromoterFlow/Heart/dexseq/merge/raw_promoter_counts.rds"))
proactiv_counts <- clean_mat(readRDS("/mnt/citadel2/research/shared/AltPromoterFlow/Heart/proactiv/merge/raw_promoter_counts.rds"))

# Keep no intron and no internal
cage_counts     <- cage_counts[rownames(cage_counts)     %in% valid_ids, , drop = FALSE]
salmon_counts   <- salmon_counts[rownames(salmon_counts) %in% valid_ids, , drop = FALSE]
dexseq_counts   <- dexseq_counts[rownames(dexseq_counts) %in% valid_ids, , drop = FALSE]
proactiv_counts <- proactiv_counts[rownames(proactiv_counts) %in% valid_ids, , drop = FALSE]

methods <- list(
  Salmon   = salmon_counts,
  DEXSeq   = dexseq_counts
)

all_df <- do.call(
  rbind,
  lapply(cell_lines, function(cl) {
    do.call(rbind, lapply(names(methods), function(m) {
      message("▶ ", cl, "  vs CAGE  —  ", m)
      make_df_one_method_one_cl(cage_counts, methods[[m]], cl, m)
    }))
  })
)

cnt_ok <- with(all_df, cage_norm > 10 & other_norm > 10)
set2   <- subset(all_df, cnt_ok)

for (cl in cell_lines) {
  for (m in names(methods)) {
    df <- subset(set2, cell_line == cl & method == m)
    if (nrow(df) == 0) {
        message("Skipping ", cl, " - ", m, " due to empty data.")
        next
    }
    ttl <- sprintf("%s ‒ %s (no intron, counts ≥10)", cl, m)
    of  <- sprintf("scatter_nointron_%s_%s.pdf", tolower(cl), tolower(m))
    plot_density_scatter(df, ttl, of, use_density = TRUE)
  }
}


#-----------------------------------------------------------------------------------------------------------
#Plot scatterplot of all promoter counts of Salmon and DEXseq vs CAGE counts
cage_counts     <- clean_mat(readRDS("/mnt/citadel2/research/shared/AltPromoterFlow/Heart/cage/merge/raw_promoter_counts.rds"))
salmon_counts   <- clean_mat(readRDS("/mnt/citadel2/research/shared/AltPromoterFlow/Heart/salmon/merge/raw_promoter_counts.rds"))
dexseq_counts   <- clean_mat(readRDS("/mnt/citadel2/research/shared/AltPromoterFlow/Heart/dexseq/merge/raw_promoter_counts.rds"))
proactiv_counts <- clean_mat(readRDS("/mnt/citadel2/research/shared/AltPromoterFlow/Heart/proactiv/merge/raw_promoter_counts.rds"))


methods <- list(
  Salmon   = salmon_counts,
  DEXSeq   = dexseq_counts,
  proActiv = proactiv_counts
)

all_df <- do.call(
  rbind,
  lapply(cell_lines, function(cl) {
    do.call(rbind, lapply(names(methods), function(m) {
      message("▶ ", cl, "  vs CAGE  —  ", m)
      make_df_one_method_one_cl(cage_counts, methods[[m]], cl, m)
    }))
  })
)

cnt_ok <- with(all_df, cage_norm > 10 & other_norm > 10)
set2   <- subset(all_df, cnt_ok)

for (cl in cell_lines) {
  for (m in names(methods)) {
    df <- subset(set2, cell_line == cl & method == m)
    
    if (nrow(df) == 0) next
    ttl <- sprintf("%s ‒ %s (counts ≥10)", cl, m)
    of  <- sprintf("scatter_%s_%s.pdf", tolower(cl), tolower(m))
    plot_density_scatter(df, ttl, of)
  }
}


#---------------------------------------------------------------------------
#Plot scatterplot of all promoter counts of Salmon and DEXseq vs CAGE counts

promoter_anno <- readRDS("/mnt/citadel2/research/shared/AltPromoterFlow/genome/organisms/hg38/Annotation/proActiv_promoter_annotation.rds")
coord <- promoter_anno@promoterCoordinates
valid_ids <- as.character(coord$promoterId[!coord$internalPromoter %in% TRUE &
                                            lengths(coord$intronId) > 0])

cage_counts     <- clean_mat(readRDS("/mnt/citadel2/research/shared/AltPromoterFlow/Heart/cage/merge/raw_promoter_counts.rds"))
salmon_counts   <- clean_mat(readRDS("/mnt/citadel2/research/shared/AltPromoterFlow/Heart/salmon/merge/raw_promoter_counts.rds"))
dexseq_counts   <- clean_mat(readRDS("/mnt/citadel2/research/shared/AltPromoterFlow/Heart/dexseq/merge/raw_promoter_counts.rds"))
proactiv_counts <- clean_mat(readRDS("/mnt/citadel2/research/shared/AltPromoterFlow/Heart/proactiv/merge/raw_promoter_counts.rds"))

# Keep with intron and no internal
cage_counts     <- cage_counts[rownames(cage_counts)     %in% valid_ids, , drop = FALSE]
salmon_counts   <- salmon_counts[rownames(salmon_counts) %in% valid_ids, , drop = FALSE]
dexseq_counts   <- dexseq_counts[rownames(dexseq_counts) %in% valid_ids, , drop = FALSE]
proactiv_counts <- proactiv_counts[rownames(proactiv_counts) %in% valid_ids, , drop = FALSE]

methods <- list(
  Salmon   = salmon_counts,
  DEXSeq   = dexseq_counts,
  proActiv = proactiv_counts
)

all_df <- do.call(
  rbind,
  lapply(cell_lines, function(cl) {
    do.call(rbind, lapply(names(methods), function(m) {
      message("▶ ", cl, "  vs CAGE  —  ", m)
      make_df_one_method_one_cl(cage_counts, methods[[m]], cl, m)
    }))
  })
)

# Filter norm > 10
cnt_ok <- with(all_df, cage_norm > 10 & other_norm > 10)
set2   <- subset(all_df, cnt_ok)

for (cl in cell_lines) {
  for (m in names(methods)) {
    df <- subset(set2, cell_line == cl & method == m)
    if (nrow(df) == 0) {
        message("Skipping ", cl, " - ", m, " due to empty data.")
        next
    }
    ttl <- sprintf("%s ‒ %s (with intron, counts ≥10)", cl, m)
    of  <- sprintf("scatter_withintron_%s_%s.pdf", tolower(cl), tolower(m))
    plot_density_scatter(df, ttl, of, use_density = TRUE)
  }
}





#--------------------------------------------------------------------------------------
#Plot scatterplot of all log2FC of Proactiv, Salmon and DEXseq vs CAGE counts
library(ggplot2)
library(MASS)
library(viridis)
library(tidyr)

dexseq <- readRDS("/mnt/citadel2/research/shared/AltPromoterFlow/test_for_sample_sheet/dexseq/differential/comparisons_/K562_vs_H1-hESC/Promoter_differential_activity_FDR0_05.rds")
salmon <- readRDS("/mnt/citadel2/research/shared/AltPromoterFlow/test_for_sample_sheet/salmon/differential/comparisons_/K562_vs_H1-hESC/Promoter_differential_activity_FDR0_05.rds")
proactiv <- readRDS("/mnt/citadel2/research/shared/AltPromoterFlow/test_for_sample_sheet/proactiv/differential/comparisons_/K562_vs_H1-hESC/Promoter_differential_activity_FDR0_05.rds")
cage <- readRDS("/mnt/citadel2/research/shared/AltPromoterFlow/test_for_sample_sheet/cage/differential/comparisons_/K562_vs_H1-hESC/Promoter_differential_activity_FDR0_05.rds")
dexseq_all <- readRDS("/mnt/citadel2/research/shared/AltPromoterFlow/test_for_sample_sheet/dexseq/differential/comparisons_/K562_vs_H1-hESC/Promoter_differential_activity.rds")
salmon_all <- readRDS("/mnt/citadel2/research/shared/AltPromoterFlow/test_for_sample_sheet/salmon/differential/comparisons_/K562_vs_H1-hESC/Promoter_differential_activity.rds")
proactiv_all <- readRDS("/mnt/citadel2/research/shared/AltPromoterFlow/test_for_sample_sheet/proactiv/differential/comparisons_/K562_vs_H1-hESC/Promoter_differential_activity.rds")
cage_all <- readRDS("/mnt/citadel2/research/shared/AltPromoterFlow/test_for_sample_sheet/cage/differential/comparisons_/K562_vs_H1-hESC/Promoter_differential_activity.rds")
out_dir <- "/mnt/citadel2/research/syidan/Projects/SnakeAltPromoterResult/K562_vs_H1-hESC"

dexseq <- readRDS("/mnt/citadel2/research/shared/AltPromoterFlow/test_for_sample_sheet/dexseq/differential/comparisons_/K562_vs_GM12878/Promoter_differential_activity_FDR0_05.rds")
salmon <- readRDS("/mnt/citadel2/research/shared/AltPromoterFlow/test_for_sample_sheet/salmon/differential/comparisons_/K562_vs_GM12878/Promoter_differential_activity_FDR0_05.rds")
proactiv <- readRDS("/mnt/citadel2/research/shared/AltPromoterFlow/test_for_sample_sheet/proactiv/differential/comparisons_/K562_vs_GM12878/Promoter_differential_activity_FDR0_05.rds")
cage <- readRDS("/mnt/citadel2/research/shared/AltPromoterFlow/test_for_sample_sheet/cage/differential/comparisons_/K562_vs_GM12878/Promoter_differential_activity_FDR0_05.rds")
dexseq_all <- readRDS("/mnt/citadel2/research/shared/AltPromoterFlow/test_for_sample_sheet/dexseq/differential/comparisons_/K562_vs_GM12878/Promoter_differential_activity.rds")
salmon_all <- readRDS("/mnt/citadel2/research/shared/AltPromoterFlow/test_for_sample_sheet/salmon/differential/comparisons_/K562_vs_GM12878/Promoter_differential_activity.rds")
proactiv_all <- readRDS("/mnt/citadel2/research/shared/AltPromoterFlow/test_for_sample_sheet/proactiv/differential/comparisons_/K562_vs_GM12878/Promoter_differential_activity.rds")
cage_all <- readRDS("/mnt/citadel2/research/shared/AltPromoterFlow/test_for_sample_sheet/cage/differential/comparisons_/K562_vs_GM12878/Promoter_differential_activity.rds")
out_dir <- "/mnt/citadel2/research/syidan/Projects/SnakeAltPromoterResult/K562_vs_GM12878"

dexseq <- readRDS("/mnt/citadel2/research/shared/AltPromoterFlow/new_test_for_sample_sheet/dexseq/differential/comparisons_/Heart_Failure_vs_Healthy/Promoter_differential_activity_FDR0_05.rds")
salmon <- readRDS("/mnt/citadel2/research/shared/AltPromoterFlow/new_test_for_sample_sheet/salmon/differential/comparisons_/Heart_Failure_vs_Healthy/Promoter_differential_activity_FDR0_05.rds")
proactiv <- readRDS("/mnt/citadel2/research/shared/AltPromoterFlow/new_test_for_sample_sheet/proactiv/differential/comparisons_/Heart_Failure_vs_Healthy/Promoter_differential_activity_FDR0_05.rds")
cage <- readRDS("/mnt/citadel2/research/shared/AltPromoterFlow/new_test_for_sample_sheet/cage/differential/comparisons_/Heart_Failure_vs_Healthy/Promoter_differential_activity_FDR0_05.rds")
dexseq_all <- readRDS("/mnt/citadel2/research/shared/AltPromoterFlow/new_test_for_sample_sheet/dexseq/differential/comparisons_/Heart_Failure_vs_Healthy/Promoter_differential_activity.rds")
salmon_all <- readRDS("/mnt/citadel2/research/shared/AltPromoterFlow/new_test_for_sample_sheet/salmon/differential/comparisons_/Heart_Failure_vs_Healthy/Promoter_differential_activity.rds")
proactiv_all <- readRDS("/mnt/citadel2/research/shared/AltPromoterFlow/new_test_for_sample_sheet/proactiv/differential/comparisons_/Heart_Failure_vs_Healthy/Promoter_differential_activity.rds")
cage_all <- readRDS("/mnt/citadel2/research/shared/AltPromoterFlow/new_test_for_sample_sheet/cage/differential/comparisons_/Heart_Failure_vs_Healthy/Promoter_differential_activity.rds")
out_dir <- "/mnt/citadel2/research/syidan/Projects/SnakeAltPromoterResult/Heart_Failure_vs_Healthy"

#gene_name = read.delim("/mnt/citadel2/research/syidan/Genomes/GRCh38/release-47-index/annotation/genes.symbol", header = FALSE, comment.char = "#")
#cage$geneSymbol <- gene_name$V2[match(cage$geneId, gene_name$V1)]
#write.table(unique(cage$geneSymbol[!is.na(cage$geneSymbol) & cage$FDR < 0.05 & cage$logFC > 0]),  file.path(out_dir,"cage_upregulated_genes.txt"),  quote = FALSE, row.names = FALSE, col.names = FALSE)
#write.table(unique(cage$geneSymbol[!is.na(cage$geneSymbol) & cage$FDR < 0.05 & cage$logFC < 0]), file.path(out_dir,"cage_downregulated_genes.txt"), quote = FALSE, row.names = FALSE, col.names = FALSE)


length(rownames(dexseq))
length(rownames(cage))
length(rownames(proactiv))
length(rownames(salmon))
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

add_method <- function(df, m){
  df$method <- m
  df[, c("promoterId", "logFC", "method")]
}
logfc_df <- rbind(
  add_method(cage,     "CAGE"),
  add_method(salmon,   "Salmon"),
  add_method(dexseq,   "DEXSeq"),
  add_method(proactiv, "proActiv")
)
logfc_df_all <- rbind(
  add_method(cage_all,     "CAGE"),
  add_method(salmon_all,   "Salmon"),
  add_method(dexseq_all,   "DEXSeq"),
  add_method(proactiv_all, "proActiv")
)


plot_density_scatter_for_logFC <- function(df, title, outfile, use_density = TRUE, n_grid = 200) {
  if (nrow(df) == 0) {
    message("No rows for ", outfile); return(invisible(NULL))
  }

  df$logC <- df$cage_log
  df$logO <- df$other_log
  df$density <- NA_real_

  if (use_density) {
    keep <- is.finite(df$logC) & is.finite(df$logO)
    dens <- with(df[keep, ], MASS::kde2d(logC, logO, n = n_grid))
    ix <- findInterval(df$logC[keep], dens$x)
    iy <- findInterval(df$logO[keep], dens$y)
    #df$density[keep] <- log10(dens$z[cbind(ix, iy)] + 1e-8)
    df$density[keep] <- dens$z[cbind(ix, iy)]
  }

  cor_labels <- do.call(rbind, lapply(split(df, df$method), function(sub) {
    x <- sub$logC; y <- sub$logO
    cp <- cor.test(x, y, method = "pearson")
    cs <- cor.test(x, y, method = "spearman")
    data.frame(
      method = unique(sub$method),
      x = min(sub$logC, na.rm = TRUE),
      y = max(sub$logO, na.rm = TRUE),
      #x = 1,
      #y = 22.5,
      label = sprintf("Pearson r = %.3f (p %.1g)\nSpearman = %.3f (p %.1g)",
                      cp$estimate, cp$p.value, cs$estimate, cs$p.value)
    )
  }))

  lims <- range(c(df$logC, df$logO), finite = TRUE)

  p <- ggplot(df, aes(x = logC, y = logO)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey50") +
    scale_x_continuous(limits = lims, breaks = scales::pretty_breaks()) +
    scale_y_continuous(limits = lims, breaks = scales::pretty_breaks()) +
    facet_wrap(~method, nrow = 1) +
    geom_text(data = cor_labels, aes(x = x, y = y, label = label),
              hjust = 0, vjust = 1, size = 3.2, inherit.aes = FALSE) +
    labs(x = "CAGE  (|log2FC|)", y = "Other method  (|log2FC|)") +
    theme_bw(base_size = 11) +
    theme(plot.title = element_text(hjust = 0.5), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.background = element_rect(fill = "white"))

  if (use_density) {
    p <- p + geom_point(aes(color = density), size = 1, alpha = 0.85) +
      scale_color_viridis_c(option = "A", name = "Density")
  } else {
    p <- p + geom_point(alpha = 0.5, size = 0.7)
  }

  ggsave(file.path(out_dir, outfile), p, width = 5, height = 4)
}

# plot boxplot and violin plot for logFC
plot_box_violin_logFC <- function(df, title, outfile, direction = NULL) {
  if (nrow(df) == 0) {
    message("No rows for ", outfile); return(invisible(NULL))
  }

  # Filter based on direction parameter
  if (!is.null(direction)) {
    if (direction == "up") {
      df <- df[df$logFC > 0, ]
    } else if (direction == "down") {
      df <- df[df$logFC < 0, ]
    } else {
      stop("Invalid direction parameter. Use 'up', 'down', or NULL.")
    }
  }

  if (nrow(df) == 0) {
    message("No rows after filtering for ", outfile); return(invisible(NULL))
  }

  p <- ggplot(df, aes(x = method, y = logFC, fill = method)) +
    geom_violin(alpha = 0.4, trim = TRUE) +
    geom_boxplot(width = 0.2, fill = "white", outlier.size = 1, outlier.alpha = 0.5) +
    labs(x = "Method", y = "log2FC", title = title) +
    scale_fill_manual(values = c("CAGE" = "#1b9e77", "Salmon" = "#d95f02", "DEXSeq" = "#7570b3", "proActiv" = "#e7298a")) +
    theme_bw(base_size = 11) +
    theme(
      plot.title = element_text(hjust = 0.5),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text = element_text(color = "black"),  # Ensure axis text is black
      legend.position = "none"  # Remove legend since x-axis labels are sufficient
    )

  ggsave(file.path(out_dir, outfile), p, width = 5, height = 4)
}

logfc_wide <- pivot_wider(logfc_df,
                          id_cols    = promoterId,
                          names_from = method,
                          values_from = logFC)
logfc_wide_all <- pivot_wider(logfc_df_all,
                          id_cols    = promoterId,
                          names_from = method,
                          values_from = logFC)


for (m in c("Salmon", "DEXSeq", "proActiv")) {
  for (sgn in c("up", "down")) {

    df_tmp <- logfc_wide[!is.na(logfc_wide[[m]]), ] #!is.na(logfc_wide$CAGE) & !is.na(logfc_wide[[m]])
    df_tmp <- logfc_wide_all[logfc_wide_all$promoterId %in% df_tmp$promoterId, ]

    if (sgn == "up") {
      df_tmp <- df_tmp[df_tmp$CAGE > 0 & df_tmp[[m]] > 0, ]
    } else {
      df_tmp <- df_tmp[df_tmp$CAGE < 0 & df_tmp[[m]] < 0, ]
    }

    if (nrow(df_tmp) == 0L) {
      message("Skip ", m, " - ", sgn, " (empty)")
      next
    }

    df_plot <- data.frame(
        promoterId = df_tmp$promoterId,
        cage_log   = as.numeric(abs(df_tmp$CAGE)),
        other_log  = as.numeric(abs(df_tmp[[m]])),
        method     = m,
        stringsAsFactors = FALSE
    )

    check_cor <- function(df_plot) {
        keep <- is.finite(df_plot$cage_log) & is.finite(df_plot$other_log)
        x <- df_plot$cage_log[keep]
        y <- df_plot$other_log[keep]
        cat("Pearson:", cor(x, y, method = "pearson"), "\n")
        cat("Spearman:", cor(x, y, method = "spearman"), "\n")
    }
    check_cor(df_plot)
    plot_density_scatter_for_logFC(
      df = df_plot,
      title = "",  
      outfile = paste0("scatter_logFC_", m, "_", sgn, ".pdf"),
      use_density = TRUE
    )
  }
}

# Generate boxplot and violin plot
plot_box_violin_logFC(
  df = logfc_df,
  title = "log2FC Distribution by Method",
  outfile = "box_violin_logFC_all_methods_down.pdf",
  direction = "down"
)
plot_box_violin_logFC(
  df = logfc_df,
  title = "log2FC Distribution by Method",
  outfile = "box_violin_logFC_all_methods_up.pdf",
  direction = "up"
)
plot_box_violin_logFC(
  df = logfc_df,
  title = "log2FC Distribution by Method",
  outfile = "box_violin_logFC_all_methods.pdf"
)



#--------------------------------------------------------------------------------------
##scatterplot of Promoter usage change
#--------------------------------------------------------------------------------------

library(SummarizedExperiment)
library(ggplot2)
library(MASS)

cage_sig <- readRDS("/mnt/citadel2/research/shared/AltPromoterFlow/3_Cell_Lines/cage/differential/Promoter_differential_activity_FDR0_05.rds")
proactiv_sig <- readRDS("/mnt/citadel2/research/shared/AltPromoterFlow/3_Cell_Lines/proactiv/differential/Promoter_differential_activity_FDR0_05.rds")
salmon_sig <- readRDS("/mnt/citadel2/research/shared/AltPromoterFlow/3_Cell_Lines/salmon/differential/Promoter_differential_activity_FDR0_05.rds")
dexseq_sig <- readRDS("/mnt/citadel2/research/shared/AltPromoterFlow/3_Cell_Lines/dexseq/differential/Promoter_differential_activity_FDR0_05.rds")
prom_ids <- unique(as.numeric(cage_sig$promoterId))

se_cage    <- readRDS("/mnt/citadel2/research/shared/AltPromoterFlow/3_Cell_Lines/cage/merge/Promoter_activity_SE.rds")
se_pro     <- readRDS("/mnt/citadel2/research/shared/AltPromoterFlow/3_Cell_Lines/proactiv/merge/Promoter_activity_SE.rds")
se_salmon  <- readRDS("/mnt/citadel2/research/shared/AltPromoterFlow/3_Cell_Lines/salmon/merge/Promoter_activity_SE.rds")
se_dexseq  <- readRDS("/mnt/citadel2/research/shared/AltPromoterFlow/3_Cell_Lines/dexseq/merge/Promoter_activity_SE.rds")
out_dir <- "/mnt/citadel2/research/syidan/Projects/SnakeAltPromoterResult/K562_vs_GM12878"

# Use SummarizedExperiment calculate Δrelative (K562 - GM12878)
compute_delta_rel <- function(se, cond1 = "K562", cond0 = "GM12878") {
  stopifnot("relativePromoterActivity" %in% names(assays(se)))
  rel <- assays(se)$relativePromoterActivity
  cond <- as.character(SummarizedExperiment::colData(se)$condition)
  # Find columns of every group
  i1 <- which(cond == cond1)
  i0 <- which(cond == cond0)
  if (length(i1) == 0 || length(i0) == 0) {
    stop("在 colData(se)$condition 里找不到 ", cond1, " 或 ", cond0)
  }
  mu1 <- rowMeans(rel[, i1, drop = FALSE], na.rm = TRUE)
  mu0 <- rowMeans(rel[, i0, drop = FALSE], na.rm = TRUE)
  data.frame(
    promoterId = as.numeric(SummarizedExperiment::rowData(se)$promoterId),
    delta_rel  = mu1 - mu0,
    stringsAsFactors = FALSE
  )
}

cage_delta   <- compute_delta_rel(se_cage)
pro_delta    <- compute_delta_rel(se_pro)
salmon_delta <- compute_delta_rel(se_salmon)
dexseq_delta <- compute_delta_rel(se_dexseq)

# Only keeping ids significant in differential activity
keep <- data.frame(promoterId = prom_ids)
cage_delta   <- merge(keep, cage_delta,   by = "promoterId", all.x = TRUE)
pro_delta    <- merge(keep, pro_delta,    by = "promoterId", all.x = TRUE,
                      suffixes = c("", "_pro"))
salmon_delta <- merge(keep, salmon_delta, by = "promoterId", all.x = TRUE,
                      suffixes = c("", "_salmon"))
dexseq_delta <- merge(keep, dexseq_delta, by = "promoterId", all.x = TRUE,
                      suffixes = c("", "_dexseq"))



split_up_down_ids <- function(sig_df) {
  list(
    up   = as.numeric(sig_df$promoterId[sig_df$logFC > 0]),
    down = as.numeric(sig_df$promoterId[sig_df$logFC < 0])
  )
}
cage_ids    <- split_up_down_ids(cage_sig)
pro_ids     <- split_up_down_ids(proactiv_sig)
salmon_ids  <- split_up_down_ids(salmon_sig)
dexseq_ids  <- split_up_down_ids(dexseq_sig)

# Scatter plot of each method versus cage
plot_vs_cage_on_ids <- function(cage_df, method_df, other_name, ids, out_pdf,
                                n_grid = 200, point_alpha = 0.8, min_n = 10, thr = 0.1) {
  df_cage  <- cage_df[cage_df$promoterId %in% ids, ]
  df_other <- method_df[method_df$promoterId %in% ids, ]
  df <- merge(df_cage, df_other, by = "promoterId", suffixes = c("_cage", "_other"))

  # Only keeping finite
  keep <- is.finite(df$delta_rel_cage) & is.finite(df$delta_rel_other)
  df <- df[keep, ]
  keep <- abs(df$delta_rel_cage) >= thr | abs(df$delta_rel_other) >= thr
  df <- df[keep, ]
  if (nrow(df) == 0) {
    message("Skip ", other_name, " (no finite points).")
    return(invisible(NULL))
  }

  x <- df$delta_rel_cage
  y <- df$delta_rel_other

  # Use points, sd to determine whether use density
  use_density <- TRUE
  if (nrow(df) < min_n ||
      stats::sd(x) == 0 || stats::sd(y) == 0 ||
      length(unique(x)) < 5 || length(unique(y)) < 5) {
    use_density <- FALSE
  }

  df$density <- NA_real_
  if (use_density) {
    # prevent NA
    hx <- max(stats::bw.nrd(x), .Machine$double.eps)
    hy <- max(stats::bw.nrd(y), .Machine$double.eps)
    kd <- try(MASS::kde2d(x, y, n = n_grid, h = c(hx, hy)), silent = TRUE)
    if (!inherits(kd, "try-error")) {
      ix <- findInterval(x, kd$x)
      iy <- findInterval(y, kd$y)
      df$density <- kd$z[cbind(ix, iy)]
    } else {
      use_density <- FALSE
      message("Fall back to plain points for ", other_name, " (kde2d failed).")
    }
  }

  # Correlation: NA but do not return error
  can_cor <- (nrow(df) >= 3) && (stats::sd(x) > 0) && (stats::sd(y) > 0)
  if (can_cor) {
    cp <- try(suppressWarnings(cor.test(x, y, method = "pearson")),  silent = TRUE)
    cs <- try(suppressWarnings(cor.test(x, y, method = "spearman")), silent = TRUE)
    pearson_r <- if (inherits(cp, "try-error")) NA_real_ else unname(cp$estimate)
    pearson_p <- if (inherits(cp, "try-error")) NA_real_ else cp$p.value
    spear_r   <- if (inherits(cs, "try-error")) NA_real_ else unname(cs$estimate)
    spear_p   <- if (inherits(cs, "try-error")) NA_real_ else cs$p.value
  } else {
    pearson_r <- NA_real_; pearson_p <- NA_real_
    spear_r   <- NA_real_; spear_p   <- NA_real_
  }
  ann <- sprintf("Pearson r = %s (p %s)\nSpearman = %s (p %s)",
                 ifelse(is.na(pearson_r), "NA", sprintf("%.3f", pearson_r)),
                 ifelse(is.na(pearson_p), "NA", format(pearson_p, digits = 2, scientific = TRUE)),
                 ifelse(is.na(spear_r),   "NA", sprintf("%.3f", spear_r)),
                 ifelse(is.na(spear_p),   "NA", format(spear_p, digits = 2, scientific = TRUE)))

  lims <- range(c(x, y), finite = TRUE)
  if (!is.finite(diff(lims)) || diff(lims) == 0) {
    pad <- max(1e-3, abs(lims[1]) * 0.05)
    lims <- c(lims[1] - pad, lims[2] + pad)
  }

  p <- ggplot2::ggplot(df, ggplot2::aes(x = delta_rel_cage, y = delta_rel_other)) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey50") +
    ggplot2::scale_x_continuous(limits = lims) +
    ggplot2::scale_y_continuous(limits = lims) +
    ggplot2::labs(
      x = "CAGE  (Promoter usage: K562 - GM12878)",
      y = paste0(other_name, "  (Promoter usage: K562 - GM12878)"),
      title = paste0(other_name, " vs CAGE")
    ) +
    ggplot2::annotate("text", x = lims[1], y = lims[2], hjust = 0, vjust = 1,
                      size = 3.2, label = ann) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))

  if (use_density) {
    p <- p +
      ggplot2::geom_point(ggplot2::aes(color = density), size = 1.2, alpha = point_alpha) +
      ggplot2::scale_color_viridis_c(option = "A", name = "Density")
  } else {
    p <- p + ggplot2::geom_point(alpha = 0.7, size = 1)
  }

  ggplot2::ggsave(out_pdf, p, width = 5.5, height = 4.5)
  message("Wrote: ", out_pdf)
}

# proActiv
ids_up   <- intersect(cage_ids$up,   pro_ids$up)
ids_down <- intersect(cage_ids$down, pro_ids$down)
plot_vs_cage_on_ids(cage_delta, pro_delta, "proActiv",
  ids_up,   file.path(out_dir, "scatter_deltaRel_proActiv_vs_CAGE_up.pdf"))
plot_vs_cage_on_ids(cage_delta, pro_delta, "proActiv",
  ids_down, file.path(out_dir, "scatter_deltaRel_proActiv_vs_CAGE_down.pdf"))

# Salmon
ids_up   <- intersect(cage_ids$up,   salmon_ids$up)
ids_down <- intersect(cage_ids$down, salmon_ids$down)
plot_vs_cage_on_ids(cage_delta, salmon_delta, "Salmon",
  ids_up,   file.path(out_dir, "scatter_deltaRel_Salmon_vs_CAGE_up.pdf"))
plot_vs_cage_on_ids(cage_delta, salmon_delta, "Salmon",
  ids_down, file.path(out_dir, "scatter_deltaRel_Salmon_vs_CAGE_down.pdf"))

# DEXSeq
ids_up   <- intersect(cage_ids$up,   dexseq_ids$up)
ids_down <- intersect(cage_ids$down, dexseq_ids$down)
plot_vs_cage_on_ids(cage_delta, dexseq_delta, "DEXSeq",
  ids_up,   file.path(out_dir, "scatter_deltaRel_DEXSeq_vs_CAGE_up.pdf"))
plot_vs_cage_on_ids(cage_delta, dexseq_delta, "DEXSeq",
  ids_down, file.path(out_dir, "scatter_deltaRel_DEXSeq_vs_CAGE_down.pdf"))


#-------------
# Sanity Check
#-------------

library(dplyr)

methods <- list(Salmon = salmon, DEXSeq = dexseq, proActiv = proactiv)

summary_table <- data.frame(
  Method     = character(),
  Direction  = character(),
  N_Pairs    = integer(),
  Min        = numeric(),
  Q1         = numeric(),
  Median     = numeric(),
  Mean       = numeric(),
  Q3         = numeric(),
  Max        = numeric(),
  stringsAsFactors = FALSE
)

for (method_name in names(methods)) {
  df_method <- methods[[method_name]]
  df_common <- inner_join(
    cage[, c("promoterId", "logFC")],
    df_method[, c("promoterId", "logFC")],
    by = "promoterId",
    suffix = c("_cage", "_method")
  )

  for (dir in c("up", "down")) {
    if (dir == "up") {
      df_sub <- df_common %>% filter(logFC_cage > 0, logFC_method > 0)
    } else {
      df_sub <- df_common %>% filter(logFC_cage < 0, logFC_method < 0)
    }

    if (nrow(df_sub) == 0) next

    stat <- summary(df_sub$logFC_method)
    stat_num <- as.numeric(stat)

    summary_table <- rbind(summary_table, data.frame(
    Method    = method_name,
    Direction = dir,
    N_Pairs   = nrow(df_sub),
    Min       = stat_num[1],
    Q1        = stat_num[2],
    Median    = stat_num[3],
    Mean      = stat_num[4],
    Q3        = stat_num[5],
    Max       = stat_num[6]
    ))

  }
}

print(summary_table)


    Method Direction N_Pairs         Min        Q1    Median      Mean
1   Salmon        up    1078   0.3231317  1.351549  2.664919  4.035902
2   Salmon      down    1651 -23.9154509 -8.489788 -5.099268 -5.616285
3   DEXSeq        up     583   0.6415377  2.054174  3.766849  4.314652
4   DEXSeq      down    1042 -16.4714099 -7.566906 -5.807926 -5.664206
5 proActiv        up     640   0.6317126  1.949220  3.943887  4.375141
6 proActiv      down    1066 -13.5291621 -7.126131 -5.588291 -5.375614
         Q3        Max
1  6.169641 24.4513245
2 -2.008386 -0.3412513
3  6.020754 13.4673443
4 -3.012444 -0.6648830
5  6.289326 13.2363037
6 -3.211264 -0.7706925

