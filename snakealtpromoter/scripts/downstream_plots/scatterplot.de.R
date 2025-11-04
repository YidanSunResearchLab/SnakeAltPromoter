
#/genome/.snakemake_conda/671ce29470b1e8f250516402ed04d982_

# ---------------------------
# Plot: feature counts and no intron
# ---------------------------

library(DESeq2)
library(ggplot2)
library(MASS)
library(viridis)
library(DescTools)   
library(tidyr)
library(SummarizedExperiment)
library(dplyr)


#out_dir <- "/mnt/citadel2/research/shared/SnakeAltPromoter_paper_revision/scatterplot"

##Gene annotation to get intronless promoters
load(file.path("/mnt/citadelb/publication/snakealtpromoter/Result/lengths_promoters.RData"))
promoter_anno <- readRDS("/mnt/citadel2/research/shared/AltPromoterFlow/genome/organisms/hg38/Annotation/proActiv_promoter_annotation.rds")
coord <- promoter_anno@promoterCoordinates
intronless_ids <- as.character(coord$promoterId[!coord$internalPromoter %in% TRUE &
                                            lengths(coord$intronId) == 0])
intron_ids <- as.character(coord$promoterId[!coord$internalPromoter %in% TRUE &
                                            !lengths(coord$intronId) == 0])
all_ids <- as.character(coord$promoterId[!coord$internalPromoter %in% TRUE])



#--------------------------------------------------------------------------------------
#Plot scatterplot of all log2FC of Proactiv, Salmon and DEXseq vs CAGE counts

add_method <- function(df, m){
  df$method <- m
  df[, c("promoterId", "logFC", "method")]
}


plot_density_scatter_for_logFC <- function(df, title = "", outfile, use_density = TRUE, n_grid = 200) {
  if (nrow(df) == 0) {
    message("No rows for ", outfile); return(invisible(NULL))
  }

  # --- Log values (already log2FC, no +1) ---
  df$logC <- df$cage_log
  df$logO <- df$other_log
  df$density <- NA_real_

  # --- Density ---
  if (use_density) {
    keep <- is.finite(df$logC) | is.finite(df$logO) #at least in one method
    dens <- with(df[keep, ], MASS::kde2d(logC, logO, n = n_grid))
    ix <- findInterval(df$logC[keep], dens$x)
    iy <- findInterval(df$logO[keep], dens$y)
    df$density[keep] <- dens$z[cbind(ix, iy)]  # Keep raw density (not log)
  }

  # --------------------------------------------------------------
  #  cor_labels – professional, consistent, with white box
  # --------------------------------------------------------------
  cor_labels <- do.call(rbind, lapply(split(df, df$method), function(sub) {
    x <- sub$logC
    y <- sub$logO
    ok <- is.finite(x) & is.finite(y)

    if (sum(ok) < 2) {
      return(data.frame(
        method = unique(sub$method),
        x = NA_real_, y = NA_real_, label = "N < 2",
        stringsAsFactors = FALSE
      ))
    }

    # --- Pearson ---
    cp <- try(suppressWarnings(cor.test(x[ok], y[ok], method = "pearson")), silent = TRUE)
    pearson_r <- if (inherits(cp, "try-error")) NA_real_ else unname(cp$estimate)
    pearson_p <- if (inherits(cp, "try-error")) NA_real_ else cp$p.value

    # --- Spearman ---
    cs <- try(suppressWarnings(cor.test(x[ok], y[ok], method = "spearman")), silent = TRUE)
    spear_r <- if (inherits(cs, "try-error")) NA_real_ else unname(cs$estimate)
    spear_p <- if (inherits(cs, "try-error")) NA_real_ else cs$p.value

    # --- CCC ---
    ccc_val <- tryCatch({
      CCC(x[ok], y[ok])$rho.c$est
    }, error = function(e) NA_real_)
    ccc_val <- round(ccc_val, 3)

    # --- Mean Bias ---
    mean_diff <- mean(y[ok] - x[ok])
    bias_str <- sprintf("%.3f", round(mean_diff, 3))

    # --- BA bias ---
    diff_xy <- y[ok] - x[ok]
    mean_xy <- (y[ok] + x[ok]) / 2
    ba_corr <- try(suppressWarnings(cor(mean_xy, diff_xy)), silent = TRUE)
    ba_corr <- if (inherits(ba_corr, "try-error")) NA_real_ else round(ba_corr, 3)

    # --- p-value formatting ---
    pfmt <- function(p) {
      if (is.na(p)) "NA"
      else if (p < 1e-16) "< 1e-16"
      else sprintf("%.1g", p)
    }

    # --- Professional label ---
    label <- paste0(
      "Pearson's r = ", ifelse(is.na(pearson_r), "NA", sprintf("%.3f", pearson_r)), "\n",
      "Spearman's ρ = ", ifelse(is.na(spear_r), "NA", sprintf("%.3f", spear_r)), "\n",
      "Concordance r = ", ifelse(is.na(ccc_val), "NA", ccc_val), "\n",
      "Mean bias = ", bias_str, "\n",
      "BA bias = ", ifelse(is.na(ba_corr), "NA", ba_corr)
      )

    # --- Position: 5% from top-left ---
    xlim <- range(sub$logC[ok], na.rm = TRUE)
    ylim <- range(c(sub$logC[ok], sub$logO[ok]), na.rm = TRUE)
    label_x <- xlim[1] + 0.001 * diff(xlim)
    label_y <- ylim[2] - 0.001 * diff(ylim)

    data.frame(
      method = unique(sub$method),
      x = label_x,
      y = label_y,
      label = label,
      pearson_r = pearson_r,
      spearman_rho = spear_r,
      ccc_val = ccc_val,
      bias_str = bias_str,
      ba_corr = ba_corr,
      stringsAsFactors = FALSE
    )
  }))

  # --- Plot limits ---
  lims <- range(c(df$logC, df$logO), finite = TRUE)

  # --- Base plot ---
  p <- ggplot(df, aes(x = logC, y = logO)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey50") +
    scale_x_continuous(limits = lims, breaks = scales::pretty_breaks()) +
    scale_y_continuous(limits = lims, breaks = scales::pretty_breaks()) +
    facet_wrap(~method, nrow = 1) +
    labs(
      x = "CAGE (|log2FC|)",
      y = "Other method (|log2FC|)",
      title = title
    ) +
    theme_bw(base_size = 11) +
    theme(
      plot.title = element_text(hjust = 0.5),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      strip.background = element_rect(fill = "white")
    )

  # --- Points ---
  if (use_density) {
    p <- p + geom_point(aes(color = density), size = 1, alpha = 0.85) +
      scale_color_viridis_c(option = "A", name = "Density")
  } else {
    p <- p + geom_point(alpha = 0.5, size = 0.7)
  }

  # --- Label with white background (always on top) ---
  p <- p + geom_label(
    data = cor_labels,
    aes(x = x, y = y, label = label),
    hjust = 0, vjust = 1,
    size = 3.2,
    fill = "white",
    color = "black",
    label.padding = unit(0.2, "lines"),
    label.r = unit(0.15, "lines"),
    alpha = 1
  ) 

  # --- Save ---
  ggsave(file.path(out_dir, outfile), p, width = 5, height = 4, dpi = 300)
  # Return the correlations
  cor_results <- cor_labels[, c("method", "pearson_r", "spearman_rho", "ccc_val", "bias_str", "ba_corr")]
  return(cor_results)
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



combine.logFC.for.each.dataset <- function(logfc_wide, logfc_wide_all, out_dir, filename) {
    dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
    if (filename=="nointron"){methods.name=c("Salmon", "DEXSeq")} else {methods.name=c("Salmon", "DEXSeq", "proActiv")}
    # Generate boxplot and violin plot
    cor_list <- list()    
    for (m in methods.name) {
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
        cor_list[[paste(m,sgn,sep="_")]] <- plot_density_scatter_for_logFC(
          df = df_plot,
          title = "",  
          outfile = paste0("scatter_logFC_",filename, "_", m, "_", sgn, ".pdf"),
          use_density = TRUE
        )
      }
    }
    cor_df <- do.call(rbind, cor_list)
    write.csv(t(cor_df), file = file.path(out_dir, sprintf("scatter_logFC_%s_correlations.csv", filename)), row.names = TRUE, col.names = FALSE)

    # Generate boxplot and violin plot
    plot_box_violin_logFC(
      df = logfc_df,
      title = "log2FC Distribution by Method",
      outfile = paste0("box_violin_logFC_",filename,"_methods_down.pdf"),
      direction = "down"
    )
    plot_box_violin_logFC(
      df = logfc_df,
      title = "log2FC Distribution by Method",
      outfile = paste0("box_violin_logFC_",filename,"_methods_up.pdf"),
      direction = "up"
    )
    plot_box_violin_logFC(
      df = logfc_df,
      title = "log2FC Distribution by Method",
      outfile = paste0("box_violin_logFC_",filename,"_methods_all.pdf")
    )
}



# dexseq <- readRDS("/mnt/citadel2/research/shared/AltPromoterFlow/test_for_sample_sheet/dexseq/differential/comparisons_/K562_vs_H1-hESC/Promoter_differential_activity_FDR0_05.rds")
# salmon <- readRDS("/mnt/citadel2/research/shared/AltPromoterFlow/test_for_sample_sheet/salmon/differential/comparisons_/K562_vs_H1-hESC/Promoter_differential_activity_FDR0_05.rds")
# proactiv <- readRDS("/mnt/citadel2/research/shared/AltPromoterFlow/test_for_sample_sheet/proactiv/differential/comparisons_/K562_vs_H1-hESC/Promoter_differential_activity_FDR0_05.rds")
# cage <- readRDS("/mnt/citadel2/research/shared/AltPromoterFlow/test_for_sample_sheet/cage/differential/comparisons_/K562_vs_H1-hESC/Promoter_differential_activity_FDR0_05.rds")
# dexseq_all <- readRDS("/mnt/citadel2/research/shared/AltPromoterFlow/test_for_sample_sheet/dexseq/differential/comparisons_/K562_vs_H1-hESC/Promoter_differential_activity.rds")
# salmon_all <- readRDS("/mnt/citadel2/research/shared/AltPromoterFlow/test_for_sample_sheet/salmon/differential/comparisons_/K562_vs_H1-hESC/Promoter_differential_activity.rds")
# proactiv_all <- readRDS("/mnt/citadel2/research/shared/AltPromoterFlow/test_for_sample_sheet/proactiv/differential/comparisons_/K562_vs_H1-hESC/Promoter_differential_activity.rds")
# cage_all <- readRDS("/mnt/citadel2/research/shared/AltPromoterFlow/test_for_sample_sheet/cage/differential/comparisons_/K562_vs_H1-hESC/Promoter_differential_activity.rds")
# #out_dir <- "/mnt/citadelb/publication/snakealtpromoter/Result/K562_vs_H1-hESC"
# out_dir <- "/mnt/citadel2/research/shared/SnakeAltPromoter_paper_revision/scatterplot/K562_vs_H"

#gene_name = read.delim("/mnt/citadel2/research/syidan/Genomes/GRCh38/release-47-index/annotation/genes.symbol", header = FALSE, comment.char = "#")
#cage$geneSymbol <- gene_name$V2[match(cage$geneId, gene_name$V1)]
#write.table(unique(cage$geneSymbol[!is.na(cage$geneSymbol) & cage$FDR < 0.05 & cage$logFC > 0]),  file.path(out_dir,"cage_upregulated_genes.txt"),  quote = FALSE, row.names = FALSE, col.names = FALSE)
#write.table(unique(cage$geneSymbol[!is.na(cage$geneSymbol) & cage$FDR < 0.05 & cage$logFC < 0]), file.path(out_dir,"cage_downregulated_genes.txt"), quote = FALSE, row.names = FALSE, col.names = FALSE)

#--------------------------------------------------------
###For heart healthy and failure
#--------------------------------------------------------
if(TRUE){
  out_dir <- "/mnt/citadelb/publication/snakealtpromoter/Result/Heart_Failure_vs_Healthy"
  # dexseq <- readRDS("/mnt/citadel2/research/shared/AltPromoterFlow/new_test_for_sample_sheet/dexseq/differential/comparisons_/Heart_Failure_vs_Healthy/Promoter_differential_activity_FDR0_05.rds")
  # salmon <- readRDS("/mnt/citadel2/research/shared/AltPromoterFlow/new_test_for_sample_sheet/salmon/differential/comparisons_/Heart_Failure_vs_Healthy/Promoter_differential_activity_FDR0_05.rds")
  # proactiv <- readRDS("/mnt/citadel2/research/shared/AltPromoterFlow/new_test_for_sample_sheet/proactiv/differential/comparisons_/Heart_Failure_vs_Healthy/Promoter_differential_activity_FDR0_05.rds")
  # cage <- readRDS("/mnt/citadel2/research/shared/AltPromoterFlow/new_test_for_sample_sheet/cage/differential/comparisons_/Heart_Failure_vs_Healthy/Promoter_differential_activity_FDR0_05.rds")
  # dexseq_all <- readRDS("/mnt/citadel2/research/shared/AltPromoterFlow/new_test_for_sample_sheet/dexseq/differential/comparisons_/Heart_Failure_vs_Healthy/Promoter_differential_activity.rds")
  # salmon_all <- readRDS("/mnt/citadel2/research/shared/AltPromoterFlow/new_test_for_sample_sheet/salmon/differential/comparisons_/Heart_Failure_vs_Healthy/Promoter_differential_activity.rds")
  # proactiv_all <- readRDS("/mnt/citadel2/research/shared/AltPromoterFlow/new_test_for_sample_sheet/proactiv/differential/comparisons_/Heart_Failure_vs_Healthy/Promoter_differential_activity.rds")
  # cage_all <- readRDS("/mnt/citadel2/research/shared/AltPromoterFlow/new_test_for_sample_sheet/cage/differential/comparisons_/Heart_Failure_vs_Healthy/Promoter_differential_activity.rds")
  dexseq <- readRDS("/mnt/citadelb/publication/snakealtpromoter/Processed/heart_AltPromoter/RNAseq/dexseq/differential/comparisons_/Failure_vs_Healthy/Promoter_differential_activity_FDR0_05.rds")
  salmon <- readRDS("/mnt/citadelb/publication/snakealtpromoter/Processed/heart_AltPromoter/RNAseq/salmon/differential/comparisons_/Failure_vs_Healthy/Promoter_differential_activity_FDR0_05.rds")
  proactiv <- readRDS("/mnt/citadelb/publication/snakealtpromoter/Processed/heart_AltPromoter/RNAseq/proactiv/differential/comparisons_/Failure_vs_Healthy/Promoter_differential_activity_FDR0_05.rds")
  cage <- readRDS("/mnt/citadelb/publication/snakealtpromoter/Processed/heart_AltPromoter/CAGE/cage/differential/comparisons_/Failure_vs_Healthy/Promoter_differential_activity_FDR0_05.rds")
  dexseq_all <- readRDS("/mnt/citadelb/publication/snakealtpromoter/Processed/heart_AltPromoter/RNAseq/dexseq/differential/comparisons_/Failure_vs_Healthy/Promoter_differential_activity.rds")
  salmon_all <- readRDS("/mnt/citadelb/publication/snakealtpromoter/Processed/heart_AltPromoter/RNAseq/salmon/differential/comparisons_/Failure_vs_Healthy/Promoter_differential_activity.rds")
  proactiv_all <- readRDS("/mnt/citadelb/publication/snakealtpromoter/Processed/heart_AltPromoter/RNAseq/proactiv/differential/comparisons_/Failure_vs_Healthy/Promoter_differential_activity.rds")
  cage_all <- readRDS("/mnt/citadelb/publication/snakealtpromoter/Processed/heart_AltPromoter/CAGE/cage/differential/comparisons_/Failure_vs_Healthy/Promoter_differential_activity.rds")
  logfc_df     <- rbind(add_method(cage, "CAGE"), add_method(salmon, "Salmon"), add_method(dexseq, "DEXSeq"), add_method(proactiv, "proActiv"))
  logfc_df_all <- rbind(add_method(cage_all, "CAGE"), add_method(salmon_all, "Salmon"), add_method(dexseq_all, "DEXSeq"), add_method(proactiv_all, "proActiv"))

  ##Intronless
  logfc_wide <- pivot_wider(logfc_df[logfc_df$promoterId %in% intronless_ids,],id_cols    = promoterId,names_from = method,values_from = logFC)
  logfc_wide_all <- pivot_wider(logfc_df_all[logfc_df_all$promoterId %in% intronless_ids,], id_cols    = promoterId, names_from = method, values_from = logFC)
  combine.logFC.for.each.dataset(logfc_wide, logfc_wide_all, out_dir, filename="nointron")

  ##With Intron
  logfc_wide <- pivot_wider(logfc_df[logfc_df$promoterId %in% intron_ids,],id_cols    = promoterId,names_from = method,values_from = logFC)
  logfc_wide_all <- pivot_wider(logfc_df_all[logfc_df_all$promoterId %in% intron_ids,], id_cols    = promoterId, names_from = method, values_from = logFC)
  combine.logFC.for.each.dataset(logfc_wide, logfc_wide_all, out_dir, filename="withintron")

}

#--------------------------------------------------------
###For Brain female and male
#--------------------------------------------------------
if(TRUE){
  out_dir <- "/mnt/citadelb/publication/snakealtpromoter/Result/Brain"
  dexseq <- readRDS("/mnt/citadelb/publication/snakealtpromoter/Processed/brain_AltPromoter/RNAseq_test/dexseq/differential/comparisons_/male_vs_female/Promoter_differential_activity_FDR0_05.rds")
  salmon <- readRDS("/mnt/citadelb/publication/snakealtpromoter/Processed/brain_AltPromoter/RNAseq_test/salmon/differential/comparisons_/male_vs_female/Promoter_differential_activity_FDR0_05.rds")
  proactiv <- readRDS("/mnt/citadelb/publication/snakealtpromoter/Processed/brain_AltPromoter/RNAseq_test/proactiv/differential/comparisons_/male_vs_female/Promoter_differential_activity_FDR0_05.rds")
  cage <- readRDS("/mnt/citadelb/publication/snakealtpromoter/Processed/brain_AltPromoter/CAGE_test/cage/differential/comparisons_/male_vs_female/Promoter_differential_activity_FDR0_05.rds")
  dexseq_all <- readRDS("/mnt/citadelb/publication/snakealtpromoter/Processed/brain_AltPromoter/RNAseq_test/dexseq/differential/comparisons_/male_vs_female/Promoter_differential_activity.rds")
  salmon_all <- readRDS("/mnt/citadelb/publication/snakealtpromoter/Processed/brain_AltPromoter/RNAseq_test/salmon/differential/comparisons_/male_vs_female/Promoter_differential_activity.rds")
  proactiv_all <- readRDS("/mnt/citadelb/publication/snakealtpromoter/Processed/brain_AltPromoter/RNAseq_test/proactiv/differential/comparisons_/male_vs_female/Promoter_differential_activity.rds")
  cage_all <- readRDS("/mnt/citadelb/publication/snakealtpromoter/Processed/brain_AltPromoter/CAGE_test/cage/differential/comparisons_/male_vs_female/Promoter_differential_activity.rds")
  logfc_df     <- rbind(add_method(cage, "CAGE"), add_method(salmon, "Salmon"), add_method(dexseq, "DEXSeq"), add_method(proactiv, "proActiv"))
  logfc_df_all <- rbind(add_method(cage_all, "CAGE"), add_method(salmon_all, "Salmon"), add_method(dexseq_all, "DEXSeq"), add_method(proactiv_all, "proActiv"))

  ##Intronless
  logfc_wide <- pivot_wider(logfc_df[logfc_df$promoterId %in% intronless_ids,],id_cols    = promoterId,names_from = method,values_from = logFC)
  logfc_wide_all <- pivot_wider(logfc_df_all[logfc_df_all$promoterId %in% intronless_ids,], id_cols    = promoterId, names_from = method, values_from = logFC)
  combine.logFC.for.each.dataset(logfc_wide, logfc_wide_all, out_dir, filename="nointron")

  ##With Intron
  logfc_wide <- pivot_wider(logfc_df[logfc_df$promoterId %in% intron_ids,],id_cols    = promoterId,names_from = method,values_from = logFC)
  logfc_wide_all <- pivot_wider(logfc_df_all[logfc_df_all$promoterId %in% intron_ids,], id_cols    = promoterId, names_from = method, values_from = logFC)
  combine.logFC.for.each.dataset(logfc_wide, logfc_wide_all, out_dir, filename="withintron")

}


#--------------------------------------------------------
###For cell lines
#--------------------------------------------------------
if(TRUE){
  out_dir <- "/mnt/citadelb/publication/snakealtpromoter/Result/K562_vs_GM12878"
  # dexseq <- readRDS("/mnt/citadel2/research/shared/AltPromoterFlow/test_for_sample_sheet/dexseq/differential/comparisons_/K562_vs_GM12878/Promoter_differential_activity_FDR0_05.rds")
  # salmon <- readRDS("/mnt/citadel2/research/shared/AltPromoterFlow/test_for_sample_sheet/salmon/differential/comparisons_/K562_vs_GM12878/Promoter_differential_activity_FDR0_05.rds")
  # proactiv <- readRDS("/mnt/citadel2/research/shared/AltPromoterFlow/test_for_sample_sheet/proactiv/differential/comparisons_/K562_vs_GM12878/Promoter_differential_activity_FDR0_05.rds")
  # cage <- readRDS("/mnt/citadel2/research/shared/AltPromoterFlow/test_for_sample_sheet/cage/differential/comparisons_/K562_vs_GM12878/Promoter_differential_activity_FDR0_05.rds")
  # dexseq_all <- readRDS("/mnt/citadel2/research/shared/AltPromoterFlow/test_for_sample_sheet/dexseq/differential/comparisons_/K562_vs_GM12878/Promoter_differential_activity.rds")
  # salmon_all <- readRDS("/mnt/citadel2/research/shared/AltPromoterFlow/test_for_sample_sheet/salmon/differential/comparisons_/K562_vs_GM12878/Promoter_differential_activity.rds")
  # proactiv_all <- readRDS("/mnt/citadel2/research/shared/AltPromoterFlow/test_for_sample_sheet/proactiv/differential/comparisons_/K562_vs_GM12878/Promoter_differential_activity.rds")
  # cage_all <- readRDS("/mnt/citadel2/research/shared/AltPromoterFlow/test_for_sample_sheet/cage/differential/comparisons_/K562_vs_GM12878/Promoter_differential_activity.rds")
  dexseq <- readRDS("/mnt/citadelb/publication/snakealtpromoter/Processed/cellline_AltPromoter/RNAseq/dexseq/differential/comparisons_/K562_vs_GM12878/Promoter_differential_activity_FDR0_05.rds")
  salmon <- readRDS("/mnt/citadelb/publication/snakealtpromoter/Processed/cellline_AltPromoter/RNAseq/salmon/differential/comparisons_/K562_vs_GM12878/Promoter_differential_activity_FDR0_05.rds")
  proactiv <- readRDS("/mnt/citadelb/publication/snakealtpromoter/Processed/cellline_AltPromoter/RNAseq/proactiv/differential/comparisons_/K562_vs_GM12878/Promoter_differential_activity_FDR0_05.rds")
  cage <- readRDS("/mnt/citadelb/publication/snakealtpromoter/Processed/cellline_AltPromoter/CAGE/cage/differential/comparisons_/K562_vs_GM12878/Promoter_differential_activity_FDR0_05.rds")
  dexseq_all <- readRDS("/mnt/citadelb/publication/snakealtpromoter/Processed/cellline_AltPromoter/RNAseq/dexseq/differential/comparisons_/K562_vs_GM12878/Promoter_differential_activity.rds")
  salmon_all <- readRDS("/mnt/citadelb/publication/snakealtpromoter/Processed/cellline_AltPromoter/RNAseq/salmon/differential/comparisons_/K562_vs_GM12878/Promoter_differential_activity.rds")
  proactiv_all <- readRDS("/mnt/citadelb/publication/snakealtpromoter/Processed/cellline_AltPromoter/RNAseq/proactiv/differential/comparisons_/K562_vs_GM12878/Promoter_differential_activity.rds")
  cage_all <- readRDS("/mnt/citadelb/publication/snakealtpromoter/Processed/cellline_AltPromoter/CAGE/cage/differential/comparisons_/K562_vs_GM12878/Promoter_differential_activity.rds")
  logfc_df     <- rbind(add_method(cage, "CAGE"), add_method(salmon, "Salmon"), add_method(dexseq, "DEXSeq"), add_method(proactiv, "proActiv"))
  logfc_df_all <- rbind(add_method(cage_all, "CAGE"), add_method(salmon_all, "Salmon"), add_method(dexseq_all, "DEXSeq"), add_method(proactiv_all, "proActiv"))

  ##All
  logfc_wide <- pivot_wider(logfc_df,id_cols    = promoterId,names_from = method,values_from = logFC)
  logfc_wide_all <- pivot_wider(logfc_df_all, id_cols    = promoterId, names_from = method, values_from = logFC)
  combine.logFC.for.each.dataset(logfc_wide, logfc_wide_all, out_dir, filename="all")

  ##Intronless
  logfc_wide <- pivot_wider(logfc_df[logfc_df$promoterId %in% intronless_ids,],id_cols    = promoterId,names_from = method,values_from = logFC)
  logfc_wide_all <- pivot_wider(logfc_df_all[logfc_df_all$promoterId %in% intronless_ids,], id_cols    = promoterId, names_from = method, values_from = logFC)
  combine.logFC.for.each.dataset(logfc_wide, logfc_wide_all, out_dir, filename="nointron")

  ##With Intron
  logfc_wide <- pivot_wider(logfc_df[logfc_df$promoterId %in% intron_ids,],id_cols    = promoterId,names_from = method,values_from = logFC)
  logfc_wide_all <- pivot_wider(logfc_df_all[logfc_df_all$promoterId %in% intron_ids,], id_cols    = promoterId, names_from = method, values_from = logFC)
  combine.logFC.for.each.dataset(logfc_wide, logfc_wide_all, out_dir, filename="withintron")

}




#--------------------------------------------------------------------------------------
##scatterplot of Promoter usage change
#--------------------------------------------------------------------------------------
# Use SummarizedExperiment calculate Δrelative (K562 - GM12878)
compute_delta_rel <- function(se, cond0 = "GM12878", cond1 = "K562") {
  stopifnot("relativePromoterActivity" %in% names(assays(se)))
  rel <- assays(se)$relativePromoterActivity
  cond <- as.character(SummarizedExperiment::colData(se)$condition)
  # Find columns of every group
  i1 <- which(cond == cond1)
  i0 <- which(cond == cond0)
  # mu1 <- rowMeans(rel[, i1, drop = FALSE], na.rm = TRUE)
  # mu0 <- rowMeans(rel[, i0, drop = FALSE], na.rm = TRUE)
  mu1 <- rowMeans(rel[, grep(cond0, colnames(rel)), drop = FALSE], na.rm = TRUE)
  mu0 <- rowMeans(rel[, grep(cond1, colnames(rel)), drop = FALSE], na.rm = TRUE)
  data.frame(
    promoterId = as.numeric(SummarizedExperiment::rowData(se)$promoterId),
    delta_rel  = mu1 - mu0,
    stringsAsFactors = FALSE
  )
}

split_up_down_ids <- function(sig_df) {
  list(
    up   = as.numeric(sig_df$promoterId[sig_df$logFC > 0]),
    down = as.numeric(sig_df$promoterId[sig_df$logFC < 0])
  )
}



# Scatter plot of each method versus cage
plot_vs_cage_on_ids <- function(cage_df, method_df, other_name, ids, out_pdf,
                                n_grid = 200, point_alpha = 0.8, min_n = 10, thr = 0.1) {
  # --- Merge and filter ---
  df_cage  <- cage_df[cage_df$promoterId %in% ids, ]
  df_other <- method_df[method_df$promoterId %in% ids, ]
  df <- merge(df_cage, df_other, by = "promoterId", suffixes = c("_cage", "_other"))

  # Keep only finite + threshold
  keep <- is.finite(df$delta_rel_cage) & is.finite(df$delta_rel_other)
  df <- df[keep, ]
  keep <- abs(df$delta_rel_cage) >= thr | abs(df$delta_rel_other) >= thr
  df <- df[keep, ]
  if (nrow(df) == 0) {
    message("Skip ", other_name, " (no finite points).")
    return(invisible(NULL))
  }

  x <- df$delta_rel_cage
  y <- c(df$delta_rel_other)

  # --- Decide on density ---
  use_density <- nrow(df) >= min_n &&
                 stats::sd(x) > 0 && stats::sd(y) > 0 &&
                 length(unique(x)) >= 5 && length(unique(y)) >= 5

  df$density <- NA_real_
  if (use_density) {
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

  # --- Correlation stats ---
  can_cor <- nrow(df) >= 3 && stats::sd(x) > 0 && stats::sd(y) > 0
  if (can_cor) {
    # Pearson
    cp <- try(cor.test(x, y, method = "pearson"), silent = TRUE)
    pearson_r <- if (inherits(cp, "try-error")) NA_real_ else unname(cp$estimate)
    pearson_p <- if (inherits(cp, "try-error")) NA_real_ else cp$p.value

    # Spearman
    cs <- try(suppressWarnings(cor.test(x, y, method = "spearman")), silent = TRUE)
    spear_r <- if (inherits(cs, "try-error")) NA_real_ else unname(cs$estimate)
    spear_p <- if (inherits(cs, "try-error")) NA_real_ else cs$p.value

    # CCC
    ccc_val <- tryCatch(CCC(x, y)$rho.c$est, error = function(e) NA_real_)
    ccc_val <- round(ccc_val, 3)

    # Mean bias
    bias_str <- round(mean(y - x), 3)

    # BA bias
    diff_xy <- y - x
    mean_xy <- (y + x) / 2
    ba_corr <- try(suppressWarnings(cor(mean_xy, diff_xy)), silent = TRUE)
    ba_corr <- if (inherits(ba_corr, "try-error")) NA_real_ else round(ba_corr, 3)
  } else {
    pearson_r <- pearson_p <- spear_r <- spear_p <- ccc_val <- bias_str <- ba_corr <- NA_real_
  }

  # --- p-value formatting ---
  pfmt <- function(p) {
    if (is.na(p)) "NA"
    else if (p < 1e-16) "< 1e-16"
    else sprintf("%.1g", p)
  }

  # --- Professional label ---
  label <- paste0(
    "Pearson's r = ", ifelse(is.na(pearson_r), "NA", sprintf("%.3f", pearson_r)),
    "\n",
    "Spearman's ρ = ", ifelse(is.na(spear_r), "NA", sprintf("%.3f", spear_r)),
    "\n",
    "Concordance r = ", ifelse(is.na(ccc_val), "NA", ccc_val), "\n",
    "Mean bias = ", bias_str, "\n",
    "BA bias = ", ifelse(is.na(ba_corr), "NA", ba_corr)
  )

  # --- Label position: 1% from top-left ---
  xlim <- range(x, na.rm = TRUE)
  ylim <- range(c(x,y), na.rm = TRUE)
  label_x <- xlim[1] + 0.001 * diff(xlim)
  label_y <- ylim[2] - 0.001 * diff(ylim)

  cor_labels <- data.frame(
    method = other_name,
    x = label_x, y = label_y, label = label,
      pearson_r = pearson_r,
      spearman_rho = spear_r,
      ccc_val = ccc_val,
      bias_str = bias_str,
      ba_corr = ba_corr,
    stringsAsFactors = FALSE
  )

  # --- Plot limits ---
  lims <- range(c(x, y), finite = TRUE)
  if (diff(lims) == 0) {
    pad <- max(1e-3, abs(lims[1]) * 0.05)
    lims <- c(lims[1] - pad, lims[2] + pad)
  }

  # --- Base plot ---
  p <- ggplot(df, aes(x = delta_rel_cage, y = delta_rel_other)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey50") +
    scale_x_continuous(limits = lims, breaks = scales::pretty_breaks()) +
    scale_y_continuous(limits = lims, breaks = scales::pretty_breaks()) +
    labs(
      x = "CAGE (Promoter usage: K562 - GM12878)",
      y = paste0(other_name, " (Promoter usage: K562 - GM12878)"),
      title = paste0(other_name, " vs CAGE")
    ) +
    theme_bw(base_size = 11) +
    theme(
      plot.title = element_text(hjust = 0.5),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    )

  # --- Points ---
  if (use_density) {
    p <- p + geom_point(aes(color = density), size = 1.2, alpha = point_alpha) +
      scale_color_viridis_c(option = "A", name = "Density")
  } else {
    p <- p + geom_point(alpha = 0.7, size = 1)
  }

  # --- White label (on top) ---
  p <- p + geom_label(
    data = cor_labels,
    aes(x = x, y = y, label = label),
    hjust = 0, vjust = 1,
    size = 3.2,
    fill = "white", color = "black",
    label.padding = unit(0.2, "lines"),
    label.r = unit(0.15, "lines"),
    alpha = 1
  ) +
    coord_cartesian(expand = FALSE)

  # --- Save ---
  ggsave(out_pdf, p, width = 5.5, height = 4.5, dpi = 300)
  # Return the correlations
  cor_results <- cor_labels[, c("method","pearson_r", "spearman_rho", "ccc_val", "bias_str", "ba_corr")]
  return(cor_results)
  message("Wrote: ", out_pdf)
}




#--------------------------------------------------------
###For cell lines
#--------------------------------------------------------
out_dir <- "/mnt/citadelb/publication/snakealtpromoter/Result/K562_vs_GM12878"
# cage_sig <- readRDS("/mnt/citadel2/research/shared/AltPromoterFlow/3_Cell_Lines/cage/differential/Promoter_differential_activity_FDR0_05.rds")
# proactiv_sig <- readRDS("/mnt/citadel2/research/shared/AltPromoterFlow/3_Cell_Lines/proactiv/differential/Promoter_differential_activity_FDR0_05.rds")
# salmon_sig <- readRDS("/mnt/citadel2/research/shared/AltPromoterFlow/3_Cell_Lines/salmon/differential/Promoter_differential_activity_FDR0_05.rds")
# dexseq_sig <- readRDS("/mnt/citadel2/research/shared/AltPromoterFlow/3_Cell_Lines/dexseq/differential/Promoter_differential_activity_FDR0_05.rds")
# prom_ids <- unique(as.numeric(cage_sig$promoterId))

# se_cage    <- readRDS("/mnt/citadel2/research/shared/AltPromoterFlow/3_Cell_Lines/cage/merge/Promoter_activity_SE.rds")
# se_pro     <- readRDS("/mnt/citadel2/research/shared/AltPromoterFlow/3_Cell_Lines/proactiv/merge/Promoter_activity_SE.rds")
# se_salmon  <- readRDS("/mnt/citadel2/research/shared/AltPromoterFlow/3_Cell_Lines/salmon/merge/Promoter_activity_SE.rds")
# se_dexseq  <- readRDS("/mnt/citadel2/research/shared/AltPromoterFlow/3_Cell_Lines/dexseq/merge/Promoter_activity_SE.rds")

cage_sig <- readRDS("/mnt/citadelb/publication/snakealtpromoter/Processed/cellline_AltPromoter/CAGE/cage/differential/Promoter_differential_activity_FDR0_05.rds")
proactiv_sig <- readRDS("/mnt/citadelb/publication/snakealtpromoter/Processed/cellline_AltPromoter/RNAseq/proactiv/differential/Promoter_differential_activity_FDR0_05.rds")
salmon_sig <- readRDS("/mnt/citadelb/publication/snakealtpromoter/Processed/cellline_AltPromoter/RNAseq/salmon/differential/Promoter_differential_activity_FDR0_05.rds")
dexseq_sig <- readRDS("/mnt/citadelb/publication/snakealtpromoter/Processed/cellline_AltPromoter/RNAseq/dexseq/differential/Promoter_differential_activity_FDR0_05.rds")
prom_ids <- unique(as.numeric(cage_sig$promoterId))

se_cage    <- readRDS("/mnt/citadelb/publication/snakealtpromoter/Processed/cellline_AltPromoter/CAGE/cage/promoter_classification_total/Promoter_activity_SE.rds")
se_pro     <- readRDS("/mnt/citadelb/publication/snakealtpromoter/Processed/cellline_AltPromoter/RNAseq/proactiv/promoter_classification_total/Promoter_activity_SE.rds")
se_salmon  <- readRDS("/mnt/citadelb/publication/snakealtpromoter/Processed/cellline_AltPromoter/RNAseq/salmon/promoter_classification_total/Promoter_activity_SE.rds")
se_dexseq  <- readRDS("/mnt/citadelb/publication/snakealtpromoter/Processed/cellline_AltPromoter/RNAseq/dexseq/promoter_classification_total/Promoter_activity_SE.rds")

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

cage_ids    <- split_up_down_ids(cage_sig)
pro_ids     <- split_up_down_ids(proactiv_sig)
salmon_ids  <- split_up_down_ids(salmon_sig)
dexseq_ids  <- split_up_down_ids(dexseq_sig)

################
##With Introns
################
# proActiv
ids_up   <- intersect(cage_ids$up,   pro_ids$up)
ids_down <- intersect(cage_ids$down, pro_ids$down)
ids_up <- ids_up[ids_up %in% intron_ids]
ids_down <- ids_down[ids_down %in% intron_ids]
cor_list1 <- plot_vs_cage_on_ids(cage_delta, pro_delta, "proActiv",
  ids_up,   file.path(out_dir, "scatter_deltaRel_withintron_proActiv_vs_CAGE_up.pdf"))
cor_list2 <- plot_vs_cage_on_ids(cage_delta, pro_delta, "proActiv",
  ids_down, file.path(out_dir, "scatter_deltaRel_withintron_proActiv_vs_CAGE_down.pdf"))

# Salmon
ids_up   <- intersect(cage_ids$up,   salmon_ids$up)
ids_down <- intersect(cage_ids$down, salmon_ids$down)
ids_up <- ids_up[ids_up %in% intron_ids]
ids_down <- ids_down[ids_down %in% intron_ids]
cor_list3 <- plot_vs_cage_on_ids(cage_delta, salmon_delta, "Salmon",
  ids_up,   file.path(out_dir, "scatter_deltaRel_withintron_Salmon_vs_CAGE_up.pdf"))
cor_list4 <- plot_vs_cage_on_ids(cage_delta, salmon_delta, "Salmon",
  ids_down, file.path(out_dir, "scatter_deltaRel_withintron_Salmon_vs_CAGE_down.pdf"))

# DEXSeq
ids_up   <- intersect(cage_ids$up,   dexseq_ids$up)
ids_down <- intersect(cage_ids$down, dexseq_ids$down)
ids_up <- ids_up[ids_up %in% intron_ids]
ids_down <- ids_down[ids_down %in% intron_ids]
cor_list5 <- plot_vs_cage_on_ids(cage_delta, dexseq_delta, "DEXSeq",
  ids_up,   file.path(out_dir, "scatter_deltaRel_withintron_DEXSeq_vs_CAGE_up.pdf"))
cor_list6 <- plot_vs_cage_on_ids(cage_delta, dexseq_delta, "DEXSeq",
  ids_down, file.path(out_dir, "scatter_deltaRel_withintron_DEXSeq_vs_CAGE_down.pdf"))

##export csv files
cor_df_up <- rbind(cor_list1, cor_list3, cor_list5)
cor_df_down <- rbind(cor_list2, cor_list4, cor_list6)
write.csv(t(cor_df_up), file = file.path(out_dir, "scatter_deltaRel_withintron_correlations_up.csv"), row.names = TRUE, col.names = FALSE)
write.csv(t(cor_df_down), file = file.path(out_dir, "scatter_deltaRel_withintron_correlations_down.csv"), row.names = TRUE, col.names = FALSE)


################
##No Introns
################
# proActiv
ids_up   <- intersect(cage_ids$up,   pro_ids$up)
ids_down <- intersect(cage_ids$down, pro_ids$down)
ids_up <- ids_up[ids_up %in% intronless_ids]
ids_down <- ids_down[ids_down %in% intronless_ids]
plot_vs_cage_on_ids(cage_delta, pro_delta, "proActiv",
  ids_up,   file.path(out_dir, "scatter_deltaRel_nointron_proActiv_vs_CAGE_up.pdf"))
plot_vs_cage_on_ids(cage_delta, pro_delta, "proActiv",
  ids_down, file.path(out_dir, "scatter_deltaRel_nointron_proActiv_vs_CAGE_down.pdf"))

# Salmon
ids_up   <- intersect(cage_ids$up,   salmon_ids$up)
ids_down <- intersect(cage_ids$down, salmon_ids$down)
ids_up <- ids_up[ids_up %in% intronless_ids]
ids_down <- ids_down[ids_down %in% intronless_ids]
plot_vs_cage_on_ids(cage_delta, salmon_delta, "Salmon",
  ids_up,   file.path(out_dir, "scatter_deltaRel_nointron_Salmon_vs_CAGE_up.pdf"))
plot_vs_cage_on_ids(cage_delta, salmon_delta, "Salmon",
  ids_down, file.path(out_dir, "scatter_deltaRel_nointron_Salmon_vs_CAGE_down.pdf"))

# DEXSeq
ids_up   <- intersect(cage_ids$up,   dexseq_ids$up)
ids_down <- intersect(cage_ids$down, dexseq_ids$down)
ids_up <- ids_up[ids_up %in% intronless_ids]
ids_down <- ids_down[ids_down %in% intronless_ids]
plot_vs_cage_on_ids(cage_delta, dexseq_delta, "DEXSeq",
  ids_up,   file.path(out_dir, "scatter_deltaRel_nointron_DEXSeq_vs_CAGE_up.pdf"))
plot_vs_cage_on_ids(cage_delta, dexseq_delta, "DEXSeq",
  ids_down, file.path(out_dir, "scatter_deltaRel_nointron_DEXSeq_vs_CAGE_down.pdf"))



