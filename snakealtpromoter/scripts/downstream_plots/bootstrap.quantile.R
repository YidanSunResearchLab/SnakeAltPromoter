# ---------------------------
# Plot: feature counts and no intron
# ---------------------------

library(DESeq2)
library(ggplot2)
library(MASS)
library(viridis)
library(DescTools)   
library(data.table)

##Gene annotation to get intronless promoters
coord <- promoter_anno@promoterCoordinates
intronless_ids <- as.character(coord$promoterId[!coord$internalPromoter %in% TRUE &
                                            lengths(coord$intronId) == 0])
intron_ids <- as.character(coord$promoterId[!coord$internalPromoter %in% TRUE &
                                            !lengths(coord$intronId) == 0])
all_ids <- as.character(coord$promoterId[!coord$internalPromoter %in% TRUE])


# define consistent plotting colors for methods
cols_method <- c(
  "Salmon"   = "#3C78D8",
  "DEXSeq"   = "#D64B4B",
  "proActiv" = "#7A5AA6"
)

theme_rev <- function() {
  # define publication-style theme for revision figures
  theme_classic(base_size = 14, base_family = "Helvetica") +
    theme(
      axis.title = element_text(size = 15, color = "black"),
      axis.text  = element_text(size = 13, color = "black"),
      strip.text = element_text(size = 14, color = "black"),
      legend.title = element_blank(),
      legend.text = element_text(size = 13, color = "black"),
      plot.title = element_text(size = 15, color = "black", hjust = 0.5)
    )
}


#--------------------------------------------------------
#Functions to Plot scatterplot of promoter counts of Proactiv Salmon and DEXseq vs CAGE counts (Fig.3)
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


# --------------------------------------------------------------------
# plot_density_scatter trimming
# --------------------------------------------------------------------
plot_density_scatter <- function(
    df, title = "", outfile, out_dir,
    use_density = TRUE, n_grid = 200,
    trim_prop = 0,
    diff_mode = c("lowess_resid"),
    pseudocount = 1
) {
  diff_mode <- match.arg(diff_mode)
  
  if (nrow(df) == 0) {
    message("No rows for ", outfile)
    return(invisible(NULL))
  }
  
  df$logC <- log10(df$cage_norm + 1)
  df$logO <- log10(df$other_norm + 1)
  df$density <- NA_real_
  
  if (use_density) {
    keep <- is.finite(df$logC) & is.finite(df$logO)
    dens <- with(df[keep, ], MASS::kde2d(logC, logO, n = n_grid))
    ix <- findInterval(df$logC[keep], dens$x)
    iy <- findInterval(df$logO[keep], dens$y)
    df$density[keep] <- log10(dens$z[cbind(ix, iy)] + 1e-8)
  }
  
  cor_labels <- do.call(rbind, lapply(split(df, df$method), function(sub) {
    x_raw <- sub$cage_norm
    y_raw <- sub$other_norm
    ok    <- is.finite(x_raw) & is.finite(y_raw)
    
    if (sum(ok) < 2) {
      return(data.frame(
        method = unique(sub$method),
        x = NA_real_, y = NA_real_, label = "N < 2",
        pearson_r = NA_real_, spearman_rho = NA_real_, ccc_val = NA_real_,
        bias_str = NA_character_, ba_corr = NA_real_,
        stringsAsFactors = FALSE
      ))
    }
    
    xv <- x_raw[ok]
    yv <- y_raw[ok]
    lx <- xv
    ly <- yv
    
    if (trim_prop > 0) {
      if (diff_mode == "abs_log2") {
        diffs <- abs(ly - lx)
      } else {
        fit  <- stats::lowess(lx, ly, f = 0.5)
        yhat <- approx(fit$x, fit$y, xout = lx, rule = 2)$y
        diffs <- abs(ly - yhat)
      }
      cutoff <- stats::quantile(diffs, 1 - trim_prop, na.rm = TRUE)
      keep_trim <- diffs <= cutoff
      lx_t <- lx[keep_trim]
      ly_t <- ly[keep_trim]
    } else {
      keep_trim <- rep(TRUE, length(lx))
      lx_t <- lx
      ly_t <- ly
      cutoff <- NA_real_
    }
    
    if (sum(keep_trim) < 2) {
      pearson_r <- spearman_rho <- ccc_val <- ba_corr <- NA_real_
      bias_str <- "NA"
    } else {
      pearson_r    <- round(cor(lx_t, ly_t, method = "pearson"), 3)
      spearman_rho <- round(cor(lx_t, ly_t, method = "spearman"), 3)
      ccc_val <- tryCatch({
        DescTools::CCC(lx_t, ly_t)$rho.c$est
      }, error = function(e) NA_real_)
      if (is.finite(ccc_val)) ccc_val <- round(ccc_val, 3)
      
      diff_xy   <- (ly_t - lx_t)
      mean_xy   <- (ly_t + lx_t) / 2
      mean_diff <- mean(diff_xy)
      bias_str  <- sprintf("%.3f", round(mean_diff, 3))
      ba_corr <- suppressWarnings(cor(mean_xy, diff_xy, use = "complete.obs"))
      ba_corr <- round(ifelse(is.na(ba_corr), NA, ba_corr), 3)
    }
    
    label <- paste0(
      "Pearson's r = ", pearson_r, "\n",
      "Spearman's rho = ", spearman_rho, "\n",
      "Concordance r = ", ifelse(is.na(ccc_val), "NA", ccc_val), "\n",
      "Mean bias = ", bias_str, "\n",
      "BA-corr = ", ifelse(is.na(ba_corr), "NA", ba_corr)
    )
    
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
      spearman_rho = spearman_rho,
      ccc_val = ccc_val,
      bias_str = bias_str,
      ba_corr = ba_corr,
      stringsAsFactors = FALSE
    )
  }))
  
  lims <- range(c(df$logC, df$logO), finite = TRUE)
  
  p <- ggplot(df, aes(x = logC, y = logO)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey50") +
    scale_x_continuous(limits = lims) +
    scale_y_continuous(limits = lims) +
    facet_wrap(~method, nrow = 1) +
    labs(x = "CAGE (log10 normalized)", y = "RNA-seq based method (log10 normalized)") +
    theme_bw(base_size = 11) +
    theme(plot.title = element_blank())
  
  if (use_density) {
    p <- p +
      geom_point(aes(color = density), size = 1, alpha = 0.85) +
      scale_color_viridis_c(option = "A", name = "Density")
  } else {
    p <- p + geom_point(alpha = 0.5, size = 0.7)
  }
  
  p <- p + geom_label(
    data = cor_labels,
    aes(x = x, y = y, label = label),
    hjust = 0, vjust = 1, size = 3.2,
    fill = "white", color = "black",
    label.padding = unit(0.1, "lines"),
    label.r = unit(0.15, "lines"),
    fontface = "plain", alpha = 1
  )
  
  ggsave(file.path(out_dir, outfile), p, width = 5, height = 4)
  
  cor_results <- cor_labels[, c("method", "pearson_r", "spearman_rho", "ccc_val", "bias_str", "ba_corr")]
  cor_results
}

# ============================================================
# RPKM
# ============================================================

# RPKM conversion given counts matrix and length vector
to_rpkm_by_len <- function(counts_mat, len_bp_named) {
  counts_mat <- as.matrix(counts_mat)
  ids <- intersect(rownames(counts_mat), names(len_bp_named))
  if (length(ids) == 0L) stop("RPKM: no overlapping promoter IDs with length vector.")
  x <- counts_mat[ids, , drop = FALSE]
  L <- len_bp_named[ids]; L[!is.finite(L) | L <= 0] <- 1
  lib <- colSums(x);      lib[!is.finite(lib) | lib <= 0] <- 1
  t(t(x) / lib) * 1e9 / L
}


# ============================================================
# BOOTSTRAP + STRATIFIED
# ============================================================

calc_corr_metrics <- function(x, y) {
  # keep only finite paired values
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]
  y <- y[ok]
  
  # return NA metrics if too few observations remain
  if (length(x) < 2) {
    return(data.frame(
      pearson_r = NA_real_,
      spearman_rho = NA_real_,
      ccc_val = NA_real_,
      ba_corr = NA_real_
    ))
  }
  
  # compute pairwise difference for BA analysis
  diff_xy <- y - x
  
  # compute pairwise mean for BA analysis
  mean_xy <- (y + x) / 2
  
  # compute correlation between mean and difference
  ba_corr <- suppressWarnings(cor(mean_xy, diff_xy, use = "complete.obs"))
  
  # return all correlation / agreement metrics in one row
  data.frame(
    pearson_r = suppressWarnings(cor(x, y, method = "pearson")),
    spearman_rho = suppressWarnings(cor(x, y, method = "spearman")),
    ccc_val = tryCatch(DescTools::CCC(x, y)$rho.c$est, error = function(e) NA_real_),
    ba_corr = ba_corr
  )
}

bootstrap_corr_long <- function(df, n_boot = 1000, seed = 1) {
  # set seed for reproducible bootstrap resampling
  set.seed(seed)
  
  # preallocate list to store bootstrap results
  out <- vector("list", n_boot)
  
  # iterate over bootstrap replicates
  for (b in seq_len(n_boot)) {
    # sample row indices with replacement
    idx <- sample(seq_len(nrow(df)), size = nrow(df), replace = TRUE)
    
    # calculate metrics on bootstrap-resampled data
    met <- calc_corr_metrics(df$cage_norm[idx], df$other_norm[idx])
    
    # store bootstrap result together with method and cell line labels
    out[[b]] <- data.frame(
      bootstrap_id = b,
      method = unique(df$method),
      cell_line = unique(df$cell_line),
      pearson_r = met$pearson_r,
      spearman_rho = met$spearman_rho,
      ccc_val = met$ccc_val,
      ba_corr = met$ba_corr
    )
  }
  
  # combine all bootstrap results into one data.frame
  do.call(rbind, out)
}

summarize_quantile_counts <- function(shared_df, out_dir, prefix) {
  # create output directory if needed
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  # convert to data.table for grouped summaries
  shared_df <- as.data.table(shared_df)
  
  # summarize CAGE expression distribution within each cell line and quantile
  cage_summary <- shared_df[
    ,
    .(
      n_promoters = uniqueN(promoterId),
      min = min(cage_norm, na.rm = TRUE),
      q1 = quantile(cage_norm, 0.25, na.rm = TRUE),
      median = median(cage_norm, na.rm = TRUE),
      mean = mean(cage_norm, na.rm = TRUE),
      q3 = quantile(cage_norm, 0.75, na.rm = TRUE),
      max = max(cage_norm, na.rm = TRUE)
    ),
    by = .(cell_line, quantile)
  ]
  
  # save per-quantile CAGE summary table
  fwrite(
    cage_summary,
    file.path(out_dir, paste0(prefix, "_quantile_cage_counts_summary.tsv")),
    sep = "\t"
  )
  
  # print CAGE summary for quick inspection in console
  print(cage_summary)
  
  # summarize method-specific expression distribution within each cell line and quantile
  other_summary <- shared_df[
    ,
    .(
      n_promoters = uniqueN(promoterId),
      min = min(other_norm, na.rm = TRUE),
      q1 = quantile(other_norm, 0.25, na.rm = TRUE),
      median = median(other_norm, na.rm = TRUE),
      mean = mean(other_norm, na.rm = TRUE),
      q3 = quantile(other_norm, 0.75, na.rm = TRUE),
      max = max(other_norm, na.rm = TRUE)
    ),
    by = .(cell_line, method, quantile)
  ]
  
  # save per-method per-quantile summary table
  fwrite(
    other_summary,
    file.path(out_dir, paste0(prefix, "_quantile_other_counts_summary.tsv")),
    sep = "\t"
  )
  
  # visualize CAGE expression distribution across shared quantiles
  p <- ggplot(shared_df, aes(x = quantile, y = cage_norm)) +
    geom_boxplot() +
    facet_wrap(~cell_line) +
    scale_y_log10() +
    theme_bw(base_size = 12) +
    labs(
      x = "CAGE expression quantile",
      y = "CAGE normalized counts (log10)"
    )
  
  # save CAGE quantile boxplot
  ggsave(
    file.path(out_dir, paste0(prefix, "_quantile_cage_counts_boxplot.pdf")),
    p,
    width = 6,
    height = 4
  )
}

make_shared_cage_quantiles <- function(set2, n_quantiles = 4) {
  # initialize list to store per-cell-line shared promoter data
  out_list <- list()
  
  # index for sequential list filling
  idx <- 1
  
  # loop over cell lines independently
  for (cl in unique(set2$cell_line)) {
    # subset data for current cell line
    df_cl <- set2[set2$cell_line == cl, , drop = FALSE]
    
    # split promoter IDs by method within the current cell line
    promoter_sets <- split(df_cl$promoterId, df_cl$method)
    
    # skip if fewer than two methods are available
    if (length(promoter_sets) < 2) {
      message("Skipping ", cl, ": fewer than 2 methods available.")
      next
    }
    
    # keep only promoters shared across all methods
    common_ids <- Reduce(intersect, promoter_sets)
    
    # skip if too few shared promoters to define requested quantiles
    if (length(common_ids) < n_quantiles) {
      message("Skipping shared quantiles for ", cl, ": too few common promoters.")
      next
    }
    
    # extract one unique CAGE value per shared promoter
    cage_ref <- unique(df_cl[df_cl$promoterId %in% common_ids, c("promoterId", "cage_norm")])
    
    # retain only finite CAGE values for quantile assignment
    cage_ref <- cage_ref[is.finite(cage_ref$cage_norm), , drop = FALSE]
    
    # skip if too few finite values remain
    if (nrow(cage_ref) < n_quantiles) {
      message("Skipping shared quantiles for ", cl, ": too few finite CAGE values.")
      next
    }
    
    # define quantile cut points from shared CAGE values
    qs <- quantile(
      cage_ref$cage_norm,
      probs = seq(0, 1, length.out = n_quantiles + 1),
      na.rm = TRUE
    )
    
    # force full coverage at both ends
    qs[1] <- -Inf
    qs[length(qs)] <- Inf
    
    # assign each shared promoter to a CAGE-based quantile bin
    cage_ref$quantile <- cut(
      cage_ref$cage_norm,
      breaks = qs,
      include.lowest = TRUE,
      labels = paste0("Q", seq_len(n_quantiles))
    )
    
    # print quick summary of shared promoter counts per quantile
    cat("Cell line:", cl, " common promoters:", length(common_ids), "\n")
    print(table(cage_ref$quantile))
    
    # merge shared quantile labels back to all method rows for shared promoters
    df_cl2 <- merge(
      df_cl[df_cl$promoterId %in% common_ids, , drop = FALSE],
      cage_ref[, c("promoterId", "quantile")],
      by = "promoterId",
      all.x = TRUE
    )
    
    # store current cell-line result
    out_list[[idx]] <- df_cl2
    idx <- idx + 1
  }
  
  # return NULL if no valid cell line produced output
  if (length(out_list) == 0) return(NULL)
  
  # combine all cell-line results into one data.frame
  do.call(rbind, out_list)
}

shared_stratified_corr <- function(set2, n_quantiles = 4) {
  # construct shared promoter dataset with CAGE-defined quantiles
  shared_df <- make_shared_cage_quantiles(set2, n_quantiles = n_quantiles)
  
  # return NULL if no valid shared dataset is available
  if (is.null(shared_df) || nrow(shared_df) == 0) return(NULL)
  
  # split by cell line, method, and quantile for stratified metric calculation
  split_df <- split(
    shared_df,
    list(shared_df$cell_line, shared_df$method, shared_df$quantile),
    drop = TRUE
  )
  
  # calculate metrics within each stratified subset
  out <- lapply(split_df, function(sub) {
    # skip subsets with too few rows
    if (nrow(sub) < 2) return(NULL)
    
    # compute metrics comparing method values to CAGE values
    met <- calc_corr_metrics(sub$cage_norm, sub$other_norm)
    
    # return one summary row for this stratum
    data.frame(
      cell_line = unique(sub$cell_line),
      method = unique(sub$method),
      quantile = unique(as.character(sub$quantile)),
      n_promoters = nrow(sub),
      pearson_r = met$pearson_r,
      spearman_rho = met$spearman_rho,
      ccc_val = met$ccc_val,
      ba_corr = met$ba_corr
    )
  })
  
  # remove empty subsets
  out <- out[!vapply(out, is.null, logical(1))]
  
  # return NULL if all subsets were empty
  if (length(out) == 0) return(NULL)
  
  # combine all stratified results
  do.call(rbind, out)
}


plot_bootstrap_boxplot <- function(boot_long, outfile, plot_title = NULL) {
  # enforce consistent facet and legend ordering
  boot_long$cell_line <- factor(boot_long$cell_line,
                               levels = c("Healthy", "Failure"))
  boot_long$metric <- factor(boot_long$metric, levels = c("Pearson", "Spearman", "CCC", "BA-corr"))
  boot_long$method <- factor(boot_long$method, levels = c("proActiv","Salmon", "DEXSeq"))
  
  # draw bootstrap distribution of metrics across methods and cell lines
  p <- ggplot(boot_long, aes(x = method, y = value, fill = method)) +
    geom_boxplot(width = 0.7, outlier.size = 0.3, alpha = 0.85) +
    scale_fill_manual(values = cols_method) +
    facet_grid(cell_line ~ metric, scales = "fixed") +
    labs(
      title = plot_title,
      x = NULL,
      y = "Bootstrap correlation"
    ) +
    theme_rev() +
    theme(
      legend.position = "none",
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank()
    )
  
  # save bootstrap boxplot figure
  ggsave(outfile, p, width = 9.5, height = 5.2)
}

plot_stratified_lineplot <- function(strat_long, outfile, plot_title = NULL) {
  # enforce consistent facet, line, and x-axis ordering
  strat_long$cell_line <- factor(strat_long$cell_line,
                               levels = c("Healthy", "Failure"))
  strat_long$metric <- factor(strat_long$metric, levels = c("Pearson", "Spearman", "CCC", "BA-corr"))
  strat_long$method <- factor(strat_long$method, levels = c("proActiv","Salmon", "DEXSeq"))
  strat_long$quantile <- factor(strat_long$quantile, levels = c("Q1", "Q2", "Q3", "Q4"))
  
  # draw per-quantile trend lines for each method and metric
  p <- ggplot(strat_long, aes(x = quantile, y = value, color = method, group = method)) +
    geom_line(linewidth = 1.2, alpha = 0.7) +
    geom_point(size = 2.6) +
    scale_color_manual(values = cols_method) +
    facet_grid(cell_line ~ metric, scales = "fixed") +
    labs(
      title = plot_title,
      x = "CAGE expression quantile",
      y = "Correlation"
    ) +
    theme_rev() +
    theme(
      legend.position = "none"
    )
  
  # save quantile-stratified lineplot
  ggsave(outfile, p, width = 9.5, height = 5.2)
}

run_bootstrap_and_stratified <- function(set2, out_dir, prefix, n_boot = 1000, seed = 1) {
  # ensure output directory exists
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  # bootstrap correlation metrics separately within each cell line and method
  boot_list <- lapply(
    split(set2, list(set2$cell_line, set2$method), drop = TRUE),
    function(df_sub) {
      # skip subsets that are too small for resampling-based correlation analysis
      if (nrow(df_sub) < 2) return(NULL)
      
      # run bootstrap metric estimation on the current subset
      bootstrap_corr_long(df_sub, n_boot = n_boot, seed = seed)
    }
  )
  
  # drop empty bootstrap subsets
  boot_list <- boot_list[!vapply(boot_list, is.null, logical(1))]
  
  # stop early if no bootstrap results were generated
  if (length(boot_list) == 0) {
    message("No bootstrap results for ", prefix)
    return(invisible(NULL))
  }
  
  # combine all bootstrap results
  boot_df <- do.call(rbind, boot_list)
  
  # reshape bootstrap results to long format for plotting
  boot_long <- reshape(
    boot_df,
    varying = c("pearson_r", "spearman_rho", "ccc_val", "ba_corr"),
    times = c("Pearson", "Spearman", "CCC", "BA-corr"),
    v.names = "value",
    timevar = "metric",
    direction = "long"
  )
  
  # clean row names after reshape
  rownames(boot_long) <- NULL
  
  # save wide bootstrap table
  write.csv(
    boot_df,
    file = file.path(out_dir, paste0(prefix, "_bootstrap_metrics_wide.csv")),
    row.names = FALSE
  )
  
  # save long bootstrap table
  write.csv(
    boot_long,
    file = file.path(out_dir, paste0(prefix, "_bootstrap_metrics_long.csv")),
    row.names = FALSE
  )
  
  # generate bootstrap boxplot
  plot_bootstrap_boxplot(
    boot_long = boot_long,
    outfile = file.path(out_dir, paste0(prefix, "_bootstrap_boxplot.pdf")),
    plot_title = NULL
  )
  
  # create shared-promoter quantile dataset based on CAGE expression
  shared_df <- make_shared_cage_quantiles(set2, n_quantiles = 4)
  
  # if no valid shared dataset exists, return bootstrap results only
  if (is.null(shared_df) || nrow(shared_df) == 0) {
    message("No shared quantile dataset for ", prefix)
    return(invisible(list(
      bootstrap = boot_df,
      stratified = NULL
    )))
  }
  
  # save descriptive summaries for shared CAGE quantiles
  summarize_quantile_counts(shared_df, out_dir, prefix)
  
  # compute correlation metrics within each shared CAGE quantile
  strat_df <- shared_stratified_corr(set2, n_quantiles = 4)
  
  # if no stratified results were produced, return bootstrap results only
  if (is.null(strat_df) || nrow(strat_df) == 0) {
    message("No stratified results for ", prefix)
    return(invisible(list(
      bootstrap = boot_df,
      stratified = NULL
    )))
  }
  
  # reshape stratified results to long format for plotting
  strat_long <- reshape(
    strat_df,
    varying = c("pearson_r", "spearman_rho", "ccc_val", "ba_corr"),
    times = c("Pearson", "Spearman", "CCC", "BA-corr"),
    v.names = "value",
    timevar = "metric",
    direction = "long"
  )
  
  # clean row names after reshape
  rownames(strat_long) <- NULL
  
  # save wide stratified table
  write.csv(
    strat_df,
    file = file.path(out_dir, paste0(prefix, "_stratified_metrics_wide.csv")),
    row.names = FALSE
  )
  
  # save long stratified table
  write.csv(
    strat_long,
    file = file.path(out_dir, paste0(prefix, "_stratified_metrics_long.csv")),
    row.names = FALSE
  )
  
  # generate stratified lineplot
  plot_stratified_lineplot(
    strat_long = strat_long,
    outfile = file.path(out_dir, paste0(prefix, "_stratified_lineplot.pdf")),
    plot_title = NULL
  )
  
  # invisibly return both bootstrap and stratified results for downstream use
  invisible(list(
    bootstrap = boot_df,
    stratified = strat_df
  ))
}
# ============================================================
# MAIN ANALYSIS FUNCTION
# ============================================================

combine.counts.for.each.dataset <- function(
    cage_counts_ori,
    salmon_counts_ori,
    dexseq_counts_ori,
    proactiv_counts_ori,
    cell_lines,
    out_dir,
    filename = "counts",
    counts_cutoff = 2,
    trim_prop = 0
) {
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  cat("combine.counts.for.each.dataset trim_prop:", trim_prop, "\n")
  
  # 1. ALL PROMOTERS
  cage_counts     <- cage_counts_ori
  salmon_counts   <- salmon_counts_ori
  dexseq_counts   <- dexseq_counts_ori
  proactiv_counts <- proactiv_counts_ori
  
  methods <- list(
    Salmon   = salmon_counts,
    DEXSeq   = dexseq_counts,
    proActiv = proactiv_counts
  )
  
  all_df <- do.call(
    rbind,
    lapply(cell_lines, function(cl) {
      do.call(rbind, lapply(names(methods), function(m) {
        message("▶ ", cl, " vs CAGE — ", m)
        make_df_one_method_one_cl(cage_counts, methods[[m]], cl, m)
      }))
    })
  )
  
  cnt_ok <- with(all_df, cage_norm > counts_cutoff & other_norm > counts_cutoff)
  set2   <- subset(all_df, cnt_ok)
  
  for (cl in cell_lines) {
    cor_list <- list()
    
    for (m in names(methods)) {
      df <- subset(set2, cell_line == cl & method == m)
      if (nrow(df) == 0) next
      
      ttl <- sprintf("%s - %s (counts %s)", cl, m, counts_cutoff)
      of  <- sprintf("scatter_%s_all_%s_%s.pdf", filename, tolower(cl), tolower(m))
      
      cor_list[[m]] <- plot_density_scatter(
        df = df,
        title = ttl,
        outfile = of,
        out_dir = out_dir,
        use_density = TRUE,
        trim_prop = trim_prop
      )
    }
    
    cor_df <- do.call(rbind, cor_list)
    cor_df$cell_line <- cl
    
    write.csv(
      t(cor_df),
      file = file.path(out_dir, sprintf("scatter_%s_all_correlations_%s.csv", filename, tolower(cl))),
      row.names = TRUE,
      col.names = FALSE
    )
  }
  
  # robustness
  run_bootstrap_and_stratified(
    set2 = set2,
    out_dir = file.path(out_dir, "robustness"),
    prefix = paste0(filename, "_all"),
    n_boot = 1000,
    seed = 1
  )
  
  # 2. INTRONLESS
  cage_counts     <- cage_counts_ori[rownames(cage_counts_ori) %in% intronless_ids, , drop = FALSE]
  salmon_counts   <- salmon_counts_ori[rownames(salmon_counts_ori) %in% intronless_ids, , drop = FALSE]
  dexseq_counts   <- dexseq_counts_ori[rownames(dexseq_counts_ori) %in% intronless_ids, , drop = FALSE]
  proactiv_counts <- proactiv_counts_ori[rownames(proactiv_counts_ori) %in% intronless_ids, , drop = FALSE]
  
  methods <- list(
    Salmon = salmon_counts,
    DEXSeq = dexseq_counts
  )
  
  all_df <- do.call(
    rbind,
    lapply(cell_lines, function(cl) {
      do.call(rbind, lapply(names(methods), function(m) {
        message("▶ ", cl, " vs CAGE — ", m)
        make_df_one_method_one_cl(cage_counts, methods[[m]], cl, m)
      }))
    })
  )
  
  cnt_ok <- with(all_df, cage_norm > counts_cutoff & other_norm > counts_cutoff)
  set2   <- subset(all_df, cnt_ok)
  
  for (cl in cell_lines) {
    cor_list <- list()
    
    for (m in names(methods)) {
      df <- subset(set2, cell_line == cl & method == m)
      if (nrow(df) == 0) {
        message("Skipping ", cl, " - ", m, " due to empty data.")
        next
      }
      
      ttl <- sprintf("%s - %s (no intron, counts %s)", cl, m, counts_cutoff)
      of  <- sprintf("scatter_%s_nointron_%s_%s.pdf", filename, tolower(cl), tolower(m))
      
      cor_list[[m]] <- plot_density_scatter(
        df = df,
        title = ttl,
        outfile = of,
        out_dir = out_dir,
        use_density = TRUE,
        trim_prop = trim_prop
      )
    }
    
    cor_df <- do.call(rbind, cor_list)
    cor_df$cell_line <- cl
    
    write.csv(
      t(cor_df),
      file = file.path(out_dir, sprintf("scatter_%s_nointron_correlations_%s.csv", filename, tolower(cl))),
      row.names = TRUE,
      col.names = FALSE
    )
  }
  
  run_bootstrap_and_stratified(
    set2 = set2,
    out_dir = file.path(out_dir, "robustness"),
    prefix = paste0(filename, "_nointron"),
    n_boot = 1000,
    seed = 1
  )


  # 3. WITH INTRON
  cage_counts     <- cage_counts_ori[rownames(cage_counts_ori) %in% intron_ids, , drop = FALSE]
  salmon_counts   <- salmon_counts_ori[rownames(salmon_counts_ori) %in% intron_ids, , drop = FALSE]
  dexseq_counts   <- dexseq_counts_ori[rownames(dexseq_counts_ori) %in% intron_ids, , drop = FALSE]
  proactiv_counts <- proactiv_counts_ori[rownames(proactiv_counts_ori) %in% intron_ids, , drop = FALSE]
  
  methods <- list(
    Salmon   = salmon_counts,
    DEXSeq   = dexseq_counts,
    proActiv = proactiv_counts
  )
  
  all_df <- do.call(
    rbind,
    lapply(cell_lines, function(cl) {
      do.call(rbind, lapply(names(methods), function(m) {
        message("▶ ", cl, " vs CAGE — ", m)
        make_df_one_method_one_cl(cage_counts, methods[[m]], cl, m)
      }))
    })
  )
  
  cnt_ok <- with(all_df, cage_norm > counts_cutoff & other_norm > counts_cutoff)
  set2   <- subset(all_df, cnt_ok)
  
  for (cl in cell_lines) {
    cor_list <- list()
    
    for (m in names(methods)) {
      df <- subset(set2, cell_line == cl & method == m)
      if (nrow(df) == 0) {
        message("Skipping ", cl, " - ", m, " due to empty data.")
        next
      }
      
      ttl <- sprintf("%s - %s (with intron, counts >= %s)", cl, m, counts_cutoff)
      of  <- sprintf("scatter_%s_withintron_%s_%s.pdf", filename, tolower(cl), tolower(m))
      
      cor_list[[m]] <- plot_density_scatter(
        df = df,
        title = ttl,
        outfile = of,
        out_dir = out_dir,
        use_density = TRUE,
        trim_prop = trim_prop
      )
    }
    
    cor_df <- do.call(rbind, cor_list)
    cor_df$cell_line <- cl
    
    write.csv(
      t(cor_df),
      file = file.path(out_dir, sprintf("scatter_%s_withintron_correlations_%s.csv", filename, tolower(cl))),
      row.names = TRUE,
      col.names = FALSE
    )
  }

  run_bootstrap_and_stratified(
    set2 = set2,
    out_dir = file.path(out_dir, "robustness"),
    prefix = paste0(filename, "_withintron"),
    n_boot = 1000,
    seed = 1
  )
}
