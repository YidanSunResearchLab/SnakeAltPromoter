
# ---------------------------
# Plot: feature counts and no intron
# ---------------------------

library(DESeq2)
library(ggplot2)
library(MASS)
library(viridis)
library(DescTools)   


##Gene annotation to get intronless promoters
coord <- promoter_anno@promoterCoordinates
intronless_ids <- as.character(coord$promoterId[!coord$internalPromoter %in% TRUE &
                                            lengths(coord$intronId) == 0])
intron_ids <- as.character(coord$promoterId[!coord$internalPromoter %in% TRUE &
                                            !lengths(coord$intronId) == 0])
all_ids <- as.character(coord$promoterId[!coord$internalPromoter %in% TRUE])



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
  df, title = "", outfile, use_density = TRUE, n_grid = 200,
  trim_prop = 0,                      # fraction to trim
  diff_mode = c("lowess_resid"), #,"abs_log2"
  pseudocount = 1
) {
  diff_mode <- match.arg(diff_mode)

  if (nrow(df) == 0) {
    message("No rows for ", outfile); return(invisible(NULL))
  }

  # Coordinates for plotting (always show all points)
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

  # ----------------------------
  # Correlations per method WITH TRIMMING
  # ----------------------------
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

    # vectors for trimming/correlation
    xv <- x_raw[ok]; yv <- y_raw[ok]
    lx <- xv
    ly <- yv

    # difference metric
    print(paste0("plot_density_scatter trimming", trim_prop))
    if (trim_prop > 0) {
      if (diff_mode == "abs_log2") {
        diffs <- abs(ly - lx)
      } else { # lowess_resid on log2-scale
        fit  <- stats::lowess(lx, ly, f = 0.5)
        yhat <- approx(fit$x, fit$y, xout = lx, rule = 2)$y
        diffs <- abs(ly - yhat)
      }
      cutoff <- stats::quantile(diffs, 1 - trim_prop, na.rm = TRUE)
      keep_trim <- diffs <= cutoff
      lx_t <- lx[keep_trim]; ly_t <- ly[keep_trim]
    } else {
      keep_trim <- rep(TRUE, length(lx))
      lx_t <- lx; ly_t <- ly
      cutoff <- NA_real_
    }

    # correlations on trimmed set (log2 scale)
    if (sum(keep_trim) < 2) {
      pearson_r <- spearman_rho <- ccc_val <- ba_corr <- NA_real_
      bias_str <- "NA"
    } else {
      pearson_r    <- round(cor(lx_t, ly_t, method = "pearson"), 3)
      spearman_rho <- round(cor(lx_t, ly_t, method = "spearman"), 3)
      ccc_val <- tryCatch({
        CCC(lx_t, ly_t)$rho.c$est
      }, error = function(e) NA_real_)
      if (is.finite(ccc_val)) ccc_val <- round(ccc_val, 3)

      diff_xy   <- (ly_t - lx_t)
      mean_xy   <- (ly_t + lx_t) / 2
      mean_diff <- mean(diff_xy)
      bias_str  <- sprintf("%.3f", round(mean_diff, 3))
      ba_corr <- suppressWarnings(cor(mean_xy, diff_xy, use = "complete.obs"))
      ba_corr <- round(ifelse(is.na(ba_corr), NA, ba_corr), 3)
    }

    # Label text (note: explicitly state trimming)
    label <- paste0(
      "Pearson's r = ", pearson_r, "\n",
      "Spearman's rho = ", spearman_rho, "\n",
      "Concordance r = ", ifelse(is.na(ccc_val), "NA", ccc_val), "\n",
      "Mean bias = ", bias_str, "\n",
      "BA bias = ", ifelse(is.na(ba_corr), "NA", ba_corr)
    )

    # Position label
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

  # Coordinate range for plotting (all points)
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
    p <- p + geom_point(aes(color = density), size = 1, alpha = 0.85) +
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

  # Return correlations + kept counts
  cor_results <- cor_labels[, c("method","pearson_r","spearman_rho","ccc_val","bias_str","ba_corr")]
  return(cor_results)
}


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


combine.counts.for.each.dataset <- function(cage_counts_ori, salmon_counts_ori, dexseq_counts_ori, proactiv_counts_ori, cell_lines, out_dir,filename="counts", counts_cutoff=2, trim_prop=0) {
    dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
    print(paste0("combine.counts.for.each.dataset trim_prop: ",trim_prop))
    #-----------------------------------------------------------------------------------------------------------
    #Plot scatterplot of all promoter counts of Salmon and DEXseq vs CAGE counts
    cage_counts = cage_counts_ori
    salmon_counts = salmon_counts_ori
    dexseq_counts = dexseq_counts_ori
    proactiv_counts = proactiv_counts_ori

    # Combine different methods
    methods <- list(
      Salmon   = salmon_counts_ori,
      DEXSeq   = dexseq_counts_ori,
      proActiv = proactiv_counts_ori
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

    ##Remove low counts
    cnt_ok <- with(all_df, cage_norm > counts_cutoff & other_norm > counts_cutoff)
    set2   <- subset(all_df, cnt_ok)

    ##Plot scatterplot
    for (cl in cell_lines) {
      cor_list <- list()
      for (m in names(methods)) {
        df <- subset(set2, cell_line == cl & method == m)
        
        if (nrow(df) == 0) next
        ttl <- sprintf("%s ‒ %s (counts  %s)", cl, m, counts_cutoff)
        of  <- sprintf("scatter_%s_all_%s_%s.pdf", filename, tolower(cl), tolower(m))
        cor_list[[m]] <- plot_density_scatter(df, ttl, of, use_density = TRUE, trim_prop=trim_prop)
      }
      cor_df <- do.call(rbind, cor_list)
      cor_df$cell_line <- cl  # Add cell line column
      write.csv(t(cor_df), file = file.path(out_dir, sprintf("scatter_%s_all_correlations_%s.csv", filename, tolower(cl))), row.names = TRUE, col.names = FALSE)
    }


    #--------------------------------------------------------
    #Plot scatterplot of intronless promoter counts of Proactiv Salmon and DEXseq vs CAGE counts (Fig.3)
    # Keep only intronless and remove internal
    cage_counts     <- cage_counts_ori[rownames(cage_counts_ori)     %in% intronless_ids, , drop = FALSE]
    salmon_counts   <- salmon_counts_ori[rownames(salmon_counts_ori) %in% intronless_ids, , drop = FALSE]
    dexseq_counts   <- dexseq_counts_ori[rownames(dexseq_counts_ori) %in% intronless_ids, , drop = FALSE]
    proactiv_counts <- proactiv_counts_ori[rownames(proactiv_counts_ori) %in% intronless_ids, , drop = FALSE]

    # Combine different methods
    methods <- list(Salmon   = salmon_counts,DEXSeq   = dexseq_counts)
    all_df <- do.call(
      rbind,
      lapply(cell_lines, function(cl) {
        do.call(rbind, lapply(names(methods), function(m) {
          message("▶ ", cl, "  vs CAGE  —  ", m)
          make_df_one_method_one_cl(cage_counts, methods[[m]], cl, m)
        }))
      })
    )

    ##Remove low counts
    cnt_ok <- with(all_df, cage_norm > counts_cutoff & other_norm > counts_cutoff)
    set2   <- subset(all_df, cnt_ok)

    ##Plot scatterplot
    for (cl in cell_lines) {
      cor_list <- list()
      for (m in names(methods)) {
        df <- subset(set2, cell_line == cl & method == m)
        if (nrow(df) == 0) {
            message("Skipping ", cl, " - ", m, " due to empty data.")
            next
        }
        ttl <- sprintf("%s ‒ %s (no intron, counts %s)", cl, m, counts_cutoff)
        of  <- sprintf("scatter_%s_nointron_%s_%s.pdf", filename, tolower(cl), tolower(m))
        cor_list[[m]] <- plot_density_scatter(df, ttl, of, use_density = TRUE, trim_prop=trim_prop)
      }
      cor_df <- do.call(rbind, cor_list)
      cor_df$cell_line <- cl  # Add cell line column
      write.csv(t(cor_df), file = file.path(out_dir, sprintf("scatter_%s_nointron_correlations_%s.csv", filename, tolower(cl))), row.names = TRUE, col.names = FALSE)
    }

    #---------------------------------------------------------------------------
    #Plot scatterplot of all promoter counts of Salmon and DEXseq vs CAGE counts
    # Keep with intron and no internal
    cage_counts     <- cage_counts_ori[rownames(cage_counts_ori)     %in% intron_ids, , drop = FALSE]
    salmon_counts   <- salmon_counts_ori[rownames(salmon_counts_ori) %in% intron_ids, , drop = FALSE]
    dexseq_counts   <- dexseq_counts_ori[rownames(dexseq_counts_ori) %in% intron_ids, , drop = FALSE]
    proactiv_counts <- proactiv_counts_ori[rownames(proactiv_counts_ori) %in% intron_ids, , drop = FALSE]

    # Combine different methods
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

    ##Remove low counts
    cnt_ok <- with(all_df, cage_norm > counts_cutoff & other_norm > counts_cutoff)
    set2   <- subset(all_df, cnt_ok)

    ##Plot scatterplot
    for (cl in cell_lines) {
      cor_list <- list()
      for (m in names(methods)) {
        df <- subset(set2, cell_line == cl & method == m)
        if (nrow(df) == 0) {
            message("Skipping ", cl, " - ", m, " due to empty data.")
            next
        }
        ttl <- sprintf("%s ‒ %s (with intron, counts ≥ %s)", cl, m, counts_cutoff)
        of  <- sprintf("scatter_%s_withintron_%s_%s.pdf", filename, tolower(cl), tolower(m))
        cor_list[[m]] <- plot_density_scatter(df, ttl, of, use_density = TRUE, trim_prop=trim_prop)
      }
      cor_df <- do.call(rbind, cor_list)
      cor_df$cell_line <- cl  # Add cell line column
      write.csv(t(cor_df), file = file.path(out_dir, sprintf("scatter_%s_withintron_correlations_%s.csv", filename, tolower(cl))), row.names = TRUE, col.names = FALSE)
    }
}


summarize_correlations <- function(base_dir, settings, cell_lines, out_tsv,
                                  setting_col = "setting") {

  all_res <- list()

  for (s in settings) {
    for (cl in cell_lines) {

      f <- file.path(base_dir, s,
                     paste0("scatter_rpkm_withintron_correlations_", cl, ".csv"))
      if (!file.exists(f)) next

      dt <- fread(f, header = TRUE)
      setnames(dt, 1, "metric")
      dt <- dt[metric != "method"]

      # ---- normalize metric names: ccc_val -> ccc ----
      dt[metric == "ccc_val", metric := "ccc"]

      # ---- keep metrics (now ccc is available) ----
      keep_metrics <- c("pearson_r","spearman_rho","ccc","bias_str","ba_corr")
      dt <- dt[metric %in% keep_metrics]

      dt_long <- melt(
        dt,
        id.vars = "metric",
        variable.name = "analysis_method",
        value.name = "value"
      )

      dt_wide <- dcast(
        dt_long,
        analysis_method ~ metric,
        value.var = "value"
      )

      # ensure numeric columns (some files have quotes/spaces)
      for (col in keep_metrics) {
        if (!col %in% names(dt_wide)) dt_wide[, (col) := NA_real_]
        dt_wide[, (col) := as.numeric(gsub("[[:space:]\"]+", "", as.character(get(col))))]
      }

      dt_wide[, (setting_col) := s]
      dt_wide[, cell_line := toupper(cl)]

      all_res[[length(all_res) + 1]] <- dt_wide
    }
  }

  final_dt <- rbindlist(all_res, use.names = TRUE, fill = TRUE)

  setcolorder(final_dt, c(setting_col, "cell_line", "analysis_method",
                         "pearson_r","spearman_rho","ccc","bias_str","ba_corr"))

  fwrite(final_dt, out_tsv, sep = "\t", quote = FALSE)
  cat("Written:", out_tsv, "\n")
  invisible(final_dt)
}