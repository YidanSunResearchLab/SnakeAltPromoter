#Run at the start
#Load results from scatterplot.count.rpkm.ran and summarize correlations
library(data.table)
down_base <- "/home/yuqing/shared/SnakeAltPromoter_0228/downstream_fixed/downsample"
depths <- c("10M","20M","30M","40M","80M")
cell_lines <- c("healthy","failure")

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


summarize_correlations(
  base_dir = down_base,
  settings = depths,
  cell_lines = cell_lines,
  out_tsv = file.path(down_base, "downsample_correlations_summary.tsv"),
  setting_col = "downsample"
)

len_base <- "/home/yuqing/shared/SnakeAltPromoter_0228/downstream_fixed/readlensim"
lens <- c("50","75","100","150")
cell_lines <- c("healthy","failure")

summarize_correlations(
  base_dir = len_base,
  settings = lens,
  cell_lines = cell_lines,
  out_tsv = file.path(len_base, "readlensim_correlations_summary.tsv"),
  setting_col = "readlensim"
)


#Read summary tsv
base_dir <- "/home/yuqing/shared/SnakeAltPromoter_0228/downstream_fixed"
down_tsv <- file.path(base_dir, "downsample/downsample_correlations_summary.tsv")
len_tsv  <- file.path(base_dir, "readlensim/readlensim_correlations_summary.tsv")
#Output directory for sequence depth and read length
out_dir_depth    <- file.path(base_dir, "downsample_plots_transposed")
out_dir_len      <- file.path(base_dir, "readlensim_plots_transposed")
#Output directory for major and minor promoter classification
out_dir_depth_mm <- file.path(base_dir, "downsample_major_minor_plots_transposed")
out_dir_len_mm   <- file.path(base_dir, "readlen_major_minor_plots_transposed")
out_dir_thr      <- file.path(base_dir, "threshold_major_minor_plots_transposed")