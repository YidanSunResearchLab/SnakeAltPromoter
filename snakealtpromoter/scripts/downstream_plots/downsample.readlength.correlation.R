library(data.table)
library(ggplot2)
library(scales)

# ============================================================
# GLOBAL STYLE
# ============================================================

cols_method <- c(
  "Salmon"   = "#3C78D8",
  "DEXSeq"   = "#D64B4B",
  "proActiv" = "#7A5AA6"
)

theme_paper_panel <- function() {
  theme_classic(base_size = 20, base_family = "Helvetica") +
    theme(
      plot.title = element_text(size = 20, hjust = 0.5, face = "plain", color = "black"),
      axis.title = element_text(size = 20, color = "black"),
      axis.text  = element_text(size = 20, color = "black"),
      legend.title = element_blank(),
      legend.text = element_text(size = 20, color = "black")
    )
}

# ============================================================
# PATHS
# ============================================================

dir.create(out_dir_depth,    showWarnings = FALSE, recursive = TRUE)
dir.create(out_dir_len,      showWarnings = FALSE, recursive = TRUE)
dir.create(out_dir_depth_mm, showWarnings = FALSE, recursive = TRUE)
dir.create(out_dir_len_mm,   showWarnings = FALSE, recursive = TRUE)
dir.create(out_dir_thr,      showWarnings = FALSE, recursive = TRUE)

# ============================================================
# SETTINGS
# ============================================================

methods <- c("proActiv", "Salmon", "DEXSeq")
cell_lines <- c("HEALTHY", "FAILURE")
groups <- c("major", "minor")

metrics <- c("pearson_r", "spearman_rho", "ccc", "ba_corr", "bias_str")

metric_map <- c(
  pearson_r    = "Pearson",
  spearman_rho = "Spearman",
  ccc          = "CCC",
  ba_corr      = "BA-corr",
  bias_str     = "Bias-str"
)

depth_levels <- c("10M", "20M", "30M", "40M", "80M")
len_levels   <- c(50, 75, 100, 150)
thr_levels   <- c("0.1", "0.25", "0.5")

# ============================================================
# HELPERS
# ============================================================

pct_to_num <- function(x) {
  as.numeric(gsub("%", "", as.character(x)))
}

standardize_method <- function(x) {
  x <- trimws(as.character(x))
  x[x %in% c("proactiv", "ProActiv", "Proactiv", "proActiv")] <- "proActiv"
  x[x %in% c("DEXSeq", "DexSeq", "dexseq")] <- "DEXSeq"
  x[x %in% c("Salmon", "salmon")] <- "Salmon"
  x
}

to_long <- function(dt, xcol) {
  long <- melt(
    dt,
    id.vars = c(xcol, "cell_line", "analysis_method"),
    measure.vars = metrics,
    variable.name = "metric",
    value.name = "value"
  )

  long[, metric := as.character(metric)]
  long[, metric := metric_map[metric]]
  long[, analysis_method := standardize_method(analysis_method)]
  long
}

read_overlap_series <- function(base_dir, settings, suffix = "_major_minor") {

  read_one_overlap <- function(setting) {
    fpath <- file.path(base_dir, paste0(setting, suffix), "overlap_summary.tsv")
    if (!file.exists(fpath)) stop("Missing file: ", fpath)

    dt <- fread(fpath)
    dt[, precision := pct_to_num(Percent_in_Set1)]
    dt[, recall    := pct_to_num(Percent_in_Set2)]
    dt[, method := standardize_method(Set1)]
    dt[, group_label := ifelse(grepl("^major", Group), "major",
                        ifelse(grepl("^minor", Group), "minor", Group))]
    dt[, setting := setting]
    dt
  }

  rbindlist(lapply(settings, read_one_overlap), fill = TRUE)
}

read_threshold_data <- function() {
  files <- c(
    file.path(base_dir, "venn_major_minor_0.1/overlap_summary.tsv"),
    file.path(base_dir, "venn_major_minor/overlap_summary.tsv"),
    file.path(base_dir, "venn_major_minor_0.5/overlap_summary.tsv")
  )

  read_one <- function(fpath, thr_label) {
    if (!file.exists(fpath)) stop("Missing file: ", fpath)

    dt <- fread(fpath)
    dt[, precision := pct_to_num(Percent_in_Set1)]
    dt[, recall    := pct_to_num(Percent_in_Set2)]
    dt[, method := standardize_method(Set1)]
    dt[, group_label := ifelse(grepl("^major", Group), "major",
                        ifelse(grepl("^minor", Group), "minor", Group))]
    dt[, threshold := thr_label]
    dt
  }

  out <- rbindlist(list(
    read_one(files[1], "0.1"),
    read_one(files[2], "0.25"),
    read_one(files[3], "0.5")
  ), fill = TRUE)

  out[, threshold := factor(threshold, levels = thr_levels)]
  out
}

# ============================================================
# PANEL FUNCTIONS
# ============================================================

plot_corr_metric_panel <- function(long_dt, xcol, cell_name, metric_name, xlab, ylab = NULL) {
  sub <- long_dt[
    cell_line == cell_name &
    metric == metric_name
  ]

  sub[, analysis_method := factor(analysis_method, levels = methods)]

  ggplot(sub, aes_string(x = xcol, y = "value", color = "analysis_method", group = "analysis_method")) +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.4, color = "grey70") +
    geom_line(linewidth = 1.2, alpha = 0.7) +
    geom_point(size = 2.6) +
    scale_color_manual(values = cols_method, drop = FALSE) +
    labs(
      title = NULL,
      x = xlab,
      y = ylab,
      color = NULL
    ) +
    theme_paper_panel() +
    theme(
      legend.position = "none"
    )
}

plot_overlap_metric_panel <- function(dt, xcol, grp, metric_name, xlab, ylab = NULL) {
  sub <- copy(dt[group_label == grp])

  long <- melt(
    sub,
    id.vars = c(xcol, "method", "group_label"),
    measure.vars = c("precision", "recall"),
    variable.name = "metric",
    value.name = "value"
  )

  metric_map2 <- c(
    precision = "Precision",
    recall    = "Recall"
  )

  long[, metric := metric_map2[as.character(metric)]]
  long <- long[metric == metric_name]
  long[, method := factor(method, levels = methods)]

  ggplot(long, aes_string(x = xcol, y = "value", color = "method", group = "method")) +
    geom_line(linewidth = 1.2, alpha = 0.7) +
    geom_point(size = 2.6) +
    scale_color_manual(values = cols_method, drop = FALSE) +
    labs(
      title = NULL,
      x = xlab,
      y = ylab,
      color = NULL
    ) +
    theme_paper_panel() +
    theme(
      legend.position = "none"
    )
}

# ============================================================
# LOAD DEPTH / READLEN CORRELATION DATA
# ============================================================

if (!file.exists(down_tsv)) stop("Missing: ", down_tsv)
if (!file.exists(len_tsv))  stop("Missing: ", len_tsv)

down_dt <- fread(down_tsv)
len_dt  <- fread(len_tsv)

down_dt[, cell_line := toupper(cell_line)]
len_dt[,  cell_line := toupper(cell_line)]

down_dt[, analysis_method := trimws(analysis_method)]
len_dt[,  analysis_method := trimws(analysis_method)]

down_dt[, downsample := factor(downsample, levels = depth_levels)]

len_dt[, readlensim := as.integer(readlensim)]
len_dt[, readlensim := factor(readlensim, levels = len_levels)]

down_long <- to_long(down_dt, "downsample")
len_long  <- to_long(len_dt,  "readlensim")

# ============================================================
# DEPTH CORRELATION / BIAS FIGURES
# ============================================================

metric_labels <- c(
  "Pearson"  = "Pearson",
  "Spearman" = "Spearman",
  "CCC"      = "CCC",
  "BA-corr"  = "BA-corr",
  "Bias-str" = "Mean bias"
)

for (cl in cell_lines) {
  for (met in c("Pearson", "Spearman", "CCC", "BA-corr", "Bias-str")) {

    p <- plot_corr_metric_panel(
      down_long,
      "downsample",
      cl,
      met,
      "Downsample depth",
      metric_labels[[met]]
    )

    ggsave(
      file.path(out_dir_depth, paste0("Depth_", cl, "_", met, ".pdf")),
      p, width = 5.8, height = 4.8
    )
    ggsave(
      file.path(out_dir_depth, paste0("Depth_", cl, "_", met, ".png")),
      p, width = 5.8, height = 4.8, dpi = 300
    )
  }
}

# ============================================================
# READ LENGTH CORRELATION / BIAS FIGURES
# ============================================================

for (cl in cell_lines) {
  for (met in c("Pearson", "Spearman", "CCC", "BA-corr", "Bias-str")) {

    p <- plot_corr_metric_panel(
      len_long,
      "readlensim",
      cl,
      met,
      "Simulated read length (bp)",
      metric_labels[[met]]
    )

    ggsave(
      file.path(out_dir_len, paste0("Readlen_", cl, "_", met, ".pdf")),
      p, width = 5.8, height = 4.8
    )
    ggsave(
      file.path(out_dir_len, paste0("Readlen_", cl, "_", met, ".png")),
      p, width = 5.8, height = 4.8, dpi = 300
    )
  }
}

cat("All transposed figures generated.\n")