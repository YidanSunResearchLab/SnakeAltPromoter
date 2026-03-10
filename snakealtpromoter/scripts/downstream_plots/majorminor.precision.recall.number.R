library(data.table)
library(ggplot2)
library(scales)

# ============================================================
# GLOBAL STYLE
# ============================================================

theme_paper_panel <- function() {
  
  base <- 20
  
  theme_classic(base_size = base, base_family = "Helvetica") +
    theme(
      plot.title = element_text(size = base + 2, hjust = 0.5, face = "plain", color = "black"),
      axis.title = element_text(size = base + 2, color = "black"),
      axis.text  = element_text(size = base, color = "black"),
      legend.position = "none"
    )
}

# ============================================================
# COLORS
# ============================================================

cols_method <- c(
  "Salmon"   = "#3C78D8",  # blue
  "DEXSeq"   = "#D64B4B",  # red
  "proActiv" = "#7A5AA6"   # purple
)

# ============================================================
# PATHS
# ============================================================

dir.create(out_dir_depth_mm, showWarnings = FALSE, recursive = TRUE)
dir.create(out_dir_len_mm,   showWarnings = FALSE, recursive = TRUE)
dir.create(out_dir_thr,      showWarnings = FALSE, recursive = TRUE)

# ============================================================
# SETTINGS
# ============================================================

methods <- c("Salmon", "DEXSeq", "proActiv")
groups  <- c("major", "minor")
metrics_mm <- c("Number", "Precision", "Recall")

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
  file_map <- list(
    "0.1"  = file.path(base_dir, "venn_major_minor_0.1/overlap_summary.tsv"),
    "0.25" = file.path(base_dir, "venn_major_minor/overlap_summary.tsv"),
    "0.5"  = file.path(base_dir, "venn_major_minor_0.5/overlap_summary.tsv")
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
  
  out <- rbindlist(
    lapply(names(file_map), function(thr) read_one(file_map[[thr]], thr)),
    fill = TRUE
  )
  
  out[, threshold := factor(threshold, levels = thr_levels)]
  out
}

add_number_ratio <- function(dt, xcol, baseline_level) {
  dt <- copy(dt)
  
  base_dt <- dt[get(xcol) == baseline_level, .(
    method,
    group_label,
    baseline_size = Size_Set1
  )]
  
  dt <- merge(dt, base_dt, by = c("method", "group_label"), all.x = TRUE)
  dt[, number := 100 * Size_Set1 / baseline_size]
  dt
}

make_mm_long <- function(dt, xcol, xlab_text) {
  dt2 <- copy(dt)
  
  long <- melt(
    dt2,
    id.vars = c(xcol, "method", "group_label"),
    measure.vars = c("number", "precision", "recall"),
    variable.name = "metric",
    value.name = "value"
  )
  
  metric_map2 <- c(
    number    = "Number",
    precision = "Precision",
    recall    = "Recall"
  )
  
  long[, metric := metric_map2[as.character(metric)]]
  long[, xlab := xlab_text]
  setnames(long, old = xcol, new = "x_value")
  long
}

# ============================================================
# PANEL FUNCTION
# ============================================================

plot_mm_panel <- function(long_dt, grp, metric_name, xlab_text, ylab_text = NULL) {
  
  sub <- long_dt[group_label == grp & metric == metric_name]
  
  ggplot(sub, aes(x = x_value, y = value, color = method, group = method)) +
    geom_line(linewidth = 1.4, alpha = 0.7) +
    geom_point(size = 3.0, alpha = 0.7) +
    scale_color_manual(values = cols_method, drop = FALSE) +
    labs(
      title = NULL,
      x = xlab_text,
      y = ylab_text,
      color = NULL
    ) +
    theme_paper_panel()
}

# ============================================================
# DEPTH MAJOR/MINOR
# baseline = 80M
# ============================================================

depth_mm_dt <- read_overlap_series(
  base_dir = base_dir,
  settings = depth_levels
)
depth_mm_dt[, method := standardize_method(method)]
depth_mm_dt[, setting := factor(setting, levels = depth_levels)]
depth_mm_dt <- add_number_ratio(depth_mm_dt, "setting", "80M")

depth_long <- make_mm_long(
  dt = depth_mm_dt,
  xcol = "setting",
  xlab_text = "Downsample depth"
)
depth_long[, x_value := factor(x_value, levels = depth_levels)]
depth_long[, method := factor(method, levels = methods)]

for (g in groups) {
  for (met in metrics_mm) {
    
    p <- plot_mm_panel(
      long_dt = depth_long,
      grp = g,
      metric_name = met,
      xlab_text = "Downsample depth",
      ylab_text = met
    )
    
    ggsave(
      file.path(out_dir_depth_mm, paste0("Depth_", g, "_", met, ".pdf")),
      p, width = 6.2, height = 5.2
    )
    ggsave(
      file.path(out_dir_depth_mm, paste0("Depth_", g, "_", met, ".png")),
      p, width = 6.2, height = 5.2, dpi = 300
    )
  }
}

# ============================================================
# READ LENGTH MAJOR/MINOR
# baseline = 150
# ============================================================

len_mm_dt <- read_overlap_series(
  base_dir = base_dir,
  settings = as.character(len_levels)
)
len_mm_dt[, method := standardize_method(method)]
len_mm_dt[, setting := as.integer(setting)]
len_mm_dt <- add_number_ratio(len_mm_dt, "setting", 150)
len_mm_dt[, setting := factor(setting, levels = len_levels)]

len_long <- make_mm_long(
  dt = len_mm_dt,
  xcol = "setting",
  xlab_text = "Simulated read length (bp)"
)
len_long[, x_value := factor(x_value, levels = len_levels)]
len_long[, method := factor(method, levels = methods)]

for (g in groups) {
  for (met in metrics_mm) {
    
    p <- plot_mm_panel(
      long_dt = len_long,
      grp = g,
      metric_name = met,
      xlab_text = "Simulated read length (bp)",
      ylab_text = met
    )
    
    ggsave(
      file.path(out_dir_len_mm, paste0("Readlen_", g, "_", met, ".pdf")),
      p, width = 6.2, height = 5.2
    )
    ggsave(
      file.path(out_dir_len_mm, paste0("Readlen_", g, "_", met, ".png")),
      p, width = 6.2, height = 5.2, dpi = 300
    )
  }
}

# ============================================================
# THRESHOLD MAJOR/MINOR
# baseline = 0.1
# ============================================================

thr_dt <- read_threshold_data()
thr_dt[, method := standardize_method(method)]
thr_dt <- add_number_ratio(thr_dt, "threshold", "0.1")

thr_long <- make_mm_long(
  dt = thr_dt,
  xcol = "threshold",
  xlab_text = "Threshold"
)
thr_long[, x_value := factor(x_value, levels = thr_levels)]
thr_long[, method := factor(method, levels = methods)]

for (g in groups) {
  for (met in metrics_mm) {
    
    p <- plot_mm_panel(
      long_dt = thr_long,
      grp = g,
      metric_name = met,
      xlab_text = "Threshold",
      ylab_text = met
    )
    
    ggsave(
      file.path(out_dir_thr, paste0("Threshold_", g, "_", met, ".pdf")),
      p, width = 6.2, height = 5.2
    )
    ggsave(
      file.path(out_dir_thr, paste0("Threshold_", g, "_", met, ".png")),
      p, width = 6.2, height = 5.2, dpi = 300
    )
  }
}

cat("All major/minor figures generated.\n")

library(ggplot2)
library(cowplot)

cols_method <- c(
  "proActiv" = "#7A5AA6",
  "Salmon"   = "#3C78D8",
  "DEXSeq"   = "#D64B4B"
)

df <- data.frame(
  x = 1:3,
  y = 1,
  method = factor(
    c("proActiv","Salmon","DEXSeq"),
    levels=c("proActiv","Salmon","DEXSeq")
  )
)


p <- ggplot(df, aes(x, y, color = method)) +
  geom_line(linewidth = 1.8) +
  geom_point(size = 3.8) +
  scale_color_manual(values = cols_method) +
  theme_void(base_family = "Helvetica") +
  theme(
    legend.position = "right",
    legend.text = element_text(size = 16, color = "black", family="Helvetica"),
    legend.title = element_blank()
  )

# legend
legend <- cowplot::get_legend(p)

# legend
ggsave(
  file.path(out_dir_depth_mm, "method_legend_vertical.pdf"),
  legend,
  width = 2,
  height = 3,
  bg = "white"
)

ggsave(
  file.path(out_dir_depth_mm, "method_legend_vertical.png"),
  legend,
  width = 2,
  height = 3,
  dpi = 300,
  bg = "white"
)
