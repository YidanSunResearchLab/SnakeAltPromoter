# ---------------------------
# Load Required Libraries
# ---------------------------
library(ggvenn)
library(ggplot2)
library(gridExtra)
library(grid)

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
fill_colors <- c("#ffeda0", "#a6bddb")

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
    scale_fill_manual(values = c(Major = "#fec44f", Minor = "#fff7bc")) +
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
#--------------------------------------------------------------------------------------------------------

# ---------------------------
# Helper function of logFC IQR plots
# ---------------------------
plot_log2fc_iqr <- function(data_list, condition_name, value_col, file_prefix) {
  directions <- c("up" = paste0(file_prefix, "_up.pdf"),
                  "down" = paste0(file_prefix, "_down.pdf"))
  
  for (dir in names(directions)) {
    stats <- do.call(rbind, lapply(names(data_list), function(method) {
      vec <- data_list[[method]][[value_col]]
      vec <- vec[is.finite(vec)]
      if (dir == "up")   vec <- vec[vec > 0]
      if (dir == "down") vec <- vec[vec < 0]
      if (length(vec) < 3) return(NULL)
      data.frame(
        method = method,
        median = median(vec),
        q25 = quantile(vec, 0.25),
        q75 = quantile(vec, 0.75),
        ymin = max(min(vec), quantile(vec, 0.25) - 1.5 * IQR(vec)),  # lower whisker
        ymax = min(max(vec), quantile(vec, 0.75) + 1.5 * IQR(vec))   # upper whisker
      )
    }))
    
    if (nrow(stats) > 0) {
      stats$xval <- as.numeric(factor(stats$method, levels = unique(stats$method))) 
      
      p <- ggplot(stats, aes(x = xval)) +
        geom_rect(aes(xmin = xval - 0.25, xmax = xval + 0.25,
                      ymin = q25, ymax = q75), fill = "grey80", color = "black") +
        geom_segment(aes(x = xval - 0.25, xend = xval + 0.25,
                         y = median, yend = median), color = "black", size = 0.6) +
        geom_segment(aes(x = xval, xend = xval, y = ymin, yend = q25), linetype = "solid") +
        geom_segment(aes(x = xval, xend = xval, y = q75, yend = ymax), linetype = "solid") +
        scale_x_continuous(breaks = stats$xval, labels = stats$method) + 
        ylab("log2 FC") + xlab("") +
        ggtitle(paste0(condition_name, " (", dir, "regulated)")) +
        theme_bw() +
        theme(plot.title = element_text(hjust = 0.5))
      
      ggsave(file.path(out_dir, directions[[dir]]), p, width = 6, height = 4)
    }
  }
}

# ---------------------------
# Prepare data list
# ---------------------------

print(colnames(dexseq_promoter))
print(length(dexseq_promoter$geneId[dexseq_promoter$pvalue_usage < 0.05]))
print(length(dexseq_promoter$geneId[dexseq_promoter$padj_usage < 0.05]))
print(length(dexseq_promoter$promoterId[dexseq_promoter$pvalue_usage < 0.05]))
print(length(dexseq_promoter$promoterId[dexseq_promoter$padj_usage < 0.05]))


promoter_padj <- list(
  CAGE     = CAGE_promoter[as.numeric(CAGE_promoter$padj_usage) < 0.05, ],
  Salmon   = salmon_promoter[as.numeric(salmon_promoter$padj_usage) < 0.05, ],
  proActiv = proactiv_promoter[as.numeric(proactiv_promoter$padj_usage) < 0.05, ],
  DEXSeq   = dexseq_promoter[as.numeric(dexseq_promoter$padj_usage) < 0.05, ]
)

promoter_pval <- list(
  CAGE     = CAGE_promoter[as.numeric(CAGE_promoter$pvalue_usage) < 0.05, ],
  Salmon   = salmon_promoter[as.numeric(salmon_promoter$pvalue_usage) < 0.05, ],
  proActiv = proactiv_promoter[as.numeric(proactiv_promoter$pvalue_usage) < 0.05, ],
  DEXSeq   = dexseq_promoter[as.numeric(dexseq_promoter$pvalue_usage) < 0.05, ]
)

print(colnames(dexseq_gene))
print(length(dexseq_gene$geneId[dexseq_gene$FDR < 0.05]))


gene_fdr <- list(
  CAGE     = CAGE_gene[as.numeric(CAGE_gene$FDR) < 0.05, ],
  Salmon   = salmon_gene[as.numeric(salmon_gene$FDR) < 0.05, ],
  proActiv = proactiv_gene[as.numeric(proactiv_gene$FDR) < 0.05, ],
  DEXSeq   = dexseq_gene[as.numeric(dexseq_gene$FDR) < 0.05, ]
)

# ---------------------------
# Generate all 6 plots
# ---------------------------
plot_log2fc_iqr(promoter_padj, "Promoter log2FC (padj < 0.05)", "log2FC_usage", "log2FC_IQR_padj")
plot_log2fc_iqr(promoter_pval, "Promoter log2FC (pvalue < 0.05)", "log2FC_usage", "log2FC_IQR_pval")
plot_log2fc_iqr(gene_fdr,      "Gene log2FC (FDR < 0.05)",        "logFC",        "log2FC_gene_IQR_FDR")




#--------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------


# ---------------------------
# Barplot: Single stacked bar per method (Up/Down segments)
# ---------------------------

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
    scale_fill_manual(values = c("Up" = "#ffeda0", "Down" = "#a6bddb")) +
    labs(x = NULL, y = "Significant count", title = title) +
    theme_bw() +
    theme(legend.position = "right",
          plot.title = element_text(hjust = 0.5, face = "bold"))

  ggsave(file.path(out_dir, outfile), p, width = 6, height = 4.5)
}








# ---------------------------
# Prepare data lists
# ---------------------------
promoter_list <- list(
  CAGE     = CAGE_promoter,
  Salmon   = salmon_promoter,
  proActiv = proactiv_promoter,
  DEXSeq   = dexseq_promoter
)

gene_list <- list(
  CAGE     = CAGE_gene,
  Salmon   = salmon_gene,
  proActiv = proactiv_gene,
  DEXSeq   = dexseq_gene
)

# ---------------------------
# Count Up/Down
# ---------------------------
promoter_padj_df <- count_up_down_by_method(promoter_list, "log2FC_usage", "padj_usage", 0.05)
promoter_pval_df <- count_up_down_by_method(promoter_list, "log2FC_usage", "pvalue_usage", 0.05)
gene_fdr_df      <- count_up_down_by_method(gene_list,      "logFC",         "FDR",          0.05)

# ---------------------------
# Plot stacked barplots
# ---------------------------
plot_stacked_updown_bar(promoter_padj_df, "Promoter log2FC (padj < 0.05)", "promoter_count_padj.pdf")
plot_stacked_updown_bar(promoter_pval_df, "Promoter log2FC (p-value < 0.05)", "promoter_count_pval.pdf")
plot_stacked_updown_bar(gene_fdr_df,      "Gene log2FC (FDR < 0.05)",          "gene_count_FDR.pdf")

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
# Promoter pie charts (padj)
# ---------------------------
threshold <- 0
cage_promoters_down <- CAGE_promoter$promoterId[CAGE_promoter$padj_usage < 0.05 & CAGE$deltaPU < threshold]
cage_promoters_up <- CAGE_promoter$promoterId[CAGE_promoter$padj_usage < 0.05 & CAGE$deltaPU > threshold]
cage_promoters <- union(cage_promoters_up, cage_promoters_down)
salmon_promoters_down <- salmon_promoter$promoterId[salmon_promoter$padj_usage < 0.05 & salmon$deltaPU < threshold]
salmon_promoters_up <- salmon_promoter$promoterId[salmon_promoter$padj_usage < 0.05 & salmon$deltaPU > threshold]
salmon_promoters <- union(salmon_promoters_up, salmon_promoters_down)
proactiv_promoters_down <- proactiv_promoter$promoterId[proactiv_promoter$padj_usage < 0.05 & proactiv$deltaPU < threshold]
proactiv_promoters_up <- proactiv_promoter$promoterId[proactiv_promoter$padj_usage < 0.05 & proactiv$deltaPU > threshold]
proactiv_promoters <- union(proactiv_promoters_up, proactiv_promoters_down)
dexseq_promoters_down <- dexseq_promoter$promoterId[dexseq_promoter$padj_usage < 0.05 & dexseq$deltaPU < threshold]
dexseq_promoters_up <- dexseq_promoter$promoterId[dexseq_promoter$padj_usage < 0.05 & dexseq$deltaPU > threshold]
dexseq_promoters <- union(dexseq_promoters_up, dexseq_promoters_down)

salmon_promoter_pie   <- get_overlap_pie(cage_promoters, salmon_promoters,   "Salmon")
dexseq_promoter_pie   <- get_overlap_pie(cage_promoters, dexseq_promoters,   "DEXSeq")
proactiv_promoter_pie <- get_overlap_pie(cage_promoters, proactiv_promoters, "proActiv")

# ---------------------------
# Gene pie charts (FDR)
# ---------------------------
cage_genes <- CAGE_gene$geneId[CAGE_gene$FDR < 0.05]
salmon_gene_ids <- salmon_gene$geneId[salmon_gene$FDR < 0.05]
dexseq_gene_ids <- dexseq_gene$geneId[dexseq_gene$FDR < 0.05]
proactiv_gene_ids <- proactiv_gene$geneId[proactiv_gene$FDR < 0.05]

salmon_gene_pie   <- get_overlap_pie(cage_genes, salmon_gene_ids,   "Salmon")
dexseq_gene_pie   <- get_overlap_pie(cage_genes, dexseq_gene_ids,   "DEXSeq")
proactiv_gene_pie <- get_overlap_pie(cage_genes, proactiv_gene_ids, "proActiv")

# ---------------------------
# Save pie plots
# ---------------------------
pdf(file.path(out_dir, "combined_promoter_overlap_pie.pdf"), width = 10, height = 4)
grid.arrange(salmon_promoter_pie, dexseq_promoter_pie, proactiv_promoter_pie,
             ncol = 3,
             top = textGrob("CAGE Overlap of Significant Promoters (padj < 0.05)", gp = gpar(fontsize = 12, fontface = "bold")))
dev.off()

pdf(file.path(out_dir, "combined_gene_overlap_pie.pdf"), width = 10, height = 4)
grid.arrange(salmon_gene_pie, dexseq_gene_pie, proactiv_gene_pie,
             ncol = 3,
             top = textGrob("CAGE Overlap of DE Genes (FDR < 0.05)", gp = gpar(fontsize = 12, fontface = "bold")))
dev.off()






# ---------------------------
# Venn diagram for padj promoter overlap
# ---------------------------
venn_pval_promoter <- list(
  CAGE = CAGE_promoter$promoterId[CAGE_promoter$pvalue_usage < 0.05],
  Salmon = salmon_promoter$promoterId[salmon_promoter$pvalue_usage < 0.05],
  DEXSeq = dexseq_promoter$promoterId[dexseq_promoter$pvalue_usage < 0.05],
  proActiv = proactiv_promoter$promoterId[proactiv_promoter$pvalue_usage < 0.05]
)

p_promoter <- ggvenn::ggvenn(venn_pval_promoter, fill_color = c("#cbcde0", "#decbe0", "#e0decb", "#cde0cb"),
                             stroke_size = 0.5, text_size = 4, show_percentage = FALSE)
ggsave(file.path(out_dir, "venn_promoter_pval.pdf"), p_promoter, width = 6, height = 6)

# ---------------------------
# Venn diagram for gene overlap (FDR)
# ---------------------------
cage_gene_ids <- cage_genes
venn_fdr_gene <- list(
  CAGE = cage_gene_ids,
  Salmon = salmon_gene$geneId[salmon_gene$FDR < 0.05],
  DEXSeq = dexseq_gene$geneId[dexseq_gene$FDR < 0.05],
  proActiv = proactiv_gene$geneId[proactiv_gene$FDR < 0.05]
)

# ---------------------------
# Individual method vs CAGE (Gene FDR)
# ---------------------------
venn_gene_salmon   <- ggvenn(list(Salmon = salmon_gene$geneId[salmon_gene$FDR < 0.05],
                                  CAGE = cage_gene_ids),
                             fill_color = c("#b3cde3", "#fbb4ae"),
                             stroke_size = 0.5, text_size = 4, show_percentage = FALSE, set_name_size = 4) +
                      labs(title = "Salmon vs CAGE (DEG)") +
                      theme(plot.title = element_text(hjust = 0.5, size = 10))

venn_gene_dexseq   <- ggvenn(list(DEXSeq = dexseq_gene$geneId[dexseq_gene$FDR < 0.05],
                                  CAGE = cage_gene_ids),
                             fill_color = c("#b3cde3", "#fbb4ae"),
                             stroke_size = 0.5, text_size = 4, show_percentage = FALSE, set_name_size = 4) +
                      labs(title = "DEXSeq vs CAGE (DEG)") +
                      theme(plot.title = element_text(hjust = 0.5, size = 10))

venn_gene_proactiv <- ggvenn(list(proActiv = proactiv_gene$geneId[proactiv_gene$FDR < 0.05],
                                  CAGE = cage_gene_ids),
                             fill_color = c("#b3cde3", "#fbb4ae"),
                             stroke_size = 0.5, text_size = 4, show_percentage = FALSE, set_name_size = 4) +
                      labs(title = "proActiv vs CAGE (DEG)") +
                      theme(plot.title = element_text(hjust = 0.5, size = 10))

# ---------------------------
# Individual method vs CAGE (Promoter p-value)
# ---------------------------
venn_promoter_salmon <- ggvenn(list(Salmon = salmon_promoter$promoterId[salmon_promoter$pvalue_usage < 0.05],
                                    CAGE = CAGE_promoter$promoterId[CAGE_promoter$pvalue_usage < 0.05]),
                               fill_color = c("#b3cde3", "#fbb4ae"),
                               stroke_size = 0.5, text_size = 4, show_percentage = FALSE, set_name_size = 4) +
                        labs(title = "Salmon vs CAGE (Promoter)") +
                        theme(plot.title = element_text(hjust = 0.5, size = 10))

venn_promoter_dexseq <- ggvenn(list(DEXSeq = dexseq_promoter$promoterId[dexseq_promoter$pvalue_usage < 0.05],
                                    CAGE = CAGE_promoter$promoterId[CAGE_promoter$pvalue_usage < 0.05]),
                               fill_color = c("#b3cde3", "#fbb4ae"),
                               stroke_size = 0.5, text_size = 4, show_percentage = FALSE, set_name_size = 4) +
                        labs(title = "DEXSeq vs CAGE (Promoter)") +
                        theme(plot.title = element_text(hjust = 0.5, size = 10))

venn_promoter_proactiv <- ggvenn(list(proActiv = proactiv_promoter$promoterId[proactiv_promoter$pvalue_usage < 0.05],
                                      CAGE = CAGE_promoter$promoterId[CAGE_promoter$pvalue_usage < 0.05]),
                                 fill_color = c("#b3cde3", "#fbb4ae"),
                                 stroke_size = 0.5, text_size = 4, show_percentage = FALSE, set_name_size = 4) +
                          labs(title = "proActiv vs CAGE (Promoter)") +
                          theme(plot.title = element_text(hjust = 0.5, size = 10))
# Combine promoter Venns
pdf(file.path(out_dir, "combined_promoter_overlap_venn.pdf"), width = 10, height = 4)
grid.arrange(
  venn_promoter_salmon, venn_promoter_dexseq, venn_promoter_proactiv,
  ncol = 3,
  top = textGrob("Venn Diagram: Promoter Overlap (p < 0.05)", gp = gpar(fontsize = 12, fontface = "bold"))
)
dev.off()

# Combine gene Venns
pdf(file.path(out_dir, "combined_gene_overlap_venn.pdf"), width = 10, height = 4)
grid.arrange(
  venn_gene_salmon, venn_gene_dexseq, venn_gene_proactiv,
  ncol = 3,
  top = textGrob("Venn Diagram: DEG Overlap (FDR < 0.05)", gp = gpar(fontsize = 12, fontface = "bold"))
)
dev.off()

# ---------------------------
# Individual method vs CAGE (Promoter padj-value)
# ---------------------------
venn_promoter_salmon <- ggvenn(list(Salmon = salmon_promoter$promoterId[salmon_promoter$padj_usage < 0.05],
                                    CAGE = CAGE_promoter$promoterId[CAGE_promoter$padj_usage < 0.05]),
                               fill_color = c("#b3cde3", "#fbb4ae"),
                               stroke_size = 0.5, text_size = 4, show_percentage = FALSE) +
                        labs(title = "Salmon vs CAGE (Promoter)") +
                        theme(plot.title = element_text(hjust = 0.5, size = 10))

venn_promoter_dexseq <- ggvenn(list(DEXSeq = dexseq_promoter$promoterId[dexseq_promoter$padj_usage < 0.05],
                                    CAGE = CAGE_promoter$promoterId[CAGE_promoter$padj_usage < 0.05]),
                               fill_color = c("#b3cde3", "#fbb4ae"),
                               stroke_size = 0.5, text_size = 4, show_percentage = FALSE) +
                        labs(title = "DEXSeq vs CAGE (Promoter)") +
                        theme(plot.title = element_text(hjust = 0.5, size = 10))

venn_promoter_proactiv <- ggvenn(list(proActiv = proactiv_promoter$promoterId[proactiv_promoter$padj_usage < 0.05],
                                      CAGE = CAGE_promoter$promoterId[CAGE_promoter$padj_usage < 0.05]),
                                 fill_color = c("#b3cde3", "#fbb4ae"),
                                 stroke_size = 0.5, text_size = 4, show_percentage = FALSE) +
                          labs(title = "proActiv vs CAGE (Promoter)") +
                          theme(plot.title = element_text(hjust = 0.5, size = 10))
# Combine promoter Venns
pdf(file.path(out_dir, "padj_combined_promoter_overlap_venn.pdf"), width = 10, height = 4)
grid.arrange(
  venn_promoter_salmon, venn_promoter_dexseq, venn_promoter_proactiv,
  ncol = 3,
  top = textGrob("Venn Diagram: Promoter Overlap (padj < 0.05)", gp = gpar(fontsize = 12, fontface = "bold"))
)
dev.off()
