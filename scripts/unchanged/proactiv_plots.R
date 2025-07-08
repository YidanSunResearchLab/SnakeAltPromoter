#!/usr/bin/env Rscript

library(proActiv)
library(ggplot2)
library(reshape2)
library(ggpubr)

# Parse arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 6) {
  stop("Usage: Rscript proactiv_plots.R <output_dir> <promoter_rds> <abs_activity> <gene_expr> <category> <norm_method>")
}
output_dir <- args[1]
promoter_rds <- args[2]
abs_activity_file <- args[3]
gene_expr_file <- args[4]
category_file <- args[5]
norm_method <- args[6]

# Validate inputs
if (!file.exists(promoter_rds)) stop("Promoter RDS file does not exist")
if (!file.exists(abs_activity_file)) stop("Absolute activity file does not exist")
if (!file.exists(gene_expr_file)) stop("Gene expression file does not exist")
if (!file.exists(category_file)) stop("Category file does not exist")

# Create output directory
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Load data
promoterAnnotationData <- readRDS(promoter_rds)
promoterCoordinates(promoterAnnotationData)[is.na(promoterCoordinates(promoterAnnotationData)$internalPromoter)]$internalPromoter <- FALSE
promoterCoordinates(promoterAnnotationData) <- promoterCoordinates(promoterAnnotationData)[!promoterCoordinates(promoterAnnotationData)$internalPromoter]
absolutePromoterActivity <- readRDS(abs_activity_file)
geneExpression <- readRDS(gene_expr_file)
absolutePromoterActivityCategory <- readRDS(category_file)

# Theme for plots
theme_classic_add <- theme_classic() + theme(
  legend.position = "none",
  strip.background = element_blank(),
  strip.text = element_blank(),
  plot.title = element_text(hjust = 0.5, family = "Helvetica"),
  axis.text = element_text(colour = "black", family = "Helvetica"),
  axis.title = element_text(colour = "black", family = "Helvetica"),
  text = element_text(colour = "black", family = "Helvetica")
)

# Plot 1: Promoter category percentage (Fig. 1c)
category_counts <- as.data.frame(table(absolutePromoterActivityCategory$PromoterCategory))
category_counts$Var1 <- factor(category_counts$Var1, levels = c("Major", "Minor", "Inactive"))
colnames(category_counts)[1] <- "Type"
p1 <- ggplot(category_counts, aes(x = 1, y = Freq, fill = Type, width = 0.5)) +
  geom_bar(position = "stack", stat = "identity") +
  labs(y = "Numbers of promoters") +
  scale_fill_manual(values = c("#ff1e56", "#ffac41", "#323232")) +
  scale_color_manual(values = c("#ff1e56", "#ffac41", "#323232")) +
  theme_classic_add +
  theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(), legend.position = "right")
ggsave(p1, file = file.path(output_dir, "promoter_activity_category_percentage.pdf"), width = 2.7, height = 2.9, units = "in", dpi = 300)

# Plot 2: Promoter activity comparison (Fig. 1b)
minor_genes <- absolutePromoterActivityCategory$geneId[absolutePromoterActivityCategory$PromoterCategory == "Minor"]
unique_plot <- absolutePromoterActivityCategory[absolutePromoterActivityCategory$geneId %in% minor_genes, ]
unique_plot$PromoterCategory <- factor(unique_plot$PromoterCategory, levels = c("Major", "Minor", "Inactive"))
p2 <- ggplot(unique_plot, aes(x = PromoterCategory, y = AverageActivity, color = PromoterCategory, fill = PromoterCategory)) +
  geom_boxplot(alpha = 1, outlier.shape = NA, notch = TRUE, width = 0.5) +
  scale_fill_manual(values = c("#ff1e56", "#ffac41", "#323232")) +
  scale_color_manual(values = rep("black", 3)) +
  labs(x = "", y = "Average promoter activity") +
  stat_compare_means(label = "p.signif", method = "wilcox.test", ref.group = "Major", na.rm = TRUE, label.y = max(unique_plot$AverageActivity, na.rm = TRUE) - 1) +
  theme_classic_add + theme(plot.title = element_text(hjust = 0.5), legend.position = "none")
ggsave(p2, file = file.path(output_dir, "promoter_activity_category_comparison.pdf"), width = 2.5, height = 2.5, units = "in", dpi = 300)

# Plot 3: Promoter position category (Fig. 1d)
gene_split <- split(absolutePromoterActivityCategory$PromoterCategory, absolutePromoterActivityCategory$geneId)
active_promoter_position <- function(z, number) {
  if (any(z %in% "Major")) {
    position <- if (length(z) >= number) z[number] else "None"
  } else {
    position <- "None"
  }
  return(position)
}
positions <- list(
  pos1 = table(sapply(gene_split, function(x) active_promoter_position(x, 1))),
  pos2 = table(sapply(gene_split, function(x) active_promoter_position(x, 2))),
  pos3 = table(sapply(gene_split, function(x) active_promoter_position(x, 3))),
  pos4 = table(sapply(gene_split, function(x) active_promoter_position(x, 4))),
  pos5 = table(unlist(sapply(gene_split, function(x) {
    if (length(x) >= 5) x[5:min(length(x), 100)] else "None"
  })))
)
position_plot <- do.call(cbind, positions)
colnames(position_plot) <- 1:5
plot_data <- as.data.frame(position_plot)[c("Major", "Minor"), ]
plot_data_melt <- melt(t(plot_data))
plot_width <- apply(plot_data, 2, sum) / sum(plot_data)
p3 <- ggplot(plot_data_melt, aes(x = Var1, y = value, fill = Var2)) +
  geom_bar(position = "fill", stat = "identity", width = log10(plot_width * 100), colour = rep("black", 5)) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.2)) +
  scale_fill_manual(values = c("#ff1e56", "#ffac41")) +
  theme_classic_add +
  labs(y = "Proportion of promoter types", x = "Promoter position (5' to 3')") +
  theme(legend.title = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())
ggsave(p3, file = file.path(output_dir, "promoter_activity_position_category.pdf"), width = 3, height = 3, units = "in", dpi = 300)

# Plot 4: Gene expression correlation (Fig. 1e)
geneExpression$AverageExpression <- rowMeans(geneExpression[, setdiff(colnames(geneExpression), c("geneId", "promoterId"))], na.rm = TRUE)
major_promoters <- absolutePromoterActivityCategory[absolutePromoterActivityCategory$PromoterCategory == "Major", ]
merged_data <- merge(major_promoters, geneExpression[, c("geneId", "AverageExpression")], by = "geneId")
p4 <- ggplot(merged_data, aes(x = AverageExpression, y = AverageActivity)) +
  geom_point(size = 0.05, colour = "black", alpha = 0.9) +
  xlab("Average gene expression") +
  ylab("Average promoter activity") +
  theme_classic_add + xlim(0, 15) + ylim(0, 15)
ggsave(p4, file = file.path(output_dir, "promoter_activity_geneexpression_correlation.pdf"), width = 3, height = 3, units = "in", dpi = 300)

# Plot 5: Promoter category percentage (genewise) (Fig. 1c alternative)
category_counts_genewise <- as.data.frame(table(absolutePromoterActivityCategory$PromoterCategory))
category_counts_genewise$Var1 <- factor(category_counts_genewise$Var1, levels = c("Major", "Minor", "Inactive"))
category_counts_genewise$Freq <- c(
  length(gene_split) - length(unlist(sapply(gene_split, function(x) grep("Major", x)))) - length(unlist(sapply(gene_split, function(x) grep("Minor", x)))),
  length(unlist(sapply(gene_split, function(x) grep("Major", x)))),
  length(unlist(sapply(gene_split, function(x) grep("Minor", x))))
)
colnames(category_counts_genewise)[1] <- "Type"
p5 <- ggplot(category_counts_genewise[c(2, 3, 1), ], aes(x = 1, y = Freq, fill = Type, width = 0.5)) +
  geom_bar(position = "stack", stat = "identity") +
  labs(y = "Numbers of promoters") +
  scale_fill_manual(values = c("#ff1e56", "#ffac41", "#323232")) +
  scale_color_manual(values = c("#ff1e56", "#ffac41", "#323232")) +
  theme_classic_add +
  theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(), legend.position = "right")
ggsave(p5, file = file.path(output_dir, "promoter_activity_category_percentage_genewise.pdf"), width = 2.7, height = 3, units = "in", dpi = 300)

# Plot 6: Single/multiple promoter categories (Fig. 1f)
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
category_plot <- as.data.frame(table(sapply(gene_split, active_promoter_classification)))
colnames(category_plot) <- c("Category", "Freq")
category_plot$Category <- c("Inactive", "Multiple active promoters\n(Multi promoter genes)", 
                           "Single active promoters\n(Multi promoter genes)", 
                           "Single active promoters\n(Single promoter genes)")
p6 <- ggplot(category_plot[-1, ], aes(x = 1, y = Freq / sum(Freq) * 100, fill = Category)) +
  geom_bar(position = "stack", stat = "identity") +
  labs(y = "% of active promoters") + coord_flip() +
  scale_fill_manual(values = c("#ff6464", "#db3056", "#851d41")) +
  scale_color_manual(values = rep("black", 3)) +
  theme_classic_add +
  theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        legend.position = "top", legend.direction = "vertical", legend.title = element_blank())
ggsave(p6, file = file.path(output_dir, "promoter_activity_single_multiple_category.pdf"), width = 2.5, height = 2.5, units = "in", dpi = 300)

# Plot 7: Promoter number histogram (all) (Fig. S1a)
promoter_id_mapping <- promoterIdMapping(promoterAnnotationData)
gene_id_table <- as.data.frame(table(unique(promoter_id_mapping[, c("promoterId", "geneId")])$geneId))
gene_id_table[gene_id_table$Freq > 11, "Freq"] <- 11
plot_data <- as.data.frame(table(gene_id_table$Freq))
p7 <- ggplot(plot_data, aes(x = Var1, y = Freq)) +
  geom_bar(stat = "identity", color = "black", fill = "white") +
  geom_text(aes(label = Freq), vjust = -0.2, size = 2.5) +
  scale_x_discrete(labels = c(1:10, ">=11")) +
  labs(x = "Promoter number", y = "Number of genes") +
  ggtitle("Promoter number distribution") + theme_classic_add + theme(plot.title = element_text(hjust = 0.5))
ggsave(p7, file = file.path(output_dir, "promoter_activity_number_hist_all.pdf"), width = 3.5, height = 3, units = "in", dpi = 300)

# Plot 8: Promoter number histogram (without single promoter)
plot_data2 <- plot_data[-1, ]
p8 <- ggplot(plot_data2, aes(x = Var1, y = Freq)) +
  geom_bar(stat = "identity", color = "black", fill = "white") +
  geom_text(aes(label = Freq), vjust = -0.2, size = 2.5) +
  scale_x_discrete(labels = c(2:10, ">=11")) +
  labs(x = "Promoter number", y = "Number of genes") +
  ggtitle("Promoter number distribution") + theme_classic() + theme(plot.title = element_text(hjust = 0.5))
ggsave(p8, file = file.path(output_dir, "promoter_activity_number_hist_without1.pdf"), width = 3.5, height = 3, units = "in", dpi = 300)