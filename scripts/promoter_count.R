df <- read.table("/mnt/citadel2/research/shared/AltPromoterFlow/data/processed/GSE147236/GSE147236_RNA_Seq_Differential_Genes.txt", header = TRUE, sep = "\t")
gene_ids <- df$GeneID[df$FDR_Pvalue<0.05]
gene_ids_up <- df$GeneID[df$FDR_Pvalue<0.05 & df$logFC<0] 
gene_ids_down <- df$GeneID[df$FDR_Pvalue<0.05 & df$logFC>0]


salmon <- readRDS("/mnt/citadel2/research/shared/AltPromoterFlow/RNA_HEART/salmon/compare/DESeq2_promoter_allResults_keepall.rds")
salmon_down <- unique(salmon$geneId[salmon$logFC < 0 & salmon$FDR < 0.05])
salmon_up <- unique(salmon$geneId[salmon$logFC > 0 & salmon$FDR < 0.05])
total_salmon <- union(salmon_down, salmon_up)

CAGE <- readRDS("/mnt/citadel2/research/shared/AltPromoterFlow/CAGE_HEART/salmon/compare/DESeq2_promoter_allResults_keepall.rds")
CAGE_down <- unique(CAGE$geneId[CAGE$logFC < 0 & CAGE$FDR < 0.05])
CAGE_up <- unique(CAGE$geneId[CAGE$logFC > 0 & CAGE$FDR < 0.05])
total_CAGE <- union(CAGE_down, CAGE_up)

out_dir <- "/mnt/citadel2/research/shared/AltPromoterFlow/RNA_HEART/comparison"

#----------------------------------
## Venn diagram of pval less than 0.05 and both T FALSE
#----------------------------------
library(ggvenn)
library(ggplot2)
library(gridExtra)
library(grid)

# Set fill colors
fill_colors <- c("#b3cde3", "#fbb4ae") 

# Create list objects
venn_total <- list(Paper = gene_ids, CAGE = total_CAGE)
venn_up    <- list(Paper = gene_ids_up, CAGE = CAGE_up)
venn_down  <- list(Paper = gene_ids_down, CAGE = CAGE_down)

# Define a plot function (no set label, no percentage)
plot_venn <- function(data_list, title_text) {
  ggvenn(
    data_list,
    fill_color = fill_colors,
    stroke_size = 0.5,
    text_size = 5,
    show_percentage = FALSE,
    set_name_size = 4
  ) +
    labs(title = title_text) +
    theme(plot.title = element_text(size = 10, hjust = 0.5))
}

# Generate plots
p_up    <- plot_venn(venn_up,    "Upregulated Promoters")
p_down  <- plot_venn(venn_down,  "Downregulated Promoters")
p_total <- plot_venn(venn_total, "All Differential Promoters")

# Save combined Venn diagram to one PDF
pdf(file.path(out_dir, "paper_vs_CAGE_differential_promoters_combined.pdf"), width = 10, height = 4)
grid.arrange(
  p_up, p_down, p_total, ncol = 3,
  top = textGrob("Paper vs CAGE (Venn Diagrams for Promoters)",
                 gp = gpar(fontsize = 10, fontface = "bold"))
)
dev.off()

cat("Saved combined Venn diagram to:", file.path(out_dir, "paper_vs_CAGE_differential_promoters_combined.pdf"), "\n")



#----------------------------------
## Venn diagram of pval less than 0.05 and both T FALSE
#----------------------------------
library(ggvenn)
library(ggplot2)
library(gridExtra)
library(grid)

# Set fill colors
fill_colors <- c("#b3cde3", "#fbb4ae") 

# Create list objects
venn_total <- list(Paper = gene_ids, RNAseq = total_salmon)
venn_up    <- list(Paper = gene_ids_up, RNAseq = salmon_up)
venn_down  <- list(Paper = gene_ids_down, RNAseq = salmon_down)

# Define a plot function (no set label, no percentage)
plot_venn <- function(data_list, title_text) {
  ggvenn(
    data_list,
    fill_color = fill_colors,
    stroke_size = 0.5,
    text_size = 5,
    show_percentage = FALSE,
    set_name_size = 4
  ) +
    labs(title = title_text) +
    theme(plot.title = element_text(size = 10, hjust = 0.5))
}

# Generate plots
p_up    <- plot_venn(venn_up,    "Upregulated Promoters")
p_down  <- plot_venn(venn_down,  "Downregulated Promoters")
p_total <- plot_venn(venn_total, "All Differential Promoters")

# Save combined Venn diagram to one PDF
pdf(file.path(out_dir, "paper_vs_RNAseq_differential_promoters_combined.pdf"), width = 10, height = 4)
grid.arrange(
  p_up, p_down, p_total, ncol = 3,
  top = textGrob("Paper vs RNAseq (Venn Diagrams for Promoters)",
                 gp = gpar(fontsize = 10, fontface = "bold"))
)
dev.off()

cat("Saved combined Venn diagram to:", file.path(out_dir, "paper_vs_RNAseq_differential_promoters_combined.pdf"), "\n")


