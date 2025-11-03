# activate environment: /mnt/citadel2/research/shared/AltPromoterFlow/genome/.snakemake_conda/3ea3da478f24d1f8a438361e1dd371ca_


#Libraries
library(ggvenn)
library(ggplot2)
library(ggVennDiagram)



build.overlap.function = function(CAGE_minor, CAGE_major,
                                  salmon_minor, salmon_major,
                                  dexseq_minor, dexseq_major,
                                  proactiv_minor, proactiv_major,
                                  out_dir){
    dir.create(out_dir)
    #Find common
    cage_salmon <- intersect(CAGE_major, salmon_major)
    cage_proactiv <- intersect(CAGE_major, proactiv_major)
    cage_dexseq <- intersect(CAGE_major, dexseq_major)
    shared_all <- Reduce(intersect, list(salmon_major, proactiv_major, dexseq_major, CAGE_major))
    write.table(cage_salmon, file = file.path(out_dir, "cage_salmon_major.txt"), quote = FALSE, row.names = FALSE, col.names = FALSE)
    write.table(cage_proactiv, file = file.path(out_dir, "cage_proactiv_major.txt"), quote = FALSE, row.names = FALSE, col.names = FALSE)
    write.table(cage_dexseq, file = file.path(out_dir, "cage_dexseq_major.txt"), quote = FALSE, row.names = FALSE, col.names = FALSE)
    write.table(shared_all, file = file.path(out_dir, "cage_all3_major.txt"), quote = FALSE, row.names = FALSE, col.names = FALSE)


    #Major list
    major_list_s <- list(
      Salmon = salmon_major,
      CAGE = CAGE_major
    )

    #Major list
    major_list_d <- list(
      DEXSeq = dexseq_major,
      CAGE = CAGE_major
    )

    #Major list
    major_list_p <- list(
      proActiv = proactiv_major,
      CAGE = CAGE_major
    )


    #Find common
    cage_salmon <- intersect(CAGE_minor, salmon_minor)
    cage_proactiv <- intersect(CAGE_minor, proactiv_minor)
    cage_dexseq <- intersect(CAGE_minor, dexseq_minor)
    shared_all <- Reduce(intersect, list(salmon_minor, proactiv_minor, dexseq_minor, CAGE_minor))
    write.table(cage_salmon, file = file.path(out_dir, "cage_salmon_minor.txt"), quote = FALSE, row.names = FALSE, col.names = FALSE)
    write.table(cage_proactiv, file = file.path(out_dir, "cage_proactiv_minor.txt"), quote = FALSE, row.names = FALSE, col.names = FALSE)
    write.table(cage_dexseq, file = file.path(out_dir, "cage_dexseq_minor.txt"), quote = FALSE, row.names = FALSE, col.names = FALSE)
    write.table(shared_all, file = file.path(out_dir, "cage_all3_minor.txt"), quote = FALSE, row.names = FALSE, col.names = FALSE)

    #Minor list
    minor_list_s <- list(
      Salmon = salmon_minor,
      CAGE = CAGE_minor
    )

    minor_list_d <- list(
      DEXSeq = dexseq_minor,
      CAGE = CAGE_minor
    )

    minor_list_p <- list(
      proActiv = proactiv_minor,
      CAGE = CAGE_minor
    )

    library(VennDiagram)
    save_venn_diagram <- function(lst, filename, out_dir, w = 5, h = 5) {
      stopifnot(length(lst) == 2)
      lst <- lapply(lst, as.character)

      # 数值计算
      area1 <- length(lst[[1]])
      area2 <- length(lst[[2]])
      cross <- length(intersect(lst[[1]], lst[[2]]))

      # 输出路径
      out_file <- file.path(out_dir, filename)
      pdf(out_file, width = w, height = h, family = "Times")  # 尝试 Times New Roman
      grid::grid.newpage()
      VennDiagram::draw.pairwise.venn(
        area1 = area1,
        area2 = area2,
        cross.area = cross,
        category = names(lst),
        fill = c("#e95280","#FFDBB6"),
        alpha = rep(0.5, 2),
        lwd = 1.2,

        cex = 2.2,              # ⬅️ 区域数字字体更大（默认是 1.5）
        cat.cex = 1.6,          # 类别标签大小
        cat.pos = c(-10, 10),
        cat.dist = 0.035,
        cat.fontfamily = "serif",  # serif字体组，兼容Times
        fontfamily = "serif"       # 区域数字字体也用 serif
      )
      dev.off()
    }




    # --- Major plots ---
    save_venn_diagram(major_list_s, "sc_major_promoter_venn_VennDiagram.pdf", out_dir)
    save_venn_diagram(major_list_d, "dc_major_promoter_venn_VennDiagram.pdf", out_dir)
    save_venn_diagram(major_list_p, "pc_major_promoter_venn_VennDiagram.pdf", out_dir)

    # --- Minor plots ---
    save_venn_diagram(minor_list_s, "sc_minor_promoter_venn_VennDiagram.pdf", out_dir)
    save_venn_diagram(minor_list_d, "dc_minor_promoter_venn_VennDiagram.pdf", out_dir)
    save_venn_diagram(minor_list_p, "pc_minor_promoter_venn_VennDiagram.pdf", out_dir)



    compute_overlap_table <- function(lst) {
      stopifnot(length(lst) == 2)

      set1 <- names(lst)[1]
      set2 <- names(lst)[2]
      s1 <- lst[[1]]
      s2 <- lst[[2]]

      overlap <- length(intersect(s1, s2))
      len1 <- length(s1)
      len2 <- length(s2)

      data.frame(
        Set1 = set1,
        Set2 = set2,
        Overlap = overlap,
        Size_Set1 = len1,
        Size_Set2 = len2,
        Percent_in_Set1 = sprintf("%.1f%%", 100 * overlap / len1),
        Percent_in_Set2 = sprintf("%.1f%%", 100 * overlap / len2),
        stringsAsFactors = FALSE
      )
    }


    all_tables <- dplyr::bind_rows(
      compute_overlap_table(major_list_s) |> dplyr::mutate(Group = "major_s"),
      compute_overlap_table(major_list_d) |> dplyr::mutate(Group = "major_d"),
      compute_overlap_table(major_list_p) |> dplyr::mutate(Group = "major_p"),
      compute_overlap_table(minor_list_s) |> dplyr::mutate(Group = "minor_s"),
      compute_overlap_table(minor_list_d) |> dplyr::mutate(Group = "minor_d"),
      compute_overlap_table(minor_list_p) |> dplyr::mutate(Group = "minor_p")
    )

    readr::write_tsv(all_tables, file.path(out_dir, "overlap_summary.tsv"))


}



#Load files
CAGE_minor <- readRDS("/mnt/citadelb/publication/snakealtpromoter/Processed/heart_AltPromoter/CAGE/cage/promoter_classification_total/Minor_promoterId_overall.rds")
CAGE_major <- readRDS("/mnt/citadelb/publication/snakealtpromoter/Processed/heart_AltPromoter/CAGE/cage/promoter_classification_total/Major_promoterId_overall.rds")
salmon_minor <- readRDS("/mnt/citadelb/publication/snakealtpromoter/Processed/heart_AltPromoter/RNAseq/salmon/promoter_classification_total/Minor_promoterId_overall.rds")
salmon_major <- readRDS("/mnt/citadelb/publication/snakealtpromoter/Processed/heart_AltPromoter/RNAseq/salmon/promoter_classification_total/Major_promoterId_overall.rds")
dexseq_minor <- readRDS("/mnt/citadelb/publication/snakealtpromoter/Processed/heart_AltPromoter/RNAseq/dexseq/promoter_classification_total/Minor_promoterId_overall.rds")
dexseq_major <- readRDS("/mnt/citadelb/publication/snakealtpromoter/Processed/heart_AltPromoter/RNAseq/dexseq/promoter_classification_total/Major_promoterId_overall.rds")
proactiv_major <- readRDS("/mnt/citadelb/publication/snakealtpromoter/Processed/heart_AltPromoter/RNAseq/proactiv/promoter_classification_total/Major_promoterId_overall.rds")
proactiv_minor <- readRDS("/mnt/citadelb/publication/snakealtpromoter/Processed/heart_AltPromoter/RNAseq/proactiv/promoter_classification_total/Minor_promoterId_overall.rds")
out_dir <- "/mnt/citadelb/publication/snakealtpromoter/Result/Heart_Failure_vs_Healthy"
build.overlap.function(CAGE_minor, CAGE_major,salmon_minor, salmon_major,dexseq_minor, dexseq_major,proactiv_minor, proactiv_major,out_dir)



#Load files
CAGE_minor <- readRDS("/mnt/citadelb/publication/snakealtpromoter/Processed/cellline_AltPromoter/CAGE/cage/promoter_classification_total/Minor_promoterId_overall.rds")
CAGE_major <- readRDS("/mnt/citadelb/publication/snakealtpromoter/Processed/cellline_AltPromoter/CAGE/cage/promoter_classification_total/Major_promoterId_overall.rds")
salmon_minor <- readRDS("/mnt/citadelb/publication/snakealtpromoter/Processed/cellline_AltPromoter/RNAseq/salmon/promoter_classification_total/Minor_promoterId_overall.rds")
salmon_major <- readRDS("/mnt/citadelb/publication/snakealtpromoter/Processed/cellline_AltPromoter/RNAseq/salmon/promoter_classification_total/Major_promoterId_overall.rds")
dexseq_minor <- readRDS("/mnt/citadelb/publication/snakealtpromoter/Processed/cellline_AltPromoter/RNAseq/dexseq/promoter_classification_total/Minor_promoterId_overall.rds")
dexseq_major <- readRDS("/mnt/citadelb/publication/snakealtpromoter/Processed/cellline_AltPromoter/RNAseq/dexseq/promoter_classification_total/Major_promoterId_overall.rds")
proactiv_major <- readRDS("/mnt/citadelb/publication/snakealtpromoter/Processed/cellline_AltPromoter/RNAseq/proactiv/promoter_classification_total/Major_promoterId_overall.rds")
proactiv_minor <- readRDS("/mnt/citadelb/publication/snakealtpromoter/Processed/cellline_AltPromoter/RNAseq/proactiv/promoter_classification_total/Minor_promoterId_overall.rds")
out_dir <- "/mnt/citadelb/publication/snakealtpromoter/Result/K562_vs_GM12878"
build.overlap.function(CAGE_minor, CAGE_major,salmon_minor, salmon_major,dexseq_minor, dexseq_major,proactiv_minor, proactiv_major,out_dir)



#Load files
CAGE_minor <- readRDS("/mnt/citadelb/publication/snakealtpromoter/Processed/cellline_AltPromoter/CAGE/cage/promoter_classification_condition_wise/Minor_promoterId_K562.rds")
CAGE_major <- readRDS("/mnt/citadelb/publication/snakealtpromoter/Processed/cellline_AltPromoter/CAGE/cage/promoter_classification_condition_wise/Major_promoterId_K562.rds")
salmon_minor <- readRDS("/mnt/citadelb/publication/snakealtpromoter/Processed/cellline_AltPromoter/RNAseq/salmon/promoter_classification_condition_wise/Minor_promoterId_K562.rds")
salmon_major <- readRDS("/mnt/citadelb/publication/snakealtpromoter/Processed/cellline_AltPromoter/RNAseq/salmon/promoter_classification_condition_wise/Major_promoterId_K562.rds")
dexseq_minor <- readRDS("/mnt/citadelb/publication/snakealtpromoter/Processed/cellline_AltPromoter/RNAseq/dexseq/promoter_classification_condition_wise/Minor_promoterId_K562.rds")
dexseq_major <- readRDS("/mnt/citadelb/publication/snakealtpromoter/Processed/cellline_AltPromoter/RNAseq/dexseq/promoter_classification_condition_wise/Major_promoterId_K562.rds")
proactiv_major <- readRDS("/mnt/citadelb/publication/snakealtpromoter/Processed/cellline_AltPromoter/RNAseq/proactiv/promoter_classification_condition_wise/Major_promoterId_K562.rds")
proactiv_minor <- readRDS("/mnt/citadelb/publication/snakealtpromoter/Processed/cellline_AltPromoter/RNAseq/proactiv/promoter_classification_condition_wise/Minor_promoterId_K562.rds")
out_dir <- "/mnt/citadelb/publication/snakealtpromoter/Result/K562_vs_GM12878/K562"
build.overlap.function(CAGE_minor, CAGE_major,salmon_minor, salmon_major,dexseq_minor, dexseq_major,proactiv_minor, proactiv_major,out_dir)


#Load files
CAGE_minor <- readRDS("/mnt/citadelb/publication/snakealtpromoter/Processed/cellline_AltPromoter/CAGE/cage/promoter_classification_condition_wise/Minor_promoterId_GM12878.rds")
CAGE_major <- readRDS("/mnt/citadelb/publication/snakealtpromoter/Processed/cellline_AltPromoter/CAGE/cage/promoter_classification_condition_wise/Major_promoterId_GM12878.rds")
salmon_minor <- readRDS("/mnt/citadelb/publication/snakealtpromoter/Processed/cellline_AltPromoter/RNAseq/salmon/promoter_classification_condition_wise/Minor_promoterId_GM12878.rds")
salmon_major <- readRDS("/mnt/citadelb/publication/snakealtpromoter/Processed/cellline_AltPromoter/RNAseq/salmon/promoter_classification_condition_wise/Major_promoterId_GM12878.rds")
dexseq_minor <- readRDS("/mnt/citadelb/publication/snakealtpromoter/Processed/cellline_AltPromoter/RNAseq/dexseq/promoter_classification_condition_wise/Minor_promoterId_GM12878.rds")
dexseq_major <- readRDS("/mnt/citadelb/publication/snakealtpromoter/Processed/cellline_AltPromoter/RNAseq/dexseq/promoter_classification_condition_wise/Major_promoterId_GM12878.rds")
proactiv_major <- readRDS("/mnt/citadelb/publication/snakealtpromoter/Processed/cellline_AltPromoter/RNAseq/proactiv/promoter_classification_condition_wise/Major_promoterId_GM12878.rds")
proactiv_minor <- readRDS("/mnt/citadelb/publication/snakealtpromoter/Processed/cellline_AltPromoter/RNAseq/proactiv/promoter_classification_condition_wise/Minor_promoterId_GM12878.rds")
out_dir <- "/mnt/citadelb/publication/snakealtpromoter/Result/K562_vs_GM12878/GM12878"
build.overlap.function(CAGE_minor, CAGE_major,salmon_minor, salmon_major,dexseq_minor, dexseq_major,proactiv_minor, proactiv_major,out_dir)



#/home/yuqing/shared/SnakeAltPromoter_paper_revision/processed_brain/dexseq/promoter_classification_total/Major_promoterId_overall.rds
#Load files
CAGE_minor <- readRDS("/mnt/citadelb/publication/snakealtpromoter/Processed/brain_AltPromoter/CAGE_test/cage/promoter_classification_total/Minor_promoterId_overall.rds")
CAGE_major <- readRDS("/mnt/citadelb/publication/snakealtpromoter/Processed/brain_AltPromoter/CAGE_test/cage/promoter_classification_total/Major_promoterId_overall.rds")
salmon_minor <- readRDS("/mnt/citadelb/publication/snakealtpromoter/Processed/brain_AltPromoter/RNAseq_test/salmon/promoter_classification_total/Minor_promoterId_overall.rds")
salmon_major <- readRDS("/mnt/citadelb/publication/snakealtpromoter/Processed/brain_AltPromoter/RNAseq_test/salmon/promoter_classification_total/Major_promoterId_overall.rds")
dexseq_minor <- readRDS("/mnt/citadelb/publication/snakealtpromoter/Processed/brain_AltPromoter/RNAseq_test/dexseq/promoter_classification_total/Minor_promoterId_overall.rds")
dexseq_major <- readRDS("/mnt/citadelb/publication/snakealtpromoter/Processed/brain_AltPromoter/RNAseq_test/dexseq/promoter_classification_total/Major_promoterId_overall.rds")
proactiv_major <- readRDS("/mnt/citadelb/publication/snakealtpromoter/Processed/brain_AltPromoter/RNAseq_test/proactiv/promoter_classification_total/Major_promoterId_overall.rds")
proactiv_minor <- readRDS("/mnt/citadelb/publication/snakealtpromoter/Processed/brain_AltPromoter/RNAseq_test/proactiv/promoter_classification_total/Minor_promoterId_overall.rds")
out_dir <- "/mnt/citadelb/publication/snakealtpromoter/Result/Brain"
build.overlap.function(CAGE_minor, CAGE_major,salmon_minor, salmon_major,dexseq_minor, dexseq_major,proactiv_minor, proactiv_major,out_dir)

#/home/yuqing/shared/SnakeAltPromoter_paper_revision/processed_brain/dexseq/promoter_classification_condition_wise/Major_promoterId_female.rds
CAGE_minor <- readRDS("/mnt/citadelb/publication/snakealtpromoter/Processed/brain_AltPromoter/CAGE_test/cage/promoter_classification_condition_wise/Minor_promoterId_female.rds")
CAGE_major <- readRDS("/mnt/citadelb/publication/snakealtpromoter/Processed/brain_AltPromoter/CAGE_test/cage/promoter_classification_condition_wise/Major_promoterId_female.rds")
salmon_minor <- readRDS("/mnt/citadelb/publication/snakealtpromoter/Processed/brain_AltPromoter/RNAseq_test/salmon/promoter_classification_condition_wise/Minor_promoterId_female.rds")
salmon_major <- readRDS("/mnt/citadelb/publication/snakealtpromoter/Processed/brain_AltPromoter/RNAseq_test/salmon/promoter_classification_condition_wise/Major_promoterId_female.rds")
dexseq_minor <- readRDS("/mnt/citadelb/publication/snakealtpromoter/Processed/brain_AltPromoter/RNAseq_test/dexseq/promoter_classification_condition_wise/Minor_promoterId_female.rds")
dexseq_major <- readRDS("/mnt/citadelb/publication/snakealtpromoter/Processed/brain_AltPromoter/RNAseq_test/dexseq/promoter_classification_condition_wise/Major_promoterId_female.rds")
proactiv_major <- readRDS("/mnt/citadelb/publication/snakealtpromoter/Processed/brain_AltPromoter/RNAseq_test/proactiv/promoter_classification_condition_wise/Major_promoterId_female.rds")
proactiv_minor <- readRDS("/mnt/citadelb/publication/snakealtpromoter/Processed/brain_AltPromoter/RNAseq_test/proactiv/promoter_classification_condition_wise/Minor_promoterId_female.rds")
out_dir <- "/mnt/citadelb/publication/snakealtpromoter/Result/Brain/female"
build.overlap.function(CAGE_minor, CAGE_major,salmon_minor, salmon_major,dexseq_minor, dexseq_major,proactiv_minor, proactiv_major,out_dir)

#/home/yuqing/shared/SnakeAltPromoter_paper_revision/processed_brain/dexseq/promoter_classification_condition_wise/Major_promoterId_female.rds
CAGE_minor <- readRDS("/mnt/citadelb/publication/snakealtpromoter/Processed/brain_AltPromoter/CAGE_test/cage/promoter_classification_condition_wise/Minor_promoterId_male.rds")
CAGE_major <- readRDS("/mnt/citadelb/publication/snakealtpromoter/Processed/brain_AltPromoter/CAGE_test/cage/promoter_classification_condition_wise/Major_promoterId_male.rds")
salmon_minor <- readRDS("/mnt/citadelb/publication/snakealtpromoter/Processed/brain_AltPromoter/RNAseq_test/salmon/promoter_classification_condition_wise/Minor_promoterId_male.rds")
salmon_major <- readRDS("/mnt/citadelb/publication/snakealtpromoter/Processed/brain_AltPromoter/RNAseq_test/salmon/promoter_classification_condition_wise/Major_promoterId_male.rds")
dexseq_minor <- readRDS("/mnt/citadelb/publication/snakealtpromoter/Processed/brain_AltPromoter/RNAseq_test/dexseq/promoter_classification_condition_wise/Minor_promoterId_male.rds")
dexseq_major <- readRDS("/mnt/citadelb/publication/snakealtpromoter/Processed/brain_AltPromoter/RNAseq_test/dexseq/promoter_classification_condition_wise/Major_promoterId_male.rds")
proactiv_major <- readRDS("/mnt/citadelb/publication/snakealtpromoter/Processed/brain_AltPromoter/RNAseq_test/proactiv/promoter_classification_condition_wise/Major_promoterId_male.rds")
proactiv_minor <- readRDS("/mnt/citadelb/publication/snakealtpromoter/Processed/brain_AltPromoter/RNAseq_test/proactiv/promoter_classification_condition_wise/Minor_promoterId_male.rds")
out_dir <- "/mnt/citadelb/publication/snakealtpromoter/Result/Brain/male"
build.overlap.function(CAGE_minor, CAGE_major,salmon_minor, salmon_major,dexseq_minor, dexseq_major,proactiv_minor, proactiv_major,out_dir)







# library(ggplot2)
# library(gridExtra)

# save_overlap_table <- function(lst, filename_prefix, out_dir) {
#   stopifnot(length(lst) == 2)
#   lst <- lapply(lst, as.character)

#   set1 <- names(lst)[1]
#   set2 <- names(lst)[2]
#   s1 <- lst[[1]]
#   s2 <- lst[[2]]
#   inter <- intersect(s1, s2)

#   overlap <- length(inter)
#   len1 <- length(s1)
#   len2 <- length(s2)
#   percent1 <- round(overlap / len1 * 100, 2)
#   percent2 <- round(overlap / len2 * 100, 2)

#   df <- data.frame(
#     Group = filename_prefix,
#     Set1 = set1,
#     Set2 = set2,
#     Overlap = overlap,
#     Size_Set1 = len1,
#     Size_Set2 = len2,
#     Percent_in_Set1 = paste0(percent1, "%"),
#     Percent_in_Set2 = paste0(percent2, "%")
#   )

#   out_tsv <- file.path(out_dir, paste0(filename_prefix, "_overlap_summary.tsv"))
#   write.table(df, file = out_tsv, sep = "\t", row.names = FALSE, quote = FALSE)

#   out_pdf <- file.path(out_dir, paste0(filename_prefix, "_overlap_summary.pdf"))
#   pdf(out_pdf, width = 6, height = 2, family = "Times")
#   gridExtra::grid.table(df, rows = NULL, theme = ttheme_default(
#     core = list(fg_params = list(fontfamily = "Times", fontsize = 11)),
#     colhead = list(fg_params = list(fontfamily = "Times", fontsize = 12, fontface = "bold"))
#   ))
#   dev.off()
# }



# anno <- readRDS("/mnt/citadel2/research/shared/AltPromoterFlow/genome/organisms/hg38/Annotation/proActiv_promoter_annotation.rds")
# pc <- anno@promoterCoordinates
# internal_map <- as.data.frame(mcols(pc))[, c("promoterId", "internalPromoter")]

# CAGE_minor     <- readRDS("/mnt/citadelb/publication/snakealtpromoter/Processed/heart_AltPromoter/RNAseq/cage/promoter_classification_total/Minor_promoterId_overall.rds")
# CAGE_major     <- readRDS("/mnt/citadelb/publication/snakealtpromoter/Processed/heart_AltPromoter/RNAseq/cage/promoter_classification_total/Major_promoterId_overall.rds")
# salmon_minor   <- readRDS("/mnt/citadelb/publication/snakealtpromoter/Processed/heart_AltPromoter/RNAseq/salmon/promoter_classification_total/Minor_promoterId_overall.rds")
# salmon_major   <- readRDS("/mnt/citadelb/publication/snakealtpromoter/Processed/heart_AltPromoter/RNAseq/salmon/promoter_classification_total/Major_promoterId_overall.rds")
# dexseq_minor   <- readRDS("/mnt/citadelb/publication/snakealtpromoter/Processed/heart_AltPromoter/RNAseq/dexseq/promoter_classification_total/Minor_promoterId_overall.rds")
# dexseq_major   <- readRDS("/mnt/citadelb/publication/snakealtpromoter/Processed/heart_AltPromoter/RNAseq/dexseq/promoter_classification_total/Major_promoterId_overall.rds")
# proactiv_minor <- readRDS("/mnt/citadelb/publication/snakealtpromoter/Processed/heart_AltPromoter/RNAseq/proactiv/promoter_classification_total/Minor_promoterId_overall.rds")
# proactiv_major <- readRDS("/mnt/citadelb/publication/snakealtpromoter/Processed/heart_AltPromoter/RNAseq/proactiv/promoter_classification_total/Major_promoterId_overall.rds")

# groups <- list(
#   CAGE_major     = CAGE_major,
#   CAGE_minor     = CAGE_minor,
#   Salmon_major   = salmon_major,
#   Salmon_minor   = salmon_minor,
#   DEXSeq_major   = dexseq_major,
#   DEXSeq_minor   = dexseq_minor,
#   proActiv_major = proactiv_major,
#   proActiv_minor = proactiv_minor
# )

# check_internal_status <- function(ids, group_name) {
#   df <- merge(
#     data.frame(promoterId = ids),
#     internal_map,
#     by = "promoterId",
#     all.x = TRUE
#   )
#   total <- nrow(df)
#   na_count    <- sum(is.na(df$internalPromoter))
#   true_count  <- sum(df$internalPromoter %in% TRUE, na.rm = TRUE)

#   cat(sprintf("%-15s | total: %5d | NA: %4d (%.1f%%) | TRUE: %4d (%.1f%%)\n",
#               group_name, total, na_count, 100 * na_count / total,
#               true_count, 100 * true_count / total))

#   write.table(df$promoterId[is.na(df$internalPromoter)],
#               file = paste0(group_name, "_internalNA.txt"),
#               row.names = FALSE, col.names = FALSE, quote = FALSE)

#   write.table(df$promoterId[df$internalPromoter %in% TRUE],
#               file = paste0(group_name, "_internalTRUE.txt"),
#               row.names = FALSE, col.names = FALSE, quote = FALSE)
# }

# cat("Group           | Total |  NA   |  TRUE\n")
# cat("---------------------------------------------\n")
# invisible(lapply(names(groups), function(g) {
#   check_internal_status(groups[[g]], g)
# }))














# > dexseq  <- readRDS("/mnt/citadelb/publication/snakealtpromoter/Processed/heart_AltPromoter/RNAseq/dexseq/promoter_classification_total/Summary_classified.rds")
# > proactiv <- readRDS("/mnt/citadelb/publication/snakealtpromoter/Processed/heart_AltPromoter/RNAseq/proactiv/promoter_classification_total/Summary_classified.rds")
# > salmon <- readRDS("/mnt/citadelb/publication/snakealtpromoter/Processed/heart_AltPromoter/RNAseq/salmon/promoter_classification_total/Summary_classified.rds")
# > cage <- readRDS("/mnt/citadelb/publication/snakealtpromoter/Processed/heart_AltPromoter/RNAseq/cage/promoter_classification_total/Summary_classified.rds")

# methods <- list(
#   CAGE     = cage,
#   Salmon   = salmon,
#   DEXSeq   = dexseq,
#   proActiv = proactiv
# )

# inactive_thresholds <- c(0.1, 0.13, 0.15, 0.2, 0.23)
# major_thresholds    <- c(0.27, 0.3, 0.35)

# result <- lapply(names(methods), function(method_name) {
#   df <- methods[[method_name]]

#   inactive_sub <- df[df$overall.class == "Inactive", ]
#   major_sub    <- df[df$overall.class == "Major",    ]


#   inact_counts <- sapply(inactive_thresholds, function(th) {
#     sum(inactive_sub$overall.mean > th, na.rm = TRUE)
#   })
#   major_counts <- sapply(major_thresholds, function(th) {
#     sum(major_sub$overall.mean < th, na.rm = TRUE)
#   })


#   inact_below_025 <- sum(inactive_sub$overall.mean < 0.25, na.rm = TRUE)
#   major_above_025 <- sum(major_sub$overall.mean > 0.25, na.rm = TRUE)

#   data.frame(
#     Method = method_name,
#     t(inact_counts),
#     t(major_counts),
#     Inactive_lt_0.25 = inact_below_025,
#     Major_gt_0.25    = major_above_025,
#     row.names = NULL
#   )
# })

# final_df <- do.call(rbind, result)
# colnames(final_df) <- c(
#   "Method",
#   paste0("Inactive> ", inactive_thresholds),
#   paste0("Major< ", major_thresholds),
#   "Inactive< 0.25",
#   "Major> 0.25"
# )


# # ----------------------------------
# final_df

#     Method Inactive> 0.1 Inactive> 0.13 Inactive> 0.15 Inactive> 0.2
# 1     CAGE          8151           7160           3679          2662
# 2   Salmon         12986           7422           6866          2917
# 3   DEXSeq         13716           8539           6880          2667
# 4 proActiv          6062           3375           3023          1848
#   Inactive> 0.23 Major< 0.27 Major< 0.3 Major< 0.35 Inactive< 0.25 Major> 0.25
# 1            831         181        729        1421          92558       20508
# 2           1762         757       1689        3064          58628       40997
# 3           1105         529       1661        3006          63684       38508
# 4            484         100        412         903          96491       16422