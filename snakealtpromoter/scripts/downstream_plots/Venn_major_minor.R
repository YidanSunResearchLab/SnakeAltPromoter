# Load packages required for Venn plotting
library(ggvenn)
library(ggplot2)
library(ggVennDiagram)
library(VennDiagram)

build.overlap.function = function(CAGE_minor, CAGE_major,
                                  salmon_minor, salmon_major,
                                  dexseq_minor, dexseq_major,
                                  proactiv_minor, proactiv_major,
                                  out_dir){
    dir.create(out_dir)

    # Major overlaps
    cage_salmon <- intersect(CAGE_major, salmon_major)
    cage_proactiv <- intersect(CAGE_major, proactiv_major)
    cage_dexseq <- intersect(CAGE_major, dexseq_major)
    shared_all <- Reduce(intersect, list(salmon_major, proactiv_major, dexseq_major, CAGE_major))

    write.table(cage_salmon, file = file.path(out_dir, "cage_salmon_major.txt"), quote = FALSE, row.names = FALSE, col.names = FALSE)
    write.table(cage_proactiv, file = file.path(out_dir, "cage_proactiv_major.txt"), quote = FALSE, row.names = FALSE, col.names = FALSE)
    write.table(cage_dexseq, file = file.path(out_dir, "cage_dexseq_major.txt"), quote = FALSE, row.names = FALSE, col.names = FALSE)
    write.table(shared_all, file = file.path(out_dir, "cage_all3_major.txt"), quote = FALSE, row.names = FALSE, col.names = FALSE)

    major_list_s <- list(
      Salmon = salmon_major,
      CAGE = CAGE_major
    )

    major_list_d <- list(
      DEXSeq = dexseq_major,
      CAGE = CAGE_major
    )

    major_list_p <- list(
      proActiv = proactiv_major,
      CAGE = CAGE_major
    )

    # Minor overlaps
    cage_salmon <- intersect(CAGE_minor, salmon_minor)
    cage_proactiv <- intersect(CAGE_minor, proactiv_minor)
    cage_dexseq <- intersect(CAGE_minor, dexseq_minor)
    shared_all <- Reduce(intersect, list(salmon_minor, proactiv_minor, dexseq_minor, CAGE_minor))

    write.table(cage_salmon, file = file.path(out_dir, "cage_salmon_minor.txt"), quote = FALSE, row.names = FALSE, col.names = FALSE)
    write.table(cage_proactiv, file = file.path(out_dir, "cage_proactiv_minor.txt"), quote = FALSE, row.names = FALSE, col.names = FALSE)
    write.table(cage_dexseq, file = file.path(out_dir, "cage_dexseq_minor.txt"), quote = FALSE, row.names = FALSE, col.names = FALSE)
    write.table(shared_all, file = file.path(out_dir, "cage_all3_minor.txt"), quote = FALSE, row.names = FALSE, col.names = FALSE)

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

      # Pairwise Venn for two sets
      area1 <- length(lst[[1]])
      area2 <- length(lst[[2]])
      cross <- length(intersect(lst[[1]], lst[[2]]))

      out_file <- file.path(out_dir, filename)
      pdf(out_file, width = w, height = h, family = "Times")
      grid::grid.newpage()
      VennDiagram::draw.pairwise.venn(
        area1 = area1,
        area2 = area2,
        cross.area = cross,
        category = names(lst),
        fill = c("#e95280","#FFDBB6"),
        alpha = rep(0.5, 2),
        lwd = 1.2,
        cex = 2.2,
        cat.cex = 1.6,
        cat.pos = c(-10, 10),
        cat.dist = 0.035,
        cat.fontfamily = "serif",
        fontfamily = "serif"
      )
      dev.off()
    }

    # Major plots
    save_venn_diagram(major_list_s, "sc_major_promoter_venn_VennDiagram.pdf", out_dir)
    save_venn_diagram(major_list_d, "dc_major_promoter_venn_VennDiagram.pdf", out_dir)
    save_venn_diagram(major_list_p, "pc_major_promoter_venn_VennDiagram.pdf", out_dir)

    # Minor plots
    save_venn_diagram(minor_list_s, "sc_minor_promoter_venn_VennDiagram.pdf", out_dir)
    save_venn_diagram(minor_list_d, "dc_minor_promoter_venn_VennDiagram.pdf", out_dir)
    save_venn_diagram(minor_list_p, "pc_minor_promoter_venn_VennDiagram.pdf", out_dir)

    compute_overlap_table <- function(lst) {
      stopifnot(length(lst) == 2)

      set1 <- names(lst)[1]
      set2 <- names(lst)[2]
      s1 <- lst[[1]]
      s2 <- lst[[2]]

      # Summarize overlap size and percentages for one pair
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

    # Collect all pairwise overlap summaries
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