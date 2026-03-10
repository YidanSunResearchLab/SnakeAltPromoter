#Run at the start
# load result from RPKM_length.R
load(file.path("/mnt/citadelb/publication/snakealtpromoter/Result/lengths_promoters.RData"))
# load promoter annotation file
promoter_anno <- readRDS("/mnt/citadel2/research/shared/AltPromoterFlow/genome/organisms/hg38/Annotation/proActiv_promoter_annotation.rds")

# Run at the end
if(TRUE){
  out_dir <- "/mnt/citadelb/publication/snakealtpromoter/Result/Heart_Failure_vs_Healthy"
  dexseq <- readRDS("/mnt/citadel2/research/shared/SnakeAltPromoter_0228/normal6/dexseq/differential/comparisons_/Failure_vs_Healthy/Promoter_differential_activity_FDR0_05.rds")
  salmon <- readRDS("/mnt/citadel2/research/shared/SnakeAltPromoter_0228/normal6/RNAseq/salmon/differential/comparisons_/Failure_vs_Healthy/Promoter_differential_activity_FDR0_05.rds")
  proactiv <- readRDS("/mnt/citadel2/research/shared/SnakeAltPromoter_0228/normal6/RNAseq/proactiv/differential/comparisons_/Failure_vs_Healthy/Promoter_differential_activity_FDR0_05.rds")
  cage <- readRDS("/mnt/citadel2/research/shared/SnakeAltPromoter_0228/normal6/CAGE/cage/differential/comparisons_/Failure_vs_Healthy/Promoter_differential_activity_FDR0_05.rds")
  dexseq_all <- readRDS("/mnt/citadel2/research/shared/SnakeAltPromoter_0228/normal6/RNAseq/dexseq/differential/comparisons_/Failure_vs_Healthy/Promoter_differential_activity.rds")
  salmon_all <- readRDS("/mnt/citadel2/research/shared/SnakeAltPromoter_0228/normal6/RNAseq/salmon/differential/comparisons_/Failure_vs_Healthy/Promoter_differential_activity.rds")
  proactiv_all <- readRDS("/mnt/citadel2/research/shared/SnakeAltPromoter_0228/normal6/RNAseq/proactiv/differential/comparisons_/Failure_vs_Healthy/Promoter_differential_activity.rds")
  cage_all <- readRDS("/mnt/citadel2/research/shared/SnakeAltPromoter_0228/normal6/CAGE/cage/differential/comparisons_/Failure_vs_Healthy/Promoter_differential_activity.rds")
  logfc_df     <- rbind(add_method(cage, "CAGE"), add_method(salmon, "Salmon"), add_method(dexseq, "DEXSeq"), add_method(proactiv, "proActiv"))
  logfc_df_all <- rbind(add_method(cage_all, "CAGE"), add_method(salmon_all, "Salmon"), add_method(dexseq_all, "DEXSeq"), add_method(proactiv_all, "proActiv"))

  ##Intronless
  logfc_wide <- pivot_wider(logfc_df[logfc_df$promoterId %in% intronless_ids,],id_cols    = promoterId,names_from = method,values_from = logFC)
  logfc_wide_all <- pivot_wider(logfc_df_all[logfc_df_all$promoterId %in% intronless_ids,], id_cols    = promoterId, names_from = method, values_from = logFC)
  combine.logFC.for.each.dataset(logfc_wide, logfc_wide_all, out_dir, filename="nointron")

  ##With Intron
  logfc_wide <- pivot_wider(logfc_df[logfc_df$promoterId %in% intron_ids,],id_cols    = promoterId,names_from = method,values_from = logFC)
  logfc_wide_all <- pivot_wider(logfc_df_all[logfc_df_all$promoterId %in% intron_ids,], id_cols    = promoterId, names_from = method, values_from = logFC)
  combine.logFC.for.each.dataset(logfc_wide, logfc_wide_all, out_dir, filename="withintron")

}

