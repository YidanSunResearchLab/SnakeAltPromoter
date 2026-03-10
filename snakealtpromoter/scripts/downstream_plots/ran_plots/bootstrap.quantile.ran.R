# conda activate /mnt/citadel2/research/shared/SnakeAltPromoter_0228/genome/.snakemake_conda/22f9e1ae900777c8c256cf2a87359557_

#Run at the start
#load from RPKM_length result
load(file.path("/mnt/citadelb/publication/snakealtpromoter/Result/lengths_promoters.RData"))
#load promoter annotation file
promoter_anno <- readRDS("/mnt/citadel2/research/shared/AltPromoterFlow/genome/organisms/hg38/Annotation/proActiv_promoter_annotation.rds")
cell_lines <- c("Healthy", "Failure")
out_dir <- "/mnt/citadel2/research/shared/SnakeAltPromoter_0228/downstream_fixed/0309bootstrapwcounts"


#Run after clean_mat
cage_counts_ori <- clean_mat(
  readRDS("/mnt/citadel2/research/shared/SnakeAltPromoter_0228/normal6/cage/counts_merged/cage_counts.rds")
)
salmon_counts_ori <- clean_mat(
  readRDS("/mnt/citadel2/research/shared/SnakeAltPromoter_0228/normal6/salmon/counts_merged/merged_promoter_counts.rds")
)
dexseq_counts_ori <- clean_mat(
  readRDS("/mnt/citadel2/research/shared/SnakeAltPromoter_0228/normal6/dexseq/counts_merged/merged_promoter_counts.rds")
)
proactiv_counts_ori <- clean_mat(
  readRDS("/mnt/citadel2/research/shared/SnakeAltPromoter_0228/normal6/proactiv/counts_merged/proactiv_raw_counts.rds")
)


# Run after the previous block
#names of samples
colnames(cage_counts_ori) <- c(
  "CAGE_HEALTHY_1", "CAGE_HEALTHY_2", "CAGE_HEALTHY_3",
  "CAGE_FAILURE_1", "CAGE_FAILURE_2", "CAGE_FAILURE_3", "CAGE_FAILURE_4"
)
colnames(salmon_counts_ori) <- c(
  "RNA_HEALTHY_1", "RNA_HEALTHY_2", "RNA_HEALTHY_3",
  "RNA_FAILURE_1", "RNA_FAILURE_2", "RNA_FAILURE_3", "RNA_FAILURE_4"
)
colnames(dexseq_counts_ori) <- c(
  "RNA_HEALTHY_1", "RNA_HEALTHY_2", "RNA_HEALTHY_3",
  "RNA_FAILURE_1", "RNA_FAILURE_2", "RNA_FAILURE_3", "RNA_FAILURE_4"
)
colnames(proactiv_counts_ori) <- c(
  "RNA_HEALTHY_1", "RNA_HEALTHY_2", "RNA_HEALTHY_3",
  "RNA_FAILURE_1", "RNA_FAILURE_2", "RNA_FAILURE_3", "RNA_FAILURE_4"
)


#Run after to_rpkm_by_len
#use rpkm for Salmon and DEXSeq
salmon_rpkm <- to_rpkm_by_len(salmon_counts_ori, salmon_len_bp)
dexseq_rpkm <- to_rpkm_by_len(dexseq_counts_ori, dex_len_bp)


#Run to get the results
# counts
combine.counts.for.each.dataset(
  cage_counts_ori = cage_counts_ori,
  salmon_counts_ori = salmon_counts_ori,
  dexseq_counts_ori = dexseq_counts_ori,
  proactiv_counts_ori = proactiv_counts_ori,
  cell_lines = cell_lines,
  out_dir = out_dir,
  filename = "counts",
  counts_cutoff = 10,
  trim_prop = 0
)

# rpkm
 combine.counts.for.each.dataset(
   cage_counts_ori = cage_counts_ori,
   salmon_counts_ori = salmon_rpkm,
   dexseq_counts_ori = dexseq_rpkm,
   proactiv_counts_ori = proactiv_counts_ori,
   cell_lines = cell_lines,
   out_dir = out_dir,
   filename = "rpkm",
   counts_cutoff = 5,
   trim_prop = 0
 )

cat("Heart-only scatter + bootstrap + stratified analysis finished.\n")