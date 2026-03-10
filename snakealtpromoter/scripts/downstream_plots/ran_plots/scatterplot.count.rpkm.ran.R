# conda activate /mnt/citadel2/research/shared/SnakeAltPromoter_0228/genome/.snakemake_conda/22f9e1ae900777c8c256cf2a87359557_

#Run in the beginning
# Load result from RPKM_length.R
load(file.path("/mnt/citadelb/publication/snakealtpromoter/Result/lengths_promoters.RData"))
# Load promoter annotation file
promoter_anno <- readRDS("/mnt/citadel2/research/shared/AltPromoterFlow/genome/organisms/hg38/Annotation/proActiv_promoter_annotation.rds")

#Run at the end
if(TRUE){
  cell_lines <- c("Healthy", "Failure")
  out_dir <- "/mnt/citadel2/research/shared/SnakeAltPromoter_0228/downstream/0309newscatter/update"
  cage_counts_ori     <- clean_mat(readRDS("/mnt/citadel2/research/shared/SnakeAltPromoter_0228/normal6/cage/counts_merged/cage_counts.rds"))
  salmon_counts_ori   <- clean_mat(readRDS("/mnt/citadel2/research/shared/SnakeAltPromoter_0228/normal6/salmon/counts_merged/merged_promoter_counts.rds"))
  dexseq_counts_ori   <- clean_mat(readRDS("/mnt/citadel2/research/shared/SnakeAltPromoter_0228/normal6/dexseq/counts_merged/merged_promoter_counts.rds"))
  proactiv_counts_ori <- clean_mat(readRDS("/mnt/citadel2/research/shared/SnakeAltPromoter_0228/normal6/proactiv/counts_merged/proactiv_raw_counts.rds"))
  # sample names
  colnames(cage_counts_ori) <- c("CAGE_HEALTHY_1", "CAGE_HEALTHY_2", "CAGE_HEALTHY_3", "CAGE_FAILURE_1", "CAGE_FAILURE_2", "CAGE_FAILURE_3", "CAGE_FAILURE_4")
  colnames(salmon_counts_ori) <- c("RNA_HEALTHY_1", "RNA_HEALTHY_2", "RNA_HEALTHY_3","RNA_FAILURE_1", "RNA_FAILURE_2", "RNA_FAILURE_3", "RNA_FAILURE_4")
  colnames(dexseq_counts_ori) <- c("RNA_HEALTHY_1", "RNA_HEALTHY_2", "RNA_HEALTHY_3","RNA_FAILURE_1", "RNA_FAILURE_2", "RNA_FAILURE_3", "RNA_FAILURE_4")
  colnames(proactiv_counts_ori) <- c("RNA_HEALTHY_1", "RNA_HEALTHY_2", "RNA_HEALTHY_3","RNA_FAILURE_1", "RNA_FAILURE_2", "RNA_FAILURE_3", "RNA_FAILURE_4")
  salmon_rpkm <- to_rpkm_by_len(salmon_counts_ori, salmon_len_bp)
  dexseq_rpkm <- to_rpkm_by_len(dexseq_counts_ori, dex_len_bp)
  combine.counts.for.each.dataset(cage_counts_ori, salmon_counts_ori, dexseq_counts_ori, proactiv_counts_ori, cell_lines, out_dir, filename="counts", counts_cutoff=10, trim_prop=0) 
  combine.counts.for.each.dataset(cage_counts_ori, salmon_rpkm, dexseq_rpkm, proactiv_counts_ori, cell_lines, out_dir, filename="rpkm", counts_cutoff=5, trim_prop=0) 
}


if(TRUE){
  cell_lines <- c("Healthy", "Failure")
  out_dir <- "/mnt/citadel2/research/shared/SnakeAltPromoter_0228/downstream/normal4"
  cage_counts_ori     <- clean_mat(readRDS("/mnt/citadel2/research/shared/SnakeAltPromoter_0228/normal6/cage/counts_merged/cage_counts.rds"))
  salmon_counts_ori   <- clean_mat(readRDS("/mnt/citadel2/research/shared/SnakeAltPromoter_0228/normal6/salmon/counts_merged/merged_promoter_counts.rds"))
  dexseq_counts_ori   <- clean_mat(readRDS("/mnt/citadel2/research/shared/SnakeAltPromoter_0228/normal6/dexseq/counts_merged/merged_promoter_counts.rds"))
  proactiv_counts_ori <- clean_mat(readRDS("/mnt/citadel2/research/shared/SnakeAltPromoter_0228/normal6/proactiv/counts_merged/proactiv_raw_counts.rds"))
  colnames(cage_counts_ori) <- c("CAGE_HEALTHY_1", "CAGE_HEALTHY_2", "CAGE_HEALTHY_3", "CAGE_FAILURE_1", "CAGE_FAILURE_2", "CAGE_FAILURE_3", "CAGE_FAILURE_4")
  colnames(salmon_counts_ori) <- c("RNA_HEALTHY_1", "RNA_HEALTHY_2", "RNA_HEALTHY_3","RNA_FAILURE_1", "RNA_FAILURE_2", "RNA_FAILURE_3", "RNA_FAILURE_4")
  colnames(dexseq_counts_ori) <- c("RNA_HEALTHY_1", "RNA_HEALTHY_2", "RNA_HEALTHY_3","RNA_FAILURE_1", "RNA_FAILURE_2", "RNA_FAILURE_3", "RNA_FAILURE_4")
  colnames(proactiv_counts_ori) <- c("RNA_HEALTHY_1", "RNA_HEALTHY_2", "RNA_HEALTHY_3","RNA_FAILURE_1", "RNA_FAILURE_2", "RNA_FAILURE_3", "RNA_FAILURE_4")
  salmon_rpkm <- to_rpkm_by_len(salmon_counts_ori, salmon_len_bp)
  dexseq_rpkm <- to_rpkm_by_len(dexseq_counts_ori, dex_len_bp)
  combine.counts.for.each.dataset(cage_counts_ori, salmon_counts_ori, dexseq_counts_ori, proactiv_counts_ori, cell_lines, out_dir, filename="counts", counts_cutoff=10, trim_prop=0) 
  combine.counts.for.each.dataset(cage_counts_ori, salmon_rpkm, dexseq_rpkm, proactiv_counts_ori, cell_lines, out_dir, filename="rpkm", counts_cutoff=5, trim_prop=0) 

}

if(TRUE){
  cell_lines <- c("Healthy", "Failure")
  out_dir <- "/mnt/citadel2/research/shared/SnakeAltPromoter_0228/downstream_fixed/downsample/10M"
  cage_counts_ori     <- clean_mat(readRDS("/mnt/citadel2/research/shared/SnakeAltPromoter_0228/normal6/cage/counts_merged/cage_counts.rds"))
  salmon_counts_ori   <- clean_mat(readRDS("/mnt/citadel2/research/shared/SnakeAltPromoter_0228/downsample/10M/salmon/counts_merged/merged_promoter_counts.rds"))
  dexseq_counts_ori   <- clean_mat(readRDS("/mnt/citadel2/research/shared/SnakeAltPromoter_0228/downsample/10M/dexseq/counts_merged/merged_promoter_counts.rds"))
  proactiv_counts_ori <- clean_mat(readRDS("/mnt/citadel2/research/shared/SnakeAltPromoter_0228/downsample/10M/proactiv/counts_merged/proactiv_raw_counts.rds"))
  colnames(cage_counts_ori) <- c("CAGE_HEALTHY_1", "CAGE_HEALTHY_2", "CAGE_HEALTHY_3", "CAGE_FAILURE_1", "CAGE_FAILURE_2", "CAGE_FAILURE_3", "CAGE_FAILURE_4")
  colnames(salmon_counts_ori) <- c("RNA_HEALTHY_1", "RNA_HEALTHY_2", "RNA_HEALTHY_3","RNA_FAILURE_1", "RNA_FAILURE_2", "RNA_FAILURE_3", "RNA_FAILURE_4")
  colnames(dexseq_counts_ori) <- c("RNA_HEALTHY_1", "RNA_HEALTHY_2", "RNA_HEALTHY_3","RNA_FAILURE_1", "RNA_FAILURE_2", "RNA_FAILURE_3", "RNA_FAILURE_4")
  colnames(proactiv_counts_ori) <- c("RNA_HEALTHY_1", "RNA_HEALTHY_2", "RNA_HEALTHY_3","RNA_FAILURE_1", "RNA_FAILURE_2", "RNA_FAILURE_3", "RNA_FAILURE_4")
  salmon_rpkm <- to_rpkm_by_len(salmon_counts_ori, salmon_len_bp)
  dexseq_rpkm <- to_rpkm_by_len(dexseq_counts_ori, dex_len_bp)
  combine.counts.for.each.dataset(cage_counts_ori, salmon_counts_ori, dexseq_counts_ori, proactiv_counts_ori, cell_lines, out_dir, filename="counts", counts_cutoff=10, trim_prop=0) 
  combine.counts.for.each.dataset(cage_counts_ori, salmon_rpkm, dexseq_rpkm, proactiv_counts_ori, cell_lines, out_dir, filename="rpkm", counts_cutoff=5, trim_prop=0) 

}

if(TRUE){
  cell_lines <- c("Healthy", "Failure")
  out_dir <- "/mnt/citadel2/research/shared/SnakeAltPromoter_0228/downstream_fixed/downsample/20M"
  cage_counts_ori     <- clean_mat(readRDS("/mnt/citadel2/research/shared/SnakeAltPromoter_0228/normal6/cage/counts_merged/cage_counts.rds"))
  salmon_counts_ori   <- clean_mat(readRDS("/mnt/citadel2/research/shared/SnakeAltPromoter_0228/downsample/20M/salmon/counts_merged/merged_promoter_counts.rds"))
  dexseq_counts_ori   <- clean_mat(readRDS("/mnt/citadel2/research/shared/SnakeAltPromoter_0228/downsample/20M/dexseq/counts_merged/merged_promoter_counts.rds"))
  proactiv_counts_ori <- clean_mat(readRDS("/mnt/citadel2/research/shared/SnakeAltPromoter_0228/downsample/20M/proactiv/counts_merged/proactiv_raw_counts.rds"))
  colnames(cage_counts_ori) <- c("CAGE_HEALTHY_1", "CAGE_HEALTHY_2", "CAGE_HEALTHY_3", "CAGE_FAILURE_1", "CAGE_FAILURE_2", "CAGE_FAILURE_3", "CAGE_FAILURE_4")
  colnames(salmon_counts_ori) <- c("RNA_HEALTHY_1", "RNA_HEALTHY_2", "RNA_HEALTHY_3","RNA_FAILURE_1", "RNA_FAILURE_2", "RNA_FAILURE_3", "RNA_FAILURE_4")
  colnames(dexseq_counts_ori) <- c("RNA_HEALTHY_1", "RNA_HEALTHY_2", "RNA_HEALTHY_3","RNA_FAILURE_1", "RNA_FAILURE_2", "RNA_FAILURE_3", "RNA_FAILURE_4")
  colnames(proactiv_counts_ori) <- c("RNA_HEALTHY_1", "RNA_HEALTHY_2", "RNA_HEALTHY_3","RNA_FAILURE_1", "RNA_FAILURE_2", "RNA_FAILURE_3", "RNA_FAILURE_4")
  salmon_rpkm <- to_rpkm_by_len(salmon_counts_ori, salmon_len_bp)
  dexseq_rpkm <- to_rpkm_by_len(dexseq_counts_ori, dex_len_bp)
  combine.counts.for.each.dataset(cage_counts_ori, salmon_counts_ori, dexseq_counts_ori, proactiv_counts_ori, cell_lines, out_dir, filename="counts", counts_cutoff=10, trim_prop=0) 
  combine.counts.for.each.dataset(cage_counts_ori, salmon_rpkm, dexseq_rpkm, proactiv_counts_ori, cell_lines, out_dir, filename="rpkm", counts_cutoff=5, trim_prop=0) 

}


if(TRUE){
  cell_lines <- c("Healthy", "Failure")
  out_dir <- "/mnt/citadel2/research/shared/SnakeAltPromoter_0228/downstream_fixed/downsample/30M"
  cage_counts_ori     <- clean_mat(readRDS("/mnt/citadel2/research/shared/SnakeAltPromoter_0228/normal6/cage/counts_merged/cage_counts.rds"))
  salmon_counts_ori   <- clean_mat(readRDS("/mnt/citadel2/research/shared/SnakeAltPromoter_0228/downsample/30M/salmon/counts_merged/merged_promoter_counts.rds"))
  dexseq_counts_ori   <- clean_mat(readRDS("/mnt/citadel2/research/shared/SnakeAltPromoter_0228/downsample/30M/dexseq/counts_merged/merged_promoter_counts.rds"))
  proactiv_counts_ori <- clean_mat(readRDS("/mnt/citadel2/research/shared/SnakeAltPromoter_0228/downsample/30M/proactiv/counts_merged/proactiv_raw_counts.rds"))
  colnames(cage_counts_ori) <- c("CAGE_HEALTHY_1", "CAGE_HEALTHY_2", "CAGE_HEALTHY_3", "CAGE_FAILURE_1", "CAGE_FAILURE_2", "CAGE_FAILURE_3", "CAGE_FAILURE_4")
  colnames(salmon_counts_ori) <- c("RNA_HEALTHY_1", "RNA_HEALTHY_2", "RNA_HEALTHY_3","RNA_FAILURE_1", "RNA_FAILURE_2", "RNA_FAILURE_3", "RNA_FAILURE_4")
  colnames(dexseq_counts_ori) <- c("RNA_HEALTHY_1", "RNA_HEALTHY_2", "RNA_HEALTHY_3","RNA_FAILURE_1", "RNA_FAILURE_2", "RNA_FAILURE_3", "RNA_FAILURE_4")
  colnames(proactiv_counts_ori) <- c("RNA_HEALTHY_1", "RNA_HEALTHY_2", "RNA_HEALTHY_3","RNA_FAILURE_1", "RNA_FAILURE_2", "RNA_FAILURE_3", "RNA_FAILURE_4")
  salmon_rpkm <- to_rpkm_by_len(salmon_counts_ori, salmon_len_bp)
  dexseq_rpkm <- to_rpkm_by_len(dexseq_counts_ori, dex_len_bp)
  combine.counts.for.each.dataset(cage_counts_ori, salmon_counts_ori, dexseq_counts_ori, proactiv_counts_ori, cell_lines, out_dir, filename="counts", counts_cutoff=10, trim_prop=0) 
  combine.counts.for.each.dataset(cage_counts_ori, salmon_rpkm, dexseq_rpkm, proactiv_counts_ori, cell_lines, out_dir, filename="rpkm", counts_cutoff=5, trim_prop=0) 

}

if(TRUE){
  cell_lines <- c("Healthy", "Failure")
  out_dir <- "/mnt/citadel2/research/shared/SnakeAltPromoter_0228/downstream_fixed/downsample/40M"
  cage_counts_ori     <- clean_mat(readRDS("/mnt/citadel2/research/shared/SnakeAltPromoter_0228/normal6/cage/counts_merged/cage_counts.rds"))
  salmon_counts_ori   <- clean_mat(readRDS("/mnt/citadel2/research/shared/SnakeAltPromoter_0228/downsample/40M/salmon/counts_merged/merged_promoter_counts.rds"))
  dexseq_counts_ori   <- clean_mat(readRDS("/mnt/citadel2/research/shared/SnakeAltPromoter_0228/downsample/40M/dexseq/counts_merged/merged_promoter_counts.rds"))
  proactiv_counts_ori <- clean_mat(readRDS("/mnt/citadel2/research/shared/SnakeAltPromoter_0228/downsample/40M/proactiv/counts_merged/proactiv_raw_counts.rds"))
  colnames(cage_counts_ori) <- c("CAGE_HEALTHY_1", "CAGE_HEALTHY_2", "CAGE_HEALTHY_3", "CAGE_FAILURE_1", "CAGE_FAILURE_2", "CAGE_FAILURE_3", "CAGE_FAILURE_4")
  colnames(salmon_counts_ori) <- c("RNA_HEALTHY_1", "RNA_HEALTHY_2", "RNA_HEALTHY_3","RNA_FAILURE_1", "RNA_FAILURE_2", "RNA_FAILURE_3", "RNA_FAILURE_4")
  colnames(dexseq_counts_ori) <- c("RNA_HEALTHY_1", "RNA_HEALTHY_2", "RNA_HEALTHY_3","RNA_FAILURE_1", "RNA_FAILURE_2", "RNA_FAILURE_3", "RNA_FAILURE_4")
  colnames(proactiv_counts_ori) <- c("RNA_HEALTHY_1", "RNA_HEALTHY_2", "RNA_HEALTHY_3","RNA_FAILURE_1", "RNA_FAILURE_2", "RNA_FAILURE_3", "RNA_FAILURE_4")
  salmon_rpkm <- to_rpkm_by_len(salmon_counts_ori, salmon_len_bp)
  dexseq_rpkm <- to_rpkm_by_len(dexseq_counts_ori, dex_len_bp)
  combine.counts.for.each.dataset(cage_counts_ori, salmon_counts_ori, dexseq_counts_ori, proactiv_counts_ori, cell_lines, out_dir, filename="counts", counts_cutoff=10, trim_prop=0) 
  combine.counts.for.each.dataset(cage_counts_ori, salmon_rpkm, dexseq_rpkm, proactiv_counts_ori, cell_lines, out_dir, filename="rpkm", counts_cutoff=5, trim_prop=0) 

}

if(TRUE){
  cell_lines <- c("Healthy", "Failure")
  out_dir <- "/mnt/citadel2/research/shared/SnakeAltPromoter_0228/downstream_fixed/downsample/80M"
  cage_counts_ori     <- clean_mat(readRDS("/mnt/citadel2/research/shared/SnakeAltPromoter_0228/normal6/cage/counts_merged/cage_counts.rds"))
  salmon_counts_ori   <- clean_mat(readRDS("/mnt/citadel2/research/shared/SnakeAltPromoter_0228/downsample/80M/salmon/counts_merged/merged_promoter_counts.rds"))
  dexseq_counts_ori   <- clean_mat(readRDS("/mnt/citadel2/research/shared/SnakeAltPromoter_0228/downsample/80M/dexseq/counts_merged/merged_promoter_counts.rds"))
  proactiv_counts_ori <- clean_mat(readRDS("/mnt/citadel2/research/shared/SnakeAltPromoter_0228/downsample/80M/proactiv/counts_merged/proactiv_raw_counts.rds"))
  colnames(cage_counts_ori) <- c("CAGE_HEALTHY_1", "CAGE_HEALTHY_2", "CAGE_HEALTHY_3", "CAGE_FAILURE_1", "CAGE_FAILURE_2", "CAGE_FAILURE_3", "CAGE_FAILURE_4")
  colnames(salmon_counts_ori) <- c("RNA_HEALTHY_1", "RNA_HEALTHY_2", "RNA_HEALTHY_3","RNA_FAILURE_1", "RNA_FAILURE_2", "RNA_FAILURE_3", "RNA_FAILURE_4")
  colnames(dexseq_counts_ori) <- c("RNA_HEALTHY_1", "RNA_HEALTHY_2", "RNA_HEALTHY_3","RNA_FAILURE_1", "RNA_FAILURE_2", "RNA_FAILURE_3", "RNA_FAILURE_4")
  colnames(proactiv_counts_ori) <- c("RNA_HEALTHY_1", "RNA_HEALTHY_2", "RNA_HEALTHY_3","RNA_FAILURE_1", "RNA_FAILURE_2", "RNA_FAILURE_3", "RNA_FAILURE_4")
  salmon_rpkm <- to_rpkm_by_len(salmon_counts_ori, salmon_len_bp)
  dexseq_rpkm <- to_rpkm_by_len(dexseq_counts_ori, dex_len_bp)
  combine.counts.for.each.dataset(cage_counts_ori, salmon_counts_ori, dexseq_counts_ori, proactiv_counts_ori, cell_lines, out_dir, filename="counts", counts_cutoff=10, trim_prop=0) 
  combine.counts.for.each.dataset(cage_counts_ori, salmon_rpkm, dexseq_rpkm, proactiv_counts_ori, cell_lines, out_dir, filename="rpkm", counts_cutoff=5, trim_prop=0) 

}

if(TRUE){
  cell_lines <- c("Healthy", "Failure")
  out_dir <- "/mnt/citadel2/research/shared/SnakeAltPromoter_0228/downstream_fixed/readlensim/50"
  cage_counts_ori     <- clean_mat(readRDS("/mnt/citadel2/research/shared/SnakeAltPromoter_0228/normal6/cage/counts_merged/cage_counts.rds"))
  salmon_counts_ori   <- clean_mat(readRDS("/mnt/citadel2/research/shared/SnakeAltPromoter_0228/readlensim/50/salmon/counts_merged/merged_promoter_counts.rds"))
  dexseq_counts_ori   <- clean_mat(readRDS("/mnt/citadel2/research/shared/SnakeAltPromoter_0228/readlensim/50/dexseq/counts_merged/merged_promoter_counts.rds"))
  proactiv_counts_ori <- clean_mat(readRDS("/mnt/citadel2/research/shared/SnakeAltPromoter_0228/readlensim/50/proactiv/counts_merged/proactiv_raw_counts.rds"))
  colnames(cage_counts_ori) <- c("CAGE_HEALTHY_1", "CAGE_HEALTHY_2", "CAGE_HEALTHY_3", "CAGE_FAILURE_1", "CAGE_FAILURE_2", "CAGE_FAILURE_3", "CAGE_FAILURE_4")
  colnames(salmon_counts_ori) <- c("RNA_HEALTHY_1", "RNA_HEALTHY_2", "RNA_HEALTHY_3","RNA_FAILURE_1", "RNA_FAILURE_2", "RNA_FAILURE_3", "RNA_FAILURE_4")
  colnames(dexseq_counts_ori) <- c("RNA_HEALTHY_1", "RNA_HEALTHY_2", "RNA_HEALTHY_3","RNA_FAILURE_1", "RNA_FAILURE_2", "RNA_FAILURE_3", "RNA_FAILURE_4")
  colnames(proactiv_counts_ori) <- c("RNA_HEALTHY_1", "RNA_HEALTHY_2", "RNA_HEALTHY_3","RNA_FAILURE_1", "RNA_FAILURE_2", "RNA_FAILURE_3", "RNA_FAILURE_4")
  salmon_rpkm <- to_rpkm_by_len(salmon_counts_ori, salmon_len_bp)
  dexseq_rpkm <- to_rpkm_by_len(dexseq_counts_ori, dex_len_bp)
  combine.counts.for.each.dataset(cage_counts_ori, salmon_counts_ori, dexseq_counts_ori, proactiv_counts_ori, cell_lines, out_dir, filename="counts", counts_cutoff=10, trim_prop=0) 
  combine.counts.for.each.dataset(cage_counts_ori, salmon_rpkm, dexseq_rpkm, proactiv_counts_ori, cell_lines, out_dir, filename="rpkm", counts_cutoff=5, trim_prop=0) 

}

if(TRUE){
  cell_lines <- c("Healthy", "Failure")
  out_dir <- "/mnt/citadel2/research/shared/SnakeAltPromoter_0228/downstream_fixed/readlensim/75"
  cage_counts_ori     <- clean_mat(readRDS("/mnt/citadel2/research/shared/SnakeAltPromoter_0228/normal6/cage/counts_merged/cage_counts.rds"))
  salmon_counts_ori   <- clean_mat(readRDS("/mnt/citadel2/research/shared/SnakeAltPromoter_0228/readlensim/75/salmon/counts_merged/merged_promoter_counts.rds"))
  dexseq_counts_ori   <- clean_mat(readRDS("/mnt/citadel2/research/shared/SnakeAltPromoter_0228/readlensim/75/dexseq/counts_merged/merged_promoter_counts.rds"))
  proactiv_counts_ori <- clean_mat(readRDS("/mnt/citadel2/research/shared/SnakeAltPromoter_0228/readlensim/75/proactiv/counts_merged/proactiv_raw_counts.rds"))
  colnames(cage_counts_ori) <- c("CAGE_HEALTHY_1", "CAGE_HEALTHY_2", "CAGE_HEALTHY_3", "CAGE_FAILURE_1", "CAGE_FAILURE_2", "CAGE_FAILURE_3", "CAGE_FAILURE_4")
  colnames(salmon_counts_ori) <- c("RNA_HEALTHY_1", "RNA_HEALTHY_2", "RNA_HEALTHY_3","RNA_FAILURE_1", "RNA_FAILURE_2", "RNA_FAILURE_3", "RNA_FAILURE_4")
  colnames(dexseq_counts_ori) <- c("RNA_HEALTHY_1", "RNA_HEALTHY_2", "RNA_HEALTHY_3","RNA_FAILURE_1", "RNA_FAILURE_2", "RNA_FAILURE_3", "RNA_FAILURE_4")
  colnames(proactiv_counts_ori) <- c("RNA_HEALTHY_1", "RNA_HEALTHY_2", "RNA_HEALTHY_3","RNA_FAILURE_1", "RNA_FAILURE_2", "RNA_FAILURE_3", "RNA_FAILURE_4")
  salmon_rpkm <- to_rpkm_by_len(salmon_counts_ori, salmon_len_bp)
  dexseq_rpkm <- to_rpkm_by_len(dexseq_counts_ori, dex_len_bp)
  combine.counts.for.each.dataset(cage_counts_ori, salmon_counts_ori, dexseq_counts_ori, proactiv_counts_ori, cell_lines, out_dir, filename="counts", counts_cutoff=10, trim_prop=0) 
  combine.counts.for.each.dataset(cage_counts_ori, salmon_rpkm, dexseq_rpkm, proactiv_counts_ori, cell_lines, out_dir, filename="rpkm", counts_cutoff=5, trim_prop=0) 

}

if(TRUE){
  cell_lines <- c("Healthy", "Failure")
  out_dir <- "/mnt/citadel2/research/shared/SnakeAltPromoter_0228/downstream_fixed/readlensim/100"
  cage_counts_ori     <- clean_mat(readRDS("/mnt/citadel2/research/shared/SnakeAltPromoter_0228/normal6/cage/counts_merged/cage_counts.rds"))
  salmon_counts_ori   <- clean_mat(readRDS("/mnt/citadel2/research/shared/SnakeAltPromoter_0228/readlensim/100/salmon/counts_merged/merged_promoter_counts.rds"))
  dexseq_counts_ori   <- clean_mat(readRDS("/mnt/citadel2/research/shared/SnakeAltPromoter_0228/readlensim/100/dexseq/counts_merged/merged_promoter_counts.rds"))
  proactiv_counts_ori <- clean_mat(readRDS("/mnt/citadel2/research/shared/SnakeAltPromoter_0228/readlensim/100/proactiv/counts_merged/proactiv_raw_counts.rds"))
  colnames(cage_counts_ori) <- c("CAGE_HEALTHY_1", "CAGE_HEALTHY_2", "CAGE_HEALTHY_3", "CAGE_FAILURE_1", "CAGE_FAILURE_2", "CAGE_FAILURE_3", "CAGE_FAILURE_4")
  colnames(salmon_counts_ori) <- c("RNA_HEALTHY_1", "RNA_HEALTHY_2", "RNA_HEALTHY_3","RNA_FAILURE_1", "RNA_FAILURE_2", "RNA_FAILURE_3", "RNA_FAILURE_4")
  colnames(dexseq_counts_ori) <- c("RNA_HEALTHY_1", "RNA_HEALTHY_2", "RNA_HEALTHY_3","RNA_FAILURE_1", "RNA_FAILURE_2", "RNA_FAILURE_3", "RNA_FAILURE_4")
  colnames(proactiv_counts_ori) <- c("RNA_HEALTHY_1", "RNA_HEALTHY_2", "RNA_HEALTHY_3","RNA_FAILURE_1", "RNA_FAILURE_2", "RNA_FAILURE_3", "RNA_FAILURE_4")
  salmon_rpkm <- to_rpkm_by_len(salmon_counts_ori, salmon_len_bp)
  dexseq_rpkm <- to_rpkm_by_len(dexseq_counts_ori, dex_len_bp)
  combine.counts.for.each.dataset(cage_counts_ori, salmon_counts_ori, dexseq_counts_ori, proactiv_counts_ori, cell_lines, out_dir, filename="counts", counts_cutoff=10, trim_prop=0) 
  combine.counts.for.each.dataset(cage_counts_ori, salmon_rpkm, dexseq_rpkm, proactiv_counts_ori, cell_lines, out_dir, filename="rpkm", counts_cutoff=5, trim_prop=0) 

}

if(TRUE){
  cell_lines <- c("Healthy", "Failure")
  out_dir <- "/mnt/citadel2/research/shared/SnakeAltPromoter_0228/downstream_fixed/readlensim/150"
  cage_counts_ori     <- clean_mat(readRDS("/mnt/citadel2/research/shared/SnakeAltPromoter_0228/normal6/cage/counts_merged/cage_counts.rds"))
  salmon_counts_ori   <- clean_mat(readRDS("/mnt/citadel2/research/shared/SnakeAltPromoter_0228/readlensim/150/salmon/counts_merged/merged_promoter_counts.rds"))
  dexseq_counts_ori   <- clean_mat(readRDS("/mnt/citadel2/research/shared/SnakeAltPromoter_0228/readlensim/150/dexseq/counts_merged/merged_promoter_counts.rds"))
  proactiv_counts_ori <- clean_mat(readRDS("/mnt/citadel2/research/shared/SnakeAltPromoter_0228/readlensim/150/proactiv/counts_merged/proactiv_raw_counts.rds"))
  colnames(cage_counts_ori) <- c("CAGE_HEALTHY_1", "CAGE_HEALTHY_2", "CAGE_HEALTHY_3", "CAGE_FAILURE_1", "CAGE_FAILURE_2", "CAGE_FAILURE_3", "CAGE_FAILURE_4")
  colnames(salmon_counts_ori) <- c("RNA_HEALTHY_1", "RNA_HEALTHY_2", "RNA_HEALTHY_3","RNA_FAILURE_1", "RNA_FAILURE_2", "RNA_FAILURE_3", "RNA_FAILURE_4")
  colnames(dexseq_counts_ori) <- c("RNA_HEALTHY_1", "RNA_HEALTHY_2", "RNA_HEALTHY_3","RNA_FAILURE_1", "RNA_FAILURE_2", "RNA_FAILURE_3", "RNA_FAILURE_4")
  colnames(proactiv_counts_ori) <- c("RNA_HEALTHY_1", "RNA_HEALTHY_2", "RNA_HEALTHY_3","RNA_FAILURE_1", "RNA_FAILURE_2", "RNA_FAILURE_3", "RNA_FAILURE_4")
  salmon_rpkm <- to_rpkm_by_len(salmon_counts_ori, salmon_len_bp)
  dexseq_rpkm <- to_rpkm_by_len(dexseq_counts_ori, dex_len_bp)
  combine.counts.for.each.dataset(cage_counts_ori, salmon_counts_ori, dexseq_counts_ori, proactiv_counts_ori, cell_lines, out_dir, filename="counts", counts_cutoff=10, trim_prop=0) 
  combine.counts.for.each.dataset(cage_counts_ori, salmon_rpkm, dexseq_rpkm, proactiv_counts_ori, cell_lines, out_dir, filename="rpkm", counts_cutoff=5, trim_prop=0) 

}

if(TRUE){
  cell_lines <- c("Healthy", "Failure")
  out_dir <- "/mnt/citadel2/research/shared/SnakeAltPromoter_0228/downstream/threshold/0.1"
  cage_counts_ori     <- clean_mat(readRDS("/mnt/citadel2/research/shared/SnakeAltPromoter_0228/normal6/cage/counts_merged/cage_counts.rds"))
  salmon_counts_ori   <- clean_mat(readRDS("/mnt/citadel2/research/shared/SnakeAltPromoter_0228/threshold0.1/salmon/counts_merged/merged_promoter_counts.rds"))
  dexseq_counts_ori   <- clean_mat(readRDS("/mnt/citadel2/research/shared/SnakeAltPromoter_0228/threshold0.1/dexseq/counts_merged/merged_promoter_counts.rds"))
  proactiv_counts_ori <- clean_mat(readRDS("/mnt/citadel2/research/shared/SnakeAltPromoter_0228/threshold0.1/proactiv/counts_merged/proactiv_raw_counts.rds"))
  colnames(cage_counts_ori) <- c("CAGE_HEALTHY_1", "CAGE_HEALTHY_2", "CAGE_HEALTHY_3", "CAGE_FAILURE_1", "CAGE_FAILURE_2", "CAGE_FAILURE_3", "CAGE_FAILURE_4")
  colnames(salmon_counts_ori) <- c("RNA_HEALTHY_1", "RNA_HEALTHY_2", "RNA_HEALTHY_3","RNA_FAILURE_1", "RNA_FAILURE_2", "RNA_FAILURE_3", "RNA_FAILURE_4")
  colnames(dexseq_counts_ori) <- c("RNA_HEALTHY_1", "RNA_HEALTHY_2", "RNA_HEALTHY_3","RNA_FAILURE_1", "RNA_FAILURE_2", "RNA_FAILURE_3", "RNA_FAILURE_4")
  colnames(proactiv_counts_ori) <- c("RNA_HEALTHY_1", "RNA_HEALTHY_2", "RNA_HEALTHY_3","RNA_FAILURE_1", "RNA_FAILURE_2", "RNA_FAILURE_3", "RNA_FAILURE_4")
  salmon_rpkm <- to_rpkm_by_len(salmon_counts_ori, salmon_len_bp)
  dexseq_rpkm <- to_rpkm_by_len(dexseq_counts_ori, dex_len_bp)
  combine.counts.for.each.dataset(cage_counts_ori, salmon_counts_ori, dexseq_counts_ori, proactiv_counts_ori, cell_lines, out_dir, filename="counts", counts_cutoff=10, trim_prop=0) 
  combine.counts.for.each.dataset(cage_counts_ori, salmon_rpkm, dexseq_rpkm, proactiv_counts_ori, cell_lines, out_dir, filename="rpkm", counts_cutoff=5, trim_prop=0) 

}

if(TRUE){
  cell_lines <- c("Healthy", "Failure")
  out_dir <- "/mnt/citadel2/research/shared/SnakeAltPromoter_0228/downstream/threshold/0.5"
  cage_counts_ori     <- clean_mat(readRDS("/mnt/citadel2/research/shared/SnakeAltPromoter_0228/normal6/cage/counts_merged/cage_counts.rds"))
  salmon_counts_ori   <- clean_mat(readRDS("/mnt/citadel2/research/shared/SnakeAltPromoter_0228/threshold0.5/salmon/counts_merged/merged_promoter_counts.rds"))
  dexseq_counts_ori   <- clean_mat(readRDS("/mnt/citadel2/research/shared/SnakeAltPromoter_0228/threshold0.5/dexseq/counts_merged/merged_promoter_counts.rds"))
  proactiv_counts_ori <- clean_mat(readRDS("/mnt/citadel2/research/shared/SnakeAltPromoter_0228/threshold0.5/proactiv/counts_merged/proactiv_raw_counts.rds"))
  colnames(cage_counts_ori) <- c("CAGE_HEALTHY_1", "CAGE_HEALTHY_2", "CAGE_HEALTHY_3", "CAGE_FAILURE_1", "CAGE_FAILURE_2", "CAGE_FAILURE_3", "CAGE_FAILURE_4")
  colnames(salmon_counts_ori) <- c("RNA_HEALTHY_1", "RNA_HEALTHY_2", "RNA_HEALTHY_3","RNA_FAILURE_1", "RNA_FAILURE_2", "RNA_FAILURE_3", "RNA_FAILURE_4")
  colnames(dexseq_counts_ori) <- c("RNA_HEALTHY_1", "RNA_HEALTHY_2", "RNA_HEALTHY_3","RNA_FAILURE_1", "RNA_FAILURE_2", "RNA_FAILURE_3", "RNA_FAILURE_4")
  colnames(proactiv_counts_ori) <- c("RNA_HEALTHY_1", "RNA_HEALTHY_2", "RNA_HEALTHY_3","RNA_FAILURE_1", "RNA_FAILURE_2", "RNA_FAILURE_3", "RNA_FAILURE_4")
  salmon_rpkm <- to_rpkm_by_len(salmon_counts_ori, salmon_len_bp)
  dexseq_rpkm <- to_rpkm_by_len(dexseq_counts_ori, dex_len_bp)
  combine.counts.for.each.dataset(cage_counts_ori, salmon_counts_ori, dexseq_counts_ori, proactiv_counts_ori, cell_lines, out_dir, filename="counts", counts_cutoff=10, trim_prop=0) 
  combine.counts.for.each.dataset(cage_counts_ori, salmon_rpkm, dexseq_rpkm, proactiv_counts_ori, cell_lines, out_dir, filename="rpkm", counts_cutoff=5, trim_prop=0) 

}



#--------------------------------------------------------
###For heart healthy and failure
#--------------------------------------------------------
if(TRUE){
  cell_lines <- c("Healthy", "Failure")
  out_dir <- "/mnt/citadelb/publication/snakealtpromoter/Result/Heart_Failure_vs_Healthy"
  # cage_counts_ori     <- clean_mat(readRDS("/mnt/citadel2/research/shared/AltPromoterFlow/Heart/cage/merge/raw_promoter_counts.rds"))
  # salmon_counts_ori   <- clean_mat(readRDS("/mnt/citadel2/research/shared/AltPromoterFlow/Heart/salmon/merge/raw_promoter_counts.rds"))
  # dexseq_counts_ori   <- clean_mat(readRDS("/mnt/citadel2/research/shared/AltPromoterFlow/Heart/dexseq/merge/raw_promoter_counts.rds"))
  # proactiv_counts_ori <- clean_mat(readRDS("/mnt/citadel2/research/shared/AltPromoterFlow/Heart/proactiv/merge/raw_promoter_counts.rds"))
  cage_counts_ori     <- clean_mat(readRDS("/mnt/citadelb/publication/snakealtpromoter/Processed/heart_AltPromoter/CAGE/cage/counts_merged/cage_counts.rds"))
  salmon_counts_ori   <- clean_mat(readRDS("/mnt/citadelb/publication/snakealtpromoter/Processed/heart_AltPromoter/RNAseq/salmon/counts_merged/merged_promoter_counts.rds"))
  dexseq_counts_ori   <- clean_mat(readRDS("/mnt/citadelb/publication/snakealtpromoter/Processed/heart_AltPromoter/RNAseq/dexseq/counts_merged/merged_promoter_counts.rds"))
  proactiv_counts_ori <- clean_mat(readRDS("/mnt/citadelb/publication/snakealtpromoter/Processed/heart_AltPromoter/RNAseq/proactiv/counts_merged/proactiv_raw_counts.rds"))
  colnames(cage_counts_ori) <- c("CAGE_HEALTHY_1", "CAGE_HEALTHY_2", "CAGE_HEALTHY_3", "CAGE_FAILURE_1", "CAGE_FAILURE_2", "CAGE_FAILURE_3", "CAGE_FAILURE_4")
  colnames(salmon_counts_ori) <- c("RNA_HEALTHY_1", "RNA_HEALTHY_2", "RNA_HEALTHY_3","RNA_FAILURE_1", "RNA_FAILURE_2", "RNA_FAILURE_3", "RNA_FAILURE_4")
  colnames(dexseq_counts_ori) <- c("RNA_HEALTHY_1", "RNA_HEALTHY_2", "RNA_HEALTHY_3","RNA_FAILURE_1", "RNA_FAILURE_2", "RNA_FAILURE_3", "RNA_FAILURE_4")
  colnames(proactiv_counts_ori) <- c("RNA_HEALTHY_1", "RNA_HEALTHY_2", "RNA_HEALTHY_3","RNA_FAILURE_1", "RNA_FAILURE_2", "RNA_FAILURE_3", "RNA_FAILURE_4")
  salmon_rpkm <- to_rpkm_by_len(salmon_counts_ori, salmon_len_bp)
  dexseq_rpkm <- to_rpkm_by_len(dexseq_counts_ori, dex_len_bp)
  combine.counts.for.each.dataset(cage_counts_ori, salmon_counts_ori, dexseq_counts_ori, proactiv_counts_ori, cell_lines, out_dir, filename="counts", counts_cutoff=10, trim_prop=0) 
  combine.counts.for.each.dataset(cage_counts_ori, salmon_rpkm, dexseq_rpkm, proactiv_counts_ori, cell_lines, out_dir, filename="rpkm", counts_cutoff=5, trim_prop=0) 

}

#--------------------------------------------------------
###For new tissue brain
#--------------------------------------------------------
if(TRUE){
  cell_lines <- c("BRAIN")
  out_dir <- "/mnt/citadelb/publication/snakealtpromoter/Result/Brain"
  cage_counts_ori <- clean_mat(readRDS("/mnt/citadelb/publication/snakealtpromoter/Processed/brain_AltPromoter/CAGE_test/cage/counts_merged/cage_counts.rds"))  
  salmon_counts_ori <- clean_mat(readRDS("/mnt/citadelb/publication/snakealtpromoter/Processed/brain_AltPromoter/RNAseq_test/salmon/counts_merged/merged_promoter_counts.rds"))  
  dexseq_counts_ori <- clean_mat(readRDS("/mnt/citadelb/publication/snakealtpromoter/Processed/brain_AltPromoter/RNAseq_test/dexseq/counts_merged/merged_promoter_counts.rds"))  
  proactiv_counts_ori <- clean_mat(readRDS("/mnt/citadelb/publication/snakealtpromoter/Processed/brain_AltPromoter/RNAseq_test/proactiv/counts_merged/proactiv_raw_counts.rds"))
  # cage_counts_ori <- clean_mat(readRDS("/mnt/citadel2/research/shared/SnakeAltPromoter_paper_revision/processed_brain/cage/counts_merged/cage_counts.rds"))  
  # salmon_counts_ori <- clean_mat(readRDS("/mnt/citadel2/research/shared/SnakeAltPromoter_paper_revision/processed_brain/salmon/counts_merged/merged_promoter_counts.rds"))  
  # dexseq_counts_ori <- clean_mat(readRDS("/mnt/citadel2/research/shared/SnakeAltPromoter_paper_revision/processed_brain/dexseq/counts_merged/merged_promoter_counts.rds"))  
  # proactiv_counts_ori <- clean_mat(readRDS("/mnt/citadel2/research/shared/SnakeAltPromoter_paper_revision/processed_brain/proactiv/counts_merged/proactiv_raw_counts.rds"))
  colnames(cage_counts_ori) <- c("CAGE_BRAIN_FEMALE_1","CAGE_BRAIN_FEMALE_2","CAGE_BRAIN_MALE_1","CAGE_BRAIN_MALE_2")
  colnames(salmon_counts_ori) <- c("RNA_BRAIN_FEMALE_1","RNA_BRAIN_FEMALE_2","RNA_BRAIN_MALE_1","RNA_BRAIN_MALE_2")
  colnames(dexseq_counts_ori) <- c("RNA_BRAIN_FEMALE_1","RNA_BRAIN_FEMALE_2","RNA_BRAIN_MALE_1","RNA_BRAIN_MALE_2")
  colnames(proactiv_counts_ori) <- c("RNA_BRAIN_FEMALE_1","RNA_BRAIN_FEMALE_2","RNA_BRAIN_MALE_1","RNA_BRAIN_MALE_2")
  # cage_counts_ori = cage_counts_ori[,c(3),drop=FALSE]
  # salmon_counts_ori = salmon_counts_ori[,c(3),drop=FALSE]
  # dexseq_counts_ori = dexseq_counts_ori[,c(3),drop=FALSE]
  # proactiv_counts_ori = proactiv_counts_ori[,c(3),drop=FALSE]
  salmon_rpkm <- to_rpkm_by_len(salmon_counts_ori, salmon_len_bp)
  dexseq_rpkm <- to_rpkm_by_len(dexseq_counts_ori, dex_len_bp)
  combine.counts.for.each.dataset(cage_counts_ori, salmon_counts_ori, dexseq_counts_ori, proactiv_counts_ori, cell_lines, out_dir, filename="counts", counts_cutoff=10, trim_prop=0) 
  combine.counts.for.each.dataset(cage_counts_ori, salmon_rpkm, dexseq_rpkm, proactiv_counts_ori, cell_lines, out_dir, filename="rpkm", counts_cutoff=5, trim_prop=0) 
}

#--------------------------------------------------------
###For new tissue testis
#--------------------------------------------------------
if(FALSE){
  cell_lines <- c("TESTIS")
  out_dir <- "/mnt/citadelb/publication/snakealtpromoter/Result/Testis"
  cage_counts_ori <- clean_mat(readRDS("/mnt/citadel2/research/shared/SnakeAltPromoter_paper_revision/processed_cage/cage/counts_merged/cage_counts.rds"))  
  salmon_counts_ori <- clean_mat(readRDS("/mnt/citadel2/research/shared/SnakeAltPromoter_paper_revision/processed/salmon/counts_merged/merged_promoter_counts.rds"))  
  dexseq_counts_ori <- clean_mat(readRDS("/mnt/citadel2/research/shared/SnakeAltPromoter_paper_revision/processed/dexseq/counts_merged/merged_promoter_counts.rds"))  
  proactiv_counts_ori <- clean_mat(readRDS("/mnt/citadel2/research/shared/SnakeAltPromoter_paper_revision/processed/proactiv/counts_merged/proactiv_raw_counts.rds"))
  colnames(cage_counts_ori) <- c("CAGE_TESTIS_1","CAGE_TESTIS_2")
  colnames(salmon_counts_ori) <- c("RNA_TESTIS_1","RNA_TESTIS_2")
  colnames(dexseq_counts_ori) <- c("RNA_TESTIS_1","RNA_TESTIS_2")
  colnames(proactiv_counts_ori) <- c("RNA_TESTIS_1","RNA_TESTIS_2")
  salmon_rpkm <- to_rpkm_by_len(salmon_counts_ori, salmon_len_bp)
  dexseq_rpkm <- to_rpkm_by_len(dexseq_counts_ori, dex_len_bp)
  combine.counts.for.each.dataset(cage_counts_ori, salmon_counts_ori, dexseq_counts_ori, proactiv_counts_ori, cell_lines, out_dir, filename="counts", counts_cutoff=10, trim_prop=0) 
  combine.counts.for.each.dataset(cage_counts_ori, salmon_rpkm, dexseq_rpkm, proactiv_counts_ori, cell_lines, out_dir, filename="rpkm", counts_cutoff=5, trim_prop=0) 
}


#--------------------------------------------------------
###For cell lines
#--------------------------------------------------------
if(TRUE){
  cell_lines <- c("K562", "GM12878")
  out_dir <- "/mnt/citadelb/publication/snakealtpromoter/Result/K562_vs_GM12878"
  # cage_counts_ori <- clean_mat(readRDS("/mnt/citadel2/research/shared/AltPromoterFlow/test_for_sample_sheet/cage/counts_merged/cage_counts.rds"))  
  # salmon_counts_ori <- clean_mat(readRDS("/mnt/citadel2/research/shared/AltPromoterFlow/test_for_sample_sheet/salmon/counts_merged/merged_promoter_counts.rds"))  
  # dexseq_counts_ori <- clean_mat(readRDS("/mnt/citadel2/research/shared/AltPromoterFlow/test_for_sample_sheet/dexseq/counts_merged/merged_promoter_counts.rds"))  
  # proactiv_counts_ori <- clean_mat(readRDS("/mnt/citadel2/research/shared/AltPromoterFlow/test_for_sample_sheet/proactiv/counts_merged/proactiv_raw_counts.rds"))
  cage_counts_ori <- clean_mat(readRDS("/mnt/citadelb/publication/snakealtpromoter/Processed/cellline_AltPromoter/CAGE/cage/counts_merged/cage_counts.rds"))  
  salmon_counts_ori <- clean_mat(readRDS("/mnt/citadelb/publication/snakealtpromoter/Processed/cellline_AltPromoter/RNAseq/salmon/counts_merged/merged_promoter_counts.rds"))  
  dexseq_counts_ori <- clean_mat(readRDS("/mnt/citadelb/publication/snakealtpromoter/Processed/cellline_AltPromoter/RNAseq/dexseq/counts_merged/merged_promoter_counts.rds"))  
  proactiv_counts_ori <- clean_mat(readRDS("/mnt/citadelb/publication/snakealtpromoter/Processed/cellline_AltPromoter/RNAseq/proactiv/counts_merged/proactiv_raw_counts.rds"))
  colnames(cage_counts_ori) <- c("CAGE_GM12878_1","CAGE_GM12878_2","CAGE_K562_1","CAGE_K562_2") #"CAGE_H1-hESC_1","CAGE_H1-hESC_2",
  colnames(salmon_counts_ori) <- c("RNA_GM12878_1","RNA_GM12878_2","RNA_GM12878_3","RNA_K562_1","RNA_K562_2","RNA_K562_3") #"RNA_H1-hESC_1","RNA_H1-hESC_2","RNA_H1-hESC_3",
  colnames(dexseq_counts_ori) <- c("RNA_GM12878_1","RNA_GM12878_2","RNA_GM12878_3","RNA_K562_1","RNA_K562_2","RNA_K562_3") #"RNA_H1-hESC_1","RNA_H1-hESC_2","RNA_H1-hESC_3",
  colnames(proactiv_counts_ori) <- c("RNA_GM12878_1","RNA_GM12878_2","RNA_GM12878_3","RNA_K562_1","RNA_K562_2","RNA_K562_3") #"RNA_H1-hESC_1","RNA_H1-hESC_2","RNA_H1-hESC_3",
  salmon_rpkm <- to_rpkm_by_len(salmon_counts_ori, salmon_len_bp)
  dexseq_rpkm <- to_rpkm_by_len(dexseq_counts_ori, dex_len_bp)
  combine.counts.for.each.dataset(cage_counts_ori, salmon_counts_ori, dexseq_counts_ori, proactiv_counts_ori, cell_lines, out_dir, filename="counts", counts_cutoff=10) 
  combine.counts.for.each.dataset(cage_counts_ori, salmon_rpkm, dexseq_rpkm, proactiv_counts_ori, cell_lines, out_dir, filename="rpkm", counts_cutoff=5) 
}








# plot_density_scatter <- function(df, title = "", outfile, use_density = TRUE, n_grid = 200) {
#   if (nrow(df) == 0) {
#     message("No rows for ", outfile); return(invisible(NULL))
#   }

#   # Use log10 to calculate coordinates
#   df$logC <- log10(df$cage_norm + 1)
#   df$logO <- log10(df$other_norm + 1)
#   df$density <- NA_real_

#   # Use log 10 to calculate density
#   if (use_density) {
#     keep <- is.finite(df$logC) & is.finite(df$logO)
#     dens <- with(df[keep, ], MASS::kde2d(logC, logO, n = n_grid))
#     ix <- findInterval(df$logC[keep], dens$x)
#     iy <- findInterval(df$logO[keep], dens$y)
#     df$density[keep] <- log10(dens$z[cbind(ix, iy)] + 1e-8)
#   }


#   # --------------------------------------------------------------
#   #  cor_labels – drop-in replacement
#   # --------------------------------------------------------------
#   cor_labels <- do.call(rbind, lapply(split(df, df$method), function(sub) {
#     x_raw <- sub$cage_norm
#     y_raw <- sub$other_norm
#     ok    <- is.finite(x_raw) & is.finite(y_raw)

#     # ---- safety: need at least 2 points ---------------------------------
#     if (sum(ok) < 2) {
#       return(data.frame(
#         method = unique(sub$method),
#         x = NA_real_, y = NA_real_, label = "N < 2",
#         stringsAsFactors = FALSE
#       ))
#     }

#     x <- x_raw[ok]
#     y <- y_raw[ok]

    
#     # ---- Pearson --------------------------------------------------------
#     pearson_r <- round(cor(x, y), 3)

#     # ---- Spearman -------------------------------------------------------
#     spearman_rho <- round(cor(x, y, method = "spearman"), 3)

#     # ---- Concordance Correlation Coefficient (CCC) ----------------------
#     ccc_val <- tryCatch({
#       CCC(x, y)$rho.c$est
#     }, error = function(e) NA_real_)
#     ccc_val <- round(ccc_val, 3)

#     # ---- Bland-Altman ---------------------------------------------------
#     diff_xy   <- y - x
#     mean_xy   <- (y + x) / 2
#     mean_diff <- mean(diff_xy)
#     bias_str  <- sprintf("%.3f", round(mean_diff, 3))

#     ba_corr <- suppressWarnings(
#       cor(mean_xy, diff_xy, use = "complete.obs")
#     )
#     ba_corr <- round(ifelse(is.na(ba_corr), NA, ba_corr), 3)

#     # ---- Log2 transformed correlations ----------------------------------
#     #pearson_r <- round(cor(log2(x + 1), log2(y + 1), method = "pearson"), 3)
#     #spearman_rho <- round(cor(log2(x + 1), log2(y + 1), method = "spearman"), 3)
#     #ccc_val <- tryCatch({
#     #  CCC(log2(x + 1), log2(y + 1))$rho.c$est
#     #}, error = function(e) NA_real_)

#     # ---- Build the label ------------------------------------------------
#     label <- paste0(
#       "Pearson's r = ", pearson_r, "\n",
#       "Spearman's rho = ", spearman_rho, "\n",
#       "Concordance r = ", ifelse(is.na(ccc_val), "NA", ccc_val), "\n",
#       "Mean bias = ", bias_str, "\n",
#       "BA bias = ", ifelse(is.na(ba_corr), "NA", ba_corr)
#     )

#     # --- Position: 5% from top-left ---
#     xlim <- range(sub$logC[ok], na.rm = TRUE)
#     ylim <- range(c(sub$logC[ok], sub$logO[ok]), na.rm = TRUE)
#     label_x <- xlim[1] + 0.001 * diff(xlim)
#     label_y <- ylim[2] - 0.001 * diff(ylim)

#     data.frame(
#       method = unique(sub$method),
#       x = label_x,
#       y = label_y,
#       label = label,
#       pearson_r = pearson_r,
#       spearman_rho = spearman_rho,
#       ccc_val = ccc_val,
#       bias_str = bias_str,
#       ba_corr = ba_corr,
#       stringsAsFactors = FALSE
#     )
#   }))

#   # Coordinate range
#   lims <- range(c(df$logC, df$logO), finite = TRUE)

#   # Plot
#   p <- ggplot(df, aes(x = logC, y = logO)) +
#     geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey50") +
#     scale_x_continuous(limits = lims) +
#     scale_y_continuous(limits = lims) +
#     facet_wrap(~method, nrow = 1) +
#     labs(x = "CAGE (log10 normalized)", y = "Other method (log10 normalized)") +
#     theme_bw(base_size = 11) +
#     theme(plot.title = element_blank())

#   if (use_density) {
#     p <- p + geom_point(aes(color = density), size = 1, alpha = 0.85) +
#       scale_color_viridis_c(option = "A", name = "Density")
#   } else {
#     p <- p + geom_point(alpha = 0.5, size = 0.7)
#   }

#   p <- p + geom_label(
#     data = cor_labels,
#     aes(x = x, y = y, label = label),
#     hjust = 0,           # left-aligned
#     vjust = 1,           # top-aligned
#     size = 3.2,
#     fill = "white",      # ← white background
#     color = "black",
#     label.padding = unit(0.1, "lines"),  # soft padding
#     label.r = unit(0.15, "lines"),       # slight rounded corners
#     fontface = "plain",
#     alpha = 1            # fully opaque
#   )
#   ggsave(file.path(out_dir, outfile), p, width = 5, height = 4)
#   # Return the correlations
#   cor_results <- cor_labels[, c("method", "pearson_r", "spearman_rho", "ccc_val", "bias_str", "ba_corr")]
#   return(cor_results)
# }








# #-------------
# # Sanity Check
# #-------------

# library(dplyr)

# methods <- list(Salmon = salmon, DEXSeq = dexseq, proActiv = proactiv)

# summary_table <- data.frame(
#   Method     = character(),
#   Direction  = character(),
#   N_Pairs    = integer(),
#   Min        = numeric(),
#   Q1         = numeric(),
#   Median     = numeric(),
#   Mean       = numeric(),
#   Q3         = numeric(),
#   Max        = numeric(),
#   stringsAsFactors = FALSE
# )

# for (method_name in names(methods)) {
#   df_method <- methods[[method_name]]
#   df_common <- inner_join(
#     cage[, c("promoterId", "logFC")],
#     df_method[, c("promoterId", "logFC")],
#     by = "promoterId",
#     suffix = c("_cage", "_method")
#   )

#   for (dir in c("up", "down")) {
#     if (dir == "up") {
#       df_sub <- df_common %>% filter(logFC_cage > 0, logFC_method > 0)
#     } else {
#       df_sub <- df_common %>% filter(logFC_cage < 0, logFC_method < 0)
#     }

#     if (nrow(df_sub) == 0) next

#     stat <- summary(df_sub$logFC_method)
#     stat_num <- as.numeric(stat)

#     summary_table <- rbind(summary_table, data.frame(
#     Method    = method_name,
#     Direction = dir,
#     N_Pairs   = nrow(df_sub),
#     Min       = stat_num[1],
#     Q1        = stat_num[2],
#     Median    = stat_num[3],
#     Mean      = stat_num[4],
#     Q3        = stat_num[5],
#     Max       = stat_num[6]
#     ))

#   }
# }

# print(summary_table)


#     Method Direction N_Pairs         Min        Q1    Median      Mean
# 1   Salmon        up    1078   0.3231317  1.351549  2.664919  4.035902
# 2   Salmon      down    1651 -23.9154509 -8.489788 -5.099268 -5.616285
# 3   DEXSeq        up     583   0.6415377  2.054174  3.766849  4.314652
# 4   DEXSeq      down    1042 -16.4714099 -7.566906 -5.807926 -5.664206
# 5 proActiv        up     640   0.6317126  1.949220  3.943887  4.375141
# 6 proActiv      down    1066 -13.5291621 -7.126131 -5.588291 -5.375614
#          Q3        Max
# 1  6.169641 24.4513245
# 2 -2.008386 -0.3412513
# 3  6.020754 13.4673443
# 4 -3.012444 -0.6648830
# 5  6.289326 13.2363037
# 6 -3.211264 -0.7706925

