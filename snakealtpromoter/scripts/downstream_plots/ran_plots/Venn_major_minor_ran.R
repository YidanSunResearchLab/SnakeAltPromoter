#activate environment: /mnt/citadel2/research/shared/AltPromoterFlow/genome/.snakemake_conda/3ea3da478f24d1f8a438361e1dd371ca_

#Run at the end
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

#Load files
CAGE_minor <- readRDS("/mnt/citadel2/research/shared/SnakeAltPromoter_0228/normal6/cage/promoter_classification_total/Minor_promoterId_overall.rds")
CAGE_major <- readRDS("/mnt/citadel2/research/shared/SnakeAltPromoter_0228/normal6/cage/promoter_classification_total/Major_promoterId_overall.rds")
salmon_minor <- readRDS("/mnt/citadel2/research/shared/SnakeAltPromoter_0228/normal6/salmon/promoter_classification_total/Minor_promoterId_overall.rds")
salmon_major <- readRDS("/mnt/citadel2/research/shared/SnakeAltPromoter_0228/normal6/salmon/promoter_classification_total/Major_promoterId_overall.rds")
dexseq_minor <- readRDS("/mnt/citadel2/research/shared/SnakeAltPromoter_0228/normal6/dexseq/promoter_classification_total/Minor_promoterId_overall.rds")
dexseq_major <- readRDS("/mnt/citadel2/research/shared/SnakeAltPromoter_0228/normal6/dexseq/promoter_classification_total/Major_promoterId_overall.rds")
proactiv_major <- readRDS("/mnt/citadel2/research/shared/SnakeAltPromoter_0228/normal6/proactiv/promoter_classification_total/Major_promoterId_overall.rds")
proactiv_minor <- readRDS("/mnt/citadel2/research/shared/SnakeAltPromoter_0228/normal6/proactiv/promoter_classification_total/Minor_promoterId_overall.rds")
out_dir <- "/mnt/citadel2/research/shared/SnakeAltPromoter_0228/downstream_fixed/venn_major_minor"
build.overlap.function(CAGE_minor, CAGE_major,salmon_minor, salmon_major,dexseq_minor, dexseq_major,proactiv_minor, proactiv_major,out_dir)

# -------- after normal6 load --------
normal_CAGE_minor <- CAGE_minor
normal_CAGE_major <- CAGE_major
normal_salmon_minor <- salmon_minor
normal_salmon_major <- salmon_major
normal_dexseq_minor <- dexseq_minor
normal_dexseq_major <- dexseq_major
normal_proactiv_minor <- proactiv_minor
normal_proactiv_major <- proactiv_major

#Load files
CAGE_minor <- readRDS("/mnt/citadel2/research/shared/SnakeAltPromoter_0228/normal6/cage/promoter_classification_total/Minor_promoterId_overall.rds")
CAGE_major <- readRDS("/mnt/citadel2/research/shared/SnakeAltPromoter_0228/normal6/cage/promoter_classification_total/Major_promoterId_overall.rds")
salmon_minor <- readRDS("/mnt/citadel2/research/shared/SnakeAltPromoter_0228/threshold0.5/salmon/promoter_classification_total/Minor_promoterId_overall.rds")
salmon_major <- readRDS("/mnt/citadel2/research/shared/SnakeAltPromoter_0228/threshold0.5/salmon/promoter_classification_total/Major_promoterId_overall.rds")
dexseq_minor <- readRDS("/mnt/citadel2/research/shared/SnakeAltPromoter_0228/threshold0.5/dexseq/promoter_classification_total/Minor_promoterId_overall.rds")
dexseq_major <- readRDS("/mnt/citadel2/research/shared/SnakeAltPromoter_0228/threshold0.5/dexseq/promoter_classification_total/Major_promoterId_overall.rds")
proactiv_major <- readRDS("/mnt/citadel2/research/shared/SnakeAltPromoter_0228/threshold0.5/proactiv/promoter_classification_total/Major_promoterId_overall.rds")
proactiv_minor <- readRDS("/mnt/citadel2/research/shared/SnakeAltPromoter_0228/threshold0.5/proactiv/promoter_classification_total/Minor_promoterId_overall.rds")
out_dir <- "/mnt/citadel2/research/shared/SnakeAltPromoter_0228/downstream_fixed/venn_major_minor_0.5"
build.overlap.function(CAGE_minor, CAGE_major,salmon_minor, salmon_major,dexseq_minor, dexseq_major,proactiv_minor, proactiv_major,out_dir)


#Load files
CAGE_minor <- readRDS("/mnt/citadel2/research/shared/SnakeAltPromoter_0228/normal6/cage/promoter_classification_total/Minor_promoterId_overall.rds")
CAGE_major <- readRDS("/mnt/citadel2/research/shared/SnakeAltPromoter_0228/normal6/cage/promoter_classification_total/Major_promoterId_overall.rds")
salmon_minor <- readRDS("/mnt/citadel2/research/shared/SnakeAltPromoter_0228/threshold0.1/salmon/promoter_classification_total/Minor_promoterId_overall.rds")
salmon_major <- readRDS("/mnt/citadel2/research/shared/SnakeAltPromoter_0228/threshold0.1/salmon/promoter_classification_total/Major_promoterId_overall.rds")
dexseq_minor <- readRDS("/mnt/citadel2/research/shared/SnakeAltPromoter_0228/threshold0.1/dexseq/promoter_classification_total/Minor_promoterId_overall.rds")
dexseq_major <- readRDS("/mnt/citadel2/research/shared/SnakeAltPromoter_0228/threshold0.1/dexseq/promoter_classification_total/Major_promoterId_overall.rds")
proactiv_major <- readRDS("/mnt/citadel2/research/shared/SnakeAltPromoter_0228/threshold0.1/proactiv/promoter_classification_total/Major_promoterId_overall.rds")
proactiv_minor <- readRDS("/mnt/citadel2/research/shared/SnakeAltPromoter_0228/threshold0.1/proactiv/promoter_classification_total/Minor_promoterId_overall.rds")
out_dir <- "/mnt/citadel2/research/shared/SnakeAltPromoter_0228/downstream_fixed/venn_major_minor_0.1"
build.overlap.function(CAGE_minor, CAGE_major,salmon_minor, salmon_major,dexseq_minor, dexseq_major,proactiv_minor, proactiv_major,out_dir)


#Load files
CAGE_minor <- readRDS("/mnt/citadel2/research/shared/SnakeAltPromoter_0228/normal6/cage/promoter_classification_total/Minor_promoterId_overall.rds")
CAGE_major <- readRDS("/mnt/citadel2/research/shared/SnakeAltPromoter_0228/normal6/cage/promoter_classification_total/Major_promoterId_overall.rds")
salmon_minor <- readRDS("/mnt/citadel2/research/shared/SnakeAltPromoter_0228/downsample/10M/salmon/promoter_classification_total/Minor_promoterId_overall.rds")
salmon_major <- readRDS("/mnt/citadel2/research/shared/SnakeAltPromoter_0228/downsample/10M/salmon/promoter_classification_total/Major_promoterId_overall.rds")
dexseq_minor <- readRDS("/mnt/citadel2/research/shared/SnakeAltPromoter_0228/downsample/10M/dexseq/promoter_classification_total/Minor_promoterId_overall.rds")
dexseq_major <- readRDS("/mnt/citadel2/research/shared/SnakeAltPromoter_0228/downsample/10M/dexseq/promoter_classification_total/Major_promoterId_overall.rds")
proactiv_major <- readRDS("/mnt/citadel2/research/shared/SnakeAltPromoter_0228/downsample/10M/proactiv/promoter_classification_total/Major_promoterId_overall.rds")
proactiv_minor <- readRDS("/mnt/citadel2/research/shared/SnakeAltPromoter_0228/downsample/10M/proactiv/promoter_classification_total/Minor_promoterId_overall.rds")
out_dir <- "/mnt/citadel2/research/shared/SnakeAltPromoter_0228/downstream_fixed/10M_major_minor"
build.overlap.function(CAGE_minor, CAGE_major,salmon_minor, salmon_major,dexseq_minor, dexseq_major,proactiv_minor, proactiv_major,out_dir)


CAGE_minor <- readRDS("/mnt/citadel2/research/shared/SnakeAltPromoter_0228/normal6/cage/promoter_classification_total/Minor_promoterId_overall.rds")
CAGE_major <- readRDS("/mnt/citadel2/research/shared/SnakeAltPromoter_0228/normal6/cage/promoter_classification_total/Major_promoterId_overall.rds")
salmon_minor <- readRDS("/mnt/citadel2/research/shared/SnakeAltPromoter_0228/downsample/20M/salmon/promoter_classification_total/Minor_promoterId_overall.rds")
salmon_major <- readRDS("/mnt/citadel2/research/shared/SnakeAltPromoter_0228/downsample/20M/salmon/promoter_classification_total/Major_promoterId_overall.rds")
dexseq_minor <- readRDS("/mnt/citadel2/research/shared/SnakeAltPromoter_0228/downsample/20M/dexseq/promoter_classification_total/Minor_promoterId_overall.rds")
dexseq_major <- readRDS("/mnt/citadel2/research/shared/SnakeAltPromoter_0228/downsample/20M/dexseq/promoter_classification_total/Major_promoterId_overall.rds")
proactiv_major <- readRDS("/mnt/citadel2/research/shared/SnakeAltPromoter_0228/downsample/20M/proactiv/promoter_classification_total/Major_promoterId_overall.rds")
proactiv_minor <- readRDS("/mnt/citadel2/research/shared/SnakeAltPromoter_0228/downsample/20M/proactiv/promoter_classification_total/Minor_promoterId_overall.rds")
out_dir <- "/mnt/citadel2/research/shared/SnakeAltPromoter_0228/downstream_fixed/20M_major_minor"
build.overlap.function(CAGE_minor, CAGE_major,salmon_minor, salmon_major,dexseq_minor, dexseq_major,proactiv_minor, proactiv_major,out_dir)

CAGE_minor <- readRDS("/mnt/citadel2/research/shared/SnakeAltPromoter_0228/normal6/cage/promoter_classification_total/Minor_promoterId_overall.rds")
CAGE_major <- readRDS("/mnt/citadel2/research/shared/SnakeAltPromoter_0228/normal6/cage/promoter_classification_total/Major_promoterId_overall.rds")
salmon_minor <- readRDS("/mnt/citadel2/research/shared/SnakeAltPromoter_0228/downsample/30M/salmon/promoter_classification_total/Minor_promoterId_overall.rds")
salmon_major <- readRDS("/mnt/citadel2/research/shared/SnakeAltPromoter_0228/downsample/30M/salmon/promoter_classification_total/Major_promoterId_overall.rds")
dexseq_minor <- readRDS("/mnt/citadel2/research/shared/SnakeAltPromoter_0228/downsample/30M/dexseq/promoter_classification_total/Minor_promoterId_overall.rds")
dexseq_major <- readRDS("/mnt/citadel2/research/shared/SnakeAltPromoter_0228/downsample/30M/dexseq/promoter_classification_total/Major_promoterId_overall.rds")
proactiv_major <- readRDS("/mnt/citadel2/research/shared/SnakeAltPromoter_0228/downsample/30M/proactiv/promoter_classification_total/Major_promoterId_overall.rds")
proactiv_minor <- readRDS("/mnt/citadel2/research/shared/SnakeAltPromoter_0228/downsample/30M/proactiv/promoter_classification_total/Minor_promoterId_overall.rds")
out_dir <- "/mnt/citadel2/research/shared/SnakeAltPromoter_0228/downstream_fixed/30M_major_minor"
build.overlap.function(CAGE_minor, CAGE_major,salmon_minor, salmon_major,dexseq_minor, dexseq_major,proactiv_minor, proactiv_major,out_dir)

CAGE_minor <- readRDS("/mnt/citadel2/research/shared/SnakeAltPromoter_0228/normal6/cage/promoter_classification_total/Minor_promoterId_overall.rds")
CAGE_major <- readRDS("/mnt/citadel2/research/shared/SnakeAltPromoter_0228/normal6/cage/promoter_classification_total/Major_promoterId_overall.rds")
salmon_minor <- readRDS("/mnt/citadel2/research/shared/SnakeAltPromoter_0228/downsample/40M/salmon/promoter_classification_total/Minor_promoterId_overall.rds")
salmon_major <- readRDS("/mnt/citadel2/research/shared/SnakeAltPromoter_0228/downsample/40M/salmon/promoter_classification_total/Major_promoterId_overall.rds")
dexseq_minor <- readRDS("/mnt/citadel2/research/shared/SnakeAltPromoter_0228/downsample/40M/dexseq/promoter_classification_total/Minor_promoterId_overall.rds")
dexseq_major <- readRDS("/mnt/citadel2/research/shared/SnakeAltPromoter_0228/downsample/40M/dexseq/promoter_classification_total/Major_promoterId_overall.rds")
proactiv_major <- readRDS("/mnt/citadel2/research/shared/SnakeAltPromoter_0228/downsample/40M/proactiv/promoter_classification_total/Major_promoterId_overall.rds")
proactiv_minor <- readRDS("/mnt/citadel2/research/shared/SnakeAltPromoter_0228/downsample/40M/proactiv/promoter_classification_total/Minor_promoterId_overall.rds")
out_dir <- "/mnt/citadel2/research/shared/SnakeAltPromoter_0228/downstream_fixed/40M_major_minor"
build.overlap.function(CAGE_minor, CAGE_major,salmon_minor, salmon_major,dexseq_minor, dexseq_major,proactiv_minor, proactiv_major,out_dir)

CAGE_minor <- readRDS("/mnt/citadel2/research/shared/SnakeAltPromoter_0228/normal6/cage/promoter_classification_total/Minor_promoterId_overall.rds")
CAGE_major <- readRDS("/mnt/citadel2/research/shared/SnakeAltPromoter_0228/normal6/cage/promoter_classification_total/Major_promoterId_overall.rds")
salmon_minor <- readRDS("/mnt/citadel2/research/shared/SnakeAltPromoter_0228/downsample/80M/salmon/promoter_classification_total/Minor_promoterId_overall.rds")
salmon_major <- readRDS("/mnt/citadel2/research/shared/SnakeAltPromoter_0228/downsample/80M/salmon/promoter_classification_total/Major_promoterId_overall.rds")
dexseq_minor <- readRDS("/mnt/citadel2/research/shared/SnakeAltPromoter_0228/downsample/80M/dexseq/promoter_classification_total/Minor_promoterId_overall.rds")
dexseq_major <- readRDS("/mnt/citadel2/research/shared/SnakeAltPromoter_0228/downsample/80M/dexseq/promoter_classification_total/Major_promoterId_overall.rds")
proactiv_major <- readRDS("/mnt/citadel2/research/shared/SnakeAltPromoter_0228/downsample/80M/proactiv/promoter_classification_total/Major_promoterId_overall.rds")
proactiv_minor <- readRDS("/mnt/citadel2/research/shared/SnakeAltPromoter_0228/downsample/80M/proactiv/promoter_classification_total/Minor_promoterId_overall.rds")
out_dir <- "/mnt/citadel2/research/shared/SnakeAltPromoter_0228/downstream_fixed/80M_major_minor"
build.overlap.function(CAGE_minor, CAGE_major,salmon_minor, salmon_major,dexseq_minor, dexseq_major,proactiv_minor, proactiv_major,out_dir)


#Load files
CAGE_minor <- readRDS("/mnt/citadel2/research/shared/SnakeAltPromoter_0228/normal6/cage/promoter_classification_total/Minor_promoterId_overall.rds")
CAGE_major <- readRDS("/mnt/citadel2/research/shared/SnakeAltPromoter_0228/normal6/cage/promoter_classification_total/Major_promoterId_overall.rds")
salmon_minor <- readRDS("/mnt/citadel2/research/shared/SnakeAltPromoter_0228/readlensim/50/salmon/promoter_classification_total/Minor_promoterId_overall.rds")
salmon_major <- readRDS("/mnt/citadel2/research/shared/SnakeAltPromoter_0228/readlensim/50/salmon/promoter_classification_total/Major_promoterId_overall.rds")
dexseq_minor <- readRDS("/mnt/citadel2/research/shared/SnakeAltPromoter_0228/readlensim/50/dexseq/promoter_classification_total/Minor_promoterId_overall.rds")
dexseq_major <- readRDS("/mnt/citadel2/research/shared/SnakeAltPromoter_0228/readlensim/50/dexseq/promoter_classification_total/Major_promoterId_overall.rds")
proactiv_major <- readRDS("/mnt/citadel2/research/shared/SnakeAltPromoter_0228/readlensim/50/proactiv/promoter_classification_total/Major_promoterId_overall.rds")
proactiv_minor <- readRDS("/mnt/citadel2/research/shared/SnakeAltPromoter_0228/readlensim/50/proactiv/promoter_classification_total/Minor_promoterId_overall.rds")
out_dir <- "/mnt/citadel2/research/shared/SnakeAltPromoter_0228/downstream_fixed/50_major_minor"
build.overlap.function(CAGE_minor, CAGE_major,salmon_minor, salmon_major,dexseq_minor, dexseq_major,proactiv_minor, proactiv_major,out_dir)

CAGE_minor <- readRDS("/mnt/citadel2/research/shared/SnakeAltPromoter_0228/normal6/cage/promoter_classification_total/Minor_promoterId_overall.rds")
CAGE_major <- readRDS("/mnt/citadel2/research/shared/SnakeAltPromoter_0228/normal6/cage/promoter_classification_total/Major_promoterId_overall.rds")
salmon_minor <- readRDS("/mnt/citadel2/research/shared/SnakeAltPromoter_0228/readlensim/100/salmon/promoter_classification_total/Minor_promoterId_overall.rds")
salmon_major <- readRDS("/mnt/citadel2/research/shared/SnakeAltPromoter_0228/readlensim/100/salmon/promoter_classification_total/Major_promoterId_overall.rds")
dexseq_minor <- readRDS("/mnt/citadel2/research/shared/SnakeAltPromoter_0228/readlensim/100/dexseq/promoter_classification_total/Minor_promoterId_overall.rds")
dexseq_major <- readRDS("/mnt/citadel2/research/shared/SnakeAltPromoter_0228/readlensim/100/dexseq/promoter_classification_total/Major_promoterId_overall.rds")
proactiv_major <- readRDS("/mnt/citadel2/research/shared/SnakeAltPromoter_0228/readlensim/100/proactiv/promoter_classification_total/Major_promoterId_overall.rds")
proactiv_minor <- readRDS("/mnt/citadel2/research/shared/SnakeAltPromoter_0228/readlensim/100/proactiv/promoter_classification_total/Minor_promoterId_overall.rds")
out_dir <- "/mnt/citadel2/research/shared/SnakeAltPromoter_0228/downstream_fixed/100_major_minor"
build.overlap.function(CAGE_minor, CAGE_major,salmon_minor, salmon_major,dexseq_minor, dexseq_major,proactiv_minor, proactiv_major,out_dir)

CAGE_minor <- readRDS("/mnt/citadel2/research/shared/SnakeAltPromoter_0228/normal6/cage/promoter_classification_total/Minor_promoterId_overall.rds")
CAGE_major <- readRDS("/mnt/citadel2/research/shared/SnakeAltPromoter_0228/normal6/cage/promoter_classification_total/Major_promoterId_overall.rds")
salmon_minor <- readRDS("/mnt/citadel2/research/shared/SnakeAltPromoter_0228/readlensim/150/salmon/promoter_classification_total/Minor_promoterId_overall.rds")
salmon_major <- readRDS("/mnt/citadel2/research/shared/SnakeAltPromoter_0228/readlensim/150/salmon/promoter_classification_total/Major_promoterId_overall.rds")
dexseq_minor <- readRDS("/mnt/citadel2/research/shared/SnakeAltPromoter_0228/readlensim/150/dexseq/promoter_classification_total/Minor_promoterId_overall.rds")
dexseq_major <- readRDS("/mnt/citadel2/research/shared/SnakeAltPromoter_0228/readlensim/150/dexseq/promoter_classification_total/Major_promoterId_overall.rds")
proactiv_major <- readRDS("/mnt/citadel2/research/shared/SnakeAltPromoter_0228/readlensim/150/proactiv/promoter_classification_total/Major_promoterId_overall.rds")
proactiv_minor <- readRDS("/mnt/citadel2/research/shared/SnakeAltPromoter_0228/readlensim/150/proactiv/promoter_classification_total/Minor_promoterId_overall.rds")
out_dir <- "/mnt/citadel2/research/shared/SnakeAltPromoter_0228/downstream_fixed/150_major_minor"
build.overlap.function(CAGE_minor, CAGE_major,salmon_minor, salmon_major,dexseq_minor, dexseq_major,proactiv_minor, proactiv_major,out_dir)

CAGE_minor <- readRDS("/mnt/citadel2/research/shared/SnakeAltPromoter_0228/normal6/cage/promoter_classification_total/Minor_promoterId_overall.rds")
CAGE_major <- readRDS("/mnt/citadel2/research/shared/SnakeAltPromoter_0228/normal6/cage/promoter_classification_total/Major_promoterId_overall.rds")
salmon_minor <- readRDS("/mnt/citadel2/research/shared/SnakeAltPromoter_0228/readlensim/75/salmon/promoter_classification_total/Minor_promoterId_overall.rds")
salmon_major <- readRDS("/mnt/citadel2/research/shared/SnakeAltPromoter_0228/readlensim/75/salmon/promoter_classification_total/Major_promoterId_overall.rds")
dexseq_minor <- readRDS("/mnt/citadel2/research/shared/SnakeAltPromoter_0228/readlensim/75/dexseq/promoter_classification_total/Minor_promoterId_overall.rds")
dexseq_major <- readRDS("/mnt/citadel2/research/shared/SnakeAltPromoter_0228/readlensim/75/dexseq/promoter_classification_total/Major_promoterId_overall.rds")
proactiv_major <- readRDS("/mnt/citadel2/research/shared/SnakeAltPromoter_0228/readlensim/75/proactiv/promoter_classification_total/Major_promoterId_overall.rds")
proactiv_minor <- readRDS("/mnt/citadel2/research/shared/SnakeAltPromoter_0228/readlensim/75/proactiv/promoter_classification_total/Minor_promoterId_overall.rds")
out_dir <- "/mnt/citadel2/research/shared/SnakeAltPromoter_0228/downstream_fixed/75_major_minor"
build.overlap.function(CAGE_minor, CAGE_major,salmon_minor, salmon_major,dexseq_minor, dexseq_major,proactiv_minor, proactiv_major,out_dir)
