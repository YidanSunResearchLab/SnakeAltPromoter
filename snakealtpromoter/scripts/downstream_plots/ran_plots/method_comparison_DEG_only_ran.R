# Run at the start
# Read in promoter-level raw count matrices
cage_counts     <- readRDS("/mnt/citadel2/research/shared/AltPromoterFlow/CAGE_HEART/salmon/merged/promoter_counts.rds")
salmon_counts   <- readRDS("/mnt/citadel2/research/shared/AltPromoterFlow/RNA_HEART/salmon/merged/promoter_counts.rds")
dexseq_counts   <- readRDS("/mnt/citadel2/research/shared/AltPromoterFlow/RNA_HEART/dexseq/merge/promoter_counts.rds")
proactiv_counts <- readRDS("/mnt/citadel2/research/shared/AltPromoterFlow/RNA_HEART/proactiv/quantify/raw_promoter_counts.rds")
out_dir <- "/mnt/citadel2/research/shared/AltPromoterFlow/RNA_HEART/comparison2"
# Read in CAGE feature counts
fc_dir   <- "/mnt/citadel2/research/shared/AltPromoterFlow/CAGE_HEART/featureCounts"
# Read in promoter annotation
anno <- readRDS("/mnt/citadel2/research/shared/AltPromoterFlow/genome/organisms/hg38/Annotation/proActiv_promoter_annotation.rds")
# Read in promoter categories
mm <- readRDS("/mnt/citadel2/research/shared/AltPromoterFlow/RNA_HEART/comparison3/CAGE_promoterId_Multipromoter.Multiactive.rds")
ms <- readRDS("/mnt/citadel2/research/shared/AltPromoterFlow/RNA_HEART/comparison3/CAGE_promoterId_Multipromoter.Singleactive.rds")
ss <- readRDS("/mnt/citadel2/research/shared/AltPromoterFlow/RNA_HEART/comparison3/CAGE_promoterId_Singlepromoter.Singleactive.rds")
# Read in major and minor promoters
major <- readRDS("/mnt/citadel2/research/shared/AltPromoterFlow/RNA_HEART/comparison3/CAGE_major_promoterId.rds")
minor <- readRDS("/mnt/citadel2/research/shared/AltPromoterFlow/RNA_HEART/comparison3/CAGE_minor_promoterId.rds")
inactive <- readRDS("/mnt/citadel2/research/shared/AltPromoterFlow/RNA_HEART/comparison3/CAGE_inactive_promoterId.rds")
intersect <- readRDS("/mnt/citadel2/research/shared/AltPromoterFlow/RNA_HEART/comparison/intersect_major_minor.rds")
