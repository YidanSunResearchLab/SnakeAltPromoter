
prewd <- "/mnt/citadel2/research/shared/AltPromoterFlow/data/test_data"
ln_dir_CAGE <- "/mnt/citadel2/research/shared/AltPromoterFlow/data/CAGE"
ln_dir <- "/mnt/citadel2/research/shared/AltPromoterFlow/data/CHIP"

sample.info.merge <- read.delim(file.path(prewd, "SraRunTable.tsv"), header = TRUE, stringsAsFactors = FALSE, sep = ",")

sample.info.merge <- subset(sample.info.merge, Organism == "Homo sapiens")

sample.info.merge$strain <- gsub("-", "", sample.info.merge$cell_line)

cols_to_clean <- c("strain", "tissue", "treatment", "LibrarySelection", "Assay.Type", "antibody")
sample.info.merge[cols_to_clean] <- lapply(sample.info.merge[cols_to_clean], function(x) {
  x[is.na(x)] <- ""
  return(x)
})

sample.info.merge$final_name <- apply(sample.info.merge, 1, function(row) {
  parts <- c(
    "Human",
    row["strain"],
    row["tissue"],
    row["treatment"],
    row["LibrarySelection"],
    row["Assay.Type"],
    row["antibody"],
    row["Run"]
  )
  paste(Filter(function(x) x != "", parts), collapse = "_")
})

sample.info.merge$final_name <- gsub(" ", "_", sample.info.merge$final_name)

for (i in 1:nrow(sample.info.merge)) {
  run_id <- sample.info.merge$Run[i]
  final_name <- sample.info.merge$final_name[i]
  library <- sample.info.merge$`LibrarySelection`[i]
  
  target_dir <- if (grepl("CAGE", library, ignore.case = TRUE)) ln_dir_CAGE else ln_dir

  from_R1 <- file.path(prewd, paste0(run_id, "_1.fastq.gz"))
  to_R1 <- file.path(target_dir, paste0(final_name, "_R1.fastq.gz"))
  if (!file.exists(to_R1)) file.symlink(from_R1, to_R1)

  from_R2 <- file.path(prewd, paste0(run_id, "_2.fastq.gz"))
  to_R2 <- file.path(target_dir, paste0(final_name, "_R2.fastq.gz"))
  if (file.exists(from_R2) && !file.exists(to_R2)) file.symlink(from_R2, to_R2)
}
