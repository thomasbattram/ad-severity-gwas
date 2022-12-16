# ----------------------------------------------------------------------
# Extract IDs to exclude
# ----------------------------------------------------------------------

## Extract individuals to exclude for phenotype data

## pkgs
library(tidyverse)


## exclusion files
# ukb_gendir <- "/mnt/storage/private/mrcieu/data/ukbiobank/genetic/variants/arrays/imputed/released/2018-09-18/data/derived"
files <- c("ancestry/data.non_white_british.qctools.txt", 
		   "ancestry/data.non_europeans.qctools.txt", 
		   "related/relateds_exclusions/data.minimal_relateds.qctools.txt", 
		   "related/relateds_exclusions/data.highly_relateds.qctools.txt",
		   "standard_exclusions/data.sex_mismatch.qctools.txt",
		   "standard_exclusions/data.putative_sex_chromosome_aneuploidy.qctools.txt", 
		   "standard_exclusions/data.het_missing_outliers.qctools.txt", 
		   "standard_exclusions/data.combined_recommended.qctools.txt")
linker_file <- "../linker.csv"

full_files <- file.path(ukb_gendir, files)

ids <- lapply(full_files, readLines)
ids <- unique(unlist(ids))
ids <- unique(gsub("\\\"", "", ids))
length(ids)

linker_dat <- read_csv(linker_file)

app_ids <- linker_dat %>%
	dplyr::filter(!ieu %in% ids) %>%
	pull(app) %>%
	as.character()

writeLines(app_ids, "data/geno-include-ids.txt")