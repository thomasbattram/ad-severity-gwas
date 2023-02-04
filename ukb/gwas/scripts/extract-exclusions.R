# ----------------------------------------------------------------------
# Extract IDs to exclude
# ----------------------------------------------------------------------

## Extract individuals to exclude for phenotype data

## pkgs
library(tidyverse)

## args
args <- commandArgs(trailingOnly = TRUE)
ukb_gendir <- args[1]
outfile <- args[2]

## manual args
# ukb_gendir <- "/mnt/storage/private/mrcieu/data/ukbiobank/genetic/variants/arrays/imputed/released/2018-09-18/data/derived"
# outfile <- "data/ukb-exclusions.tsv"

## exclusion files
files <- c(non_white_british = "ancestry/data.non_white_british.qctools.txt", 
		   non_europeans = "ancestry/data.non_europeans.qctools.txt", 
		   minimal_relateds = "related/relateds_exclusions/data.minimal_relateds.qctools.txt", 
		   highly_relateds = "related/relateds_exclusions/data.highly_relateds.qctools.txt",
		   sex_mismatch = "standard_exclusions/data.sex_mismatch.qctools.txt",
		   sex_chr_aneuploidy = "standard_exclusions/data.putative_sex_chromosome_aneuploidy.qctools.txt", 
		   het_miss_outliers = "standard_exclusions/data.het_missing_outliers.qctools.txt", 
		   combined = "standard_exclusions/data.combined_recommended.qctools.txt")

full_files <- file.path(ukb_gendir, files)
names(full_files) <- names(files)

ids <- map_dfr(1:length(full_files), function(x) {
	ids <- readLines(full_files[x])
	out <- tibble(exc_reason = names(full_files)[x], 
				  id = ids)
	return(out)
})

write.table(ids, file = outfile, col.names = T, row.names = F, quote = F, sep = "\t")