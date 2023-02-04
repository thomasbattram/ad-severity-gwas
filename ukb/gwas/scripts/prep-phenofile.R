# ------------------------------------------------------------------------------
# Prepare phenotype files for GWAS
# ------------------------------------------------------------------------------

## Prepare phenotype file for running the GWAS using bolt LMM in UKB

## Format for phenofile:
## FID IID severe-AD

## Should be whitespace delimited. The phenofile needs to be numeric - case/control = 1/0
## covariates can be categorical but can't have whitespace in the text

## pkgs
library(tidyverse) # tidy code and data

## args
args <- commandArgs(trailingOnly = TRUE)
severe_file <- args[1]
severe_def <- args[2]
linker_file <- args[3]
exclusions_file <- args[4]
exc <- args[5]
outfile <- args[6]

## manual args
# severe_file <- "../phenotype-extraction/data/case-pheno.tsv"
# severe_def <- "d4"
# linker_file <- "data/linker.csv"
# exclusions_file <- "data/ukb-exclusions.tsv"
# exc <- "sex_mismatch sex_chr_aneuploidy het_miss_outliers" ## ADD IN EXCLUSIONS FROM LIST 
# outfile <- "data/severe-ad-phenofile.txt"

## data
severe_dat <- read_tsv(severe_file)
linker <- read_csv(linker_file)
exc <- unlist(str_split(exc, " "))
exclusions <- read_tsv(exclusions_file) %>%
	dplyr::filter(exc_reason %in% exc)

# ------------------------------------------------------------------------------
# remove exclusions
# ------------------------------------------------------------------------------

severe_no_exc <- severe_dat %>%
	left_join(linker, by = c("f.eid" = "app")) %>%
	dplyr::filter(!ieu %in% exclusions$id)

nrow(severe_dat) - nrow(severe_no_exc)
## 116 people removed 

# ------------------------------------------------------------------------------
# extract severe phenotype
# ------------------------------------------------------------------------------

## check which is moderate severe and which is severe
severe_only <- severe_no_exc %>%
	dplyr::select(FID = ieu, IID = ieu, one_of(severe_def))
colnames(severe_only)[3] <- "severe"

## write out the data
write.table(severe_only, file = outfile, col.names = T, row.names = F, quote = F, sep = " ")









