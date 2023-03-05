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
library(usefunc) # own package of useful functions -- for is.binary()

## args
args <- commandArgs(trailingOnly = TRUE)
severe_file <- args[1]
severe_def <- args[2]
linker_file <- args[3]
exclusions_file <- args[4]
exc <- args[5]
withdrawl_file <- args[6]
outfile <- args[7]
plink <- args[8]

## manual args
# severe_file <- "../phenotype-extraction/data/case-pheno.tsv"
# severe_def <- "d9"
# linker_file <- "data/linker.csv"
# exclusions_file <- "data/ukb-exclusions.tsv"
# exc <- "sex_mismatch sex_chr_aneuploidy het_miss_outliers" ## ADD IN EXCLUSIONS FROM LIST 
# withdrawl_file <- "data/w15147_20210809.csv"
# outfile <- "data/mod_severe-ad-phenofile.txt"
# plink <- "TRUE"

plink <- as.logical(plink)

## data
severe_dat <- read_tsv(severe_file)
linker <- read_csv(linker_file)
exc <- unlist(str_split(exc, " "))
withdrawals <- readLines(withdrawl_file)
exclusions <- read_tsv(exclusions_file) %>%
	dplyr::filter(exc_reason %in% exc)

# ------------------------------------------------------------------------------
# remove exclusions and withdrawals
# ------------------------------------------------------------------------------

severe_no_exc <- severe_dat %>%
	dplyr::filter(!f.eid %in% withdrawals) %>%
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

## Change binary variables to being a 1/2 rather than a 0/1 because plink is stupid
if (plink) {
	bin_vars <- sapply(severe_only, is.binary)
	bin_vars <- names(bin_vars)[bin_vars]
	for (bin in bin_vars) {
		severe_only[[bin]] <- replace(severe_only[[bin]], severe_only[[bin]] == 1, 2)
		severe_only[[bin]] <- replace(severe_only[[bin]], severe_only[[bin]] == 0, 1)
	}
}

## write out the data
write.table(severe_only, file = outfile, col.names = T, row.names = F, quote = F, sep = " ")









