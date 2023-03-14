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
# exc <- "sex_mismatch sex_chr_aneuploidy het_miss_outliers non_europeans" ## ADD IN EXCLUSIONS FROM LIST 
# withdrawl_file <- "data/w15147_20210809.csv"
# outfile <- "data/mod_severe-ad-phenofile.txt"
# plink <- "FALSE"

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

## Write out exclusions list in eids for Jake
# eids_exc <- exclusions %>%
# 	left_join(linker, by = c("id" = "ieu")) %>%
# 	pull(app) %>%
# 	unique()

# eids_exc <- unique(eids_exc, withdrawals)

# writeLines(as.character(eids_exc), "data/ukb-non-eur-exclusion-list.txt")

## all eids
## Size of fam file = 488377 = size of linker file
# eids_used <- linker %>%
# 	dplyr::filter(!app %in% withdrawals, 
# 				  !ieu %in% exclusions$id) %>%
# 	pull(app)

# writeLines(as.character(eids_used), "data/ukb-non-eur-inclusion-list.txt")

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


# ## CHECKING NUMBERS
# def_exc <- c("sex_mismatch", "sex_chr_aneuploidy", "het_miss_outliers")
# other_exc <- c("minimal_relateds", "non_europeans")

# exclusions_def <- exclusions %>%
# 	dplyr::filter(exc_reason %in% def_exc)

# ## Expected exclusions
# severe_no_def_exc <- severe_dat %>%
# 	dplyr::filter(!f.eid %in% withdrawals) %>%
# 	left_join(linker, by = c("f.eid" = "app")) %>%
# 	dplyr::filter(!ieu %in% exclusions_def$id)

# table(severe_no_def_exc$d4) # severe: 605 cases, 29646 controls
# table(severe_no_def_exc$d9) # mod-severe: 6125 cases, 24126 controls

# ## Minimal related exclusions
# severe_norel <- severe_no_def_exc %>%
# 	dplyr::filter(!ieu %in% exclusions[exclusions$exc_reason == "minimal_relateds", ]$id)

# table(severe_norel$d4) # severe: 509 cases, 24754 controls
# table(severe_norel$d9) # mod-severe: 5091 cases, 20172 controls

# ## Non-European exclusions
# severe_noeur <- severe_no_def_exc %>%
# 	dplyr::filter(!ieu %in% exclusions[exclusions$exc_reason == "non_europeans", ]$id)

# table(severe_noeur$d4) # severe: 582 cases, 28353 controls
# table(severe_noeur$d9) # mod-severe: 5794 cases, 23141 controls

# ## Minimal related and Non-European exclusions
# severe_norel_eur <- severe_no_def_exc %>%
# 	dplyr::filter(!ieu %in% exclusions[exclusions$exc_reason %in% other_exc, ]$id)

# table(severe_norel_eur$d4) # severe: 490 cases, 23568 controls
# table(severe_norel_eur$d9) # mod-severe: 4788 cases, 19270 controls

