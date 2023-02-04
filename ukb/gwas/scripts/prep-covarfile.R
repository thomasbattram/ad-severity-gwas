# ------------------------------------------------------------------------------
# Prepare covariate files for GWAS
# ------------------------------------------------------------------------------

## Prepare covariate file for running the GWAS using bolt LMM in UKB

## Format for covariate file
## FID IID covar1 covar2 covar3 ...

## Should be whitespace delimited. The phenofile needs to be numeric - case/control = 1/0
## covariates can be categorical but can't have whitespace in the text

## pkgs
library(tidyverse) # tidy code and data

## args
args <- commandArgs(trailingOnly = TRUE)
pheno_file <- args[1]
pc_file <- args[2]
age_file <- args[3]
linker_file <- args[4]
outfile <- args[5]

## manual args
# pheno_file <- "data/severe-ad-phenofile.txt"
# pc_file <- "UKB_GEN_DIR/derived/principal_components/data.pca1-10.plink.txt"
# age_file <- "../phenotype-extraction/data/ukb-pheno/age-sex.txt"
# linker_file <- "data/linker.csv"
# outfile <- "data/severe-ad-covarfile.txt"

## data
pheno_dat <- read_delim(pheno_file, delim = " ")
pc_dat <- read_delim(pc_file, delim = " ", col_names = FALSE)
colnames(pc_dat) <- c("FID", "IID", paste0("PC", 1:10))
age_sex_dat <- read_delim(age_file, delim = " ")
linker <- read_csv(linker_file)

# ------------------------------------------------------------------------------
# extract covariates
# ------------------------------------------------------------------------------

## Check to see if more missing data in "age at assessment"
sum(is.na(age_sex_dat$f.21003.0.0)) # 135687
sum(is.na(age_sex_dat$f.34.0.0)) # 1
## Definitely is more missing data! So going to use YOB 

year <- as.integer(format(Sys.Date(), "%Y"))

## sort age_dat
sorted_age <- age_sex_dat %>%
	left_join(linker, by = c("f.eid" = "app")) %>%
	mutate(age = year - f.34.0.0) %>%
	dplyr::select(FID = ieu, sex = f.31.0.0, age)

covar_dat <- pc_dat %>%
	left_join(sorted_age) %>%
	dplyr::select(FID, IID, all_of(paste0("PC", 1:10)), sex, age) %>%
	dplyr::filter(FID %in% pheno_dat$FID)

write.table(covar_dat, file = outfile, col.names = T, row.names = F, quote = F, sep = " ")


