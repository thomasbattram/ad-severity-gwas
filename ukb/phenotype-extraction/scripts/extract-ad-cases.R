# ------------------------------------------------------------------------
# Extract all AD cases from UKB
# ------------------------------------------------------------------------

## Aim: Use selected UKBiobank data to define individuals as AD cases 
## NOTE: You need to have extracted specific variables from the UKB phenotype data present in the RDSF

## pkgs
library(data.table)
library(tidyverse)

## args
args <- commandArgs(trailingOnly = TRUE)
id_file <- args[1]
diag_ad_cases_file <- args[2]
icd10_all_file <- args[3]
sr_file <- args[4]
ids_outfile <- args[5]

## manual args
id_file <- "data/ukb-pheno/ids.txt"
diag_ad_cases_file <- "data/gp-diag-ad-case-eids.txt.uniq"
icd10_all_file <- "data/ukb-pheno/icd10-allcodes.txt"
sr_file <- "data/ukb-pheno/selfreport.txt"
ids_outfile <- "data/ad-ids.txt"
summary_outfile <- "data/ad-cases-summary.tsv"

## data
ids <- fread(id_file)
diag_ad <- readLines(diag_ad_cases_file)
icd10_all_dat <- fread(icd10_all_file)
sr_dat <- fread(sr_file) 

# ------------------------------------------------------------------------
# Functions to count variables within UKB data
# ------------------------------------------------------------------------

#' count number of times phenotypes occur in a given dataset by row
#' 
#' @param ids data.frame of IDs (rows = individuals)
#' @param dat data.frame - phenotype data (rows = individuals)
#' @param variables vector of variables to count presence within each individual
#' @return tibble with a count of each variable across each individual
count_pheno <- function(ids, dat, variables)
{
	dat <- cbind(ids, dat)
	count_dat <- tibble(f.eid = dat$f.eid)

	for (v in variables) {
		message("counting occurences of ", v, " for each individual")
		count_dat[[as.character(v)]] <- apply(dat, 1, function(x) {
			sum(as.character(x) == as.character(v), na.rm = TRUE)
		})
	}

	return(count_dat)
}

#' count occurences of each variable in total across the dataset
#' 
#' @param count_dat output from count_pheno()
#' @param variables variables within count_dat
#' @return tibble with count of each variable across whole dataset
get_n <- function(count_dat, variables)
{
	## get samplesize for each and total phototherapy samplesize
	all_ids <- lapply(variables, function(v) {
		count_dat[count_dat[[as.character(v)]] != 0, "f.eid", drop=T]
	})
	names(all_ids) <- as.character(variables)
	n_dat <- map_dfr(variables, function(v) {
		cv <- as.character(v)
		tibble(variable = cv, n = length(all_ids[[cv]]))
	})
	uniq_ids <- unique(unlist(all_ids))
	total <- tibble(variable = c("unique_cases", "non_unique_cases"), 
					n = c(length(uniq_ids), sum(n_dat$n))
					)

	n_dat <- bind_rows(n_dat, total)
	return(n_dat)
}

# ------------------------------------------------------------------------
# AD cases
# ------------------------------------------------------------------------

## ICD10 codes
ad_meta <- tibble(field_id = 41270, 
				  code = c("L20", "L208", "L209"), 
				  meaning = c("Atopic dermatitis", "Other atopic dermatitis", "Atopic dermatitis, unspecified"))

ad_codes <- c("L20", "L208", "L209")
ad_count <- count_pheno(ids, icd10_all_dat, ad_codes)

## Self-report
ecz_meta <- tibble(field_id = 20002, 
				   code = "1452", 
				   meaning = c("eczema/dermatitis"))

ecz_sr_code <- c(1452)
ecz_count <- count_pheno(ids, sr_dat, ecz_sr_code)

## GP diagnosis
gp_rec_dat <- tibble(f.eid = ecz_count$f.eid) %>%
			  	mutate(gp_diagnosis = ifelse(f.eid %in% diag_ad, 1, 0)) 

ad_dat <- left_join(ad_count, ecz_count) %>%
	left_join(gp_rec_dat) %>%
	mutate(ad = if_else(if_any(-f.eid, ~ . > 0), 1, 0))

ad_ids <- ad_dat %>%
	dplyr::filter(ad == 1) %>%
	pull(f.eid) %>%
	unique() %>%
	as.character()

writeLines(ad_ids, con = ids_outfile)

ad_n <- get_n(ad_dat, c(ad_codes, ecz_sr_code, "gp_diagnosis"))
all_meta <- bind_rows(ad_meta, ecz_meta)
ad_n_out <- ad_n %>%
	left_join(all_meta, by = c("variable" = "code")) %>%
	dplyr::select(variable, n, ukb_code_meaning = meaning, ukb_field_id = field_id)

write.table(ad_n_out, file = summary_outfile, row.names = F, col.names = T, sep = "\t", quote = F)

