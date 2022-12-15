# ------------------------------------------------------------------------
# Define severe and non-severe AD cases in UKB
# ------------------------------------------------------------------------

## Aim: Use selected UKBiobank data to define individuals as severe and non-severe AD cases and exclude the rest
## NOTE: You need to have extracted specific variables from the UKB phenotype data present in the RDSF

## pkgs
library(data.table)
library(tidyverse)

## args
args <- commandArgs(trailingOnly = TRUE)
id_file <- args[1]
phototherapy_file <- args[2]
treatments_file <- args[3]
icd10_file <- args[4]
icd9_file <- args[5]
icd10_all_file <- args[6]
sr_file <- args[7]
count_outfile <- args[8]
data_outfile <- args[9]

## manual args
id_file <- "data/ukb-pheno/ids.txt"
phototherapy_file <- "data/ukb-pheno/phototherapy.txt"
treatments_file <- "data/ukb-pheno/treatments.txt"
icd10_file <- "data/ukb-pheno/icd10.txt"
icd9_file <- "data/ukb-pheno/icd9.txt"
icd10_all_file <- "data/ukb-pheno/icd10-allcodes.txt"
sr_file <- "data/ukb-pheno/selfreport.txt"
summary_outfile <- "data/severity-counts.tsv" 
data_outfile <- "data/case-pheno.tsv"

gp_record_file <- "data/ukb-pheno/gp-records.txt"
gp_record <- fread(gp_record_file)

## data
ids <- fread(id_file)
pt_dat <- fread(phototherapy_file)
treat_dat <- fread(treatments_file)
icd10_dat <- fread(icd10_file)
icd9_dat <- fread(icd9_file)
icd10_all_dat <- fread(icd10_all_file)
sr_dat <- fread(sr_file)
## NOTE: currently these files = ~2Gb 

gp_scripts_ids <- fread("data/ukb-pheno/gp_scripts_ids.txt")
uniq_gp_ids <- unique(gp_scripts_ids$eid)

sum(ids$f.eid %in% uniq_gp_ids)

## To do (2022-11-07):
## 1. Finalise severe definition and implement that here
## 2. Get all eczema cases (including non-severe)
## 3. Write up meta-data table -- i.e. how codes match to variable names etc.

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
	n_dat <- map_dfr(variables, function(v) {
		tibble(variable = as.character(v), n = sum(count_dat[[as.character(v)]] != 0))
	})
	total <- tibble(variable = "total", n = sum(n_dat$n))

	n_dat <- bind_rows(n_dat, total)
	return(n_dat)
}

# ------------------------------------------------------------------------
# Phototherapy
# ------------------------------------------------------------------------

pt_meta <- tibble(field_id = 41272, 
				  code = c("S121", "S122", "S123", "S124", "S128", "S129"), 
				  meaning = c("Ultraviolet A light therapy to skin", 
				  			  "Ultraviolet B light therapy to skin", 
				  			  "Combined photochemotherapy and ultraviolet A light therapy to skin", 
				  			  "Combined photochemotherapy and ultraviolet B light therapy to skin", 
				  			  "Other specified phototherapy to skin", 
				  			  "Unspecified phototherapy to skin"))

photo_vars <- c("S121", "S122", "S123", "S124", "S128", "S129")
pt_count <- count_pheno(ids, pt_dat, photo_vars)
pt_n <- get_n(pt_count, photo_vars)

## Checking s12 - expecting 0
pt_dat <- cbind(ids, pt_dat)
s12_count <- apply(pt_dat, 1, function(x) {
	sum(x == "S12", na.rm=T)
})
s12_count

rm(pt_dat)

# ------------------------------------------------------------------------
# Systemic treatments
# ------------------------------------------------------------------------

# mtx - methotrexate (1140910036)
# methotrexate (1140869848)
# dupilumab isn't present for some reason

sys_meta <- tibble(field_id = 20003, 
				   code = c(1140909844, 1141181020, 1140909864, 1140869930, 1140925978, 1140910036, 1140869848), 
				   meaning = c("ciclosporin",
				   			   "ciclosporin product", 
				   			   "azt - azathioprine", 
				   			   "azathioprine", 
				   			   "mycophenolate", 
				   			   "mtx - methotrexate", 
				   			   "methotrexate"))

sys_treat_vars <- c(1140909844, 1141181020, 1140909864, 1140869930, 1140925978, 1140910036, 1140869848)
sys_treat_count <- count_pheno(ids, treat_dat, sys_treat_vars)
sys_n <- get_n(sys_treat_count, sys_treat_vars)

# ------------------------------------------------------------------------
# Corticosteroid treatments
# ------------------------------------------------------------------------

cort_meta <- tibble(field_id = 20003, 
					code = c(1140884696, 1141157294, 1140874896, 1140910424, 1141189464, 1140909786, 1140874790, 
			   				 1140888074, 1140888098, 1140888172),
					meaning = c("clobetasone", 
								"hydrocortisone product",
								"hydrocortisone", 
								"hc - hydrocortisone", 
								"exe-cort hydrocortisone 1% cream", 
								"beclometasone", 
								"Betamethasone", 
								"clobetasol", 
								"fluticasone", 
								"mometasone"))

cort_vars <- c(1140884696, 1141157294, 1140874896, 1140910424, 1141189464, 1140909786, 1140874790, 
			   1140888074, 1140888098, 1140888172)
cort_count <- count_pheno(ids, treat_dat, cort_vars)
cort_n <- get_n(cort_count, cort_vars)

# ------------------------------------------------------------------------
# ICD 10 codes
# ------------------------------------------------------------------------

icd10_meta <- tibble(field_id = 41202, 
					 code = c("L20", "L208", "L209"), 
					 meaning = c("Atopic dermatitis", "Other atopic dermatitis", "Atopic dermatitis, unspecified"))

icd10_codes <- c("L20", "L208", "L209")
icd10_count <- count_pheno(ids, icd10_dat, icd10_codes)
icd10_n <- get_n(icd10_count, icd10_codes)

# ------------------------------------------------------------------------
# ICD 9 codes
# ------------------------------------------------------------------------

icd9_meta <- tibble(field_id = 41203, 
					code = c("691", "6918", "69180"), 
					meaning = c("Atopic dermatitis and related conditions", 
								"Other atopic dermatitis and related conditions", 
								"Atopic dermatitis"))

icd9_codes <- c("691", "6918", "69180")
icd9_count <- count_pheno(ids, icd9_dat, icd9_codes)
icd9_n <- get_n(icd9_count, icd9_codes)

# ------------------------------------------------------------------------
# All AD cases
# ------------------------------------------------------------------------

ad_meta <- tibble(field_id = 41270, 
				  code = c("L20", "L208", "L209"), 
				  meaning = c("Atopic dermatitis", "Other atopic dermatitis", "Atopic dermatitis, unspecified"))

ad_codes <- c("L20", "L208", "L209")
ad_count <- count_pheno(ids, icd10_all_dat, ad_codes)

ecz_meta <- tibble(field_id = 20002, 
				   code = 1452, 
				   meaning = c("eczema/dermatitis"))

ecz_sr_code <- c(1452)
ecz_count <- count_pheno(ids, sr_dat, ecz_sr_code)

ad_dat <- left_join(ad_count, ecz_count) %>%
	mutate(ad = if_else(if_any(-f.eid, ~ . > 0), 1, 0))

ad_n <- get_n(ad_dat, c(ad_codes, ecz_sr_code))

ad_ids <- ad_dat %>%
	dplyr::filter(ad == 1) %>%
	pull(f.eid) %>%
	unique() %>%
	as.character()

writeLines(ad_ids, con = "data/ad-ids.txt")

# ------------------------------------------------------------------------
# Combine data and write it out
# ------------------------------------------------------------------------

## Quick combination of numbers 
comb_count <- list(pt_count, sys_treat_count, icd10_count, icd9_count) %>%
	reduce(left_join) %>%
	mutate(severe = if_else(if_any(-f.eid, ~ . > 0), 1, 0)) # asking if any are present then == 1, otherwise 0

comb_n <- list(pt_n, sys_n, icd10_n, icd9_n) %>%
	bind_rows() %>%
	dplyr::filter(variable != "total")
total <- tibble(variable = "total", n = sum(comb_n$n))
comb_n <- bind_rows(comb_n, total) ## some people will have multiple procedures/therapies
sum(comb_count$severe) ## True number of severe cases -- NOPE
## What I need is to verify what are eczema cases!

## CASES SHIT - NEED TO TIDY THIS UP SOOOOO BADLY
out_dat <- comb_count %>%
	dplyr::select(f.eid, severe) %>%
	left_join(ad_dat[, c("f.eid", "ad")]) %>%
	dplyr::filter(ad == 1)

sum(out_dat$severe)
nrow(out_dat)/100

## all count
all_count <- comb_count %>%
	dplyr::filter(f.eid %in% out_dat$f.eid) %>%
	left_join(ad_dat, by = c("f.eid" = "f.eid")) %>%
	dplyr::select(-severe, -ad)

## Need to indicate where icd10 codes are primary diagnosis 
colnames(all_count)[grepl("\\.x", colnames(all_count))] <- gsub("\\.x",
																"_primary_diagnosis", 
																colnames(all_count)[grepl("\\.x", colnames(all_count))])

colnames(all_count)[grepl("\\.y", colnames(all_count))] <- gsub("\\.y",
																"_any_diagnosis", 
																colnames(all_count)[grepl("\\.y", colnames(all_count))])

all_n <- get_n(all_count, colnames(all_count)[-1])

all_meta <- list(pt_meta, sys_meta, icd10_meta, icd9_meta, ad_meta, ecz_meta) %>%
	map(mutate_all, as.character) %>%
	bind_rows() %>%
	left_join(all_n, by = c("code" = "variable"))

l20_n <- all_n[grep("L20", all_n$variable),]
l20_n$field_id <- c(rep("41202", 3), rep("41270", 3))
l20_n$variable <- gsub("_.*", "", l20_n$variable)
l20_n <- l20_n %>%
	left_join(dplyr::select(all_meta, field_id, variable = code, meaning)) 
colnames(l20_n)[1] <- "code"

all_meta <- all_meta[!is.na(all_meta$n),] %>%
	bind_rows(l20_n)

all_meta %>% as.data.frame

write.table(all_meta, file = summary_outfile, col.names = T, row.names = F, quote = F, sep = "\t")

## Combination with corticosteroids

## WHAT WE WANT OUTPUT FOR META-DATA

## 
## AD CASES | N Severe

## LIMITED TO AD CASES
# UKB field_id | UKB code | description | N | 






