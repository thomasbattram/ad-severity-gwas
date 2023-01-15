# ------------------------------------------------------------------------
# Clean the UKB primary care data
# ------------------------------------------------------------------------

# squeue -p veryshort
# srun --job-name "InteractiveJobTest" --nodes=1 --ntasks-per-node=1 --cpus-per-task=1 --time=1:00:00 --mem=32GB --partition=test --pty bash

## At the end we want a table like:
# f.eid | corticosteroids | methotrexate | etc. | 

## pkgs
library(tidyverse)
library(data.table)
library(readxl)

## args
args <- commandArgs(trailingOnly = TRUE)
ids_file <- args[1]
ad_ids_file <- args[2]
sorted_gp_file <- args[3]
gp_scripts_file <- args[4]
treatments_file <- args[5]
sorted_sr_file <- args[6]
sr_full_outfile <- args[7]
sr_summ_outfile <- args[8]
gp_full_outfile <- args[9]
gp_summ_outfile <- args[10]

## manual args
ids_file <- "data/ukb-pheno/ids.txt"
ad_ids_file <- "data/ad-ids.txt"
sorted_gp_file <- "data/ukb-pheno/gp_scripts_ad-sorted.xlsx"
gp_scripts_file <- "data/ukb-pheno/gp_scripts_ecz_inclusive.txt"
treatments_file <- "data/ukb-pheno/treatments.txt"
sorted_sr_file <- "data/ukb-pheno/ad-sr-meds-ukb-sorted.xlsx"
sr_full_outfile <- "data/ukb-pheno/full-sr-treatment-data-clean.RData"
sr_summ_outfile <- "data/ukb-pheno/summ-sr-treatment-data-clean.tsv"
gp_full_outfile <- "data/ukb-pheno/full-gp-script-data-clean.RData"
gp_summ_outfile <- "data/ukb-pheno/summ-gp-script-data-clean.tsv"

## data
ids <- fread(ids_file)
ad_ids <- readLines(ad_ids_file)
sr_dat <- read_xlsx(sorted_sr_file)
gp_dat <- read_xlsx(sorted_gp_file, sheet = "drug_first_20", guess_max = 25000)
treat_dat <- fread(treatments_file)
gp_scripts <- fread(gp_scripts_file)
colnames(gp_scripts) <- c("f.eid", "read_2", "bnf_code", "dmd_code", "drug_name")

# ------------------------------------------------------------------------
# Functions to count variables within UKB data
# ------------------------------------------------------------------------

#' count number of times phenotypes occur in a given dataset by row
#' 
#' @param ids data.frame of IDs (rows = individuals). Set to NULL if they're already in `dat`
#' @param dat data.frame - phenotype data (rows = individuals)
#' @param variables vector of variables to count presence within each individual
#' @param ad_ids vector of IDs that are AD cases
#' @return tibble with a count of each variable across each individual
count_pheno <- function(ids, dat, variables, ad_ids = NULL)
{
	if (!is.null(ids)) {
		dat <- cbind(ids, dat)
	}
	if (!is.null(ad_ids)) {
		dat <- dat[dat$f.eid %in% ad_ids, ]
	}
	count_dat <- tibble(f.eid = dat$f.eid)

	for (v in variables) {
		message("counting occurences of ", v, " for each individual")
		count_dat[[as.character(v)]] <- apply(dat, 1, function(x) {
			sum(as.character(x) == as.character(v), na.rm = TRUE)
		})
	}

	return(count_dat)
}

#' extract all codes related to a gp prescription
#' 
#' @param dd gp_dat restricted to just rows related to one drug type (e.g. "methotrexate")
#' @param gps gp_scripts data.frame
#' @return list of read codes, bnf codes, and dmd codes related to the drug type
get_codes <- function(dd, gps)
{
	read <- unique(unlist(str_split(dd$read_2_concat, ";")))
	read <- read[read != "-"]
	bnf <- unique(unlist(str_split(dd$bnf_code_concat, ";")))
	bnf <- bnf[bnf != "-"]
	dmd <- unique(unlist(str_split(dd$dmd_code_concat, ";")))
	dmd <- dmd[dmd != "-"]
	test <- c(all(read %in% gps$read_2), 
			  all(bnf %in% gps$bnf_code), 
			  all(dmd %in% gps$dmd_code))
	if (all(test)) {
		return(list(read = read, bnf = bnf, dmd = dmd))
	} else {
		stop("Not all codes are in gp_scripts")
	}
}

#' extract IDs for all individuals who have been prescribed a specific drug using linked codes
#' 
#' @param gps gp_scripts data.frame
#' @param codes list of read codes, bnf codes, and dmd codes
#' @return vector of IDs 
get_ids <- function(gps, codes)
{
	readids <- unique(gps[gps$read_2 %in% codes$read, ]$f.eid)
	bnfids <- unique(gps[gps$bnf_code %in% codes$bnf, ]$f.eid)
	dmdids <- unique(gps[gps$dmd_code %in% codes$dmd, ]$f.eid)
	outids <- unique(c(readids, bnfids, dmdids))
	return(outids)
}


# ------------------------------------------------------------------------
# sort self-reported treatment data
# ------------------------------------------------------------------------

summary(sr_dat)
# every column with a drug name, e.g. Methotrexate, is either NA or 1. 1 indicates that code belongs to that drug category

drugs <- colnames(sr_dat)[!colnames(sr_dat) %in% c("code", "n", "drug")]

## Get all drug codes
drug_codes <- lapply(drugs, function(drug) {
	out <- sr_dat[sr_dat[[drug]] == 1, "code", drop=T]
	return(out[!is.na(out)])
})
names(drug_codes) <- drugs

ad_ids <- as.numeric(ad_ids)

## TEST
# test_count_dat <- count_pheno(ids, treat_dat, drug_codes[[1]], ad_ids = ad_ids)
# test_count_dat[[drugs[1]]] <- apply(test_count_dat, 1, function(x) {
# 	any(x == 1, na.rm = TRUE)
# })

## Count presence of each drug across each individual
full_sr_dat <- lapply(drugs, function(drug) {
	message("Extracting ", drug, " data")
	count_d <- count_pheno(ids, treat_dat, drug_codes[[drug]], ad_ids)
	count_d[[drug]] <- apply(count_d, 1, function(x) {
		any(x == 1, na.rm = TRUE)
	}) 
	return(count_d)
})
names(full_sr_dat) <- drugs

## summarise the data
summ_sr_dat <- lapply(drugs, function(drug) {
	full_sr_dat[[drug]] %>%
		dplyr::select(f.eid, all_of(drug))
})

summ_sr_tab <- reduce(summ_sr_dat, left_join)
summary(summ_sr_tab)

## Save it all boiii
save(full_sr_dat, file = sr_full_outfile)
write.table(summ_sr_tab, file = sr_summ_outfile, col.names = T, row.names = F, quote = F, sep = "\t")

## remove data to save some space
rm(list = c("full_sr_dat", "summ_sr_dat", "summ_sr_tab", "treat_dat"))

# ------------------------------------------------------------------------
# sort GP prescription data
# ------------------------------------------------------------------------

## SORT DMD CODES ##
## Problem occurred when copying numbers into excel, some got rounded and so
## need to be replaced with actual number from original data
ori_cleaned_gp_scripts <- fread("data/ukb-pheno/gp_scripts_ecz_inclusive.uniq.grouped")
ocgs <- ori_cleaned_gp_scripts %>%
	dplyr::select(dmd_correct = dmd_code_concat, drug)
all(gp_dat$drug %in% ocgs$drug)
merge_dmd <- gp_dat %>%
	dplyr::select(dmd_code_concat, drug) %>%
	left_join(ocgs)
gp_dat$dmd_code_concat <- merge_dmd$dmd_correct

## SORT BNF CODES ##
## found that some codes are in the sorted data but not in gp_scripts - looks 
## as though the 0 at the start of these codes has been removed for some reason
## updating them manually
gp_dat[gp_dat$bnf_code_concat == "8020200", "bnf_code_concat"] <- "08020200"
gp_dat[gp_dat$bnf_code_concat == "8020100", "bnf_code_concat"] <- "08020100"

## list drugs in spreadsheet
drugs <- c("Methotrexate", "Ciclosporin", "Azathioprine", "Mycophenolate", "AD antiinflammatory topicals", 
		   "Topical steroids", "Topical calcineurin inhibitors", "Topical coal tar", "Topical zinc", 
		   "Topical ichthammol", "Emollients")

## TEST
# drug_dat <- gp_dat[which(gp_dat[["Methotrexate"]] == 1),]
# read <- unique(unlist(str_split(drug_dat$read_2_concat, ";")))
# read <- read[read != "-"]
# all(read %in% gp_scripts$read_2)
# bnf <- unique(unlist(str_split(drug_dat$bnf_code_concat, ";")))
# bnf <- bnf[bnf != "-"]
# all(bnf %in% gp_scripts$bnf_code)
# bnf[which(!bnf %in% gp_scripts$bnf_code)]
# dmd <- unique(unlist(str_split(drug_dat$dmd_code_concat, ";")))
# dmd <- dmd[dmd != "-"]
# all(dmd %in% gp_scripts$dmd_code)

# readids <- unique(gp_scripts[gp_scripts$read_2 %in% read, ]$f.eid)
# bnfids <- unique(gp_scripts[gp_scripts$bnf_code %in% bnf, ]$f.eid)
# dmdids <- unique(gp_scripts[gp_scripts$dmd_code %in% dmd, ]$f.eid)
# outids <- unique(c(readids, bnfids, dmdids))

## For each drug, extract IDs of all individuals prescribed the drug
gp_outdat <- tibble(f.eid = unique(gp_scripts$f.eid))
gp_outlist <- lapply(drugs, function(drug) {
	message("Extracting ", drug, " data")
	drug_dat <- gp_dat[which(gp_dat[[drug]] == 1),]
	codes <- get_codes(drug_dat, gp_scripts)
	outids <- get_ids(gp_scripts, codes)
	gp_outdat[[drug]] <- ifelse(gp_outdat$f.eid %in% outids, TRUE, FALSE)
	return(outids)
})
names(gp_outlist) <- drugs

## Output data in tibble form
for (drug in drugs) {
	temp_ids <- gp_outlist[[drug]]
	gp_outdat[[drug]] <- ifelse(gp_outdat$f.eid %in% temp_ids, TRUE, FALSE)
}

write.table(gp_outdat, file = gp_summ_outfile, quote = F, col.names = T, row.names = F, sep = "\t")