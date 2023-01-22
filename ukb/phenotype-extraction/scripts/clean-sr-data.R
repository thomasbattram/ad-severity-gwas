# ------------------------------------------------------------------------
# Clean the UKB self-reported treatment data
# ------------------------------------------------------------------------

## Take the prescriptions highlighted by clinicians in the UKB self-reported treatment data
## and links it up to the raw UKB data so that you can see which individuals self-repot which
## treatments.

## NOTE: This script limits cases to those with GP prescription data (not all UKB participants)

## manual submission script
# srun --job-name "InteractiveJobTest" --nodes=1 --ntasks-per-node=1 --cpus-per-task=1 --time=1:00:00 --mem=32GB --partition=test --pty bash

## At the end we want a table like:
# f.eid | corticosteroids | methotrexate | etc. | 

## pkgs
library(tidyverse)
library(data.table)
library(readxl)
library(openxlsx)

## args
args <- commandArgs(trailingOnly = TRUE)
ids_file <- args[1]
ad_ids_file <- args[2]
treatments_file <- args[3]
treatments_num_file <- args[4]
gp_scripts_id_file <- args[5]
sorted_sr_file <- args[6]
sr_full_outfile <- args[7]
sr_summ_outfile <- args[8]
sr_tabs_outpre <- args[9]
sr_meta_xlsx_file <- args[10]

## manual args
ids_file <- "data/ukb-pheno/ids.txt"
ad_ids_file <- "data/ad-ids.txt"
treatments_file <- "data/ukb-pheno/treatments.txt"
treatments_num_file <- "data/ukb-pheno/treatments_num.txt"
gp_scripts_id_file <- "data/ukb-pheno/gp_scripts_ids.txt"
sorted_sr_file <- "data/ukb-pheno/ad-sr-meds-ukb-sorted.xlsx"
sr_full_outfile <- "data/ukb-pheno/full-sr-treatment-data-clean.RData"
sr_summ_outfile <- "data/ukb-pheno/summ-sr-treatment-data-clean.tsv"
sr_tabs_outpre <- "data/codelists/sr_"
sr_meta_xlsx_file <- "data/sr-script-meta.xlsx"

## data
ids <- fread(ids_file)
ad_ids <- readLines(ad_ids_file)
sr_dat <- read_xlsx(sorted_sr_file)
treat_dat <- fread(treatments_file)
treat_num_dat <- fread(treatments_num_file)
gp_script_ids <- readLines(gp_scripts_id_file)

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

# ------------------------------------------------------------------------
# Check data and extract all drugs/drug codes
# ------------------------------------------------------------------------

## Check what data looks like
str(sr_dat)
summary(sr_dat)
# every column with a drug name, e.g. Methotrexate, is either NA or 1. 1 indicates that code belongs to that drug category

# altering to keep consistent with the GP data
colnames(sr_dat)[colnames(sr_dat) == "Super-potent steroids"] <- "Super potent steroids"
drugs <- colnames(sr_dat)[!colnames(sr_dat) %in% c("code", "n", "drug")]

## Get all drug codes
drug_codes <- lapply(drugs, function(drug) {
	out <- sr_dat[sr_dat[[drug]] == 1, "code", drop=T]
	return(out[!is.na(out)])
})
names(drug_codes) <- drugs

## Limit ad_ids to just those in gp_scripts_ids
ad_ids <- ad_ids[ad_ids %in% gp_script_ids]
ad_ids <- as.numeric(ad_ids)

# ------------------------------------------------------------------------
# Count occurrence of each drug/drug combination across each individual
# ------------------------------------------------------------------------

## TEST
# test_count_dat <- count_pheno(ids, treat_dat, drug_codes[[1]], ad_ids = ad_ids)
# test_count_dat[[drugs[1]]] <- apply(test_count_dat, 1, function(x) {
# 	any(x == 1, na.rm = TRUE)
# })

## Count presence of each drug across each individual
full_sr_dat <- lapply(drugs, function(drug) {
	message("Extracting ", drug, " data")
	count_d <- count_pheno(ids, treat_dat, drug_codes[[drug]], ad_ids)
	count_d <- count_d %>%
		mutate(DRUG = if_else(if_any(-f.eid, ~ . > 0), TRUE, FALSE))
	colnames(count_d)[colnames(count_d) == "DRUG"] <- drug
	return(count_d)
})
names(full_sr_dat) <- drugs

## Look at multiple potent or multiple super potent steroids
po_sr_dat <- full_sr_dat[["Potent steroids"]] %>%
	left_join(full_sr_dat[["Super potent steroids"]]) %>%
	dplyr::select(-all_of(c("Potent steroids", "Super potent steroids"))) %>%
	mutate(count = rowSums(select(., -f.eid))) %>%
	mutate(`Potent steroids x2` = ifelse(count > 1, TRUE, FALSE))

spo_sr_dat <- full_sr_dat[["Super potent steroids"]] %>%
	dplyr::select(-all_of("Super potent steroids")) %>%
	mutate(count = rowSums(select(., -f.eid))) %>%
	mutate(`Super potent steroids x2` = ifelse(count > 1, TRUE, FALSE))

po_sr_list <- list(`Potent steroids x2` = po_sr_dat, 
				   `Super potent steroids x2` = spo_sr_dat)

full_sr_dat <- c(full_sr_dat, po_sr_list)

# ------------------------------------------------------------------------
# Summarise the data and write it out
# ------------------------------------------------------------------------

all_drugs <- c(drugs, "Potent steroids x2", "Super potent steroids x2")
summ_sr_dat <- lapply(all_drugs, function(drug) {
	full_sr_dat[[drug]] %>%
		dplyr::select(f.eid, all_of(drug))
})

summ_sr_tab <- reduce(summ_sr_dat, left_join)
summary(summ_sr_tab)


## Save it all boiii
save(full_sr_dat, file = sr_full_outfile)
write.table(summ_sr_tab, file = sr_summ_outfile, col.names = T, row.names = F, quote = F, sep = "\t")

## write out the tables of definitions and numbers
lapply(drugs, function(drug) {
	codes <- drug_codes[[drug]]
	fsd <- full_sr_dat[[drug]]
	out_tab <- map_dfr(codes, function(cod) {
		code_count <- sum(fsd[[cod]] > 0, na.rm = TRUE)
		out_tab <- sr_dat %>%
			dplyr::select(code, drug) %>%
			dplyr::filter(code == cod) %>%
			mutate(n = code_count)
		return(out_tab)
	})
	total_drug <- sum(fsd[[drug]], na.rm = TRUE)
	total_row <- tibble(code = NA, drug = toupper(drug), n = total_drug)
	out_tab <- bind_rows(out_tab, total_row)
	temp_outfile <- paste0(sr_tabs_outpre, tolower(gsub(" ", "_", drug)), ".tsv")
	write.table(out_tab, file = temp_outfile, col.names = T, row.names = F, quote = F, sep = "\t")
})

## Do the same for potent steroids x2 and super potent steroids x2
x2_drugs <- c("Potent steroids x2", "Super potent steroids x2")
lapply(x2_drugs, function(drug) {
	total_drug <- sum(full_sr_dat[[drug]][[drug]], na.rm = T)
	out_tab <- tibble(code = NA, drug = toupper(drug), n = total_drug)
	temp_outfile <- paste0(sr_tabs_outpre, tolower(gsub(" ", "_", drug)), ".tsv")
	write.table(out_tab, file = temp_outfile, col.names = T, row.names = F, quote = F, sep = "\t")		
})

## Summarise it all in a spreadsheet
wb <- createWorkbook()
lapply(all_drugs, function(drug) {
	temp_outfile <- paste0(sr_tabs_outpre, tolower(gsub(" ", "_", drug)), ".tsv")
	tab <- read_tsv(temp_outfile)
	addWorksheet(wb, sheetName = drug)
	writeData(wb, sheet = drug, x = tab,
				   colNames = TRUE, rowNames = FALSE)
})
saveWorkbook(wb, sr_meta_xlsx_file, overwrite = TRUE)

