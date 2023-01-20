# ------------------------------------------------------------------------
# Clean the UKB primary care data
# ------------------------------------------------------------------------

# squeue -p test
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
sorted_gp_file <- args[3]
gp_scripts_file <- args[4]
treatments_file <- args[5]
treatments_num_file <- args[6]
sorted_sr_file <- args[7]
sr_full_outfile <- args[8]
sr_summ_outfile <- args[9]
sr_tabs_outpre <- args[10]
sr_meta_xlsx_file <- args[11]
gp_full_outfile <- args[12]
gp_summ_outfile <- args[13]
gp_tabs_outpre <- args[14]
gp_meta_xlsx_file <- args[15]
ori_gp_scripts_file <- args[16]

## manual args
ids_file <- "data/ukb-pheno/ids.txt"
ad_ids_file <- "data/ad-ids.txt"
sorted_gp_file <- "data/ukb-pheno/gp_scripts_ad-sorted.xlsx"
gp_scripts_file <- "data/ukb-pheno/gp_scripts_ecz_inclusive.txt"
treatments_file <- "data/ukb-pheno/treatments.txt"
treatments_num_file <- "data/ukb-pheno/treatments_num.txt"
sorted_sr_file <- "data/ukb-pheno/ad-sr-meds-ukb-sorted.xlsx"
sr_full_outfile <- "data/ukb-pheno/full-sr-treatment-data-clean.RData"
sr_summ_outfile <- "data/ukb-pheno/summ-sr-treatment-data-clean.tsv"
sr_tabs_outpre <- "data/codelists/sr_"
sr_meta_xlsx_file <- "data/sr-script-meta.xlsx"
gp_full_outfile <- "data/ukb-pheno/full-gp-script-data-clean.RData"
gp_summ_outfile <- "data/ukb-pheno/summ-gp-script-data-clean.tsv"
gp_tabs_outpre <- "data/codelists/gp_"
gp_meta_xlsx_file <- "data/gp-script-meta.xlsx"
ori_gp_scripts_file <- "data/ukb-pheno/gp_scripts_ecz_inclusive.uniq.grouped"

## data
ids <- fread(ids_file)
ad_ids <- readLines(ad_ids_file)
sr_dat <- read_xlsx(sorted_sr_file)
gp_dat <- read_xlsx(sorted_gp_file, sheet = "drug_first_20", guess_max = 25000)
treat_dat <- fread(treatments_file)
treat_num_dat <- fread(treatments_num_file)
gp_scripts <- fread(gp_scripts_file)
ori_cleaned_gp_scripts <- fread(ori_gp_scripts_file)
colnames(gp_scripts) <- c("f.eid", "issue_date", "read_2", "bnf_code", "dmd_code", "drug_name")

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

# altering to keep consistent with the GP data
colnames(sr_dat)[colnames(sr_dat) == "Super-potent steroids"] <- "Super potent steroids"
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

## summarise the data
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

wb <- createWorkbook()
lapply(drugs, function(drug) {
	temp_outfile <- paste0(sr_tabs_outpre, tolower(gsub(" ", "_", drug)), ".tsv")
	tab <- read_tsv(temp_outfile)
	addWorksheet(wb, sheetName = drug)
	writeData(wb, sheet = drug, x = tab,
				   colNames = TRUE, rowNames = FALSE)
})
saveWorkbook(wb, sr_meta_xlsx_file, overwrite = TRUE)

## remove data to save some space
rm(list = c("full_sr_dat", "summ_sr_dat", "summ_sr_tab", "treat_dat"))

# ------------------------------------------------------------------------
# sort GP prescription data
# ------------------------------------------------------------------------

## SORT DMD CODES ##
## Problem occurred when copying numbers into excel, some got rounded and so
## need to be replaced with actual number from original data
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
		   "Topical ichthammol", "Emollients", "Super potent steroids", "Potent steroids", 
		   "Moderately potent steroids", "Mild steroids")

## Remove where drug_name == "[missing]"
gp_scripts_nom <- gp_scripts %>%
	dplyr::filter(drug_name != "[missing]")

gp_scripts_m <- gp_scripts %>%
	dplyr::filter(drug_name == "[missing]")

## What codes do the drugs with "[missing]" in drug_name column have?
sum(gp_scripts_m$bnf_code == "-")
sum(gp_scripts_m$dmd_code == "-")
sum(gp_scripts_m$read_2 == "-")
## Mostly read codes - probably because these codes are related to diagnoses and
## other things

nrow(gp_scripts) - nrow(gp_scripts_nom)
# removed 1762372 rows

## write out the tables of definitions and numbers
codetypes <- c("read", "bnf", "dmd") # these should match those from the "get_codes" function
names(codetypes) <- c("read_2", "bnf_code", "dmd_code") # these should match the columns from "gp_scripts"


#' Sort drug names 
#' 
#' Some codes are linked to lots of drugs, so this function extracts only
#' individuals who have drugs that were originally identified as related to AD
#' 
#' @param drug_nams vector all the drug names related to the drug of interest from gp_dat
#' @param gps gp_scripts data
#' @return sorted gp scripts data
sort_drug_names <- function(drug_nams, gps)
{
	drug_nams <- gsub("\\(!mg_dose!\\)_?", "", drug_nams)
	drug_nams <- gsub("\\(!%_conc!\\)_?", "", drug_nams)
	drug_nams <- gsub("\\(!gms_weight!\\)_?", "", drug_nams)
	drug_nams <- gsub("\\(!cm_length!\\)_?", "", drug_nams)
	drug_nams <- gsub("\\(!ml_vol!\\)_?", "", drug_nams)
	drug_nams <- gsub("\\*", "", drug_nams)
	drug_nams <- gsub("_.*", "", drug_nams)
	gps <- gps[grep(paste0(drug_nams, collapse = "|"), gps$drug_name), ]
	return(gps)
}

#' extract IDs from GP data
#' 
#' @param gp_col column from gp_scripts
#' @param code code used to extract the data
#' @param ct codetype - read, bnf, or dmd
#' @param gps gp_scripts data
#' @param drug_dat gp_dat restricted to the drug of interest
#' @param multi_presc logical. restrict to multiple prescriptions?
#' @return list of IDs and a table of counts
extract_ids <- function(gp_col, code, ct, gps, drug_dat, multi_presc = FALSE)
{
	g_dat_col <- paste0(gp_col, "_concat")
	gps <- gps[gps[[gp_col]] == code, ] %>%
		distinct()
	if (ct == "bnf" | ct == "dmd") {
		drug_nams <- drug_dat[grepl(code, drug_dat[[g_dat_col]]), "drug", drop = T]
		gps <- sort_drug_names(drug_nams, gps)
	}
	if (multi_presc) {
		ids <- unique(gps[duplicated(gps$f.eid), "f.eid", drop = T])
	} else {
		ids <- unique(gps$f.eid)
	}
	code_count <- gps %>%
		dplyr::filter(f.eid %in% ids, !duplicated(f.eid))
	drug_names <- unique(gps$drug_name)
	out_tab <- tibble(codetype = ct, code = code, 
					  drug = paste0(drug_names, collapse = " | "), 
					  n = nrow(code_count))
	out_ids <- code_count$f.eid
	return(list(ids = out_ids, tab = out_tab))
}

gp_meta_full <- lapply(drugs, function(drug) {
	message(drug)
	## Get drug data and codes
	drug_dat <- gp_dat[which(gp_dat[[drug]] == 1),]
	all_codes <- get_codes(drug_dat, gp_scripts_nom)
	## Extract the number of participants matching each code for each codetype
	outlist <- lapply(1:length(codetypes), function(x) {
		ct <- codetypes[x]
		gp_col <- names(codetypes)[x]
		g_dat_col <- paste0(gp_col, "_concat")
		codes <- all_codes[[ct]]
		c_list <- lapply(codes, function(cod) {
			extract_ids(gp_col, cod, ct, gp_scripts_nom, drug_dat)
		})
		c_tab <- map_dfr(c_list, "tab")
		codetype_ids <- unique(unlist(map(c_list, "ids")))
		## Make a row for each codetype to see whether there is much difference in cases per codetype
		total_drug <- length(codetype_ids)
		total_row <- tibble(codetype = toupper(ct), code = NA, drug = toupper(drug), n = total_drug)
		clist <- list(full = c_tab, last_row = total_row, ids = codetype_ids)
		return(clist)	
	})
	## Bind tables together, making sure that codetypes are at the bottom of the table
	out_tab <- map(outlist, "full") %>%
		bind_rows() %>%
		bind_rows(map(outlist, "last_row")) %>%
		bind_rows()
	## Get total participants prescribed the drug of interest
	outids <- unique(unlist(map(outlist, "ids")))
	total_drug <- length(outids)
	total_row <- tibble(codetype = NA, code = NA, drug = toupper(drug), n = total_drug)
	out_tab <- bind_rows(out_tab, total_row)
	temp_outfile <- paste0(gp_tabs_outpre, tolower(gsub(" ", "_", drug)), ".tsv")
	write.table(out_tab, file = temp_outfile, col.names = T, row.names = F, quote = F, sep = "\t")
	return(list(tab = out_tab, ids = outids))
})
names(gp_meta_full) <- drugs
gp_ids_list <- map(gp_meta_full, "ids")
gp_tab_list <- map(gp_meta_full, "tab")

## Get results for prescription of multiple topicals
steroids <- c("Mild steroids", "Moderately potent steroids", "Potent steroids", "Super potent steroids")
steroid_potency <- factor(steroids, levels = steroids, ordered=TRUE)
po_drugs <- c("Super potent steroids", "Potent steroids")
gp_meta_sps2 <- lapply(po_drugs, function(drug) {
	message(drug)
	## Get drug data and codes
	more_potent_drugs <- as.character(steroid_potency[steroid_potency > drug])
	all_drugs <- c(drug, more_potent_drugs)
	drug_dat <- map_dfr(all_drugs, function(dru) {
		gp_dat[which(gp_dat[[dru]] == 1), ]	
	})
	all_codes <- get_codes(drug_dat, gp_scripts_nom)
	test_gps <- gp_scripts_nom %>%
		dplyr::filter(read_2 %in% all_codes$read | 
					  bnf_code %in% all_codes$bnf |
					  dmd_code %in% all_codes$dmd)
	drug_nams <- drug_dat$drug
	test_gps2 <- sort_drug_names(drug_nams, test_gps)
	all_ids <- unique(test_gps2$f.eid)
	## This next step takes ~3 mins (when ~18,000 IDs to sort through)
	multi_gps <- map_dfr(seq_along(all_ids), function(x) {
		id <- all_ids[x]
		temp_gps <- test_gps2 %>%
			dplyr::filter(f.eid == id, !duplicated(issue_date))
		if (nrow(temp_gps) > 1) return(temp_gps)
		return(NULL)
	})
	multi_ids <- unique(multi_gps$f.eid)
	drug_out <- ifelse(length(all_drugs) == 1, paste(toupper(drug), "x2"), paste(toupper(all_drugs), collapse = " OR "))
	out_tab <- tibble(drug = drug_out, 
					  n = length(multi_ids))
	return(list(tab = out_tab, ids = multi_ids))
})
names(gp_meta_sps2) <- paste(po_drugs, "x2")
gp_po_ids_list <- map(gp_meta_sps2, "ids")
gp_po_tab_list <- map(gp_meta_sps2, "tab")

gp_tab_list <- c(gp_tab_list, gp_po_tab_list)
gp_ids_list <- c(gp_ids_list, gp_po_ids_list)
gp_ids_list[["Super potent steroids x2"]] <- gp_po_ids_list[["Super potent steroids x2"]]
gp_ids_list[["Potent steroids x2"]] <- gp_po_ids_list[["Potent steroids x2"]]

# gp_tab_list <- lapply(drugs, function(drug) {read_tsv(paste0(gp_tabs_outpre, tolower(gsub(" ", "_", drug)), ".tsv"))})
# names(gp_tab_list) <- drugs

gp_summ_outdat <- tibble(f.eid = unique(gp_scripts$f.eid))
all_drugs <- c(drugs, paste(po_drugs, "x2"))
for (drug in all_drugs) {
	temp_ids <- gp_ids_list[[drug]]
	gp_summ_outdat[[drug]] <- ifelse(gp_summ_outdat$f.eid %in% temp_ids, TRUE, FALSE)
}
summary(gp_summ_outdat)

write.table(gp_summ_outdat, file = gp_summ_outfile, col.names = T, row.names = F, quote = F, sep = "\t")

wb <- createWorkbook()
lapply(all_drugs, function(drug) {
	addWorksheet(wb, sheetName = drug)
	writeData(wb, sheet = drug, x = gp_tab_list[[drug]],
				   colNames = TRUE, rowNames = FALSE)
})
saveWorkbook(wb, gp_meta_xlsx_file, overwrite = TRUE)

# ------------------------------------------------------------------------
# Looking at individual drugs in the gp_scripts data
# ------------------------------------------------------------------------

gp_dat %>%
	dplyr::filter(grepl("13.04.03.00.00", bnf_code_concat)) %>%
	dplyr::select(total_count, drug) %>%
	arrange(desc(total_count))

gp_scripts %>%
	dplyr::filter(grepl("Mometasone 0.", drug_name)) %>%
	dplyr::filter(grepl("cream", drug_name)) 

# Betnovate 0.1% cream (GlaxoSmithKline UK Ltd)
gp_scripts %>%
	dplyr::filter(drug_name == "Betnovate 0.1% cream (GlaxoSmithKline UK Ltd)") %>%
	dplyr::filter(!duplicated(f.eid)) %>%
	nrow 
# 2782

# Mometasone 0.1% cream
gp_scripts %>%
	dplyr::filter(drug_name == "Mometasone 0.1% cream") %>%
	dplyr::filter(!duplicated(f.eid)) %>%
	nrow 
# 2013
gp_scripts %>%
	dplyr::filter(drug_name == "Fucibet cream (LEO Pharma)") %>%
	dplyr::filter(!duplicated(f.eid)) %>%
	nrow 
# 2101



