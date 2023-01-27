# ------------------------------------------------------------------------
# Clean the UKB primary care data
# ------------------------------------------------------------------------

# srun --job-name "InteractiveJobTest" --nodes=1 --ntasks-per-node=1 --cpus-per-task=1 --time=1:00:00 --mem=32GB --partition=test --pty bash
# srun --job-name "InteractiveJobTest" --account=smed001801 --partition=test --nodes=1 --ntasks-per-node=1 --cpus-per-task=1 --time=1:00:00 --mem=16GB --pty bash

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
gp_scripts_nodose_file <- args[5]
gp_full_outfile <- args[6]
gp_summ_outfile <- args[7]
gp_tabs_outpre <- args[8]
gp_meta_xlsx_file <- args[9]
ori_gp_scripts_file <- args[10]

## manual args
ids_file <- "data/ukb-pheno/ids.txt"
ad_ids_file <- "data/ad-ids.txt"
sorted_gp_file <- "data/ukb-pheno/gp_scripts_ad-sorted.xlsx"
gp_scripts_file <- "data/ukb-pheno/gp_scripts_ecz_inclusive.txt"
gp_scripts_nodose_file <- "data/gp_scripts_ecz_inclusive.nodose.txt"
gp_full_outfile <- "data/ukb-pheno/full-gp-script-data-clean.RData"
gp_summ_outfile <- "data/ukb-pheno/summ-gp-script-data-clean.tsv"
gp_tabs_outpre <- "data/codelists/gp_"
gp_meta_xlsx_file <- "data/gp-script-meta.xlsx"
ori_gp_scripts_file <- "data/ukb-pheno/gp_scripts_ecz_inclusive.uniq.grouped"

## data
ids <- fread(ids_file)
ad_ids <- readLines(ad_ids_file)
gp_dat <- read_xlsx(sorted_gp_file, sheet = "drug_first_20", guess_max = 25000)
gp_scripts <- fread(gp_scripts_file)
gp_scripts_nodose <- fread(gp_scripts_nodose_file, header = FALSE)
ori_cleaned_gp_scripts <- fread(ori_gp_scripts_file)
colnames(gp_scripts) <- c("f.eid", "issue_date", "read_2", "bnf_code", "dmd_code", "drug_name")

# ------------------------------------------------------------------------
# Clean gp scripts file to match real drug names to those given to Ravi
# ------------------------------------------------------------------------

colnames(gp_scripts_nodose) <- c("read_2", "bnf_code", "dmd_code", "drug_nodose")

## Double checking it's in the same order
all(gp_scripts_nodose$read_2 == gp_scripts$read_2)
gp_scripts_nodose[which(gp_scripts_nodose$read_2 != gp_scripts$read_2),]
gp_scripts[which(gp_scripts_nodose$read_2 != gp_scripts$read_2),]
all(gp_scripts_nodose$bnf_code == gp_scripts$bnf_code)
all(gp_scripts_nodose$dmd_code == gp_scripts$dmd_code)
## we all good! 

## bind gp scripts data together
gp_scripts <- cbind(gp_scripts, gp_scripts_nodose[,c("drug_nodose")])

## check the drug names are present in nodose too!
all(gp_dat$drug %in% gp_scripts$drug_nodose)

rm(gp_scripts_nodose)

# ------------------------------------------------------------------------
# Functions to get the codes and extract IDs from UKB data
# ------------------------------------------------------------------------

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

# ------------------------------------------------------------------------
# Manually edit some of the codes
# ------------------------------------------------------------------------

## As the drugs were manually selected from an excel spreadsheet, some of
## the codes have problems because excel is the worst

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

# ------------------------------------------------------------------------
# List drugs of interest and clean gp_scripts for data extraction
# ------------------------------------------------------------------------

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

# ------------------------------------------------------------------------
# Extract numbers prescribed each drug and write it out
# ------------------------------------------------------------------------

gp_meta_full <- lapply(drugs, function(drug) {
	message(drug)
	## Get drug data and codes
	drug_dat <- gp_dat[which(gp_dat[[drug]] == 1),]
	# all_codes <- get_codes(drug_dat, gp_scripts_nom)
	all_drug_nams <- drug_dat$drug
	stopifnot(all(all_drug_nams %in% gp_scripts_nom$drug_nodose))
	gps <- gp_scripts_nom %>%
		dplyr::filter(drug_nodose %in% all_drug_nams)
	## Extract the number of participants prescribed each drug
	outlist <- lapply(1:length(all_drug_nams), function(x) {
		dr <- all_drug_nams[x]
		gps_single <- gps %>%
			dplyr::filter(drug_nodose == dr)
		gps_ids <- unique(gps_single$f.eid)
		out_all <- gps_single %>%
			dplyr::select(-f.eid, -issue_date) %>%
			distinct() %>%
			mutate(drug = drug_name, n = length(gps_ids)) %>%
			dplyr::select(drug, read_2, bnf = bnf_code, dmd = dmd_code, n)

		out_list <- list(drug = out_all[, c("drug", "n")], 
						 codes = out_all[, c("drug", "read_2", "bnf", "dmd")], 
						 ids = gps_ids)
		return(out_list)	
	})
	drug_tab <- map_dfr(outlist, "drug") %>%
		distinct()
	codes_tab <- map_dfr(outlist, "codes") %>%
		distinct()
	## Get total participants prescribed the drug of interest
	outids <- unique(unlist(map(outlist, "ids")))
	total_drug <- length(outids)
	total_row <- tibble(drug = toupper(drug), n = total_drug)
	out_tab <- bind_rows(drug_tab, total_row)
	## write out codes
	temp_outfile <- paste0(gp_tabs_outpre, tolower(gsub(" ", "_", drug)), ".tsv")
	write.table(codes_tab, file = temp_outfile, col.names = T, row.names = F, quote = F, sep = "\t")
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
	all_drug_nams <- drug_dat$drug
	stopifnot(all(all_drug_nams %in% gp_scripts_nom$drug_nodose))
	gps <- gp_scripts_nom %>%
		dplyr::filter(drug_nodose %in% all_drug_nams)
	## Extract the number of participants prescribed each drug
	all_ids <- unique(gps$f.eid)
	## This next step takes ~3 mins (when ~18,000 IDs to sort through)
	multi_gps <- map_dfr(seq_along(all_ids), function(x) {
		## UNCOMMENT THE CODE BELOW IF YOU WANT TO TEST THE IDENTIFIED PEOPLE ARE CORRECT
		id <- all_ids[x]
		temp_gps <- gps %>%
			dplyr::filter(f.eid == id)
		# print(temp_gps)
		# Sys.sleep(5)
		uniq_date <- unique(temp_gps$issue_date)
		if (length(uniq_date) > 1) {
			# message("\n\n\nYES\n\n\n")
			return(temp_gps)
		}
		# message("\n\n\nNO\n\n\n")
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



