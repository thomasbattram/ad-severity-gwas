# ------------------------------------------------------------------------
# Define severe and non-severe AD cases in UKB
# ------------------------------------------------------------------------

## Aim: Use selected UKBiobank data to define individuals as severe and non-severe AD cases and exclude the rest
## NOTE: You need to have extracted specific variables from the UKB phenotype data present in the RDSF

## pkgs
library(data.table)
library(tidyverse)
library(openxlsx)

## args
args <- commandArgs(trailingOnly = TRUE)
id_file <- args[1]
ad_id_file <- args[2]
gp_scripts_id_file <- args[3]
phototherapy_file <- args[4]
treatments_file <- args[5]
icd10_file <- args[6]
icd9_file <- args[7]
summary_outfile <- args[8]
data_outfile <- args[9]
gp_file <- args[10]

## manual args
id_file <- "data/ukb-pheno/ids.txt"
ad_id_file <- "data/ad-ids.txt"
gp_scripts_id_file <- "data/ukb-pheno/gp_scripts_ids.txt"
phototherapy_file <- "data/ukb-pheno/phototherapy.txt"
treatments_file <- "data/ukb-pheno/summ-sr-treatment-data-clean.tsv"
icd10_file <- "data/ukb-pheno/icd10.txt"
icd9_file <- "data/ukb-pheno/icd9.txt"
summary_outfile <- "data/severity-counts.xlsx" 
data_outfile <- "data/case-pheno.tsv"
gp_file <- "data/ukb-pheno/summ-gp-script-data-clean.tsv"
full_count_outfile <- "data/full-count-out.tsv"

## data
ids <- fread(id_file)
ad_ids <- readLines(ad_id_file)
gp_ids <- readLines(gp_scripts_id_file)
pt_dat <- fread(phototherapy_file)
treat_dat <- read_tsv(treatments_file)
icd10_dat <- fread(icd10_file)
icd9_dat <- fread(icd9_file)
gp_dat <- read_tsv(gp_file)
## NOTE: currently these files = ~2Gb 

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

#' define individuals as severe
#' 
#' @param count_dat combined combination of counts for each variable of interest - comb_count
#' @param variables vector of column names within count_dat that will be used to assign people as severe/not severe
#' @param severe_only logical. keep only the ID and the severe column. Default = TRUE
#' @return tibble with an ID column and a "severe" column indicating whether an individual has severe disease
extract_severe_cases <- function(count_dat, variables, severe_only = TRUE)
{
	out <- count_dat %>%
		dplyr::select(f.eid, all_of(variables)) %>%
		mutate(severe = if_else(if_any(-f.eid, ~ . > 0), 1, 0))
	
	if (!severe_only) return(out)
	
	out %>%
		dplyr::select(f.eid, severe) %>%
		return()
}

#' get number of cases for a given definition of severe AD
#' 
#' @param count_dat combined combination of counts for each variable of interest - comb_count
#' @param systemics vector of systemic treatment names
#' @param topicals default = NULL. vector of topical steroid variables
#' @param other_treat default = NULL. list of other treatments, e.g. phototherapy treatments, ICD9 codes, ICD10 codes
#' @return tibble giving the numbers across each treatment and the unique cases 
get_definition_n <- function(count_dat, systemics, topicals = NULL, other_treat = NULL)
{
	gp_def <- paste0(c(systemics, topicals), "_gp")
	sr_def <- paste0(c(systemics, topicals), "_sr")
	all_vars <- c(gp_def, sr_def, unlist(other_treat, use.names = FALSE))
	count_def <- extract_severe_cases(count_dat, all_vars, severe_only = FALSE)
	def_n <- get_n(count_def, all_vars)
	no_cases_tab <- dplyr::filter(def_n, !grepl("cases", variable))
	## Remove gp and sr from vars
	no_cases_tab <- no_cases_tab %>%
		mutate(`Data type` = case_when(variable %in% gp_def ~ "GP prescription",
									   variable %in% sr_def ~ "Self-report")) %>%
		mutate(variable = gsub("_gp|_sr", "", variable))
	sort_tab <- map_dfr(1:nrow(no_cases_tab), function(i) {
		temp_tab <- no_cases_tab[i, ]
		v <- temp_tab$variable
		if (v %in% systemics) temp_tab[["General phenotype"]] <- "Systemics"
		if (v %in% topicals) temp_tab[["General phenotype"]] <- "Topicals"
		if (v %in% other_treat$icd9) temp_tab[["General phenotype"]] <- "Primary diagnosis - ICD9 code"
		if (v %in% other_treat$icd10) temp_tab[["General phenotype"]] <- "Primary diagnosis - ICD10 code"
		if (v %in% other_treat$photo) temp_tab[["General phenotype"]] <- "Phototherapy"
		return(temp_tab)
	})
	sort_tab <- sort_tab %>%
		arrange(`General phenotype`, variable)
	cases_tab <- def_n %>%
		dplyr::filter(grepl("cases", variable)) %>%
		mutate(`Data type` = NA, `General phenotype` = NA)
	out_tab <- bind_rows(sort_tab, cases_tab) %>%
		dplyr::select(`General phenotype`, variable, `Data type`, n)
	return(out_tab)
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
head(s12_count)
sum(s12_count)
rm(pt_dat)

# ------------------------------------------------------------------------
# Systemic treatments
# ------------------------------------------------------------------------

sr_drugs <- colnames(treat_dat)[colnames(treat_dat) != "f.eid"]
colnames(treat_dat)[colnames(treat_dat) != "f.eid"] <- paste0(sr_drugs, "_sr")

sr_treat_count <- ids %>%
	left_join(treat_dat)

sr_treat_count[sr_treat_count == TRUE] <- 1
sr_treat_count <- as_tibble(sr_treat_count)

## Sanity check
sum(treat_dat$Methotrexate_sr)
sum(sr_treat_count$Methotrexate_sr, na.rm = T)

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
# GP scripts
# ------------------------------------------------------------------------

## This is done in another script, just putting it into the same format as the others
gp_drugs <- colnames(gp_dat)[colnames(gp_dat) != "f.eid"]
colnames(gp_dat)[colnames(gp_dat) != "f.eid"] <- paste0(gp_drugs, "_gp")

gp_treat_count <- ids %>%
	left_join(gp_dat)

gp_treat_count[gp_treat_count == TRUE] <- 1
gp_treat_count <- as_tibble(gp_treat_count)

## For those with missing data, put as 0
gp_treat_count <- gp_treat_count %>%
	dplyr::filter(f.eid %in% ad_ids) 
gp_treat_count[is.na(gp_treat_count)] <- 0

## Sanity check
sum(gp_dat$Methotrexate_gp)
sum(gp_treat_count$Methotrexate_gp, na.rm = T)

# ------------------------------------------------------------------------
# Combine data and write it out
# ------------------------------------------------------------------------

## Quick combination of numbers 
comb_count <- list(pt_count, sr_treat_count, gp_treat_count, icd10_count, icd9_count) %>%
	reduce(left_join) %>%
	dplyr::filter(f.eid %in% gp_ids, f.eid %in% ad_ids)
	# mutate(severe = if_else(if_any(-f.eid, ~ . > 0), 1, 0)) # asking if any are present then == 1, otherwise 0

## Write out full count file
write.table(comb_count, file = full_count_outfile, col.names = T, row.names = F, quote = F, sep = "\t")

## Check numbers for each variety of severity definition
sys <- c("Methotrexate", "Ciclosporin", "Azathioprine", "Mycophenolate")
top1 <- c("Topical steroids", "Topical calcineurin inhibitors")
top2 <- c("Topical coal tar", "Topical zinc", "Topical ichthammol")
top3 <- c("Super potent steroids")
top4 <- c("Potent steroids")
top5 <- c("Super potent steroids x2")
top6 <- c("Potent steroids x2")

## Definition 1: systemics
def1_n <- get_definition_n(comb_count, sys)

## Definition 2: systemics or topical steroids or Topical calcineurin inhibitors
def2_n <- get_definition_n(comb_count, sys, top1)

## Definition 3: def 2 or Coal tar OR Topical zinc OR Ichthammol
def3_n <- get_definition_n(comb_count, sys, c(top1, top2))

## Definition 4: def 1 or ICD9 or ICD10 or phototherapy
other_treat <- list(photo = photo_vars, icd9 = icd9_codes, icd10 = icd10_codes)
def4_n <- get_definition_n(comb_count, sys, other_treat = other_treat)

## Definition 5: systemics or super potent topicals
def5_n <- get_definition_n(comb_count, sys, top3)

## Definition 6: systemics or super potent topicals or potent topicals
def6_n <- get_definition_n(comb_count, sys, c(top3, top4))

## Definition 7: systemics or super potent topicals x2
def7_n <- get_definition_n(comb_count, sys, top5)

## Definition 8: systemics or potent topicals x2 (or super potent topicals)
def8_n <- get_definition_n(comb_count, sys, top6)

## Definition 9: def 8 or ICD9 or ICD10 or phototherapy
def9_n <- get_definition_n(comb_count, sys, top6, other_treat = other_treat)

## Definition 10: def7 or ICD9 or ICD10 or phototherapy
def10_n <- get_definition_n(comb_count, sys, top5, other_treat = other_treat)

def_list <- list(d1 = def1_n, d2 = def2_n, d3 = def3_n, d4 = def4_n, 
				 d5 = def5_n, d6 = def6_n, d7 = def7_n, d8 = def8_n, 
				 d9 = def9_n, d10 = def10_n)

## Checking concordance between self-report and GP prescriptions
gp_def1 <- paste0(sys, "_gp")
sr_def1 <- paste0(sys, "_sr")
count_def1 <- comb_count %>%
	dplyr::select(f.eid, all_of(c(gp_def1, sr_def1))) %>%
	mutate(severe = if_else(if_any(-f.eid, ~ . > 0), 1, 0))

sum(count_def1$Methotrexate_sr == 1 & count_def1$Methotrexate_gp == 1) # 143
sum(count_def1$Ciclosporin_sr == 1 & count_def1$Ciclosporin_gp == 1) # 12
sum(count_def1$Azathioprine_sr == 1 & count_def1$Azathioprine_gp == 1) # 10
sum(count_def1$Mycophenolate_sr == 1 & count_def1$Mycophenolate_gp == 1) # 13


## QUICK AND MESSY WRITE UP - will want to alter manually later
wb <- createWorkbook()
lapply(seq_along(def_list), function(x) {
	nam <- names(def_list)[x]
	addWorksheet(wb, sheetName = nam)
	def_tab <- def_list[[nam]]
	writeData(wb, sheet = nam, x = def_tab,
				   colNames = TRUE, rowNames = FALSE)
})
saveWorkbook(wb, summary_outfile, overwrite = TRUE)

## Write out the data for analyses
severe_list <- lapply(seq_along(def_list), function(x) {
	dat <- def_list[[x]] %>%
		mutate(new_var = case_when(`Data type` == "GP prescription" ~ paste0(variable, "_gp"), 
								   `Data type` == "Self-report" ~ paste0(variable, "_sr"), 
								   is.na(`Data type`) ~ variable))
	vars <- dat$new_var[!(dat$new_var %in% c("unique_cases", "non_unique_cases"))]
	out <- extract_severe_cases(comb_count, vars, severe_only = TRUE)
	colnames(out)[colnames(out) == "severe"] <- paste0("d", x)
	return(out)
})

severe_out_tab <- reduce(severe_list, left_join)
write.table(severe_out_tab, file = data_outfile, col.names = T, row.names = F, quote = F, sep = "\t")
