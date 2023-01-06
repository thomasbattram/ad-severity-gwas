# ------------------------------------------------------------------------
# Read in the self-reported treatment data and extract only AD cases
# ------------------------------------------------------------------------

## Self-reported treatment data will also need to be looked over manually to extract the topicals

## pkgs
library(data.table)
library(openxlsx) # to write out the data as a spreadsheet
library(tidyverse)

## args
args <- commandArgs(trailingOnly = TRUE)
id_file <- args[1]
ad_id_file <- args[2]
treatments_file <- args[3]
outfile <- args[4]

## manual args
id_file <- "data/ukb-pheno/ids.txt"
ad_id_file <- "data/ad-ids.txt"
treatments_file <- "data/ukb-pheno/treatments.txt"
outfile <- "data/ukb-pheno/ad-sr-prescription-data.xlsx"

## data
ids <- fread(id_file)
ad_ids <- readLines(ad_id_file)
treat_dat <- fread(treatments_file)

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


## What we want out is a spreadsheet with these columns:
# code | n | drug | INCLUDE_FINAL | SYSTEMIC | Eczema specific treatment? | Non-eczema-specific anti-inflammatory topicals (eg topical steroids) | Emollients | Comments 
# 

## Steps:
# 1. limit the data to AD cases only
# 2. extract all unique variables (i.e. all codes within the dataset)
# 3. run the count_pheno function across all codes (may need to alter function a little)
# 4. run the get_n function across all codes (may need to alter function a little)
# 5. edit tibble so columns look like those stated above
# 6. output it as an excel spreadsheet