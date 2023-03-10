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
treatments_code_file <- args[4]
outfile <- args[5]

## manual args
id_file <- "data/ukb-pheno/ids.txt"
ad_id_file <- "data/ad-ids.txt"
treatments_file <- "data/ukb-pheno/treatments.txt"
treatments_code_file <- "data/ukb-pheno/coding4.tsv" # need to download data coding 4: https://biobank.ndph.ox.ac.uk/showcase/coding.cgi?id=4
outfile <- "data/ukb-pheno/ad-sr-meds-ukb.xlsx"

## data
ids <- fread(id_file)
ad_ids <- readLines(ad_id_file)
treat_codes <- read_tsv(treatments_code_file)
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
	# dat <- cbind(ids, dat)
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
	# uniq_ids <- unique(unlist(all_ids))
	# total <- tibble(variable = c("unique_cases", "non_unique_cases"), 
	# 				n = c(length(uniq_ids), sum(n_dat$n))
	# 				)

	# n_dat <- bind_rows(n_dat, total)
	return(n_dat)
}

# ------------------------------------------------------------------------
# Extract the self-report data
# ------------------------------------------------------------------------

## What we want out is a spreadsheet with these columns:
# code | n | drug | INCLUDE_FINAL | SYSTEMIC | Eczema specific treatment? | Non-eczema-specific anti-inflammatory topicals (eg topical steroids) | Emollients | Comments 
# 

## Limit to AD cases only
td_filt <- cbind(ids, treat_dat) %>%
	dplyr::filter(f.eid %in% ad_ids)

## Extract all unique codes
cols <- colnames(td_filt)[colnames(td_filt) != "f.eid"]
all_vals <- lapply(cols, function(col) {
	unique(td_filt[[col]])
})
val_classes <- unlist(lapply(all_vals, class))
is.na(all_vals[val_classes == "logical"]) # all "logical" classes are NA
table(val_classes) # the rest of the classes are "integer" so should be fine to put them all together and select unique values
all_vals_vec <- unique(unlist(all_vals))
all_vals_vec <- all_vals_vec[!is.na(all_vals_vec)]
length(all_vals_vec) # 2275 unique treatments reported

## run count_pheno
start_time <- proc.time()
treat_count <- count_pheno(ids, td_filt, all_vals_vec)
time_taken <- proc.time() - start_time
time_taken # ~20 mins

## run get_n
treat_n <- get_n(treat_count, all_vals_vec)

## edit tibble
treat_codes$coding <- as.character(treat_codes$coding)
all(treat_n$variable %in% treat_codes$coding) # sanity check!

out_tib <- treat_n %>%
	dplyr::rename(coding = variable) %>%
	left_join(treat_codes) %>%
	dplyr::select(code = coding, n, drug = meaning) %>%
	mutate(INCLUDE_FINAL = "", SYSTEMIC = "", `Eczema specific treatment?` = "", 
		   `Non-eczema-specific anti-inflammatory topicals (eg topical steroids)` = "",
		    Emollients = "",  Comments = "") %>%
	arrange(desc(n))

## write it out as a spreadsheet
write.xlsx(out_tib, file = outfile)

