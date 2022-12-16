#!/bin/bash

# --------------------------------------------------------------------
# Script to extract UKB variables from UKB data in RDSF space
# --------------------------------------------------------------------

## NOTE: needs to be run on login node or local dir as data in RDSF

## rdsf dir
rd='/projects/MRC-IEU/research/data/ukbiobank/phenotypic/applications/15147/released/2021-04-22/data'
wd='/user/work/tb13101/severity-gwas/data/ukb-pheno'
cd ${wd}

## get ids
awk '{print $1}' ${rd}/data.45723.tab > ids.txt

## phototherapy - field ID = 41272
## Get columns 5580-5696
cat ${rd}/data.45723.tab | awk ' BEGIN { first = 5580; last = 5696 }
{ for (i = first; i < last; i++) { printf("%s ", $i) } print $last }' > phototherapy.txt
less -S phototherapy.txt

## topical corticosteroids and systematic treatments - under same field ID = 20003
## Get columns 2325-2516
cat ${rd}/data.45723.tab | awk ' BEGIN { first = 2325; last = 2516 }
{ for (i = first; i < last; i++) { printf("%s ", $i) } print $last }' > treatments.txt

## Hospital admission for eczema - eczema needs to be primary diagnosis field ID = 41202 and 41203
## Get columns 4467-4532 (ICD 10 codes)
cat ${rd}/data.45723.tab | awk ' BEGIN { first = 4467; last = 4532 }
{ for (i = first; i < last; i++) { printf("%s ", $i) } print $last }' > icd10.txt
## Get columns 4533-4560 (ICD 9 codes)
cat ${rd}/data.45723.tab | awk ' BEGIN { first = 4533; last = 4560 }
{ for (i = first; i < last; i++) { printf("%s ", $i) } print $last }' > icd9.txt

## eczema diagnosed by ICD10 code - 41270
## Get columns 5320-5532
cat ${rd}/data.45723.tab | awk ' BEGIN { first = 5320; last = 5532 }
{ for (i = first; i < last; i++) { printf("%s ", $i) } print $last }' > icd10-allcodes.txt

## eczema/dermatitis self-reported - 20002
## Get column 2189 - 2324
cat ${rd}/data.45723.tab | awk ' BEGIN { first = 2189; last = 2324 }
{ for (i = first; i < last; i++) { printf("%s ", $i) } print $last }' > selfreport.txt

## GP prescribed data
## Get columns 6145 - 6147
cat ${rd}/data.45723.tab | awk ' BEGIN { first = 6145; last = 6147 }
{ for (i = first; i < last; i++) { printf("%s ", $i) } print $last }' > gp-records.txt

## GP scripts
awk '{print $1}' gp_scripts.txt > gp_scripts_ids.txt