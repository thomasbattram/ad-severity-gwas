#!/bin/bash

# ---------------------------------------------------------------------------
# Extract UK Biobank AD cases from GP records
# ---------------------------------------------------------------------------

## This script extracts all UK Biobank IDs that have read codes corresponding to AD as 
## diagnosed by their GP

## Probably best to run this script manually on a compute node 
## using the command below and tmux to prevent timing out
# srun --job-name "InteractiveJob" --nodes=1 --ntasks-per-node=1 --cpus-per-task=1 --time=4:00:00 --mem=16GB --pty bash

## args
# V2_LIST=$1
# CTV3_LIST=$2
# GP_FILE=$3
# OUT=$4

## manual args
V2_LIST="AD_v2.list"
CTV3_LIST="AD_ctv3.v2.uniq.list"
GP_FILE="gp_clinical.txt"
OUT="gp-diag-ad-case-eids.txt"

dos2unix ${V2_LIST}
dos2unix ${CTV3_LIST}

awk -F"\t" 'F==1{a2[$1]=$1}F==2{a3[$1]=$1}F==3{if($4 in a2 || $5 in a3){print $1}}' \
 F=1 ${V2_LIST} F=2 ${CTV3_LIST} \
 F=3 ${GP_FILE} \
 > ${OUT}

sort -u ${OUT} > ${OUT}.uniq