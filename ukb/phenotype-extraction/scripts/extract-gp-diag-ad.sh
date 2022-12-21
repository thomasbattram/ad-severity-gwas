#!/bin/bash

#SBATCH --job-name=extract-gp-dat
#SBATCH --partition=veryshort
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=4:00:00 
#SBATCH --mem=16GB

# ---------------------------------------------------------------------------
# Extract UK Biobank AD cases from GP records
# ---------------------------------------------------------------------------

## This script extracts all UK Biobank IDs that have read codes corresponding to AD as 
## diagnosed by their GP

## args
# V2_LIST=$1
# CTV3_LIST=$2
# GP_FILE=$3
# OUT=$4

## add manual args here!

dos2unix ${V2_LIST}
dos2unix ${CTV3_LIST}

awk -F"\t" 'F==1{a2[$1]=$1}F==2{a3[$1]=$1}F==3{if($4 in a2 || $5 in a3){print $1}}' \
 F=1 ${V2_LIST} F=2 ${CTV3_LIST} \
 F=3 ${GP_FILE} \
 > ${OUT}

sort -u ${OUT} > ${OUT}.uniq