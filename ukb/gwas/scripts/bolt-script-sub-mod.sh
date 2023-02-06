#!/bin/bash

#SBATCH -p mrcieu,mrcieu2,hmem
#SBATCH --job-name bolt_gwas_ad
#SBATCH -o chr_all_run.log
#SBATCH -e chr_all_run.err
#SBATCH --nodes=1 --tasks-per-node=14
#SBATCH --mem-per-cpu=4000
#SBATCH --time=5-00:00:00

SLURM_SUBMIT_DIR=""
cd $SLURM_SUBMIT_DIR

## args
# UKB_GEN_DIR=$1
# UKB_PIPELINE_DIR=$2
# PHENO_FILE=$3
# COVAR_FILE=$4
# OUTFILE_BGEN=$5
# OUTFILE=$6

## manual args
UKB_GEN_DIR="" # UK biobank genetic data directory
UKB_PIPELINE_DIR="" # UK biobank GWAS pipeline directory
PHENO_FILE="data/mod_severe-ad-phenofile.txt"
COVAR_FILE="data/severe-ad-covarfile.txt"
OUTFILE_BGEN="results/mod_severe-imputed.txt.gz"
OUTFILE="results/mod_severe-out.txt.gz"

${UKB_PIPELINE_DIR}/scripts/software/BOLT-LMM_v2.3.2/bolt \
    --bfile=${UKB_PIPELINE_DIR}/data/bolt_bfile/grm6_european_filtered_ieu \
    --bgenFile=${UKB_GEN_DIR}/dosage_bgen/data.chr0{1..9}.bgen \
    --bgenFile=${UKB_GEN_DIR}/dosage_bgen/data.chr{10..22}.bgen \
    --bgenFile=${UKB_GEN_DIR}/dosage_bgen/data.chrX.bgen \
    --sampleFile=${UKB_GEN_DIR}/dosage_bgen/data.chr1-22.sample \
    --geneticMapFile=${UKB_PIPELINE_DIR}/scripts/software/BOLT-LMM_v2.3.2/tables/genetic_map_hg19_withX.txt.gz \
    --bgenMinMAF=0.001 \
    --phenoFile=${PHENO_FILE} \
    --phenoCol=severe \
    --covarFile=${COVAR_FILE} \
    --qCovarCol=age \
    --qCovarCol=PC{1:10} \
    --covarCol=sex \
    --lmm \
    --LDscoresFile=${UKB_PIPELINE_DIR}/scripts/software/BOLT-LMM_v2.3.2/tables/LDSCORE.1000G_EUR.tab.gz \
    --LDscoresMatchBp \
    --numThreads=14 \
    --verboseStats \
    --modelSnps=${UKB_PIPELINE_DIR}/data/model_snps_for_grm/grm6_snps.prune.in \
    --statsFileBgenSnps=${OUTFILE_BGEN} \
    --covarMaxLevels=30 \
    --statsFile=${OUTFILE}

## separate whitespace-delimited file (specified with --phenoFile) with the first line
## containing column headers and subsequent lines containing records, one per individual. 
## The first two columns must be FID and IID (the PLINK identifiers of an individual). 
## Any number of columns may follow; the column containing the phenotype to analyze is 
## specified with --phenoCol. Values of -9 and NA are interpreted as missing data. 
## All other values in the column should be numeric. 
## The records in lines following the header line need not be in sorted order and need 
## not match the individuals in the genotype data (i.e., fam file)