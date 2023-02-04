#!/bin/bash

## args
UKB_GEN_DIR=$1
UKB_PIPELINE_DIR=$2
PHENO_FILE=$3
COVAR_FILE=$4
OUTFILE_BGEN=$5
OUTFILE=$6

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
