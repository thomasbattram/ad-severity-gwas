#!/bin/bash

#SBATCH -p mrcieu,mrcieu2,hmem
#SBATCH --job-name bolt_gwas_ad
#SBATCH -o chr_all_run.log
#SBATCH -e chr_all_run.err
#SBATCH --nodes=1 --tasks-per-node=14
#SBATCH --mem-per-cpu=4000
#SBATCH --time=5-00:00:00

cd $SLURM_SUBMIT_DIR

/mnt/storage/private/mrcieu/research/UKBIOBANK_GWAS_Pipeline/scripts/software/BOLT-LMM_v2.3.2/bolt \
    --bfile=/mnt/storage/private/mrcieu/research/UKBIOBANK_GWAS_Pipeline/data/bolt_bfile/grm6_european_filtered_ieu \
    --bgenFile=/mnt/storage/private/mrcieu/data/ukbiobank/genetic/variants/arrays/imputed/released/2018-09-18/data/dosage_bgen/data.chr0{1..9}.bgen \
    --bgenFile=/mnt/storage/private/mrcieu/data/ukbiobank/genetic/variants/arrays/imputed/released/2018-09-18/data/dosage_bgen/data.chr{10..22}.bgen \
    --bgenFile=/mnt/storage/private/mrcieu/data/ukbiobank/genetic/variants/arrays/imputed/released/2018-09-18/data/dosage_bgen/data.chrX.bgen \
    --sampleFile=/mnt/storage/private/mrcieu/data/ukbiobank/genetic/variants/arrays/imputed/released/2018-09-18/data/dosage_bgen/data.chr1-22.sample \
    --geneticMapFile=/mnt/storage/private/mrcieu/research/UKBIOBANK_GWAS_Pipeline/scripts/software/BOLT-LMM_v2.3.2/tables/genetic_map_hg19_withX.txt.gz \
    --bgenMinMAF=0.001 \
    --phenoFile=sample_file_prevalent_ifg_ukb_fasted.txt \
    --phenoCol=IFG \
    --covarFile=sample_file_prevalent_ifg_ukb_fasted.txt \
    --qCovarCol=assessment_age \
    --qCovarCol=PC{1:10} \
    --covarCol=sex \
    --lmm \
    --LDscoresFile=/mnt/storage/private/mrcieu/research/UKBIOBANK_GWAS_Pipeline/scripts/software/BOLT-LMM_v2.3.2/tables/LDSCORE.1000G_EUR.tab.gz \
    --LDscoresMatchBp \
    --numThreads=14 \
    --verboseStats \
    --modelSnps=/mnt/storage/private/mrcieu/research/UKBIOBANK_GWAS_Pipeline/data/model_snps_for_grm/grm6_snps.prune.in \
    --statsFileBgenSnps=IFG_fasted_imputed.txt.gz \
    --covarMaxLevels=30 \
    --statsFile=IFG_fasted_out.txt.gz


## separate whitespace-delimited file (specified with --phenoFile) with the first line
## containing column headers and subsequent lines containing records, one per individual. 
## The first two columns must be FID and IID (the PLINK identifiers of an individual). 
## Any number of columns may follow; the column containing the phenotype to analyze is 
## specified with --phenoCol. Values of -9 and NA are interpreted as missing data. 
## All other values in the column should be numeric. 
## The records in lines following the header line need not be in sorted order and need 
## not match the individuals in the genotype data (i.e., fam file)