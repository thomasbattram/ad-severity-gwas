#!/bin/bash -l

#SBATCH -p mrcieu,hmem
#SBATCH --account=smed001801
#SBATCH --job-name plink_gwas_ad
#SBATCH -o chr_%a_run.log
#SBATCH -e chr_%a_run.err
#SBATCH --ntasks=4
#SBATCH --nodes=1
#SBATCH --mem=60G
#SBATCH --time=2-00:00:00
#SBATCH --array=1

module add apps/plink/2.00

if [ "$#" -ne 5 ]; then
    echo "   ########  Five command line arguments required:"
    echo "   ########  1. Working directory"
    echo "   ########  2. Directory with UK Biobank genetic data in it"
    echo "   ########  3. Input phenotype file"
    echo "   ########  4. Input covariate file"
    echo "   ########  5. Output stem"
    exit
fi

echo "Command line arguments accepted"

WORKING_DIR=$1
UKB_GEN_DIR=$2
PHENO=$3
COVS=$4
OUTPUT=$5

## Manual arg
# WORKING_DIR=""
# UKB_GEN_DIR=""
# PHENO="data/severe-ad-phenofile.txt"
# COVS="data/severe-ad-covarfile.txt"
# OUTPUT="results/gwas-tmp/severe-gwas-out"
# SLURM_ARRAY_TASK_ID=1

cd ${WORKING_DIR}

CHR=$SLURM_ARRAY_TASK_ID

OUTPUT="${OUTPUT}-chr${CHR}"

if [[ ${CHR} -eq 1 ]]; then
	echo Command line arguments for GWAS script > ${OUTPUT}_script_commandargs
	echo 1 ${WORKING_DIR} >> ${OUTPUT}_script_commandargs
	echo 2 ${UKB_GEN_DIR} >> ${OUTPUT}_script_commandargs
	echo 3 ${PHENO} >> ${OUTPUT}_script_commandargs
	echo 4 ${COVS} >> ${OUTPUT}_script_commandargs
	echo 5 ${OUTPUT} >> ${OUTPUT}_script_commandargs
	echo >> ${OUTPUT}_script_commandargs
	echo Built in arguments for script >> ${OUTPUT}_script_commandargs
	echo Covariates ${COVS} >> ${OUTPUT}_script_commandargs
fi

if [[ $CHR -eq 23 ]]
then
	SAMPLEFILE=${UKB_GEN_DIR}/data.chrX_plink.sample
	CHR="X"
else
	SAMPLEFILE=${UKB_GEN_DIR}/data.chr1-22_plink.sample
	if [[ ${CHR} -lt 10 ]]; then
		CHR="0${CHR}"
	fi
fi

echo "PLINK processing started"

plink \
  --bgen ${UKB_GEN_DIR}/data.chr${CHR}.bgen \
  --sample $SAMPLEFILE \
  --maf 0.01 \
  --pheno ${PHENO} \
  --pheno-name severe \
  --covar ${COVS} \
  --glm hide-covar \
  --ci 0.95 \
  --out ${OUTPUT}
