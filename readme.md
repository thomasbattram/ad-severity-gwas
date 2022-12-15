# Atopic dermatitis severity GWAS

GWAS will be run across multiple cohorts and these results will be meta-analysed together. Each cohort will have a folder and within those folders, Snakemake pipelines will be setup to:

1. QC the genotype data
2. Extract and QC phenotype data
3. Run the GWAS 

Scripts that might be common to all cohorts (e.g. GWAS script) will be found in [`common-scripts`]().

