# Phenotype extraction

To submit extraction of GP diagnosed AD cases:

1. 

``` bash
SCR_DIR="/user/work/tb13101/ad-severity-gwas/ukb/pheno-extraction" # scratch directory - i.e. /user/work/...
V2_LIST=${SCR_DIR}/data/ukb-pheno/AD_v2.list
CTV3_LIST=${SCR_DIR}/data/ukb-pheno/AD_ctv3.v2.uniq.list
GP_FILE=${SCR_DIR}/data/ukb-pheno/gp_clinical.txt
OUT=${SCR_DIR}/data/gp-diag-ad-case-eids.txt
bash scripts/extract-gp-diag-ad.sh "${V2_LIST}" "${CTV3_LIST}" "${GP_FILE}" "${OUT}"
```

1. get-ukb-vars.sh
2. extract-pheno-data.R
3. simplify-gp-scripts.sh