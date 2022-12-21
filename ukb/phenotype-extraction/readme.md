# Phenotype extraction

## Script order

1. [`extract-gp-diag-ad.sh`](scripts/extract-gp-diag-ad.sh)
2. [`get-ukb-vars.sh`](scripts/get-ukb-vars.sh)
3. [`extract-ad-cases.R`](scripts/extract-ad-cases.R)
4. [`simplify-gp-scripts.sh`](scripts/simplify-gp-scripts.sh)
	+ Note: Before moving onto the next step, the GP prescription data will need to be manually sorted to extract prescriptions of systemics, corticosteroids, and emollients
5. [`extract-severe-cases.R`](scripts/extract-severe-cases.R)