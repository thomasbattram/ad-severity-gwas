# Phenotype extraction

There is a lot of looking at the data and manual editing required, so I've decided against setting up a pipeline to run the analyses and instead decided to run the scripts manually.

## Data required

* GP prescriptions
* Self-reported prescriptions
* Phototherapy treatments - Operative procedures (OPCS4) 
* Hospital Episode Statistics
* Self-reported disease

## Script order

1. [`extract-gp-diag-ad.sh`](scripts/extract-gp-diag-ad.sh)
2. [`get-ukb-vars.sh`](scripts/get-ukb-vars.sh)
3. [`extract-ad-cases.R`](scripts/extract-ad-cases.R)
4. [`simplify-gp-scripts.sh`](scripts/simplify-gp-scripts.sh)
5. [`simplify-sr-treatments.R`](scripts/simplify-sr-treatments.R)
	+ Note: Before moving onto the next step, the GP prescription data AND the self-reported prescription data will need to be manually sorted to extract prescriptions of systemics, corticosteroids, and emollients
6. [`clean-gp-data.R`](scripts/clean-gp-data.R)
7. [`clean-sr-data.R`](scripts/clean-sr-data.R)
8. [`extract-severe-cases.R`](scripts/extract-severe-cases.R)