# GWAS scripts

Requires having run the scripts in the [phenotype-extraction](../phenotype-extraction) folder.

The scripts can be executed manually or via snakemake. Here is the snakemake workflow:

0. Start a tmux session - `tmux new -s severity-gwas`
1. Activate conda env - `conda activate /user/home/tb13101/conda-envs/envs/snakemake`
2. Edit the Snakefile template and name it "Snakefile"
3. Do a dry run with `snakemake -nrp`
4. Submit pipeline as a job:

``` bash
module add tools/pandoc/2.19.2
snakemake -rp \
-j 1 \
--cluster-config bc4-cluster.json \
--cluster "sbatch \
        --account=smed001801
        --job-name={cluster.name} \
        --partition={cluster.partition} \
        --nodes={cluster.nodes} \
        --ntasks-per-node={cluster.ntask} \
        --cpus-per-task={cluster.ncpu} \
        --time={cluster.time} \
        --mem={cluster.mem} \
        --output={cluster.output} \
        --error={cluster.error}"

## snakemake -rp --cores 1
```

5. Deactivate the tmux session (CTRL+b + d)
6. Kill the tmux session when jobs are complete - `tmux kill-session -t ad-cse`

## Script order

The scripts run by the snakemake workflow are below:

1. [extract-exclusions.R](scripts/extract-exclusions.R): extract individuals who can be excluded from the GWAS
2. [prep-phenofile.R](scripts/prep-phenofile.R): prepare the phenotype data for bolt-LMM
3. [prep-covarfile.R](scripts/prep-covarfile.R): prepare the covariate data for bolt-LMM
4. [bolt-script.sh](scripts/bolt-script.sh): run the gwas in ukb using bolt-LMM


To run the plink analyses, use the code below:
``` bash
sbatch plink-script-sub.sh "/user/work/tb13101/ad-severity-gwas/ukb/gwas" \
                                     "/mnt/storage/private/mrcieu/data/ukbiobank/genetic/variants/arrays/imputed/released/2018-09-18/data/dosage_bgen" \
                                     "data/severe-ad-phenofile.txt" \
                                     "data/severe-ad-covarfile.txt" \
                                     "results/gwas-tmp/severe-gwas-out"
```


To run a test of chromosome 1, use the code below:
``` bash
sbatch plink-script-sub-chr1-only.sh "/user/work/tb13101/ad-severity-gwas/ukb/gwas" \
                                     "/mnt/storage/private/mrcieu/data/ukbiobank/genetic/variants/arrays/imputed/released/2018-09-18/data/dosage_bgen" \
                                     "data/severe-ad-phenofile.txt" \
                                     "data/severe-ad-covarfile.txt" \
                                     "results/gwas-temp/severe-gwas-out"
```