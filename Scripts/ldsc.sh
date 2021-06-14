#!/bin/sh

# file paths
SUM_STATS=/scratch1/08005/cz5959/Association_Height_50
OUTPUT=/scratch1/08005/cz5959/LD_practice
LD_SCORE=/scratch1/08005/cz5959/LD_practice/LD_scores

# set up environment
cd /scratch1/08005/cz5959/ldsc 
source activate ldsc

#### reformat summary stat
# plink
# need rsid; effect allele; non-effect allele; sample size; p-value; signed summ stat (beta or log odds, z-score, etc)
munge_sumstats.py --sumstats $SUM_STATS/linear_results_all_chrom.height.glm.linear --snp ID --a2 REF --N-col OBS_CT --out $OUTPUT/height_results --merge-alleles $LD_SCORE/w_hm3.snplist

# neale lab


#### LD Score Regression
# compute heritability and ldsc intercept
ldsc.py --h2 $OUTPUT/height_results.sumstats.gz --w-ld-chr $LD_SCORE/eur_w_ld_chr/ --ref-ld-chr $LD_SCORE/eur_w_ld_chr/ --out $OUTPUT/height_h2

