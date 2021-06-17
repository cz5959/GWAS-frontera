#!/bin/sh

# file paths
SUM_STATS=$SCRATCH/Association_Height_50
OUTPUT=$SCRATCH/LD_practice
LD_SCORE=$SCRATCH/LD_practice/LD_scores
NEALE=$SCRATCH/Neale_Lab/Standing_Height

# set up environment
cd /scratch1/08005/cz5959/ldsc 
source activate ldsc

#### reformat summary stat
# need rsid; effect allele; non-effect allele; sample size; p-value; signed summ stat (beta or log odds, z-score, etc)
# plink
munge_sumstats.py --sumstats $SUM_STATS/linear_results_all_chrom.height.glm.linear --snp ID --a2 AX --N-col OBS_CT --out $OUTPUT/height_results --merge-alleles $LD_SCORE/w_hm3.snplist
# neale lab
munge_sumstats.py --sumstats $NEALE/neale_height_tab_ldsc.txt --out $OUTPUT/neale_height --merge-alleles $LD_SCORE/w_hm3.snplist

#### LD Score Regression
# compute heritability and ldsc intercept
# plink
ldsc.py --h2 $OUTPUT/height_results.sumstats.gz --w-ld-chr $LD_SCORE/eur_w_ld_chr/ --ref-ld-chr $LD_SCORE/eur_w_ld_chr/ --out $OUTPUT/height_h2
# neale
ldsc.py --h2 $OUTPUT/neale_height.sumstats.gz --w-ld-chr $LD_SCORE/eur_w_ld_chr/ --ref-ld-chr $LD_SCORE/eur_w_ld_chr/ --out $OUTPUT/neale_height_h2