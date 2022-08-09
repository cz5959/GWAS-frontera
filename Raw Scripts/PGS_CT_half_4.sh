#!/bin/sh

while getopts p:s: flag
do
    case "${flag}" in
        p) PHENO=${OPTARG};;
        s) SET=${OPTARG};;
    esac
done
echo $PHENO; echo $SET

LD_dir=$SCRATCH/1000G/all_phase3
PGS_dir=$SCRATCH/GWAS_Results/$PHENO/PGS_$SET

### ADDITIVE ###
# clumping of base data
plink \
    --bfile $LD_dir \
    --clump-p1 1 \
    --clump-r2 0.1 \
    --clump-kb 250 \
    --clump $PGS_dir/both_sex_train.${PHENO}.glm.linear \
    --clump-snp-field ID \
    --clump-field P \
    --out $PGS_dir/both_sex_additive.${PHENO}

# extract SNP id from clumping results
awk 'NR!=1{print $3}' $PGS_dir/both_sex_additive.${PHENO}.clumped >  $PGS_dir/both_sex_additive.${PHENO}.clumped.snpid
# get snp and p-values from training set
awk '{print $3,$13}' $PGS_dir/both_sex_train.${PHENO}.glm.linear > $PGS_dir/both_sex_additive.${PHENO}.snp_pvalue

# calculate PRS
## 3:SNP    6:effective allele (A1)  12: effect size
# score is a sum across SNPs of the # of ref alleles (0,1,or 2) at that SNP multiplied by the score for that SNP
plink \
    --bfile $PGS_dir/${PHENO}_test \
    --score $PGS_dir/both_sex_train.${PHENO}.glm.linear 3 6 10 header \
    --q-score-range range_list $PGS_dir/both_sex_additive.${PHENO}.snp_pvalue \
    --extract $PGS_dir/both_sex_additive.${PHENO}.clumped.snpid \
    --out $PGS_dir/both_sex_additive_${PHENO}
