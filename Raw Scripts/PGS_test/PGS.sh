#!/bin/sh

# load post QC summary statistic
BASE=$SCRATCH/GWAS_Results/$1
mkdir -p $BASE/PRS
# load target data
TARGET=$SCRATCH/1000G

# remove duplicates in base
awk '{seen[$3]++; if(seen[$3]==1){ print}}' $BASE/female_all.height.glm.linear > $BASE/female_all_nodup.height.glm.linear
awk '{seen[$3]++; if(seen[$3]==1){ print}}' $BASE/female_all.height.glm.linear > $BASE/female_all_nodup.height.glm.linear

# clumping  for base data
plink \
    --bfile $TARGET \
    --clump-p1 1 \
    --clump-r2 0.1 \
    --clump-kb 250 \
    --clump $BASE/female_all.${1}.glm.linear \
    --clump-snp-field ID \
    --clump-field P \
    --out $BASE/PRS/female_${1}_prs

plink \
    --bfile $TARGET \
    --clump-p1 1 \
    --clump-r2 0.1 \
    --clump-kb 250 \
    --clump $BASE/male_all.${1}.glm.linear \
    --clump-snp-field ID \
    --clump-field P \
    --out $BASE/PRS/male_${1}_prs


# extract SNP id from clumping results
awk 'NR!=1{print $3}' $BASE/PRS/female_${1}_prs.clumped >  $BASE/PRS/female_${1}_prs.valid.snp
awk 'NR!=1{print $3}' $BASE/PRS/male_${1}_prs.clumped >  $BASE/PRS/male_${1}_prs.valid.snp

################### extract SNP id and p-value from results
awk '{print $3,$13}' $BASE/PRS/female_all.${1}.glm.linear > $BASE/PRS/female_${1}_SNP.pvalue
awk '{print $3,$13}' $BASE/PRS/male_all.${1}.glm.linear> $BASE/PRS/male_${1}_SNP.pvalue

# calculate PRS
## 3:SNP    6:effective allele (A1)  12: effect size
plink \
    --bfile $TARGET \
    --score $BASE/female_all.${1}.glm.linear 3 6 10 header \
    --q-score-range range_list $BASE/PRS/female_${1}_SNP.pvalue \
    --extract $BASE/PRS/female_${1}_prs.valid.snp \
    --out $BASE/PRS/EUR_female_${1}_prs
plink \
    --bfile $TARGET \
    --score $BASE/male_all.${1}.glm.linear 3 6 10 header \
    --q-score-range range_list $BASE/PRS/male_${1}_SNP.pvalue \
    --extract $BASE/PRS/male_${1}_prs.valid.snp \
    --out $BASE/PRS/EUR_male_${1}_prs


# Population Stratification
# prunning
#plink --bfile ${TARGET}.QC --indep-pairwise 200 50 0.25 --out $TARGET
# calculate the first 6 PCs
#plink --bfile ${TARGET}.QC --extract ${TARGET}.prune.in --pca 6 --out $TARGET

