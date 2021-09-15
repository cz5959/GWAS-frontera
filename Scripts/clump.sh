#!/bin/bash

FILE=$SCRATCH/GWAS_Results/$1
TARGET=$SCRATCH/1000G/EUR_all_phase3

plink \
    --bfile $TARGET \
    --clump-p1 1 \
    --clump-r2 0.1 \
    --clump-kb 250 \
    --clump $FILE/both_sex_all.${1}.glm.linear \
    --clump-snp-field ID \
    --clump-field P \
    --out $FILE/both_sex_$1

plink \
    --bfile $TARGET \
    --clump-p1 1 \
    --clump-r2 0.1 \
    --clump-kb 250 \
    --clump $FILE/female_all.${1}.glm.linear \
    --clump-snp-field ID \
    --clump-field P \
    --out $FILE/female_$1

plink \
    --bfile $TARGET \
    --clump-p1 1 \
    --clump-r2 0.1 \
    --clump-kb 250 \
    --clump $FILE/male_all.$1.glm.linear \
    --clump-snp-field ID \
    --clump-field P \
    --out $FILE/male_$1


# convert base file to plink1.9 format
# test.txt is file with . in it
#$SCRATCH/Scripts/plink2 --pfile EUR_all_phase3 --chr 1-22 --max-alleles 2 --exclude test.txt --rm-dup force-first --make-bed --out EUR_all_phase3

# change clumped results to tab delimited
sed 's/ \+/\t/g' $FILE/female_${1}.clumped | cut -f2-7 > $FILE/female_${1}_tab.clumped
sed 's/ \+/\t/g' $FILE/male_${1}.clumped | cut -f2-7 > $FILE/male_${1}_tab.clumped
