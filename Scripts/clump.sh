#!/bin/bash

while getopts p: flag
do
    case "${flag}" in
        p) PHENO=${OPTARG};;
    esac
done
echo $PHENO 

FILE=$SCRATCH/GWAS_Results/$PHENO
LD=$SCRATCH/1000G/EUR_all_phase3

declare -a arr=("both_sex" "female" "male")
for sex in "${arr[@]}"
do
plink \
    --bfile $LD \
    --clump-p1 1 \
    --clump-r2 0.1 \
    --clump-kb 250 \
    --clump $FILE/${sex}_all.${PHENO}.glm.linear \
    --clump-snp-field ID \
    --clump-field P \
    --out $FILE/${sex}_$PHENO
done





# convert base file to plink1.9 format
# test.txt is file with . in it
#$SCRATCH/Scripts/plink2 --pfile EUR_all_phase3 --chr 1-22 --max-alleles 2 --exclude test.txt --rm-dup force-first --make-bed --out EUR_all_phase3

# only keep unrelated (remove 3rd degree and closer) and europeans  #pop.txt had GBR and CEU codes 
grep -h -f pop.txt EUR_all_phase3.psam | cut -f1 - > test.txt
$SCRATCH/Scripts/plink2 --pfile EUR_all_phase3 --chr 1-22 --max-alleles 2 --keep test.txt --rm-dup exclude-all --king-cutoff 0.0442 --make-bed --out EUR_all_phase3

# change clumped results to tab delimited
sed 's/ \+/\t/g' $FILE/female_${1}.clumped | cut -f2-7 > $FILE/female_${1}_tab.clumped
sed 's/ \+/\t/g' $FILE/male_${1}.clumped | cut -f2-7 > $FILE/male_${1}_tab.clumped


#### 
PHENO="arm_fatfree_mass_L"
sex="female"
FILE=$SCRATCH/GWAS_Results/$PHENO
LD=$SCRATCH/1000G/all_phase3
plink --bfile $LD --clump-p1 5e-8 --clump-r2 0.1 --clump-kb 250 --clump $FILE/${sex}_all.${PHENO}.glm.linear --clump-snp-field ID --clump-field P --out $FILE/${sex}_$PHENO