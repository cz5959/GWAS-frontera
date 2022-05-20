#!/bin/sh

while getopts p:t:s: flag
do
    case "${flag}" in
        p) PHENO=${OPTARG};;
        t) TYPE=${OPTARG};;
        s) SET=${OPTARG};;
    esac
done
echo $PHENO; echo $TYPE; echo $SET

LD_dir=$SCRATCH/1000G/EUR_all_phase3
PGS_dir=$SCRATCH/GWAS_Results/$PHENO/PGS_$SET

### ADDITIVE ###
declare -a arr=("both_sex" "female" "male")
for sex in "${arr[@]}"
do
    # clumping of base data
    plink \
        --bfile $LD_dir \
        --clump-p1 1 \
        --clump-r2 0.1 \
        --clump-kb 250 \
        --clump $PGS_dir/${sex}_train.${PHENO}.glm.${TYPE} \
        --clump-snp-field ID \
        --clump-field P \
        --out $PGS_dir/${sex}_additive.${PHENO}
    
    # extract SNP id from clumping results
    awk 'NR!=1{print $3}' $PGS_dir/${sex}_additive.${PHENO}.clumped >  $PGS_dir/${sex}_additive.${PHENO}.clumped.snpid
    # get snp and p-values from training set
    awk '{print $3,$13}' $PGS_dir/${sex}_train.${PHENO}.glm.${TYPE} > $PGS_dir/${sex}_additive.${PHENO}.snp_pvalue
    
    # calculate PRS
    ## 3:SNP    6:effective allele (A1)  12: effect size
    # score is a sum across SNPs of the # of ref alleles (0,1,or 2) at that SNP multiplied by the score for that SNP
    plink \
        --bfile $PGS_dir/${PHENO}_test \
        --score $PGS_dir/${sex}_train.${PHENO}.glm.${TYPE} 3 6 10 header \
        --q-score-range range_list $PGS_dir/${sex}_additive.${PHENO}.snp_pvalue \
        --extract $PGS_dir/${sex}_additive.${PHENO}.clumped.snpid \
        --out $PGS_dir/${sex}_additive_${PHENO}
done

### mash ###
declare -a arr=("female" "male")
for sex in "${arr[@]}"
do
    plink \
        --bfile $LD_dir \
        --clump-p1 1 \
        --clump-r2 0.1 \
        --clump-kb 250 \
        --clump $PGS_dir/${sex}_pseudoP_pgs.${PHENO}.txt \
        --clump-snp-field ID \
        --clump-field P \
        --out $PGS_dir/${sex}_mash.${PHENO}
    
    awk 'NR!=1{print $3}' $PGS_dir/${sex}_mash.${PHENO}.clumped >  $PGS_dir/${sex}_mash.${PHENO}.clumped.snpid
    awk '{print $1,$5}' $PGS_dir/${sex}_pseudoP_pgs.${PHENO}.txt > $PGS_dir/${sex}_mash.${PHENO}.snp_pvalue

    plink \
        --bfile $PGS_dir/${PHENO}_test \
        --score $PGS_dir/${sex}_pseudoP_pgs.${PHENO}.txt 1 2 3 header \
        --q-score-range range_list $PGS_dir/${sex}_mash.${PHENO}.snp_pvalue \
        --extract $PGS_dir/${sex}_mash.${PHENO}.clumped.snpid \
        --out $PGS_dir/${sex}_mash_${PHENO}
done

# # thresholds to use
# echo "1 0 1" > range_list 
# echo "0.01 0 0.01" >> range_list
# echo "1e-5 0 1e-5" >> range_list
# echo "1e-8 0 1e-8" >> range_list
