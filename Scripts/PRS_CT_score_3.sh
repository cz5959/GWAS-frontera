#!/bin/sh


PHENO=$1
echo $1
LD_dir=$SCRATCH/1000G/EUR_all_phase3
PRS_dir=$SCRATCH/GWAS_Results/$PHENO/PRS

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
        --clump $PRS_dir/${sex}_train.${PHENO}.glm.linear \
        --clump-snp-field ID \
        --clump-field P \
        --out $PRS_dir/${sex}_additive.${PHENO}
    
    # extract SNP id from clumping results
    awk 'NR!=1{print $3}' $PRS_dir/${sex}_additive.${PHENO}.clumped >  $PRS_dir/${sex}_additive.${PHENO}.clumped.snpid
    # get snp and p-values from training set
    awk '{print $3,$13}' $PRS_dir/${sex}_train.${PHENO}.glm.linear > $PRS_dir/${sex}_additive.${PHENO}.snp_pvalue
    
    # calculate PRS
    ## 3:SNP    6:effective allele (A1)  12: effect size
    # score is a sum across SNPs of the # of ref alleles (0,1,or 2) at that SNP multiplied by the score for that SNP
    plink \
        --bfile $PRS_dir/${PHENO}_test \
        --score $PRS_dir/${sex}_train.${PHENO}.glm.linear 3 6 10 header \
        --q-score-range range_list $PRS_dir/${sex}_additive.${PHENO}.snp_pvalue \
        --extract $PRS_dir/${sex}_additive.${PHENO}.clumped.snpid \
        --out $PRS_dir/${sex}_additive_${PHENO}
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
        --clump $PRS_dir/${sex}_pseudoP.${PHENO}.txt \
        --clump-snp-field ID \
        --clump-field P \
        --out $PRS_dir/${sex}_mash.${PHENO}
    
    awk 'NR!=1{print $3}' $PRS_dir/${sex}_mash.${PHENO}.clumped >  $PRS_dir/${sex}_mash.${PHENO}.clumped.snpid
    awk '{print $1,$5}' $PRS_dir/${sex}_pseudoP.${PHENO}.txt > $PRS_dir/${sex}_mash.${PHENO}.snp_pvalue

    plink \
        --bfile $PRS_dir/${PHENO}_test \
        --score $PRS_dir/${sex}_pseudoP.${PHENO}.txt 1 2 3 header \
        --q-score-range range_list $PRS_dir/${sex}_mash.${PHENO}.snp_pvalue \
        --extract $PRS_dir/${sex}_mash.${PHENO}.clumped.snpid \
        --out $PRS_dir/${sex}_mash_${PHENO}
done

# # thresholds to use
# echo "0.001 0 0.001" > range_list 
# echo "0.05 0 0.05" >> range_list
# echo "0.1 0 0.1" >> range_list
# echo "0.2 0 0.2" >> range_list
# echo "0.3 0 0.3" >> range_list
# echo "0.4 0 0.4" >> range_list
# echo "0.5 0 0.5" >> range_list