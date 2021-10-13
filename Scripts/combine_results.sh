#!/bin/bash

while getopts p: flag
do
    case "${flag}" in
        p) PHENO=${OPTARG};;
    esac
done
echo $PHENO

FILE_PATH=$SCRATCH/GWAS_Results/$PHENO
HEADER='#CHROM\tPOS\tID\tREF\tALT\tA1\tAX\tTEST\tOBS_CT\tBETA\tSE\tT_STAT\tP\n'
echo $PHENO

declare -a arr=("both_sex" "female" "male")
for sex in "${arr[@]}"
do
    # move to directory with linear result files
    cd $FILE_PATH/$sex
    # pipeline: concatenate results vertically | remove header | remove NA rows | add header 
    cat *.linear | awk '$1 !~ /^#CHROM/' - | awk '!/'NA'/' - | { printf $HEADER; cat - ; } > $PGS_dir/${sex}_all.${PHENO}.glm.linear
done




