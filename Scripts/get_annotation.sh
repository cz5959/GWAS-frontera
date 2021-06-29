#!/bin/sh

# file paths
IO=$SCRATCH/GWAS_Results/$1
echo $1 $2
cd $SCRATCH/ensembl-vep

./vep -i $IO/${1}_both_sex_topSNPs.txt --dir_cache $SCRATCH/Annotation/cache --cache -o $IO/${1}_${2}_annotation.txt --force_overwrite --tab --symbol --fields "Uploaded_variation,Location,Allele,SYMBOL,Gene,Feature,Feature_type,Consequence" --stats_text