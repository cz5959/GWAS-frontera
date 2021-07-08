#!/bin/sh

# file paths
IO=$SCRATCH/GWAS_Results/$1
echo $1 $2
cd $SCRATCH/ensembl-vep

./vep -i $IO/${1}_both_sex_topSNPs.txt --dir_cache $SCRATCH/Annotation/cache --cache -o $IO/${1}_${2}_annotation.txt --force_overwrite \
--tab --symbol --nearest symbol \
--pick \
--fields "Uploaded_variation,Location,Allele,SYMBOL,Gene,NEAREST,Feature,Feature_type,Consequence" --stats_text

##MSG: ERROR: --nearest requires Set::IntervalTree perl module to be installed
###export PERL5LIB=$PERL5LIB:$SCRATCH/Annotation/cpanm/lib/perl5
