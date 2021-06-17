#!/bin/bash

# output path
# passed in parameter should be phenotype name
FILE_PATH=$SCRATCH/GWAS_Results
HEADER='#CHROM\tPOS\tID\tREF\tALT\tA1\tAX\tTEST\tOBS_CT\tBETA\tSE\tT_STAT\tP\n'
echo $1

# move to directory with linear result files
# pipeline: concatenate results vertically | remove header | remove NA rows | add header 
cd $FILE_PATH/${1}_both
cat *.linear | awk '$1 !~ /^#CHROM/' - | awk '!/'NA'/' - | { printf $HEADER; cat - ; } > $FILE_PATH/both_sex_all.height.glm.linear

cd $FILE_PATH/${1}_female
cat *.linear | awk '$1 !~ /^#CHROM/' - | awk '!/'NA'/' - | { printf $HEADER; cat - ; } > $FILE_PATH/female_all.height.glm.linear

cd $FILE_PATH/${1}_male
cat *.linear | awk '$1 !~ /^#CHROM/' - | awk '!/'NA'/' - | { printf $HEADER; cat - ; } > $FILE_PATH/male_all.height.glm.linear