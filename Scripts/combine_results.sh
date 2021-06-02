#!/bin/bash

# association plots - linear results 

#output path
output_path=$SCRATCH/Association_Height_50

#move to directory with linear result files
cd $SCRATCH/Association_Height_50/Results

# remove NA for chr 1
#awk '!/'NA'/' linear_results_1.height.glm.linear > ./Results/linear_results_1_2.height.glm.linear

cat *.linear | awk '$1 !~ /^#CHROM/' - | awk '!/'NA'/' - | cut -f1-12 - | { printf '#CHROM\tPOS\tID\tREF\tALT\tA1\tTEST\tOBS_CT\tBETA\tSE\tT_STAT\tP\n'; cat - ; } > $output_path/linear_results_all_chrom.height.glm.linear


#for i in {1..22}
#do
    # remove NA and header pipeline and error column  
#    awk '!/'NA'/' ./Result/linear_result_${i}.height.glm.linear | tail -n+2 - | cut -f1-12 - > ./Results_2/linear_result_${i}_2.height.glm.linear

#done

# move to results directory
#cd Results_2

# concatenate all the chromosomes into one file; add header
#cat * | { printf '#CHROM\tPOS\tID\tREF\tALT\tA1\tTEST\tOBS_CT\tBETA\tSE\tT_STAT\tP\n'; cat -; } > $output_path/linear_results_all_chrom.height.glm.linear
