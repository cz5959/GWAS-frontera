#!/bin/sh
# association - get covariates

# paths
meta_path=/corral-repl/utexas/Recombining-sex-chro/ukb/data/metadata

# go to directory with files
cd $SCRATCH/Association_Height_50

#get 10 PCs
cut -f10004-10013 $meta_path/ukb45020.txt > covar_10.txt

#get field for sex
awk -F "\t" '{print $9997}' $meta_path/ukb45020.txt > sex.txt

#get field for birth year
awk -F "\t" '{print $25}' $meta_path/ukb45020.txt > birth_year.txt

#get fields for genotype measurement batch, plate, and well
awk -F "\t" '{print $9996}' $meta_path/ukb45020.txt > geno_batch.txt
awk -F "\t" '{print $10002}' $meta_path/ukb45020.txt > geno_plate.txt
awk -F "\t" '{print $10003}' $meta_path/ukb45020.txt > geno_well.txt

#concatenate covariates: 10PCs, sex, birth year
paste -d "\t" ids.txt covar_10.txt sex.txt birth_year.txt geno_batch.txt geno_plate.txt geno_well.txt > covar.txt

#remove rows with missing values 
awk -F "\t" '{for(i=1;i<=NF;i++){if($i==""){next}}}1' covar.txt > covariates.txt



