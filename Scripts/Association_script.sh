#!/bin/sh
# association

# height, pre-computed PC, white british
# paths
input_path=$SCRATCH/QC
output_path=$SCRATCH/Association_Height_50/Results

#height linear association
#covariates include 10 PC, sex, and birth year

# move to directory with pulled meta data files
cd $SCRATCH/Association_Height_50

# PLINK2 GLM
for i in {1..22}
do 
	plink2 --memory 64000 --threads 16 --glm hide-covar --pfile $input_path/ukb_imp_chr${i}_v3_5 --keep whitebritishIDs.txt --pheno height_pheno.txt --pheno-name height --out $output_path/linear_result_${i} --covar covariates.txt --covar-col-nums 3-14,17 --covar-variance-standardize 
done

# single chromosome - 22
#plink2 --memory 64000 --threads 16 --glm hide-covar --pfile $input_path/ukb_imp_chr22_v3_5 --keep whitebritishIDs.txt --pheno height_pheno.txt --pheno-name height --out $output_path/linear_result_22 --covar covariates.txt --covar-col-nums 3-14 --covar-variance-standardize
