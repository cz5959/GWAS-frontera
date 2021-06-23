#!/bin/sh

#linear association; covariates include 10 PC, sex, and birth year

# paths and variables
INPUT=$SCRATCH/QC
OUTPUT=$SCRATCH/GWAS_Results/$1
PHENO_NAME=$1
COL_NAMES=chrom,pos,ref,alt,ax,test,nobs,beta,se,tz,p
echo $1

# move to directory with pulled metadata files
cd $SCRATCH/Phenotypes

# PLINK2 GLM
for i in {1..22}
do 
	# both sex
	plink2 --memory 64000 --threads 16 --glm no-x-sex hide-covar cols=$COL_NAMES --pfile $INPUT/ukb_imp_chr${i}_v3_7 --keep whitebritishIDs.txt --pheno pheno_${PHENO_NAME}.txt --pheno-name $PHENO_NAME --out $OUTPUT/${PHENO_NAME}_both/both_sex_${i} --covar covariates.txt --covar-col-nums 3-14 --covar-variance-standardize 
	# females
	plink2 --memory 64000 --threads 16 --glm no-x-sex hide-covar cols=$COL_NAMES --pfile $INPUT/ukb_imp_chr${i}_v3_7 --keep whitebritishIDs.txt --keep-females --pheno pheno_${PHENO_NAME}.txt --pheno-name $PHENO_NAME --out $OUTPUT/${PHENO_NAME}_female/female_${i} --covar covariates.txt --covar-col-nums 3-12,14 --covar-variance-standardize
	# males
	plink2 --memory 64000 --threads 16 --glm no-x-sex hide-covar cols=$COL_NAMES --pfile $INPUT/ukb_imp_chr${i}_v3_7 --keep whitebritishIDs.txt --keep-males --pheno pheno_${PHENO_NAME}.txt --pheno-name $PHENO_NAME --out $OUTPUT/${PHENO_NAME}_male/male_${i} --covar covariates.txt --covar-col-nums 3-12,14 --covar-variance-standardize
done