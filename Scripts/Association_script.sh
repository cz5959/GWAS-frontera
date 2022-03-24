#!/bin/sh

# phenotype command line argument
while getopts p: flag
do
    case "${flag}" in
        p) PHENO=${OPTARG};;
    esac
done
echo $PHENO 

# linear association; covariates include 10 PC, sex, and birth year
#make directories if they don't exist
cd $SCRATCH/GWAS_Results
mkdir -p $PHENO/{both_sex,female,male}
# paths and variables
QC_DIR=$SCRATCH/QC_Chr
PHENO_DIR=$SCRATCH/Phenotypes
OUTPUT_DIR=$SCRATCH/GWAS_Results/$PHENO
COL_NAMES=chrom,pos,ref,alt,ax,test,nobs,beta,se,tz,p

# move to directory with pulled metadata files
cd $SCRATCH/Phenotypes

# PLINK2 GLM
declare -a arr=("both_sex" "female" "male")
for i in {1..22}
do
	plink2 --memory 64000 --threads 16 --glm no-x-sex hide-covar cols=$COL_NAMES \
		--pfile $QC_DIR/ukb_imp_chr${i}_v3_11 \
		--pheno ${PHENO_DIR}/pheno_${PHENO}.txt --pheno-name $PHENO \
		--covar ${PHENO_DIR}/covariates.txt --covar-col-nums 3-14 --covar-variance-standardize \
		--out $OUTPUT_DIR/both_sex/both_sex_${i}
	for sex in "${arr[@]}"
	do
		# both sex
		plink2 --memory 64000 --threads 16 --glm no-x-sex hide-covar cols=$COL_NAMES \
		--pfile $QC_DIR/ukb_imp_chr${i}_v3_11 --keep-${sex}s \
		--pheno ${PHENO_DIR}/pheno_${PHENO}.txt --pheno-name $PHENO \
		--covar ${PHENO_DIR}/covariates.txt --covar-col-nums 3-14 --covar-variance-standardize \
		--out $OUTPUT_DIR/${sex}/${sex}_${i} 
	done
done

plink2 --memory 64000 --threads 16 --glm no-x-sex hide-covar cols=$COL_NAMES \
--pfile $QC_DIR/ukb_imp_chr22_v3_11 \
--pheno ${PHENO_DIR}/pheno_${PHENO}.txt --pheno-name $PHENO \
--covar ${PHENO_DIR}/covariates.txt --covar-col-nums 3-14 \
--out $OUTPUT_DIR/both_sex/both_sex_22