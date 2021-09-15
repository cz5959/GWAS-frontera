#!/bin/sh

PHENO=$1
echo $1
mkdir -p $SCRATCH/GWAS_Results/$PHENO/PRS
PRS_dir=$SCRATCH/GWAS_Results/$PHENO/PRS
mkdir -p $PRS_dir/{both_sex,female,male}
QC_dir=$SCRATCH/QC

# GWAS
COL_NAMES=chrom,pos,ref,alt,ax,test,nobs,beta,se,tz,p
cd $SCRATCH/Phenotypes
for i in {1..22}
do 
	# both sex
	plink2 --memory 64000 --threads 16 --glm no-x-sex hide-covar cols=$COL_NAMES \
    --pfile $QC_dir/ukb_imp_chr${i}_v3_prs --keep whitebritishIDs.txt \               ####### have qc file already be only white british
    --remove ${PHENO}_female_testIIDs.txt ${PHENO}_male_testIIDs.txt
    --pheno pheno_${PHENO}.txt --pheno-name $PHENO --out $PRS_dir/both_sex/both_sex_${i} \
    --covar covariates.txt --covar-col-nums 3-14 --covar-variance-standardize 
	# females
	plink2 --memory 64000 --threads 16 --glm no-x-sex hide-covar cols=$COL_NAMES \
    --pfile $QC_dir/ukb_imp_chr${i}_v3_prs --keep whitebritishIDs.txt --keep-females \
    --remove ${PHENO}_female_testIIDs.txt ${PHENO}_male_testIIDs.txt
    --pheno pheno_${PHENO}.txt --pheno-name $PHENO --out $PRS_dir/female/female_${i} \
    --covar covariates.txt --covar-col-nums 3-12,14 --covar-variance-standardize
	# males
	plink2 --memory 64000 --threads 16 --glm no-x-sex hide-covar cols=$COL_NAMES \
    --pfile $QC_dir/ukb_imp_chr${i}_v3_prs --keep whitebritishIDs.txt --keep-males \
    --remove ${PHENO}_female_testIIDs.txt ${PHENO}_male_testIIDs.txt
    --pheno pheno_${PHENO}.txt --pheno-name $PHENO --out $PRS_dir/male/male_${i} \
    --covar covariates.txt --covar-col-nums 3-12,14 --covar-variance-standardize
done

# merge test bfile
plink2 --memory 64000 --threads 16 --pmerge-list merge_pfile_chr.txt pfile \
--keep ${PHENO}_female_testIIDs.txt ${PHENO}_male_testIIDs.txt --make-bed --out ${PHENO}_test

# combine GWAS results
HEADER='#CHROM\tPOS\tID\tREF\tALT\tA1\tAX\tTEST\tOBS_CT\tBETA\tSE\tT_STAT\tP\n'
declare -a arr=("both_sex" "female" "male")
for sex in "${arr[@]}"
do
    cd $PRS_dir/$sex
    cat *.linear | awk '$1 !~ /^#CHROM/' - | awk '!/'NA'/' - | { printf $HEADER; cat - ; } > $PRS_dir/${sex}_train.${1}.glm.linear
done