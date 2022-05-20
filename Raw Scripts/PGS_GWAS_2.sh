#!/bin/sh

while getopts p:m:t:s: flag
do
    case "${flag}" in
        p) PHENO=${OPTARG};;
        m) MODE=${OPTARG};;
        t) TYPE=${OPTARG};;
        s) SET=${OPTARG};;
    esac
done
echo $PHENO; echo $MODE; echo $TYPE; echo $SET

PGS_dir=$SCRATCH/GWAS_Results/$PHENO/PGS_$SET  ####
mkdir -p $PGS_dir/{both_sex,female,male}
QC_dir=$SCRATCH/QC

# GWAS for both sex, female, and male specfic
COL_NAMES=chrom,pos,ref,alt,ax,test,nobs,beta,se,tz,p
cd $SCRATCH/Phenotypes 
declare -a arr=("female" "male")
for i in {1..22}
do
    # both sex
    plink2 --memory 64000 --threads 16 --glm no-x-sex hide-covar cols=$COL_NAMES $MODE \
    --pfile $QC_dir/ukb_imp_chr${i}_v3_11 \
    --remove $PGS_dir/${PHENO}_female_testIIDs.txt $PGS_dir/${PHENO}_male_testIIDs.txt \
    --pheno pheno_${PHENO}.txt --pheno-name $PHENO --out $PGS_dir/both_sex/both_sex_${i} \
    --covar covariates.txt --covar-col-nums 3-14 --covar-variance-standardize
    # sex-specific
    for sex in "${arr[@]}"
    do 
        plink2 --memory 64000 --threads 16 --glm no-x-sex hide-covar cols=$COL_NAMES $MODE \
        --pfile $QC_dir/ukb_imp_chr${i}_v3_11 --keep-${sex}s \
        --remove $PGS_dir/${PHENO}_female_testIIDs.txt $PGS_dir/${PHENO}_male_testIIDs.txt \
        --pheno pheno_${PHENO}.txt --pheno-name $PHENO --out $PGS_dir/${sex}/${sex}_${i} \
        --covar covariates.txt --covar-col-nums 3-14 --covar-variance-standardize 
    done
done

# create test bfile
plink2 --memory 64000 --threads 16 --pfile $QC_dir/ukb_imp_all_v3_11 --keep $PGS_dir/${PHENO}_female_testIIDs.txt $PGS_dir/${PHENO}_male_testIIDs.txt \
--make-bed --out $PGS_dir/${PHENO}_test

# combine GWAS results
HEADER='#CHROM\tPOS\tID\tREF\tALT\tA1\tAX\tTEST\tOBS_CT\tBETA\tSE\tT_STAT\tP\n'
declare -a arr=("both_sex" "female" "male")
for sex in "${arr[@]}"
do
    cd $PGS_dir/$sex
    cat *.${TYPE} | awk '$1 !~ /^#CHROM/' - | awk '!/'NA'/' - | { printf $HEADER; cat - ; } > $PGS_dir/${sex}_train.${PHENO}.glm.${TYPE}
done