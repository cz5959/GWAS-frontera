#!/bin/sh

while getopts p:s: flag
do
    case "${flag}" in
        p) PHENO=${OPTARG};;
        s) SET=${OPTARG};;
    esac
done
echo $PHENO; echo $SET

PGS_dir=$SCRATCH/GWAS_Results/$PHENO/PGS_$SET
mkdir -p $PGS_dir/both_sex
QC_dir=$SCRATCH/QC_Chr

# GWAS for both sex, female, and male specfic
COL_NAMES=chrom,pos,ref,alt,ax,test,nobs,beta,se,tz,p
cd $SCRATCH/Phenotypes 
for i in {1..22}
do
    # both sex
    plink2 --memory 64000 --threads 16 --glm no-x-sex hide-covar cols=$COL_NAMES \
    --pfile $QC_dir/ukb_imp_chr${i}_v3_11 \
    --remove $PGS_dir/${PHENO}_female_testIIDs.txt $PGS_dir/${PHENO}_male_testIIDs.txt \
    --keep $PGS_dir/${PHENO}_female_trainIIDs.txt $PGS_dir/${PHENO}_male_trainIIDs.txt \
    --pheno pheno_${PHENO}.txt --pheno-name $PHENO --out $PGS_dir/both_sex/both_sex_${i} \
    --covar covariates.txt --covar-col-nums 3-14 --covar-variance-standardize
done

# create test bfile
plink2 --memory 64000 --threads 16 --pfile $QC_dir/ukb_imp_all_v3_11 --keep $PGS_dir/${PHENO}_female_testIIDs.txt $PGS_dir/${PHENO}_male_testIIDs.txt \
--make-bed --out $PGS_dir/${PHENO}_test

125197

# combine GWAS results
HEADER='#CHROM\tPOS\tID\tREF\tALT\tA1\tAX\tTEST\tOBS_CT\tBETA\tSE\tT_STAT\tP\n'
cd $PGS_dir/both_sex
cat *.linear | awk '$1 !~ /^#CHROM/' - | awk '!/'NA'/' - | { printf $HEADER; cat - ; } > $PGS_dir/both_sex_train.${PHENO}.glm.linear


