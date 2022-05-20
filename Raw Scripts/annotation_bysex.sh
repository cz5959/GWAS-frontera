#!/bin/sh

# ls6
# get phenotype
while getopts p: flag
do
    case "${flag}" in
        p) PHENO=${OPTARG};;
    esac
done
echo $PHENO 

# clump
FILE=$SCRATCH/GWAS_Results/$PHENO
LD=$SCRATCH/1000G/all_phase3
mkdir -p $FILE/annotation
declare -a arr=("female" "male")
for sex in "${arr[@]}"
do
plink \
    --bfile $LD \
    --clump-p1 5e-8 \
    --clump-r2 0.1 \
    --clump-kb 250 \
    --clump $FILE/${sex}_all.${PHENO}.glm.linear \
    --clump-snp-field ID \
    --clump-field P \
    --out $FILE/annotation/${sex}_$PHENO

sed 's/ \+/\t/g' $FILE/annotation/${sex}_${PHENO}.clumped | cut -f4 > $FILE/annotation/${sex}_${PHENO}.clumped.ids
done

scp $FILE/female_${PHENO}.clumped.ids cz5959@frontera.tacc.utexas.edu:\$SCRATCH/Annotation
scp $FILE/male_${PHENO}.clumped.ids cz5959@frontera.tacc.utexas.edu:\$SCRATCH/Annotation

#### frontera

PHENO=height
sex=female
cd $SCRATCH/ensembl-vep
./vep -i $SCRATCH/GWAS_Results/${PHENO}/annotation/${sex}_${PHENO}.clumped.ids \
--dir_cache $SCRATCH/Annotation/cache --cache -o $SCRATCH/GWAS_Results/${PHENO}/annotation/${sex}_${PHENO}_annotation.txt \
--force_overwrite --tab --symbol --nearest symbol --pick \
--fields "Uploaded_variation,Location,Allele,SYMBOL,Gene,NEAREST,Feature,Feature_type,Consequence" --stats_text