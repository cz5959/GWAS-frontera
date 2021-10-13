#!/bin/sh

# set up environment
cd /scratch1/08005/cz5959/ldsc 
source activate ldsc

# cmd line argument - phenotype
while getopts p: flag
do
    case "${flag}" in
        p) PHENO=${OPTARG};;
    esac
done
echo $PHENO 

# file paths
LD_FILE=$SCRATCH/LD_files/$PHENO
PARTITION=$SCRATCH/LD_practice/Partitioned

BASELINE=$PARTITION/baseline/baseline.
CELLTYPES=$(cat "$PARTITION/celltypes.txt")

mkdir -p $LD_FILE/{both_sex,female,male}

declare -a arr=("both_sex" "female" "male")
for sex in "${arr[@]}"
do
    mkdir -p $LD_FILE/$sex/celltypes

    # baseline (LD phase1)
    ldsc.py	--h2 $LD_FILE/${PHENO}_${sex}.sumstats.gz \
    --ref-ld-chr $BASELINE --w-ld-chr $PARTITION/weights_hm3_no_hla/weights. \
    --overlap-annot --frqfile-chr $PARTITION/1000G_frq/1000G.mac5eur. \
    --out $LD_FILE/$sex/${PHENO}_${sex}_baseline

    # cell type group analysis (LD phase3)
    for celltype in $CELLTYPES
    do
        ldsc.py --h2 $LD_FILE/${PHENO}_both_sex.sumstats.gz \
        --ref-ld-chr ${BASELINE},$PARTITION/cell_type_groups/$celltype \
        --w-ld-chr $PARTITION/weights_hm3_no_hla/weights. --overlap-annot \
        --frqfile-chr $PARTITION/1000G_Phase3_frq/1000G.EUR.QC. \
        --out $LD_FILE/$sex/celltypes/${PHENO}_${sex}_${celltype}
    done
done

