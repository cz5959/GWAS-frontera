#!/bin/sh

# set up environment
cd /scratch/08005/cz5959/ldsc 
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
ADDITIVE=$SCRATCH/GWAS_Results/$PHENO
LD_FILE=$SCRATCH/LD_practice/$PHENO
LD_SCORE=$SCRATCH/LD_practice/LD_scores

PARTITION=$SCRATCH/LD_practice/Partitioned
BASELINE=$LD_SCORE/1000G_EUR_Phase3_baseline/baseline.
CELLTYPES=$PARTITION/1000G_Phase3_cell_type_groups

# reformat summary statistics
munge_sumstats.py --sumstats $ADDITIVE/sex_diff.${PHENO}.glm.linear \
    --snp ID --a2 AX --N-col OBS_CT --signed-sumstats BETA,0 \
    --merge-alleles $LD_SCORE/w_hm3.snplist \
    --chunksize 500000 --out $LD_FILE/${PHENO}_sexdiff

# cell type group analysis
mkdir -p $LD_FILE/celltypes
for i in {1..10}
do
    ldsc.py --h2 $LD_FILE/${PHENO}_sexdiff.sumstats.gz \
    --ref-ld-chr ${BASELINE},${CELLTYPES}/cell_type_group.${i}. \
    --w-ld-chr $LD_SCORE/weights_hm3_no_hla/weights. --overlap-annot \
    --frqfile-chr $LD_SCORE/1000G_Phase3_frq/1000G.EUR.QC. \
    --out $LD_FILE/celltypes/${PHENO}_sexdiff_${i}
done

# get list of enrichments and se    
cd $LD_FILE/celltypes
grep L2_1 ./*sexdiff*.results | cut -f1,5,6 - > ${PHENO}_sexdiff_celltype_enrichment.txt
grep "Total Observed scale h2" ./*sexdiff*.log > ${PHENO}_sexdiff_celltype_h2.txt


