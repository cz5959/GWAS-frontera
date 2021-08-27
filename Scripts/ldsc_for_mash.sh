#!/bin/sh

# set up environment
cd /scratch1/08005/cz5959/ldsc 
source activate ldsc

pheno=$1
# file paths
MASH=$SCRATCH/GWAS_Results/$pheno
FILE_PATH=$SCRATCH/LD_practice/$pheno
LD_SCORE=$SCRATCH/LD_practice/LD_scores
PARTITION=$SCRATCH/LD_practice/Partitioned
C=$PARTITION/cell_type_groups
CELL_TYPES=${C}/Adrenal_Pancreas.,${C}/Cardiovascular.,${C}/CNS.,${C}/Connective_Bone.,${C}/GI.,${C}/Immune.,${C}/Kidney.,${C}/Liver.,${C}/SkeletalMuscle.,${C}/Other.

# reformat summary statistics
munge_sumstats.py --sumstats $MASH/${pheno}_female_mash_posterior.txt --out $FILE_PATH/${pheno}_female_mash --merge-alleles $LD_SCORE/w_hm3.snplist
munge_sumstats.py --sumstats $MASH/${pheno}_male_mash_posterior.txt --out $FILE_PATH/${pheno}_male_mash --merge-alleles $LD_SCORE/w_hm3.snplist

# baseline partition
ldsc.py	--h2 $FILE_PATH/${pheno}_female_mash.sumstats.gz --ref-ld-chr $PARTITION/baseline/baseline. --w-ld-chr $PARTITION/weights_hm3_no_hla/weights. \
--overlap-annot --frqfile-chr $PARTITION/1000G_frq/1000G.mac5eur. --out $FILE_PATH/${pheno}_female_mash_baseline
ldsc.py	--h2 $FILE_PATH/${pheno}_male_mash.sumstats.gz --ref-ld-chr $PARTITION/baseline/baseline. --w-ld-chr $PARTITION/weights_hm3_no_hla/weights. \
--overlap-annot --frqfile-chr $PARTITION/1000G_frq/1000G.mac5eur. --out $FILE_PATH/${pheno}_male_mash_baseline

# cell type
ldsc.py --h2 $FILE_PATH/${pheno}_female_mash.sumstats.gz --ref-ld-chr $CELL_TYPES --w-ld-chr $PARTITION/weights_hm3_no_hla/weights. \
--overlap-annot --frqfile-chr $PARTITION/1000G_Phase3_frq/1000G.EUR.QC. --out $FILE_PATH/${pheno}_female_mash_cell_types
ldsc.py --h2 $FILE_PATH/${pheno}_male_mash.sumstats.gz --ref-ld-chr $CELL_TYPES --w-ld-chr $PARTITION/weights_hm3_no_hla/weights. \
--overlap-annot --frqfile-chr $PARTITION/1000G_Phase3_frq/1000G.EUR.QC. --out $FILE_PATH/${pheno}_male_mash_cell_types