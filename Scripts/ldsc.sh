#!/bin/sh
i=testosterone
# file paths
PLINK=$SCRATCH/GWAS_Results/$1
FILE_PATH=$SCRATCH/LD_practice/$1
LD_SCORE=$SCRATCH/LD_practice/LD_scores
NEALE=$SCRATCH/Neale_Lab/$1
PARTITION=$SCRATCH/LD_practice/Partitioned
C=$PARTITION/
CELL_TYPES=${C}/Adrenal_Pancreas.,${C}/Cardiovascular.,${C}/CNS.,${C}/Connective_Bone.,${C}/GI.,${C}/Immune.,${C}/Kidney.,${C}/Liver.,${C}/SkeletalMuscle.,${C}/Other.
TYPE=Cahoy
echo $1

# set up environment
cd /scratch1/08005/cz5959/ldsc 
conda activate ldsc

#### REFORMAT SUMM STAT
# need rsid; effect allele; non-effect allele; sample size; p-value; signed summ stat (beta or log odds, z-score, etc)
# plinkc
munge_sumstats.py --sumstats $PLINK/both_sex_all.${1}.glm.linear --snp ID --a2 AX --N-col OBS_CT --out $FILE_PATH/${1}_both_sex --merge-alleles $LD_SCORE/w_hm3.snplist
munge_sumstats.py --sumstats $PLINK/female_all.${1}.glm.linear --snp ID --a2 AX --N-col OBS_CT --out $FILE_PATH/${1}_female --merge-alleles $LD_SCORE/w_hm3.snplist
munge_sumstats.py --sumstats $PLINK/male_all.${1}.glm.linear --snp ID --a2 AX --N-col OBS_CT --out $FILE_PATH/${1}_male --merge-alleles $LD_SCORE/w_hm3.snplist
# neale lab
#munge_sumstats.py --sumstats $NEALE/neale_${1}_ldsc.txt --out $NEALE/neale_${1} --merge-alleles $LD_SCORE/w_hm3.snplist

#### LDSC REGRESSION
# compute heritability and ldsc intercept
# plink
ldsc.py --h2 $FILE_PATH/${1}_both_sex.sumstats.gz --w-ld-chr $LD_SCORE/eur_w_ld_chr/ --ref-ld-chr $LD_SCORE/eur_w_ld_chr/ --out $FILE_PATH/${1}_both_sex_h2
##ldsc.py --h2 $FILE_PATH/${1}_female.sumstats.gz --w-ld-chr $LD_SCORE/eur_w_ld_chr/ --ref-ld-chr $LD_SCORE/eur_w_ld_chr/ --out $FILE_PATH/${1}_female_h2
##ldsc.py --h2 $FILE_PATH/${1}_male.sumstats.gz --w-ld-chr $LD_SCORE/eur_w_ld_chr/ --ref-ld-chr $LD_SCORE/eur_w_ld_chr/ --out $FILE_PATH/${1}_male_h2
# neale
#ldsc.py --h2 $NEALE/neale_${1}.sumstats.gz --w-ld-chr $LD_SCORE/eur_w_ld_chr/ --ref-ld-chr $LD_SCORE/eur_w_ld_chr/ --out $NEALE/neale_${1}_h2

# ldsc summary
grep -h -A 4 "h2:" $FILE_PATH/*h2.log > ${1}_ldsc.txt

# CORRELATION
ldsc.py --rg $FILE_PATH/${1}_female.sumstats.gz,$FILE_PATH/${1}_male.sumstats.gz --w-ld-chr $LD_SCORE/eur_w_ld_chr/ --ref-ld-chr $LD_SCORE/eur_w_ld_chr/ --out $FILE_PATH/${1}_male_female
#partitioned
ldsc.py --rg $FILE_PATH/${1}_female.sumstats.gz,$FILE_PATH/${1}_male.sumstats.gz --ref-ld-chr $C --w-ld-chr $PARTITION/weights_hm3_no_hla/weights. --overlap-annot --frqfile-chr $PARTITION/1000G_frq/1000G.mac5eur. --out $FILE_PATH/${1}_part_corr

# PARTITIONED
# baseline
ldsc.py	--h2 $FILE_PATH/${1}_both_sex.sumstats.gz --ref-ld-chr $PARTITION/baseline/baseline. --w-ld-chr $PARTITION/weights_hm3_no_hla/weights. --overlap-annot --frqfile-chr $PARTITION/1000G_frq/1000G.mac5eur. --out $FILE_PATH/${1}_both_sex_baseline
ldsc.py	--h2 $FILE_PATH/${1}_female.sumstats.gz --ref-ld-chr $PARTITION/baseline/baseline. --w-ld-chr $PARTITION/weights_hm3_no_hla/weights. --overlap-annot --frqfile-chr $PARTITION/1000G_frq/1000G.mac5eur. --out $FILE_PATH/${1}_female_baseline
ldsc.py	--h2 $FILE_PATH/${1}_male.sumstats.gz --ref-ld-chr $PARTITION/baseline/baseline. --w-ld-chr $PARTITION/weights_hm3_no_hla/weights. --overlap-annot --frqfile-chr $PARTITION/1000G_frq/1000G.mac5eur. --out $FILE_PATH/${1}_male_baseline

# cell type group analysis
ldsc.py	--h2 $FILE_PATH/${1}_both_sex.sumstats.gz --ref-ld-chr $CELL_TYPES --w-ld-chr $PARTITION/weights_hm3_no_hla/weights. --overlap-annot --frqfile-chr $PARTITION/1000G_frq/1000G.mac5eur. --out $FILE_PATH/${1}_both_sex_cell_types
ldsc.py	--h2 $FILE_PATH/${1}_female.sumstats.gz --ref-ld-chr $CELL_TYPES --w-ld-chr $PARTITION/weights_hm3_no_hla/weights. --overlap-annot --frqfile-chr $PARTITION/1000G_frq/1000G.mac5eur. --out $FILE_PATH/${1}_female_cell_types
ldsc.py	--h2 $FILE_PATH/${1}_male.sumstats.gz --ref-ld-chr $CELL_TYPES --w-ld-chr $PARTITION/weights_hm3_no_hla/weights. --overlap-annot --frqfile-chr $PARTITION/1000G_frq/1000G.mac5eur. --out $FILE_PATH/${1}_male_cell_types

#partition
ldsc.py	--h2-cts $FILE_PATH/${1}_both_sex.sumstats.gz --ref-ld-chr-cts ${TYPE}.ldcts --ref-ld-chr $PARTITION/1000G_EUR_Phase3_baseline/baseline. --w-ld-chr $PARTITION/weights_hm3_no_hla/weights. --out $FILE_PATH/${1}_both_sex_$TYPE
ldsc.py	--h2-cts $FILE_PATH/${1}_female.sumstats.gz --ref-ld-chr-cts ${TYPE}.ldcts --ref-ld-chr $PARTITION/1000G_EUR_Phase3_baseline/baseline. --w-ld-chr $PARTITION/weights_hm3_no_hla/weights. --out $FILE_PATH/${1}_female_$TYPE
ldsc.py	--h2-cts $FILE_PATH/${1}_male.sumstats.gz --ref-ld-chr-cts ${TYPE}.ldcts --ref-ld-chr $PARTITION/1000G_EUR_Phase3_baseline/baseline. --w-ld-chr $PARTITION/weights_hm3_no_hla/weights. --out $FILE_PATH/${1}_male_$TYPE

echo completed
