#####################
# configuration file
#####################

# process paramerters
mem=64000
threads=16

#####################
# FILE PATHS

### R library files ###
R_LIB = "/work/08005/cz5959/ls6/R/x86_64-pc-linux-gnu-library/4.0/"

### General ###
QC_DIR="/scratch/08005/cz5959/QC"                   # quality controlled pgen sample and variant files
PHENO_DIR="/scratch/08005/cz5959/Phenotypes"        # phenotype files
GWAS_DIR="/scratch/08005/cz5959/GWAS_Results"       # summary statistics, results, intermediate tables
# mash and PGS files directories are subdirectories of $GWAS_DIR/{phenotype}/ 

### ENSEMBL annotation files ###
ENSEMBL="/scratch/08005/cz5959/ensembl-vep"         # ensembl vep github (https://github.com/Ensembl/ensembl-vep.git)
ENSEMBL_CACHE="/scratch/08005/Annotation/cache"     # homo_sapiens_CRCh38.vcf cache

### LDSC files ###
LDSC_DIR="/scratch/08005/cz5959/ldsc"               # ldsc github (https://github.com/bulik/ldsc.git)
LDSC_FILE="/scratch/08005/cz5959/LD_files"          # results from LDSC regression
LD_SCORE="/scratch/08005/cz5959/LD_scores"          # LD scores
LD_1000G="/scratch/08005/cz5959/1000G/all_phase3"   # 1000 genomes phase 3 file, Europeans, bfile

