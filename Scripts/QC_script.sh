#!/bin/bash
# QC for all chromosomes 
# Updated 
# paths
impute_path=/corral-repl/utexas/Recombining-sex-chro/ukb/data/imputation
mfi_ids_path=/scratch/08005/cz5959/QC/MFI_IDs

# chr 6 QC
	# make pgen
plink2 --memory 64000 --threads 16 --bgen $impute_path/ukb_imp_chr6_v3.bgen ref-first --sample $impute_path/ukb61666_imp_chr6_v3_s487280.sample --make-pgen --out ukb_imp_chr6_v3_1
	
	# # Info score >0.8
plink2 --memory 64000 --threads 16 --pfile ukb_imp_chr6_v3_1 --extract $mfi_ids_path/ukb_mfi_chr6_v3_IDs.txt --make-pgen --out ukb_imp_chr6_v3_2
	 
    # 	# call rate > 0.95 (missingness)
plink2 --memory 64000 --threads 16 --pfile ukb_imp_chr6_v3_2 --geno 0.05 --mind 0.05 --make-pgen --out ukb_imp_chr6_v3_3

    # 	# Alternate frequency (MAF) > 0.001 && <0.999
plink2 --memory 64000 --threads 16 --pfile ukb_imp_chr6_v3_3 --maf 0.001 --max-maf 0.999 --make-pgen --out ukb_imp_chr6_v3_4

    # 	# HWE > 1e-10
plink2 --memory 64000 --threads 16 --pfile ukb_imp_chr6_v3_4 --hwe 1e-10 --make-pgen --out ukb_imp_chr6_v3_5

for i in {1..22}
do 
	# make pgen
	plink2 --memory 64000 --threads 16 --bgen $impute_path/ukb_imp_chr${i}_v3.bgen ref-first --sample $impute_path/ukb61666_imp_chr${i}_v3_s487280.sample --make-pgen --out ukb_imp_chr${i}_v3_1
	
	# # Info score >0.8
	plink2 --memory 64000 --threads 16 --pfile ukb_imp_chr${i}_v3_1 --extract $mfi_ids_path/ukb_mfi_chr${i}_v3_IDs.txt --make-pgen --out ukb_imp_chr${i}_v3_2
	 
    # 	# call rate > 0.95 (missingness)
    # plink2 --memory 64000 --threads 16 --pfile ukb_imp_chr${i}_v3_2 --missing 
	plink2 --memory 64000 --threads 16 --pfile ukb_imp_chr${i}_v3_2 --geno 0.05 --mind 0.05 --make-pgen --out ukb_imp_chr${i}_v3_3
    # 	Rscript --no-save hist_miss.R $i
    
    # 	# Alternate frequency (MAF) > 0.001 && <0.999
    # plink2 --memory 64000 --threads 16 --pfile ukb_imp_chr${i}_v3_3 --freq --out MAF_check
	plink2 --memory 64000 --threads 16 --pfile ukb_imp_chr${i}_v3_3 --maf 0.001 --max-maf 0.999 --make-pgen --out ukb_imp_chr${i}_v3_4
    # 	Rscript --no-save MAF_check.R
    
    # 	# HWE > 1e-10
    # 	plink2 --memory 64000 --threads 16 --pfile ukb_imp_chr${i}_v3_4 --hardy 
	# plink2 --memory 64000 --threads 16 --pfile ukb_imp_chr${i}_v3_4 --hwe 1e-10 --make-pgen --out ukb_imp_chr${i}_v3_5
    # 	awk '{ if ($10 <0.000000001) print $0 }' plink2.hardy>plink2zoomhwe.hardy
    # 	Rscript --no-save hwe.R $i

	# remove duplicates and keep only snps
	#plink2 --memory 64000 --threads 16 --pfile $QC_file/ukb_imp_chr${i}_v3_5 --snps-only --rm-dup force-first --make-pgen --out $QC_file/ukb_imp_chr${i}_v3_6

    # remove indels and multiallelic snps
        awk -F '\t' 'index($3, ":") !=0 {print $3}' ukb_imp_chr${i}_v3_5.pvar  > indels_chr${i}.txt
        plink2 --pfile ukb_imp_chr${i}_v3_5 --exclude indels_chr${i}.txt --snps-only --make-pgen --out ukb_imp_chr${i}_v3_7


done



