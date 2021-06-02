#!/bin/bash
output_path=/scratch1/08005/cz5959/QC/MFI_IDs
mfi_path=/scratch1/08005/cz5959/Impute_MAF_infoscore

for i in {1..22}
do
	awk '{if($8>0.8) print $2}' $mfi_path/ukb_mfi_chr${i}_v3.txt > $output_path/ukb_mfi_chr${i}_v3_IDs.txt
	#sed -i '1i #FID\tIID' ./MFI_IDs/ukb_mfi_chr${i}_v3_IDs.txt
done
