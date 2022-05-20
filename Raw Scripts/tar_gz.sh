#!/bin/sh

while getopts t: flag
do
    case "${flag}" in
        p) TAG=${OPTARG};;
    esac
done
echo $TAG 

# mash 
cd $SCRATCH/GWAS_Results/${TAG}/mash
tar -czf $SCRATCH/mash_files/${TAG}_mash.tar.gz ./*.txt

# QC
cd $SCRATCH/QC_Chr
tar -czf $SCRATCH/QC_files/ukb_imp_chr${TAG}_v3_qc.tar.gz ./ukb_imp_chr${TAG}_v3_11.p*
