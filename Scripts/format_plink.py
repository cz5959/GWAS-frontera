#!/usr/bin/env python3

# import packages
import numpy as np
import pandas as pd
import os
import sys

def edit(sex):
    # load file
    file_name = "{1}_all.{0}.glm.linear".format(pheno, sex)
    df = pd.read_csv(file_name, sep="\t", usecols =['#CHROM','POS','ID','A1','AX','OBS_CT','BETA','P'], dtype= {'#CHROM':np.int64,'POS':np.int64,'ID':str,'A1':str,'AX':str,'OBS_CT':np.int64,'BETA':np.float64,'P':object})
    
    # convert to float
    df['ref_P']=df['P']
    print(sex)
    print(df.head(10))
    df['P'] = pd.to_numeric(results_df['P'], errors='coerce')
    NA_df = df.loc[df['P'].isna(), ['#CHROM','ID','A1','P','ref_P']]
    print("Number of SNPs can't be converted to numeric: " + str(len(NA_df)))
    if len(NA_df) > 0:
        NA_df.to_csv("{1}_{0}_NA.csv".format(pheno,sex), sep="\t", index = False)

    # drop null 
    df = df.dropna()
    df.drop(columns=['ref_P'], inplace=True)

    # create csv
    df.to_csv("{1}_all_2.{0}.glm.linear".format(pheno, sex), sep="\t", index=False)

pheno=sys.argv[1]
file_path = "/scratch1/08005/cz5959/GWAS_Results/{0}".format(pheno)
os.chdir(file_path)

edit("both_sex")
edit("female")
edit("male")

print("done")
