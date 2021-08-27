#!/usr/bin/env python3

# import packages
import numpy as np
import pandas as pd
import os
import sys

pheno=sys.argv[1]
file_path = "/scratch1/08005/cz5959/GWAS_Results/{0}".format(pheno)
os.chdir(file_path)

# load dataframes
mash_lfsr = pd.read_csv("{0}_mash_lfsr.txt".format(pheno), sep="\t", usecols=['female','male'],dtype={'female':np.float64, 'male':np.float64})
mash_pm = pd.read_csv("{0}_mash_pm.txt".format(pheno), sep="\t", usecols=['female','male'],dtype={'female':np.float64, 'male':np.float64})
df_glm = pd.read_csv("female_all.{0}.glm.linear".format(pheno), sep="\t", usecols =['ID','A1','AX','OBS_CT'],dtype= {'ID':str,'A1':str,'AX':str,'OBS_CT':np.int64,})

if len(mash_pm) != len(df_glm):
    print("lengths dont match for {0}".format(pheno))
    exit()

mash_lfsr.rename(columns={'male':'m_lfsr','female':'f_lfsr'}, inplace=True)
mash_pm.rename(columns={'male':'m_pm','female':'f_pm'}, inplace=True)

## need snp id, alleles a1 and a2(AX), ncol(OBS_CT), pvalue(lfsr), beta(pm)
# get snp id, alleles, and ncol from glm results file   # get beta and p from mash posterior file

female_df = pd.concat([df_glm, mash_pm['f_pm']], axis=1)
female_df = pd.concat([female_df, mash_lfsr['f_lfsr']], axis=1)
male_df = pd.concat([df_glm, mash_pm['m_pm']], axis=1)
male_df = pd.concat([male_df, mash_lfsr['m_lfsr']], axis=1)

male_df.rename(columns={'ID':'SNP','AX':'A2','OBS_CT':'N','m_pm':'BETA','m_lfsr':'P'}, inplace=True)
female_df.rename(columns={'ID':'SNP','AX':'A2','OBS_CT':'N','f_pm':'BETA','f_lfsr':'P'}, inplace=True)

# create csv
male_df.to_csv("{0}_male_mash_posterior.txt".format(pheno), sep="\t", index=False)
female_df.to_csv("{0}_female_mash_posterior.txt".format(pheno), sep="\t", index=False)
