#!/usr/bin/env python3

# import packages
import numpy as np
import pandas as pd
import os
import sys
import scipy.stats
import argparse

parser = argparse.ArgumentParse()
parser.add_argument("-p", "--pheno", dest = "pheno", default = "", help="phenotype")
args = parser.parse_args()
pheno = args.pheno

file_path = "/scratch1/08005/cz5959/GWAS_Results/{0}".format(pheno)
os.chdir(file_path)

# load dataframes
mash_lfsr = pd.read_csv("{0}_mash_lfsr.txt".format(pheno), sep="\t", usecols=['female','male'],dtype={'female':np.float64, 'male':np.float64})
mash_psd = pd.read_csv("{0}_mash_psd.txt".format(pheno), sep="\t", usecols=['female','male'],dtype={'female':np.float64, 'male':np.float64})
mash_pm = pd.read_csv("{0}_mash_pm.txt".format(pheno), sep="\t", usecols=['female','male'],dtype={'female':np.float64, 'male':np.float64})
df_glm = pd.read_csv("female_all.{0}.glm.linear".format(pheno), sep="\t", usecols =['ID','A1','AX','OBS_CT'],dtype= {'ID':str,'A1':str,'AX':str,'OBS_CT':np.int64})

if len(mash_pm) != len(df_glm):
    print("lengths dont match for {0}".format(pheno))
    exit()

# sample size
n = df_glm['OBS_CT'].loc[0]
# get se by  sd/sqrt(n)
mash_pse = mash_psd.apply(lambda x: x / np.sqrt(n))
# get z score
mash_z = mash_pm.div(mash_pse)
# get p value 
mash_p = mash_z.apply(lambda x: scipy.stats.norm.sf(abs(x)))

# rename
mash_pm.rename(columns={'male':'m_pm','female':'f_pm'}, inplace=True)
mash_pse.rename(columns={'male':'m_pse','female':'f_pse'}, inplace=True)
mash_z.rename(columns={'male':'m_z','female':'f_z'}, inplace=True)
mash_p.rename(columns={'male':'m_p','female':'f_p'}, inplace=True)

## need snp id, alleles a1 and a2(AX), ncol(OBS_CT), pvalue(lfsr), beta(pm)
# get snp id, alleles, and ncol from glm results file   # get beta and p from mash posterior file

female_df = pd.concat([pd.concat([pd.concat([pd.concat([df_glm, mash_pm['f_pm']], axis=1), mash_pse['f_pse']], axis=1), mash_z['f_z']], axis=1), mash_p['f_p']], axis=1)
male_df = pd.concat([pd.concat([pd.concat([pd.concat([df_glm, mash_pm['m_pm']], axis=1), mash_pse['m_pse']], axis=1), mash_z['m_z']], axis=1), mash_p['m_p']], axis=1)

male_df.rename(columns={'ID':'SNP','AX':'A2','OBS_CT':'N','m_pm':'BETA','m_pse':'SE','m_z':'Z','m_p':'P'}, inplace=True)
female_df.rename(columns={'ID':'SNP','AX':'A2','OBS_CT':'N','f_pm':'BETA','f_pse':'SE','f_z':'Z','f_p':'P'}, inplace=True)

# create csv
male_df.to_csv("{0}_male_mash_posterior.txt".format(pheno), sep="\t", index=False)
female_df.to_csv("{0}_female_mash_posterior.txt".format(pheno), sep="\t", index=False)


