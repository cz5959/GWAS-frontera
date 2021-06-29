#!/usr/bin/env python3

# import packages
import numpy as np
import pandas as pd
import os
import sys

pheno=sys.argv[1]
field_id=sys.argv[2]

# load dataframes
os.chdir("/scratch1/08005/cz5959/Neale_Lab")
variants_df = pd.read_csv("variants.tsv", sep="\t",usecols=['variant','ref','alt','rsid'],dtype={'variant':str,'ref':str,'alt':str,'rsid':str})
os.chdir("/scratch1/08005/cz5959/Neale_Lab/{0}".format(pheno))
neale_df = pd.read_csv("{0}_raw.gwas.imputed_v3.both_sexes.tsv".format(field_id), sep="\t",usecols=['variant','minor_allele','n_complete_samples','beta','pval'],dtype={'variant':str, 'minor_allele':str, 'n_complete_samples':np.int64, 'beta':np.float64, 'pval':np.float64})

# remove multialleleic and indels
variants_df = variants_df[ (variants_df['ref'].str.len() == 1) & (variants_df['alt'].str.len() == 1)]
variants_df = variants_df[~variants_df['rsid'].str.contains(":")]

# merge dataframes and remove na
neale_df_2 = pd.merge(neale_df, variants_df, how="inner", on="variant")
df.dropna(inplace=True)

# set non-effect allele and flip betas where minor allele != alt allele
neale_df_2['A2'] = np.where( neale_df_2['minor_allele'] != neale_df_2['alt'], neale_df_2['alt'], neale_df_2['ref'])
neale_df_2.loc[ (neale_df_2['minor_allele'] != neale_df_2['alt']), 'beta'] = -neale_df_2['beta']

# rename and reorder columns and save as csv
neale_df_2.rename(columns= {'rsid':'RSID', 'pval':'P', 'beta':'BETA', 'n_complete_samples':'N','minor_allele':'A1'}, inplace=True)
neale_df_2 = neale_df_2[['RSID','A1','A2','N','BETA','P']]
neale_df_2.to_csv("neale_{0}_ldsc.txt".format(pheno), index=False, sep="\t")

