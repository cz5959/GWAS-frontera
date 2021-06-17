#!/usr/bin/env python3

# import packages
import matplotlib
matplotlib.use("Agg")
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats
import os

# load in dataframes
os.chdir("/scratch1/08005/cz5959/Neale_Lab/Standing_Height")
neale_df = pd.read_csv("50_raw.gwas.imputed_v3.both_sexes.tsv", sep="\t",usecols=['variant','minor_allele','beta','pval'],dtype={'variant':str, 'minor_allele':str, 'beta':np.float64, 'pval':np.float64})
os.chdir("/scratch1/08005/cz5959/Association_Height_50")
plink_df = pd.read_csv("linear_results_all_chrom.height.glm.linear", sep="\t",usecols=['#CHROM','POS','ID','REF','ALT','A1','BETA','P'], dtype={'#CHROM':np.int64,'POS':np.int64,'ID':str,'REF':str,'ALT':str,'A1':str,'BETA':np.float64,'P':np.float64})

# concatenate chrom, position, red and alt for plink variant
plink_df['variant'] = plink_df['#CHROM'].astype(str) + ":" + plink_df['POS'].astype(str) + ":" + plink_df['REF'].astype(str) + ":" + plink_df['ALT'].astype(str) 

# check datatypes and if there are duplicates in plink
print(neale_df.dtypes)
print(plink_df.dtypes)
print(plink_df[plink_df['variant'].duplicated(keep=False)].head(10))
print( len(plink_df[plink_df['variant'].duplicated(keep=False)]['variant'].unique()) )

# join dataframes with inner merge on variant
neale_plink_df = pd.merge(neale_df, plink_df, how="inner", on="variant")
neale_plink_df.rename(columns= {'P':'plink_p', 'pval':'neale_p', 'BETA':'plink_beta', 'beta':'neale_beta'}, inplace=True)
print("length:" + str(len(neale_plink_df)))

# remove na and check joined df
neale_plink_df.dropna(inplace=True)
neale_plink_df.head(20)

# get correlation
corr = stats.pearsonr(neale_plink_df['neale_p'],neale_plink_df['plink_p'])
print("Correlation: " + str(corr[0]))

# SCATTERPLOT

# negative log 10 
neale_plink_df = neale_plink_df[(neale_plink_df['plink_p'] != 0)]
x = -np.log10(neale_plink_df['neale_p'])
y = -np.log10(neale_plink_df['plink_p'])

plt.figure(figsize=(10,10))
plt.scatter(x,y, s=0.5)
plt.xlabel("Neale -log10(p)")
plt.ylabel("plink -log10(p)")
plt.title("Scatterplot of Neale and plink -log10(p) Values")

# save figure as png
plt.tight_layout()
plt.savefig("scatter_neale_plink_pvalues.png")


# BETA
# flip neale's sign if the alt and a1 allele are mismatched; and rename
neale_plink_df.rename(columns= {'ALT':'ALT_plink', 'A1':'A1_plink', 'minor_allele':'minor_allele_neale'}, inplace=True)
neale_plink_df.loc[ (neale_plink_df['ALT_plink'] != neale_plink_df['A1_plink']), 'neale_beta'] = -neale_plink_df['neale_beta']

# beta scatterplot
x = neale_plink_df['neale_beta']
y = neale_plink_df['plink_beta']
plt.figure(figsize=(10,10))
plt.scatter(x,y, s=0.5)
plt.xlabel("Neale Betas")
plt.ylabel("Plink2 Betas")
plt.title("Scatterplot of Neale and Plink2 Beta Values")
plt.tight_layout()
plt.savefig("scatter_neale_plink_betas.png")

# betas for low p-values
low_p_betas_df = neale_plink_df.loc[ (-np.log10(neale_plink_df['plink_p']) > 10), ['variant','ID','ALT','A1','minor_allele','neale_beta','plink_beta']]
x = low_p_betas_df['neale_beta']
y = low_p_betas_df['plink_beta']
plt.figure(figsize=(10,10))
plt.scatter(x,y, s=0.5)
plt.xlabel("Neale Betas")
plt.ylabel("Plink2 Betas")
plt.title("Scatterplot of Neale and Plink2 Significant Beta Values")
plt.tight_layout()
plt.savefig("scatter_neale_plink_sig_betas.png")
