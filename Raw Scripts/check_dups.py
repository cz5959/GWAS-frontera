#!/usr/bin/env python3

# import packages
import numpy as np
import pandas as pd
import os

os.chdir("/scratch1/08005/cz5959/Association_Height_50")
plink_df = pd.read_csv("linear_results_all_chrom_dups.height.glm.linear", sep="\t",usecols=['#CHROM','ID','POS','REF','ALT','P'], dtype={'#CHROM':np.int64,'ID':str, 'POS':np.int64,'REF':str,'ALT':str,'P':np.float64})

plink_df['variant'] = plink_df['#CHROM'].astype(str) + ":" + plink_df['POS'].astype(str) + ":" + plink_df['REF'].astype(str) + ":" + plink_df['ALT'].astype(str) 

print(plink_df[plink_df['variant'].duplicated(keep=False)])

print( plink_df.loc[plink_df['ID'].duplicated(keep=False) & (plink_df['#CHROM'] == 22) , ['ID','REF','ALT','variant']]['ID'].unique() )


### CHROM 22 ###    
#check mismatch between rmdup list and before rmdup QC
plink_22 = plink_df.loc[plink_df['#CHROM'] == 22, ['ID','REF','ALT','variant']]

os.chdir("/scratch1/08005/cz5959/QC")
plink_QC5 = pd.read_csv("ukb_imp_chr22_v3_5.pvar", sep="\t",usecols=['#CHROM','ID','POS','REF','ALT'], dtype={'#CHROM':np.int64,'ID':str, 'POS':np.int64,'REF':str,'ALT':str})
plink_rmdup = pd.read_csv("ukb_imp_chr22_v3_6.rmdup.list", names=['dup_IDs'])

plink_22_dups = plink_22[plink_22['ID'].duplicated(keep=False)]
plink_QC5_dups = plink_QC5[plink_QC5['ID'].duplicated(keep=False)]

df_QC_rmdup_mismatch = plink_QC5_dups[~plink_QC5_dups['ID'].isin(plink_rmdup['dup_IDs'])]

df_22_rmdup_mismatch = plink_22_dups[~plink_22_dups['ID'].isin(plink_rmdup['dup_IDs'])]

# check QC file after removing duplicates and keeping SNPS only
plink_QC6 = pd.read_csv("ukb_imp_chr22_v3_6.pvar", sep="\t",usecols=['#CHROM','ID','POS','REF','ALT'], dtype={'#CHROM':np.int64,'ID':str, 'POS':np.int64,'REF':str,'ALT':str})
plink_QC6.loc[plink_QC6['ID'].isin(df_QC_rmdup_mismatch['ID']) ]

# check indels
plink_QC5[plink_QC5['ID'].str.contains(":")]
plink_QC6[plink_QC6['ID'].str.contains(":")]
