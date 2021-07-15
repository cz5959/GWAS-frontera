#!/usr/bin/env python3

# import packages
import matplotlib
matplotlib.use("Agg")
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
import sys

# load dataframes
pheno = sys.argv[1]
result = "cell_types"
get_names = False
print(pheno, result)

#results
os.chdir("/scratch1/08005/cz5959/LD_practice/{0}".format(pheno))
df = pd.read_csv("{0}_both_sex_{1}.results".format(pheno, result),sep="\t",usecols=['Category','Enrichment','Enrichment_std_error'],dtype={'Category':str,'Enrichment':np.float64,'Enrichment_std_error':np.float64})
df_female = pd.read_csv("{0}_female_{1}.results".format(pheno, result),sep="\t",usecols=['Category','Enrichment','Enrichment_std_error'],dtype={'Category':str,'Enrichment':np.float64,'Enrichment_std_error':np.float64})
df_male = pd.read_csv("{0}_male_{1}.results".format(pheno, result),sep="\t",usecols=['Category','Enrichment','Enrichment_std_error'],dtype={'Category':str,'Enrichment':np.float64,'Enrichment_std_error':np.float64})

# make bar plot
def make_bar(pheno, df1, df2=pd.DataFrame({'A' : []})):
    fig, ax = plt.subplots(figsize=(20,10))
    if df2.empty:
        x = np.arange(len(df1['Category']))
        df1.plot(kind='bar',x='Category',y='Enrichment',ax=ax,legend=None)
        plt.errorbar(x-0.2, df1['Enrichment'], yerr = (2*df1['Enrichment_std_error']), linestyle="none", color="darkslategray", capsize=3, label = "2 SE")
        plt.title("Partitioned Heritability for {0} - {1}".format(pheno, result.capitalize()))
        plot_name="partition_bar_{0}_{1}.png".format(pheno, result)
    else:
        df1.rename(columns={'Enrichment':'Enrichment_Female', 'Enrichment_std_error':'std_error_female'},inplace=True)
        df2.rename(columns={'Enrichment':'Enrichment_Male', 'Enrichment_std_error':'std_error_male'},inplace=True)
        df = pd.merge(df1, df2, how="outer", on="Category")
        x = np.arange(len(df['Category']))
        ax.bar(x-0.2,df['Enrichment_Female'], 0.4,color='plum',label = 'Female')
        ax.bar(x+0.2,df['Enrichment_Male'], 0.4,color='royalblue',label = 'Male')
        plt.errorbar(x-0.2, df['Enrichment_Female'], yerr = (2*df['std_error_female']), linestyle="none", color="darkslategray", capsize=5, label = "2 SE")
        plt.errorbar(x+0.2, df['Enrichment_Male'], yerr = (2*df['std_error_male']), linestyle="none", color="darkslategray", capsize=5, label = '_nolegend_')
        plt.xticks(x,df['Category'],rotation='vertical')
        plt.legend()
        plt.title("Partitioned Heritability for {0} - Male and Female - {1}".format(pheno, result.capitalize()))
        plot_name="partition_bar_{0}_{1}_both.png".format(pheno, result)
    plt.xlim(-1,x.max()+1)
    plt.ylabel("Enrichment")
    plt.xlabel("Category")
    plt.legend()
    plt.tight_layout()
    plt.savefig(plot_name)

# categories:
def categories(result)
    os.chdir("/scratch1/08005/cz5959/ldsc")
    df_names = pd.read_csv("{0}.ldcts".format(result), sep="\t", header=None, usecols=[0])
    row = pd.Series(["Control"])
    df_names2 = pd.concat([row,df_names]).reset_index(drop=True)
    os.chdir("/scratch1/08005/cz5959/LD_practice/{0}".format(pheno))
    return df_names2

if get_names:
    df_names = categories(result)
    df['Category'] = df_names.values
    df_female['Category'] = df_names.values
    df_male['Category'] = df_names.values

make_bar(pheno, df)
make_bar(pheno, df_female, df_male)
    
# CORRELATION
from scipy import stats
corr = stats.pearsonr(df['Enrichment_Female'],df['Enrichment_Male'])
print("Correlation among cell type annotations for {0}: ".format(pheno) + str(corr[0]))