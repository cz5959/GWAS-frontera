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
file_path = "/scratch1/08005/cz5959/LD_practice/{0}".format(pheno)
os.chdir(file_path)
df = pd.read_csv("{0}_both_sex_baseline.results".format(pheno),sep="\t",usecols=['Category','Enrichment'],dtype={'Category':str,'Enrichment':np.float64})
df_female = pd.read_csv("{0}_female_baseline.results".format(pheno),sep="\t",usecols=['Category','Enrichment'],dtype={'Category':str,'Enrichment':np.float64})
df_male = pd.read_csv("{0}_male_baseline.results".format(pheno),sep="\t",usecols=['Category','Enrichment'],dtype={'Category':str,'Enrichment':np.float64})

# make bar plot
def make_bar(pheno, df1,df2=pd.DataFrame({'A' : []})):
    fig, ax = plt.subplots(figsize=(20,10))
    if df2.empty:
        x = np.arange(len(df['Category']))
        df1.plot(kind='bar',x='Category',y='Enrichment',ax=ax,legend=None)
        plt.title("Partitioned Heritability for {0}".format(pheno))
        plot_name="partition_bar_{0}.png".format(pheno)
    else:
        df1.rename(columns={'Enrichment':'Enrichment_Female'},inplace=True)
        df2.rename(columns={'Enrichment':'Enrichment_Male'},inplace=True)
        df = pd.merge(df1, df2, how="outer", on="Category")
        x = np.arange(len(df['Category']))
        ax.bar(x-0.2,df['Enrichment_Female'], 0.4,color='plum',label = 'Female')
        ax.bar(x+0.2,df['Enrichment_Male'], 0.4,color='royalblue',label = 'Male')
        plt.xticks(x,df['Category'],rotation='vertical')
        plt.legend()
        plt.title("Partitioned Heritability for {0} - Male and Female".format(pheno))
        plot_name="partition_bar_{0}_both.png".format(pheno)
    plt.xlim(-1,x.max()+1)
    plt.ylabel("Enrichment")
    plt.xlabel("Category")
    plt.tight_layout()
    plt.savefig(plot_name)

make_bar(pheno, df)
make_bar(pheno, df_female, df_male)
    