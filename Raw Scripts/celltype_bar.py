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
result = "Cahoy"
file_path = "/scratch1/08005/cz5959/LD_practice/{0}".format(pheno)
os.chdir(file_path)
df = pd.read_csv("{0}_both_sex_{1}.cell_type_results.txt".format(pheno, result),sep="\t",usecols=['Name','Coefficient_P_value'],dtype={'Name':str,'Coefficient_P_value':np.float64})
df_female = pd.read_csv("{0}_female_{1}.cell_type_results.txt".format(pheno, result),sep="\t",usecols=['Name','Coefficient_P_value'],dtype={'Name':str,'Coefficient_P_value':np.float64})
df_male = pd.read_csv("{0}_male_{1}.cell_type_results.txt".format(pheno, result),sep="\t",usecols=['Name','Coefficient_P_value'],dtype={'Name':str,'Coefficient_P_value':np.float64})

# neg log p
df['NEG_LOG_P'] = -np.log10(df['Coefficient_P_value'])
df_female['NEG_LOG_P'] = -np.log10(df_female['Coefficient_P_value'])
df_male['NEG_LOG_P'] = -np.log10(df_male['Coefficient_P_value'])

# make bar plot
def make_bar(pheno, df1,df2=pd.DataFrame({'A' : []})):
    fig, ax = plt.subplots(figsize=(20,10))
    if df2.empty:
        x = np.arange(len(df1['Name']))
        df1.plot(kind='bar',x='Name',y='NEG_LOG_P',ax=ax,legend=None)
        plt.title("Cell type Analysis for {0} - {1}".format(pheno, result.capitalize()))
        plot_name="celltype_bar_{0}_{1}.png".format(pheno, result)
    else:
        df1.rename(columns={'NEG_LOG_P':'NEG_LOG_P_Female'},inplace=True)
        df2.rename(columns={'NEG_LOG_P':'NEG_LOG_P_Male'},inplace=True)
        df = pd.merge(df1, df2, how="outer", on="Name")
        x = np.arange(len(df['Name']))
        ax.bar(x-0.2,df['NEG_LOG_P_Female'], 0.4,color='plum',label = 'Female')
        ax.bar(x+0.2,df['NEG_LOG_P_Male'], 0.4,color='royalblue',label = 'Male')
        plt.xticks(x,df['Name'],rotation='vertical')
        plt.legend()
        plt.title("Cell type Analysis for {0} - Male and Female - {1}".format(pheno, result.capitalize()))
        plot_name="celltype_bar_{0}_{1}_both.png".format(pheno, result)
    plt.xlim(-1,x.max()+1)
    plt.ylabel("-LOG10 P")
    plt.xlabel("Name")
    plt.legend()
    plt.tight_layout()
    plt.savefig(plot_name)

make_bar(pheno, df)
make_bar(pheno, df_female, df_male)
    
# CORRELATION
from scipy import stats
corr = stats.pearsonr(df['NEG_LOG_P_Female'],df['NEG_LOG_P_Male'])
print("Correlation among cell type annotations for {0}: ".format(pheno) + str(corr[0]))