#!/usr/bin/env python3

### moved to JUPYTER ###

# import packages
import matplotlib
matplotlib.use("Agg")
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
import sys

# load dataframes
pheno = "height"
result = "baseline"
plot_name="partition_mash_{0}_{1}_both.pdf".format(pheno, result)

os.chdir('C:\\Users\\Carrie Zhu\\Documents\\Research\\GWAS-frontera\\LDSC\\{0}'.format(pheno))
df_male = pd.read_csv("{0}_female_mash_{1}.results".format(pheno, result),sep="\t",usecols=['Category','Enrichment','Enrichment_std_error'],dtype={'Category':str,'Enrichment':np.float64,'Enrichment_std_error':np.float64})
df_female = pd.read_csv("{0}_male_mash_{1}.results".format(pheno, result),sep="\t",usecols=['Category','Enrichment','Enrichment_std_error'],dtype={'Category':str,'Enrichment':np.float64,'Enrichment_std_error':np.float64})

df_female.rename(columns={'Enrichment':'Enrichment_Female', 'Enrichment_std_error':'std_error_female'},inplace=True)
df_male.rename(columns={'Enrichment':'Enrichment_Male', 'Enrichment_std_error':'std_error_male'},inplace=True)
df = pd.merge(df_female, df_male, how="outer", on="Category")

fig, ax = plt.subplots(figsize=(20,10))
x = np.arange(len(df['Category']))
ax.bar(x-0.2,df['Enrichment_Female'], 0.4,color='plum',label = 'Female')
ax.bar(x+0.2,df['Enrichment_Male'], 0.4,color='royalblue',label = 'Male')
plt.errorbar(x-0.2, df['Enrichment_Female'], yerr = (df['std_error_female']), linestyle="none", color="darkslategray", capsize=5, label = "2 SE")
plt.errorbar(x+0.2, df['Enrichment_Male'], yerr = (df['std_error_male']), linestyle="none", color="darkslategray", capsize=5, label = '_nolegend_')

plt.xticks(x,df['Category'],rotation='vertical',fontsize="large")
plt.yticks(fontsize="large")
plt.title("Partitioned Heritability for {0} - Male and Female - {1} - mash".format(pheno, result.capitalize()),fontsize="x-large")
plt.xlim(-1,x.max()+1)
plt.ylabel("Enrichment", fontsize="large")
plt.xlabel("Category", fontsize="large")

plt.legend()
plt.tight_layout()

plt.savefig(plot_name)