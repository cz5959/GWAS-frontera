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
result = sys.argv[2]
plot_name="partition_mash_{0}_{1}_both.png".format(pheno, result)
print(pheno, result)

os.chdir("/scratch1/08005/cz5959/LD_practice/{0}".format(pheno))
df_female = pd.read_csv("{0}_female_mash_{1}.results".format(pheno, result),sep="\t",usecols=['Category','Enrichment','Enrichment_std_error'],dtype={'Category':str,'Enrichment':np.float64,'Enrichment_std_error':np.float64})
df_male = pd.read_csv("{0}_male_mash_{1}.results".format(pheno, result),sep="\t",usecols=['Category','Enrichment','Enrichment_std_error'],dtype={'Category':str,'Enrichment':np.float64,'Enrichment_std_error':np.float64})

df1.rename(columns={'Enrichment':'Enrichment_Female', 'Enrichment_std_error':'std_error_female'},inplace=True)
df2.rename(columns={'Enrichment':'Enrichment_Male', 'Enrichment_std_error':'std_error_male'},inplace=True)
df = pd.merge(df1, df2, how="outer", on="Category")

x = np.arange(len(df['Category']))
ax.bar(x-0.2,df['Enrichment_Female'], 0.4,color='plum',label = 'Female')
ax.bar(x+0.2,df['Enrichment_Male'], 0.4,color='royalblue',label = 'Male')
plt.errorbar(x-0.2, df['Enrichment_Female'], yerr = (df['std_error_female']), linestyle="none", color="darkslategray", capsize=5, label = "2 SE")
plt.errorbar(x+0.2, df['Enrichment_Male'], yerr = (df['std_error_male']), linestyle="none", color="darkslategray", capsize=5, label = '_nolegend_')

plt.xticks(x,df['Category'],rotation='vertical',fontsize="large")
plt.yticks(fontsize="large")
plt.title("Partitioned Heritability for {0} - Male and Female - {1}".format(pheno, result.capitalize()),fontsize="large")
plt.xlim(-1,x.max()+1)
plt.ylabel("Enrichment", labelsize="large")
plt.xlabel("Category", labelsize="large")

plt.legend()
plt.tight_layout()
plt.savefig(plot_name)