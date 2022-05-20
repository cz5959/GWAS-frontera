#!/usr/bin/env python3

import matplotlib
matplotlib.use("Agg")
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sys
import os
from sklearn.linear_model import LinearRegression

pheno=sys.argv[1]
file_path = "/scratch1/08005/cz5959/GWAS_Results/{0}".format(pheno)
os.chdir(file_path)

mash_psd = pd.read_csv("{0}_mash_psd.txt".format(pheno), sep="\t", usecols=['female','male'],dtype={'female':np.float64, 'male':np.float64})
#df_mash = pd.read_csv("height_male_mash_posterior.txt", sep="\t") 
#mash_pse = df_mash['BETA'] / df_mash['Z']

df_glm = pd.read_csv("female_all.{0}.glm.linear".format(pheno), sep="\t", usecols =['ID','A1','AX','OBS_CT','SE','P'],dtype= {'ID':str,'A1':str,'AX':str,'OBS_CT':np.int64,'SE':np.float64, 'P':np.float64})

# concatenate
df = pd.concat([df_glm, mash_psd], axis=1)

# p-value thresholds
p = [0,0.05,0.1,0.2,0.3,0.4,0.5,0.75,1]
index = 0

# plot
fig, axes = plt.subplots(2,4,figsize=(14, 7))
for ax in axes.reshape(-1):
    df_sub = df.loc[(df['P'] > p[index]) & (df['P'] <= p[index+1])]
    print(len(df_sub))

    ax.scatter(df_sub['SE'], df_sub['male'], s=1)
    y_max = df_sub['SE'].max()
    ax.set_ylim([0,y_max])
    ax.set_xlim([0,y_max])
    ax.set_title('({0},{1}]'.format(p[index],p[index+1]), size='large')
    index += 1

fig.text(0.5, 0.05, 'GWAS SE', ha='center', size='large')
fig.text(0.05, 0.5, 'mash psd', va='center', rotation='vertical', size='large')
plt.suptitle('mash psd vs GWAS SE - {0} - Male'.format(pheno.upper()), size='xx-large', y=0.95)
plt.subplots_adjust(bottom=0.12, hspace=0.3)

plt.savefig("{0}_psd_SE_male.png".format(pheno))

# linear regression
x = df['SE'].values.reshape(-1,1)
y = df['male']
model = LinearRegression().fit(x, y)

print('coefficient of determination:', model.score(x, y))
print('intercept:', model.intercept_)
print('slope:', model.coef_)
