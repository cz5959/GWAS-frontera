#!/usr/bin/env python3

# import packages
import matplotlib
matplotlib.use("Agg")
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sys
import os

pheno=sys.argv[1]
field_id=sys.argv[2]

# load file to dataframe
os.chdir("/scratch1/08005/cz5959/Neale_Lab/{0}".format(pheno))
results_df = pd.read_csv("{0}_raw.gwas.imputed_v3.both_sexes.tsv".format(field_id), sep="\t", usecols =["variant","pval"])

# drop rows with any column having null/missing data
results_df = results_df.dropna()

# split variant column into chromosome and variant
split_df = results_df["variant"].str.split(":", n=3, expand = True)
results_df["chrom"] = pd.to_numeric(split_df[0], errors='coerce')
results_df["variant"] = pd.to_numeric(split_df[1], errors='coerce')
results_df = results_df.dropna()

# sort by column then position; reset index
results_df = results_df.sort_values(['chrom', 'variant'])
results_df.reset_index(inplace=True, drop=True)

#change P column to float type; if there is a value smaller than it can handle, change it to float min
if results_df.pval.dtypes == object:
    results_df['pval'] = pd.to_numeric(results_df['pval'], errors='coerce')
    results_df = results_df.dropna()

# create column with negative log p value
print(results_df[results_df['pval'] == 0])
print(len(results_df[results_df['pval'] == 0]))
results_df = results_df[(results_df['pval'] != 0)]
results_df['NEG_LOG_P'] = -np.log10(results_df['pval'])

#p values that are greater than 5
#sig_P_5 = results_df.loc[results_df['NEG_LOG_P'] > 5]
#sig_P_5.to_csv("sig_P_5.csv",index=False, sep="\t")

#make chromosome column into type category
results_df['chrom'] = results_df['chrom'].astype('category')

#index; used for x axis; assume uniform SNP distrubtion across chromosome
results_df['index'] = range(len(results_df))

#group by chromosome
grouped_df = results_df.groupby(('chrom'))

#create Manhattan plot

fig = plt.figure(figsize=(18,14))
#axes of figure - 1row,1col,1idx
ax = fig.add_subplot(111)
colors = ['#466EA6','#7251B8']
x_labels = []
x_labels_pos = []

#create subplots for each chromosome (name = #CHROM)
for num, (name, group) in enumerate(grouped_df):
    ##### plot, x is index and y is neg log p ######
    group.plot(kind='scatter', x='index', y='NEG_LOG_P',color=colors[num % len(colors)], ax=ax, s=5 marker = '.')
    #name of chr
    x_labels.append(name)
    #tick marks; middle of group
    x_labels_pos.append((group['index'].iloc[-1] - (group['index'].iloc[-1] - group['index'].iloc[0])/2))

#line
ax.plot([0,len(results_df)],[5,5])

#figure labels
ax.set_xticks(x_labels_pos)
ax.set_xticklabels(x_labels)
y_max = results_df['NEG_LOG_P'].max()
print(y_max)
ax.set_ylim([0, y_max])
ax.set_xlabel('Chromosome')
ax.set_title('Manhattan Plot for Height - Neale')

# save as png
plt.tight_layout()
plt.savefig("neale_manhattan_height.png")
