#!/usr/bin/env python3

# import packages
import matplotlib
matplotlib.use("Agg")
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import sys
import os

# set working directory
os.chdir("/scratch1/08005/cz5959/Association_Height_50")

# DATAFRAME

# load file to dataframe
results_df = pd.read_csv("linear_results_all_chrom.height.glm.linear", sep="\t", usecols =['#CHROM','POS','P','ID'])

# drop rows with any column having null/missing data
results_df = results_df.dropna()

# sort by column then position; reset index
results_df = results_df.sort_values(['#CHROM', 'POS'])
results_df.reset_index(inplace=True, drop=True)

# change P column to float type; if there is a value smaller than it can handle, change it to float min
results_df['P'] = pd.to_numeric(results_df['P'], errors='coerce')
results_df['P'] = pd.to_numeric(results_df['P'])
# min p values
#min_P = results_df[results_df['P'].isna()]
#min_P.to_csv("min_P_22.csv",index=False, sep="\t")
#results_df['P'] = results_df['P'].fillna(sys.float_info.min)

# p values that are equal to 0
#zero_P = results_df.loc[results_df['P'] == 0]
#zero_P.to_csv("zero_P_22.csv",index=False, sep="\t")

# create column with negative log p value
# ignore divide error message because fixed below
np.seterr(divide = 'ignore') 
# if p=0, make dummy -2
results_df['NEG_LOG_P']= neg_log_values = np.where(results_df['P'] != 0, -np.log10(results_df['P']), -2)

# neg log p values that are greater than 5
#sig_P_5 = results_df.loc[results_df['NEG_LOG_P'] > 5]
# create csv
#sig_P_5.to_csv("sig_P_5.csv", index=False, sep="\t")

# make chromosome column into type category
results_df['#CHROM'] = results_df['#CHROM'].astype('category')

# index; used for x axis; assume uniform SNP distrubtion across chromosome
results_df['index'] = range(len(results_df))

# group by chromosome
grouped_df = results_df.groupby(('#CHROM'))

#PLOT

fig = plt.figure(figsize=(16,6))
# axes of figure - 1row,1col,1idx
ax = fig.add_subplot(111)
colors = ['#466EA6','#7251B8']
x_labels = []
x_labels_pos = []

# create subplots for each chromosome (name = #CHROM)
for num, (name, group) in enumerate(grouped_df):
    ##### plot, x is index and y is neg log p ######
    group.plot(kind='scatter', x='index', y='NEG_LOG_P',color=colors[num % len(colors)], ax=ax, s=0.1, marker = '.')
    # name of chr
    x_labels.append(name)
    # tick marks; middle of group
    x_labels_pos.append((group['index'].iloc[-1] - (group['index'].iloc[-1] - group['index'].iloc[0])/2))

#line
ax.plot([0,len(results_df)],[5,5])

# figure labels
ax.set_xticks(x_labels_pos)
ax.set_xticklabels(x_labels)
ax.set_xlim([0, len(results_df)])
ax.set_ylim([0, 10])
ax.set_xlabel('Chromosome')
ax.set_title('Manhattan Plot')

# save as png
plt.tight_layout()
plt.savefig("manhattan_height.png")
