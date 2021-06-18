#!/usr/bin/env python3

# import packages
import matplotlib
matplotlib.use("Agg")
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import sys
import os

def manhattan(pheno="phenotype", sex="both_sex", y_max = 0):
    #set variables
    file_name = "{1}_all.{0}.glm.linear".format(pheno, sex)
    plot_name = "manhattan_{0}_{1}.png".format(pheno, sex)
    plot_title = "Manhattan Plot of {0} : {1}".format(pheno.capitalize(), sex.capitalize)

    # DATAFRAME
    # load file to dataframe
    results_df = pd.read_csv(file_name, sep="\t", usecols =['#CHROM','POS','P','ID'], dtype= {'#CHROM':np.int64,'POS':np.int64,'P':np.float64,'ID':str})

    # drop rows with any column having null/missing data ; should already be removed when combine results
    results_df = results_df.dropna()

    # sort by column then position; reset index
    results_df = results_df.sort_values(['#CHROM', 'POS'])
    results_df.reset_index(inplace=True, drop=True)

    # change P column to float type
    results_df['P'] = pd.to_numeric(results_df['P'])

    # create column with negative log p value
    print(results_df[results_df['P'] == 0])
    results_df = results_df[(results_df['P'] != 0)]
    results_df['NEG_LOG_P'] = -np.log10(results_df['P'])

    # make chromosome column into type category
    results_df['#CHROM'] = results_df['#CHROM'].astype('category')

    # index; used for x axis; assume uniform SNP distrubtion across chromosome
    results_df['index'] = range(len(results_df))

    # group by chromosome
    grouped_df = results_df.groupby(('#CHROM'))

    #PLOT

    fig = plt.figure(figsize=(18,14))
    # axes of figure - 1row,1col,1idx
    ax = fig.add_subplot(111)
    colors = ['#466EA6','#7251B8']
    x_labels = []
    x_labels_pos = []

    # create subplots for each chromosome (name = #CHROM)
    for num, (name, group) in enumerate(grouped_df):
        ##### plot, x is index and y is neg log p ######
        group.plot(kind='scatter', x='index', y='NEG_LOG_P',color=colors[num % len(colors)], ax=ax, s=5, marker = '.')
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
    ax.set_xlabel('Chromosome')
    ax.set_title(plot_title)

    #y-axis, keep as max of both_sex
    if sex == "both_sex":
        y_max = results_df['NEG_LOG_P'].max()
        print(y_max)
        ax.set_ylim([0, y_max])
        plt.tight_layout()
        plt.savefig(plot_name)
        return y_max
    else:
        ax.set_ylim([0, y_max])
        print(results_df['NEG_LOG_P'].max())
        # save as png
        plt.tight_layout()
        plt.savefig(plot_name)


# set working directory
os.chdir("/scratch1/08005/cz5959/GWAS_Results")
    
# sys.argv[1] should be phenotype name
y_max = manhattan(str(sys.argv[1]),"both_sex")
manhattan(str(sys.argv[1]),"female", y_max)
manhattan(str(sys.argv[1]),"male", y_max)

