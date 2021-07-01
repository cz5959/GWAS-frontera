#!/usr/bin/env python3

# import packages
import matplotlib
matplotlib.use("Agg")
import numpy as np
import matplotlib.pyplot as plt
from sklearn.cluster import DBSCAN
import pandas as pd
import sys
import os
import time

def manhattan(pheno="phenotype", sex="both_sex", y_max = 0):
    #set variables
    file_name = "{1}_all_2.{0}.glm.linear".format(pheno, sex)
    plot_name = "manhattan_{0}_{1}.png".format(pheno, sex)
    plot_title = "Manhattan Plot of {0} : {1}".format(pheno.capitalize(), sex.capitalize())
    #plot_title = "Manhattan Plot of : {0}".format(sex.capitalize())
    sig_snp_thresh = int(sys.argv[2]) if sex == "both_sex" else int(sys.argv[3])

    # DATAFRAME
    results_df = pd.read_csv(file_name, sep="\t", usecols =['#CHROM','POS','P','ID'], dtype= {'#CHROM':np.int64,'POS':np.int64,'P':np.float64,'ID':str})

    # drop rows with any column having null/missing data ; should already be removed when combine results
    results_df = results_df.dropna()
    print("length of df: " + str(len(results_df)))
    print(results_df.dtypes)

    # sort by column then position; reset index
    results_df = results_df.sort_values(['#CHROM', 'POS'])
    results_df.reset_index(inplace=True, drop=True)

    # change P column to float type, remove p=0, and create column with negative log p value
    results_df['P'] = pd.to_numeric(results_df['P'])
    print(results_df[results_df['P'] == 0])
    print(len(results_df[results_df['P'] == 0]))
    results_df = results_df[(results_df['P'] != 0)]
    results_df['NEG_LOG_P'] = -np.log10(results_df['P'])
    results_df['LOG10'] = -results_df['NEG_LOG_P']

    # make chromosome column into type category and create index col
    results_df['#CHROM'] = results_df['#CHROM'].astype('category')
    results_df['index'] = range(len(results_df))

    # group by chromosome
    grouped_df = results_df.groupby(('#CHROM'))

    #PLOT
    fig = plt.figure(figsize=(18,14))
    ax = fig.add_subplot(111)
    colors = ['#466EA6','#7251B8']
    x_labels = []
    x_labels_pos = []
    snp_list = pd.DataFrame(columns=['index','NEG_LOG_P','LOG10','ID'])
    above_thresh = pd.DataFrame(columns=['#CHROM','ID','NEG_LOG_P'])

    # create subplots for each chromosome (name = #CHROM)
    for num, (name, group) in enumerate(grouped_df):
        ##### plot, x is index and y is neg log p ######
        group.plot(kind='scatter', x='index', y='NEG_LOG_P',color=colors[num % len(colors)], ax=ax, s=5, marker = '.')
        # x-axis labels and ticks
        x_labels.append(name)
        x_labels_pos.append((group['index'].iloc[-1] - (group['index'].iloc[-1] - group['index'].iloc[0])/2))

        ##### annotation #######
        top_group = group.loc[group['NEG_LOG_P'] > sig_snp_thresh].copy(deep=True)
        # cluster and label top SNP
        if len(top_group) > 0:
            above_thresh = pd.concat([above_thresh,top_group['#CHROM','ID','NEG_LOG_P']])
            clustering = DBSCAN(eps=500, min_samples=0).fit(top_group[['index','NEG_LOG_P']])
            top_group.reset_index(inplace=True)
            top_group['label'] = pd.DataFrame(clustering.labels_)
            for group in top_group['label'].unique():
                y = top_group.loc[top_group['label'] == group, ['NEG_LOG_P']]
                topSNP = top_group.iloc[y.idxmax()]
                snp_list = pd.concat([snp_list,topSNP[['index','NEG_LOG_P','LOG10','ID']]])
    above_thresh.to_csv("{0}_{1}_SNPsAboveThresh.txt".format(pheno, sex), sep="\t", index=False)
    anno = get_annotation(snp_list, pheno, sex)
    for i, row in anno.iterrows():
        text = row['SYMBOL'] if row['SYMBOL'] != "-" else row['ID']
        ax.text(row['index'], row['NEG_LOG_P'], text)
    # sig line
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
        print("ymax= " + str(y_max))
        ax.set_ylim([0, y_max])
        plt.tight_layout()
        plt.savefig(plot_name)
        return y_max
    else:
        ax.set_ylim([0, y_max])
        print("ymax= " + str(results_df['NEG_LOG_P'].max()))
        # save as png
        plt.tight_layout()
        plt.savefig(plot_name)
    return results_df, anno

def manhattan_combine(pheno, female_df, male_df, female_anno, male_anno):
    female_grouped = female_df.groupby(('#CHROM'))
    male_grouped = male_df.groupby(('#CHROM'))
    plot_name = "manhattan_{0}_combined.png".format(pheno)
    plot_title = "Manhattan Plot of {0} : Combined Sex".format(pheno.capitalize())

    fig = plt.figure(figsize=(18,16))
    ax1 = fig.add_subplot(211)
    ax2 = fig.add_subplot(212, sharex = ax1)
    colors = ['#466EA6','#7251B8']
    x_labels = []
    x_labels_pos = []

    for num, (name, group) in enumerate(female_grouped):
        #ax1.axvspan(group['index'].min(),group['index'].max(), ec=None, fc=colors[num%len(colors)], alpha = 0.3)
        group.plot(kind='scatter', x='index', y='NEG_LOG_P',color=colors[num % len(colors)], ax=ax1, s=5, marker = '.')
        x_labels.append(name)
        x_labels_pos.append((group['index'].iloc[-1] - (group['index'].iloc[-1] - group['index'].iloc[0])/2))
    for num, (name, group) in enumerate(male_grouped):
        #ax2.axvspan(group['index'].min(),group['index'].max(), ec=None, fc=colors[num%len(colors)], alpha = 0.3)
        group.plot(kind='scatter', x='index', y='LOG10',color=colors[num % len(colors)], ax=ax2, s=5, marker = '.')

    # annotations
    for i, row in female_anno.iterrows():
        text = row['SYMBOL'] if row['SYMBOL'] != "-" else row['ID']
        ax1.text(row['index'], row['NEG_LOG_P'], text)
    for i, row in male_anno.iterrows():
        text = row['SYMBOL'] if row['SYMBOL'] != "-" else row['ID']
        ax2.text(row['index'], row['LOG10'], text)
    # x
    x_max = max(len(female_df), len(male_df))
    ax1.plot([0,x_max],[5,5])
    ax1.set_xticks(x_labels_pos)
    ax1.set_xticklabels(x_labels)
    ax1.xaxis.set_tick_params(labelbottom=True)
    ax1.set_xlim([0, x_max])
    ax2.set_xlabel('Chromosome')
    # y
    y_max = max(female_df['NEG_LOG_P'].max(), abs(male_df['LOG10'].min()))
    print(y_max)
    ax1.set_ylim(0, y_max)
    ax1.set_ylabel('Female - NEG LOG10 P')
    ax2.set_ylabel('Male - LOG10 P')
    ax2.set_ylim(-y_max,0)
    ax1.set_title(plot_title)

    plt.tight_layout()
    plt.savefig(plot_name)

def get_annotation(snp_list, pheno, sex):
    # create file with list of rsids
    file_name = "{0}_{1}_topSNPs.txt".format(pheno,sex)
    snp_file = open(file_name,"w")
    for i in snp_list['ID']:
        snp_file.write(str(i) + "\n")
    snp_file.close()
    #ensembl VEP
    os.chdir("/scratch1/08005/cz5959/Scripts")
    os.system("pwd") 
    cmd = "./get_annotation.sh {0} {1}".format(pheno,sex)
    os.system(cmd)
    #load and return results
    file_path = "/scratch1/08005/cz5959/GWAS_Results/{0}".format(str(sys.argv[1]))
    os.chdir(file_path)
    file_name = "{0}_{1}_annotation.txt".format(pheno,sex)
    # check if results exist
    timeout = time.time() + 60*5 #timeout in 5 min
    while os.path.isfile(file_name) == False:
        time.sleep(5)
        if time.time() > timeout:
            print("timeout")
            quit()
    # load results
    anno = pd.read_csv(file_name, sep="\t", usecols=[0,3], header=31) 
    anno.rename(columns={"#Uploaded_variation":"ID"},inplace=True)
    anno.drop_duplicates(subset=["ID"],inplace=True)
    anno = pd.merge(snp_list, anno, how="left", on="ID")
    anno.dropna(inplace=True)
    print(anno)
    print(anno.dtypes)
    return anno        
    

### MAIN ####
# set working directory
file_path = "/scratch1/08005/cz5959/GWAS_Results/{0}".format(str(sys.argv[1]))
os.chdir(file_path)
    
# sys.argv[1] should be phenotype name
y_max = manhattan(str(sys.argv[1]),"both_sex")
female_df, female_anno = manhattan(str(sys.argv[1]),"female", y_max)
male_df, male_anno = manhattan(str(sys.argv[1]),"male", y_max)
manhattan_combine(str(sys.argv[1]), female_df, male_df, female_anno, male_anno)

