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

pheno = str(sys.argv[1])
sig_snp_thresh = int(sys.argv[2])
file_path = "/scratch1/08005/cz5959/GWAS_Results/{0}".format(pheno)
os.chdir(file_path)

# load and format posteriors to dataframe
df = pd.read_csv("{0}_posteriors_strong.txt".format(pheno),sep="\t")
split_df = df.index.str.split(":", n=3, expand = True)
df["#CHROM"] = pd.to_numeric(split_df[0], errors='coerce')
df["ID"] = pd.to_numeric(split_df[1], errors='coerce')
df["POS"] = pd.to_numeric(split_df[4], errors='coerce')
df = df.dropna()
print("length of df: " + str(len(df)))

results_df = results_df.sort_values(['#CHROM', 'POS'])
results_df.reset_index(inplace=True, drop=True)
results_df['P'] = pd.to_numeric(results_df['P'])
print(results_df[results_df['P'] == 0])
print(len(results_df[results_df['P'] == 0]))
results_df = results_df[(results_df['P'] != 0)]
results_df['NEG_LOG_P'] = -np.log10(results_df['P'])
results_df['LOG10'] = -results_df['NEG_LOG_P']