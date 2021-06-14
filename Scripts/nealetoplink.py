#!/usr/bin/env python3

# CONVERT NEALE RESULTS TO PLINK FORMAT
# to perform ldsc on neale results

# import packages
import matplotlib
matplotlib.use("Agg")
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os

# load in dataframes
os.chdir("/scratch1/08005/cz5959/Neale_Lab/Standing_Height")
neale_df = pd.read_csv("50_raw.gwas.imputed_v3.both_sexes.tsv", sep="\t",usecols=['variant','beta','pval'],dtype={'variant':str, 'beta':np.float64, 'pval':np.float64})