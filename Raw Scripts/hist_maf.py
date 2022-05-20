#!/usr/bin/env python3

# import packages
import matplotlib
matplotlib.use("Agg")
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import sys

# read into dataframes
zero_df = pd.read_csv('zero_P_1000_MAF.txt', sep = "\t")
normal_df = pd.read_csv('chr22_MAF.txt', sep = "\t")

# create histogram
hist = plt.figure(figsize=(16,10))

zero = zero_df['ALT_FREQS']
normal = normal_df['ALT_FREQS']
max_ = max(zero.max(),normal.max())
plt.hist(zero, range=(0, max_), alpha = 0.5, label = 'zero', bins = 500)
plt.hist(normal, alpha = 0.5, label = 'normal', bins = 500)

plt.xlim(0,0.015)
plt.ylabel('frequency')
plt.xlabel('MAF')
plt.title('Histogram of MAF')
plt.legend()

plt.tight_layout()
plt.savefig("zoom2_hist_maf_chr22.png")
