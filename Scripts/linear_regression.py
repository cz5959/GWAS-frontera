#!/usr/bin/env python3

# import packages
import matplotlib
matplotlib.use("Agg")
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sys

# linear regression packages
from sklearn import linear_model
from sklearn.linear_model import LinearRegression
import statsmodels.api as sm
from scipy import stats

# read raw into dataframe
df = pd.read_csv('chr22_exportA.txt', sep = "\t")
#df.drop(['FID','IID','PAT','MAT','SEX'],axis=1, inplace=True)
#df = df[df['PHENOTYPE'] > 0]

# read covariates into dataframe

#plot as scatter
fig = plt.figure(figsize=(10,10))
for i in range(1,len(df.columns)):
    plt.scatter(df.iloc[:,[i]], df['PHENOTYPE'], alpha=0.2, color='red', s=0.1)
plt.xlabel("allelic dosage")
plt.ylabel("height")
plt.title("Chr22 Height Association")

# save plot as png
plt.tight_layout()
plt.savefig("linreg_height_chr22.png")

# linear regression
X = df.drop(['PHENOTYPE'], axis=1)
y = df['PHENOTYPE'].values.reshape(-1,1)
reg = LinearRegression()
reg.fit(X,y)

# summary of OLS
X1 = X
y1 = df['PHENOTYPE']
X2 = sm.add_constant(X1)
est = sm.OLS(y, X2)
est2 = est.fit()

# save summary as csv
file = open("lin_reg_summary_22.csv","w")
file.write(est2.summary().as_csv())
file.close()
