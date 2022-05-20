#!/usr/bin/env python3

import numpy as np
import pandas as pd
import sys
import os

os.chdir("/scratch1/08005/cz5959/Phenotypes")
pheno_list = ['height','testosterone','bmi','creatinine','IGF1','RBC_count',
'weight','calcium','protein_total','urea','SHBG','whole_body_fat_mass','FVC_best','HbA1c',
'albumin','arm_fatfree_mass_L','arm_fatfree_mass_R','diastolicBP_auto','systolicBP_auto','eosinophil_perc', 'lymphocyte_perc', 'hip_circ', 'waist_circ', 'waist_to_hip', 'pulse_rate', 'urate']
m_var_list = np.zeros(len(pheno_list))
f_var_list = np.zeros(len(pheno_list))
m_std_list = np.zeros(len(pheno_list))
f_std_list = np.zeros(len(pheno_list))

sex_df = pd.read_csv("sex_ids.txt", sep="\t")

for i in range(len(pheno_list)):
    pheno = pheno_list[i]
    pheno_df = pd.read_csv("pheno_{0}.txt".format(pheno), sep="\t")
    df = pd.merge(sex_df, pheno_df, how="inner", on="IID")
    male = df.loc[df['sex']==1,[pheno]]
    female = df.loc[df['sex']==0,[pheno]]
    m_var_list[i] = np.var(male)
    f_var_list[i] = np.var(female)
    m_std_list[i] = np.std(male)
    f_std_list[i] = np.std(female)

final_dict = {'pheno': pheno_list, 'm_var': m_var_list, 'm_std': m_std_list, 'f_var': f_var_list, 'f_std': f_std_list} 
final_df = pd.DataFrame(final_dict)
final_df.to_csv("pheno_variances_bysex.txt", index=False, sep="\t")

# genetic sex 21001; male=1 ; female=0