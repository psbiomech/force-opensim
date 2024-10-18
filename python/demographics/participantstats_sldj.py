# -*- coding: utf-8 -*-
"""
Participant Info: Descriptive Stats

@author: Prasanna S
"""

import os
import pandas as pd
import numpy as np
from scipy.stats import chisquare
from scipy.stats import chi2_contingency
from scipy.stats import shapiro
from scipy.stats import ttest_ind
from scipy.stats import mannwhitneyu



# Data source
srcfile = r"C:\Users\Owner\Documents\data\FORCe\inputdatabase\force_sldj_participants_for_pub.csv"

# Results
#outfile = r"C:\Users\Owner\Documents\data\FORCe\outputdatabase_sldj\csvfolder\force_sldj_participant_descriptives.csv"



# %% DATA SOURCE


# Load data
datadf0 = pd.read_csv(srcfile)


# Recode
grouprc = {"Control":0, "FORCE": 1}
sexrc = {"Male":0, "Female":1}
affsiderc = {1:0, 2:0, 3:1}     # unilateral:0
datadf1 = datadf0.replace(dict(group = grouprc, sex = sexrc, aff_side = affsiderc))



# %% DESCRIPTIVES

# For simplicity, apply to all columns, but manually select the correct
# descriptives for publication/presentation.

# Mean/SD: Continuous variables
data_mean = datadf1.groupby("group").mean()
data_sd = datadf1.groupby("group").std()

# Median/IQR:
data_median = datadf1.groupby("group").median()
data_iqr = datadf1.groupby("group").quantile(0.75) - datadf1.groupby("group").quantile(0.25)



# %% STATISTICAL ANALYSES

# Drop string columns
datadf2 = datadf1.drop(columns = ["id"])

# Chi-square: Sex
sexfreq = pd.crosstab(index = datadf2["group"], columns = datadf2["sex"]).transpose()
sex_chisq = chi2_contingency(sexfreq, correction = False)

# Chi-square: Affected side
afffreq = pd.crosstab(index = datadf2["group"], columns = datadf2["aff_side"])
aff_chisq = chisquare(afffreq)

# T-Test vs Mann-Whitney U based on Shapiro-Wilk and manual boxplot inspection.
# Run both inference tests, select which to use and report manually based on
# Shapiro-Wilks test of normality.
datadf3 = datadf2.drop(columns = ["sex", "aff_side"]) 
datadf3grps = datadf3.groupby("group") 
shapw = {}
ttest = {}
mannwu = {}                                                                                       
for c, col in enumerate(datadf3.columns):
    
    # Data
    data0 = datadf3grps.get_group(0).loc[:, col].to_frame()
    data1 = datadf3grps.get_group(1).loc[:, col].to_frame()
    
    # Shapiro-Wilks: test normality
    shapw[col] = {}
    shapw[col][0] = shapiro(data0).pvalue
    shapw[col][1] = shapiro(data1).pvalue
    
    # T-test (Welch's)
    ttest[col] = ttest_ind(data0, data1, equal_var = False, nan_policy = "omit").pvalue
    
    # Mann-Whitney U
    mannwu[col] = mannwhitneyu(data0, data1, nan_policy = "omit").pvalue
    
    
    

