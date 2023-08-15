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
srcfile = r"C:\Users\Owner\Documents\data\FORCe\inputdatabase\force_sdp_participants_for_pub.csv"

# Results
outfile = r"C:\Users\Owner\Documents\data\FORCe\outputdatabase\csvfolder\force_sdp_participant_descriptives.csv"



# %% DATA SOURCE


# Load data
datadf0 = pd.read_csv(srcfile)


# Recode
grouprc = {"Control":0, "FORCE": 1}
sexrc = {"Male":0, "Female":1}
presrc = {"Unilateral Pain":0, "Bilateral Pain":1}
datadf1 = datadf0.replace(dict(Group = grouprc, Sex = sexrc, Presentation = presrc))



# %% DESCRIPTIVES

# For simplicity, apply to all columns, but manually select the correct
# descriptives for publication/presentation.

# Mean/SD: Continuous variables
data_mean = datadf1.groupby("Group").mean()
data_sd = datadf1.groupby("Group").std()

# Median/IQR:
data_median = datadf1.groupby("Group").median()
data_iqr = datadf1.groupby("Group").quantile(0.75) - datadf1.groupby("Group").quantile(0.25)



# %% STATISTICAL ANALYSES

# Drop string columns
datadf2 = datadf1.drop(columns = ["Participant", "Sport"])

# Chi-square: Sex
sexfreq = pd.crosstab(index = datadf2["Group"], columns = datadf2["Sex"]).transpose()
sex_chisq = chi2_contingency(sexfreq, correction = False)

# Chi-square: Presentation
presfreq = pd.crosstab(index = datadf2["Group"], columns = datadf2["Presentation"])
pres_chisq = chisquare(presfreq)

# T-Test vs Mann-Whitney U based on Shapiro-Wilk and manual boxplot inspection.
# Run both inference tests, select which to use and report manually based on
# Shapiro-Wilks test of normality.
datadf3 = datadf2.drop(columns = ["Sex", "Presentation", "Average_Pain_7days_NRS"]) 
datadf3grps = datadf3.groupby("Group") 
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
    
    
    

