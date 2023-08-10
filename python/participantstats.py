# -*- coding: utf-8 -*-
"""
Participant Info: Descriptive Stats

@author: Prasanna S
"""

import os
import pandas as pd
import numpy as np
from scipy.stats import chisquare



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
sex_chisq = chisquare(sexfreq)

# Chi-square: Presentation
presfreq = pd.crosstab(index = datadf2["Group"], columns = datadf2["Presentation"]).transpose()
pres_chisq = chisquare(presfreq)

# T-Test vs Mann-Whitney U based on Shapiro-Wilk

