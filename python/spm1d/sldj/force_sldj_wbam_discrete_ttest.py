# -*- coding: utf-8 -*-
'''
FORCE SINGLE LEG DROP JUMP: MARGIN OF STABILITY T-TEST

@author: Prasanna Sritharan, August 2023
'''


import os
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib.pyplot as plt

from matplotlib import rc
rc("pdf", fonttype=42)
rc("font", **{'family':'sans-serif','sans-serif':['Arial']})


# Data file
srcpath = r"C:\Users\Owner\Documents\data\FORCe\outputdatabase_sldj\csvfolder"
srcfile = "force_sldj_results_subject_descriptives_wbam_discrete.csv"

# Output file
outpath = r"C:\Users\Owner\Documents\data\FORCe\outputdatabase_sldj\spm1d"
if not os.path.isdir(outpath): os.makedirs(outpath)
outfilename = "force_sldj_ttest_wbam_discrete_more_ctrl"


# Outliers
outliers = []
# outliers = ["FAILTCRT08", "FAILTCRT13", "FAILTCRT28",  "FAILTCRT26", "FAILTCRT31",
#             "FAILT69", "FAILT97", "FAILT146", "FAILT117"]


# %% PREPARE DATA

# Load XLS data into a dataframe
df0 = pd.read_csv(os.path.join(srcpath, srcfile))

# Select mean rows
df = df0[df0["statistic"] == "mean"]


# Data labels
variables = ["L_int_X", "L_int_Y", "L_int_Z", "L_range_X", "L_range_Y", "L_range_Z",]
subjtype = ["sym", "ctrl"]
subjtypefulllabel = ["more symptomatic", "control"]
legtype = [["more"], ["more", "less"]]

# Get data into arrays
datamat = {}
for v in variables:
    
    datamat[v] = {}
    
    # Group 0: symptomatic
    gp0data = df[~df["subject"].isin(outliers) & df["leg_type"].isin(legtype[0]) & (df["subj_type"] == subjtype[0]) & (df["variable"] == v)]
    datamat[v][0] = gp0data.loc[:, "value"].to_numpy()

    # Group 1: Controls (each limb is considered independently)
    gp1data = df[~df["subject"].isin(outliers) & df["leg_type"].isin(legtype[1]) & (df["subj_type"] == subjtype[1]) & (df["variable"] == v)]
    datamat[v][1] = gp1data.loc[:, "value"].to_numpy()


# %% RUN ANALYSES: DESCRIPTIVES, T-TESTS

# Calculate group ensemble descriptives from file
desc = {}
for v in variables:
    desc[v] = {}
    for s in range(len(subjtype)):    
        desc[v][s] = {}                      
        desc[v][s]["mean"] = np.mean(datamat[v][s], axis = 0)
        desc[v][s]["sd"] = np.std(datamat[v][s], axis = 0)
        

# Run t-tests and inference across all legs, analyses, variables and group pairs
ttest = {}
ttestinf = {}
alpha = 0.05
for v in variables:
    Y0 = datamat[v][0]
    Y1 = datamat[v][1]
    ttest[v] = stats.ttest_ind(Y0, Y1, equal_var=False)
    ttestinf[v] = ttest[v].pvalue < alpha


# %% OUTPUT TABLE


# Create output table
csvdata = []
for v in variables:
    csvrow = [v]
    for s in range(len(subjtype)):
        csvrow.append(desc[v][s]["mean"])
        csvrow.append(desc[v][s]["sd"])
    csvrow.append(ttest[v].pvalue)
    csvdata.append(csvrow)
    
# Create dataframe
headers = ["variable", "sym_mean", "sym_sd", "ctrl_mean", "ctrl_sd", "p"]
csvdf = pd.DataFrame(csvdata, columns = headers)

# Save to CSV
csvdf.to_csv(os.path.join(outpath, outfilename + ".csv"), index = False)