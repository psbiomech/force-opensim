# -*- coding: utf-8 -*-
'''
FORCE SINGLE LEG DROP JUMP: WHOLE BODY ANGULAR MOMENTUM SPM{t}

@author: Prasanna Sritharan, August 2023
'''


import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import spm1d
#import pickle as pk

# Matplotlib: write text as text not path
from matplotlib import rc
rc("pdf", fonttype=42)
rc("font", **{'family':'sans-serif','sans-serif':['Arial']})


# Data file
srcpath = r"C:\Users\Owner\Documents\data\FORCe\outputdatabase_sldj\csvfolder"
srcfile = "force_sldj_results_subject_descriptives_stability_normalised.csv"

# Output file
outpath = r"C:\Users\Owner\Documents\data\FORCe\outputdatabase_sldj\spm1d"
if not os.path.isdir(outpath): os.makedirs(outpath)
outfilename = "force_sldj_spm1dt_wbam_normalised_more_ctrl"


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
variables = ["L_X", "L_Y", "L_Z"]
subjtype = ["sym", "ctrl"]
subjtypefulllabel = ["more symptomatic", "control"]
legtype = [["more"], ["more", "less"]]

# Get data into arrays
datamat = {}
for v in ["time"] + variables:
    
    datamat[v] = {}
    
    # Group 0: symptomatic
    gp0data = df[~df["subject"].isin(outliers) & df["leg_type"].isin(legtype[0]) & (df["subj_type"] == subjtype[0]) & (df["variable"] == v)]
    datamat[v][0] = gp0data.loc[:, "t1":"t101"].to_numpy()

    # Group 1: Controls (each limb is considered independently)
    gp1data = df[~df["subject"].isin(outliers) & df["leg_type"].isin(legtype[1]) & (df["subj_type"] == subjtype[1]) & (df["variable"] == v)]
    datamat[v][1] = gp1data.loc[:, "t1":"t101"].to_numpy()

    # # Group 1: Controls (mean of both limbs)
    # # A single participant mean for each control may not be appropriate for this
    # # study as only one limb is condsidered for each symptomatic. Need to get
    # # clarification.
    # gp1data = df[~df["subject"].isin(outliers) & df["leg_type"].isin(legtype[1]) & (df["subj_type"] == subjtype[1]) & (df["variable"] == v)]
    # gp1means = gp1data.groupby("subject").mean()
    # datamat[v][1] = gp1means.loc[:, "t1":"t101"].to_numpy()




# %% RUN ANALYSES: DESCRIPTIVES, SPM{t}

# Calculate group ensemble descriptives from file
desc = {}
for v in variables:
    desc[v] = {}
    for s in range(len(subjtype)):    
        desc[v][s] = {}                      
        desc[v][s]["mean"] = np.mean(datamat[v][s], axis = 0)
        desc[v][s]["sd"] = np.std(datamat[v][s], axis = 0)
        

# Run SPM{t} and inference across all legs, analyses, variables and group pairs
spmt = {}
spmtinf = {}
for v in variables:
    Y0 = datamat[v][0]
    Y1 = datamat[v][1]
    spmt[v] = spm1d.stats.ttest2(Y0, Y1, equal_var=False)
    spmtinf[v] = spmt[v].inference(alpha = 0.05, two_tailed=True, interp=True)



# %% PLOT OUTPUT


# Plot parameters
plotheads = ["Frontal (X)", "Transverse (Y)", "Sagittal (Z)"]
plotfont = {'fontname': 'Arial'}
units = "(dimensionless)"


# Create plot area
fig = plt.figure(constrained_layout=True, figsize=(12, 5))   
fig.suptitle("Single-leg drop-jump: %s vs %s. Whole body angular momentum." % (subjtypefulllabel[0].upper(), subjtypefulllabel[1].upper()), fontsize=20)
spec = fig.add_gridspec(nrows = 2, ncols = len(variables), height_ratios = [2, 1]) 


# Create plots
x = range(101)
event0 = 48.5  # from IKID SPM script
for col in range(len(variables)):                
    
    v = variables[col]
     
    # Mean
    m0 = desc[variables[col]][0]["mean"]
    m1 = desc[variables[col]][1]["mean"]
    
    # Upper
    u0 = m0 + desc[variables[col]][0]["sd"]
    u1 = m1 + desc[variables[col]][1]["sd"]
    
    # Lower
    l0 = m0 - desc[variables[col]][0]["sd"]
    l1 = m1 - desc[variables[col]][1]["sd"]       
    
    # Plot
    ax = fig.add_subplot(spec[0, col])
    ax.set_title(plotheads[col], fontsize = 12)
    if col == 0:
        ax.set_ylabel("L " + units, fontsize = 12) 
    ax.fill_between(x, l1, u1, alpha = 0.3, linewidth = 0.0, color = "blue")
    ax.fill_between(x, l0, u0, alpha = 0.3, linewidth = 0.0, color = "red")
    ax.plot(x, m1, label = subjtypefulllabel[1], linewidth = 2.0, color = "blue") 
    ax.plot(x, m0, label = subjtypefulllabel[0], linewidth = 2.0, color = "red")
    ax.set_xlim([x[0], x[-1]])
    #ax.set_xlabel("% of stance", fontsize = 12)
    ax.axvline(x = event0, linewidth = 1.0, linestyle = ":", color = "k")
    if col == 0: ax.legend(frameon = False, loc = "lower left")

    # SPM significance shading on variable plot
    issig = [1 if abs(z) > spmtinf[v].zstar else 0 for t, z in enumerate(spmtinf[v].z)]
    issigdiff = np.diff([0] + issig + [0])
    t0s = np.where(issigdiff == 1)
    t1s = np.where(issigdiff == -1)  # Should be the same length as t0s, I hope!
    if t0s[0].tolist():
        for t in range(np.size(t0s[0])):
            ax.axvspan(t0s[0][t], t1s[0][t], alpha = 0.3, color = "grey")
    
    # SPM plot
    ax = fig.add_subplot(spec[1, col])
    ax.set_xlabel("% of stance", fontsize = 12)
    ax.axvline(x = event0, linewidth = 1.0, linestyle = ":", color = "k")
    if col == 0: ax.set_ylabel("SPM{t}", fontsize = 10) 
    spmtinf[variables[col]].plot(plot_ylabel = False)


# Save to pdf
plt.savefig(os.path.join(outpath, outfilename + ".pdf"))
plt.close(fig)


