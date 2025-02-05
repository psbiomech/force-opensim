# -*- coding: utf-8 -*-
'''
FORCE SINGLE LEG DROP JUMP: WHOLE BODY ANGULAR MOMENTUM SEGMENTAL SPM{t}

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
outfilename = "force_sldj_spm1dt_wbam_segmental_merged_normalised_more_ctrl"


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
varprefix = ["L_seg_merged_arm", "L_seg_merged_thigh", "L_seg_merged_shank", "L_seg_merged_foot"]
genericvars = ["L_seg_pelvis", "L_seg_head_torso"] + varprefix
limb = ["stance", "swing"]
subjtype = ["sym", "ctrl"]
subjtypefulllabel = ["more symptomatic", "control"]
legtype = [["more"], ["more", "less"]]
planes = ["X", "Y", "Z"]

# Get data into arrays
datamat = {}
variables = {}
for f in limb:
    
    datamat[f] = {}  
    variables[f] = ["L_seg_pelvis", "L_seg_torso"] + [v + "_" + f for v in varprefix]
    
    for v in variables[f]:
        
        datamat[f][v] = {}
        
        for p in planes:
        
            datamat[f][v][p] = {}
            
            varcomponent = v + "_" + p
            
            # Group 0: symptomatic
            gp0data = df[~df["subject"].isin(outliers) & df["leg_type"].isin(legtype[0]) & (df["subj_type"] == subjtype[0]) & (df["variable"] == varcomponent)]
            datamat[f][v][p][0] = gp0data.loc[:, "t1":"t101"].to_numpy()
        
            # Group 1: Controls (each limb is considered independently)
            gp1data = df[~df["subject"].isin(outliers) & df["leg_type"].isin(legtype[1]) & (df["subj_type"] == subjtype[1]) & (df["variable"] == varcomponent)]
            datamat[f][v][p][1] = gp1data.loc[:, "t1":"t101"].to_numpy()
    
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
for f in limb:
    desc[f] = {}
    for v in variables[f]:
        desc[f][v] = {}
        for p in planes:
            desc[f][v][p] = {}
            for s in range(len(subjtype)):    
                desc[f][v][p][s] = {}                      
                desc[f][v][p][s]["mean"] = np.mean(datamat[f][v][p][s], axis = 0)
                desc[f][v][p][s]["sd"] = np.std(datamat[f][v][p][s], axis = 0)
            

# Run SPM{t} and inference across all legs, analyses, variables and group pairs
bonferroni = 1      # 0=no, 1=yes
significance = [0.05, 0.005]  # [treat variables as independent vs Bonferroni corrected for 10 comparisons]
spmt = {}
spmtinf = {}
for f in limb:
    spmt[f] = {}
    spmtinf[f] = {}
    for v in variables[f]:
        spmt[f][v] = {}
        spmtinf[f][v] = {}
        for p in planes:
            Y0 = datamat[f][v][p][0]
            Y1 = datamat[f][v][p][1]
            spmt[f][v][p] = spm1d.stats.ttest2(Y0, Y1, equal_var=False)
            spmtinf[f][v][p] = spmt[f][v][p].inference(alpha = significance[bonferroni], two_tailed=True, interp=True)



# %% PLOT OUTPUT


# Plot parameters
planesstr = ["Frontal (X)", "Transverse (Y)", "Sagittal (Z)"]
plotfont = {'fontname': 'Arial'}
units = "(dimensionless)"
subjtypeshortlabel = ["more sym", "ctrl"]
plotvars = [["L_seg_pelvis"] + [v + "_stance" for v in varprefix],
            ["L_seg_torso"] + [v + "_swing" for v in varprefix]]
titles = [["Pelvis", "Arm (Stance limb side)", "Thigh (Stance limb)", "Shank (Stance limb)", "Foot (Stance limb)"],
          ["Torso", "Arm (Swing limb side)", "Thigh (Swing limb)", "Shank (Swing limb)", "Foot (Swing limb)"]]

# Planes: X (Frontal), Y (Transverse), Z (Sagittal)
event0 = 48.5  # max knee flexion from IKID SPM script
for pn, p in enumerate(planes): 

    # Create plot area
    fig = plt.figure(constrained_layout=True, figsize=(20, 10))   
    fig.suptitle("Single-leg drop jump. Stance limb: %s vs %s. Contributions to WBAM. %s plane." % (subjtypefulllabel[0].upper(), subjtypefulllabel[1].upper(), planesstr[pn]), fontsize=20)
    heights = [2, 1, 0.5, 2, 1]
    spec = fig.add_gridspec(nrows = 5, ncols = len(plotvars[0]), height_ratios = heights) 

    # Plot results
    for s in range(len(subjtype)): 
    
        # Create plots
        x = range(101)
        for col in range(len(plotvars[0])):        
            
            # Mean + stdev
            for r, row in enumerate([0, 3]):        
                 
                f = limb[r]
                v = plotvars[r][col]
                               
                # Mean
                m0 = desc[f][v][p][0]["mean"]
                m1 = desc[f][v][p][1]["mean"]
                
                # Upper
                u0 = m0 + desc[f][v][p][0]["sd"]
                u1 = m1 + desc[f][v][p][1]["sd"]
                
                # Lower
                l0 = m0 - desc[f][v][p][0]["sd"]
                l1 = m1 - desc[f][v][p][1]["sd"]       
                
                # Plot
                ax = fig.add_subplot(spec[row, col])
                ax.set_title(titles[r][col], fontsize = 12)
                if (row == 0) and (col == 0):
                    ax.set_ylabel("L " + units, fontsize = 12)
                elif (row == 3) and (col == 0):
                    ax.set_ylabel("L " + units, fontsize = 12)   
                ax.fill_between(x, l1, u1, alpha = 0.3, linewidth = 0.0, color = "blue")
                ax.fill_between(x, l0, u0, alpha = 0.3, linewidth = 0.0, color = "red")
                ax.plot(x, m1, label = subjtypeshortlabel[1], linewidth = 2.0, color = "blue") 
                ax.plot(x, m0, label = subjtypeshortlabel[0], linewidth = 2.0, color = "red")
                ax.set_xlim([x[0], x[-1]])
                ax.axvline(x = event0, linewidth = 1.0, linestyle = ":", color = "k")
                if (row == 0 and col == 0): ax.legend(frameon = False, loc = "lower left")
            
                # Event labels
                #for at in range(6): ax.text((eventlist[at] / 100) + eventlabeladjust[at], 0.95, eventlabels[at], transform = ax.transAxes, horizontalalignment = eventlabelalign[at], fontsize = 8)
                
                # SPM significance shading
                issig = [1 if abs(z) > spmtinf[f][v][p].zstar else 0 for t, z in enumerate(spmtinf[f][v][p].z)]
                issigdiff = np.diff([0] + issig + [0])
                t0s = np.where(issigdiff == 1)
                t1s = np.where(issigdiff == -1)  # Should be the same length as t0s, I hope!
                if t0s[0].tolist():
                    for t in range(np.size(t0s[0])):
                        ax.axvspan(t0s[0][t], t1s[0][t], alpha = 0.3, color = "grey")
                
            # SPM inference
            for r, row in enumerate([1, 4]):  
                
                f = limb[r]
                v = plotvars[r][col]
                
                # plot
                ax = fig.add_subplot(spec[row, col])
                ax.set_xlabel("% of landing", fontsize = 12)
                if col == 0: ax.set_ylabel("SPM{t}", fontsize = 10) 
                ax.axvline(x = event0, linewidth = 1.0, linestyle = ":", color = "k")
                spmtinf[f][v][p].plot(plot_ylabel = False)
                        
                        
    # save to pdf
    plt.savefig(os.path.join(outpath, outfilename + "_" + p + ".pdf"))
