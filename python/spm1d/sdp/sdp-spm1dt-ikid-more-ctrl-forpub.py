# -*- coding: utf-8 -*-
"""
STEP DOWN AND PIVOT: MORE SYMPTOMATIC VS CTRL - SPM1D T-TEST (FOR PUBLICATION)
@author: Prasanna Sritharan
"""


import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import spm1d
#import pickle as pk

from matplotlib import rc
rc("pdf", fonttype=42)
rc("font", **{'family':'sans-serif','sans-serif':['Arial']})

# Data file
srcpath = r"C:\Users\Owner\Documents\data\FORCe\outputdatabase\csvfolder"
srcfile = "force_sdp_results_subject_descriptives_normalised.csv"

# Output file
outpath = r"C:\Users\Owner\Documents\data\FORCe\outputdatabase\spm1d"
if not os.path.isdir(outpath): os.makedirs(outpath)
outfilename = "sdp-spm1dt-ikid-more-ctrl"


# %% PREPARE DATA

# Load XLS data into a dataframe
df0 = pd.read_csv(os.path.join(srcpath, srcfile))

# Remove stdev rows
df1 = df0[df0["statistic"] == "mean"]

# Data labels
analyses = ["ik", "id"]
osimvars = {}
osimvars["ik"] = ["hip_flexion", "hip_adduction", "hip_rotation", "knee_angle", "ankle_angle", "lumbar_extension", "lumbar_bending", "lumbar_rotation"]
osimvars["id"] = [k + "_moment" for k in osimvars["ik"]]
subjtype = ["sym", "ctrl"]
subjtypefulllabel = ["more symptomatic", "control"]

# Scenario labels
scenario = ["pivot", "nonpivot"]
trialcombo =[["pivot_more", "pivot"], ["pivot_less", "pivot"]]

# Get data into arrays
datamat = {}
eventmat = {}
for sn, n in enumerate(scenario):
    df = df1[df1["data_leg_role"] == n]
    datamat[n] = {}
    eventmat[n] = {}
    for a in analyses:
        datamat[n][a] = {}
        eventmat[n][a] = {}
        for v in ["time"] + osimvars[a]:
            
            datamat[n][a][v] = {}
            eventmat[n][a][v] = {}        
            
            # Group 1
            gp1data = df[(df["subj_type"] == subjtype[0]) & (df["trial_combo"] == trialcombo[sn][0]) & (df["analysis"] == a) & (df["variable"] == v)]
            datamat[n][a][v][0] = gp1data.loc[:, "t1":"t101"].to_numpy()
            eventmat[n][a][v][0] = gp1data.loc[:, "es1_PFO1":"es6_PFS4"].to_numpy()
        
            # Group 2
            gp2data = df[(df["subj_type"] == subjtype[1]) & (df["trial_combo"].str.contains(trialcombo[sn][1])) & (df["analysis"] == a) & (df["variable"] == v)]
            datamat[n][a][v][1] = gp2data.loc[:, "t1":"t101"].to_numpy()
            eventmat[n][a][v][1] = gp2data.loc[:, "es1_PFO1":"es6_PFS4"].to_numpy()



# %% EVENTS

# Event time steps and descriptives
events = {}
for n in scenario:
    events[n] = {}
    events[n]["data"] = {}
    events[n]["desc"] = {}
    descmat = np.zeros((2,6))
    for g in range(len(subjtype)):
        events[n]["data"][g] = eventmat[n]["ik"]["time"][g]
        events[n]["desc"][g] = {}
        events[n]["desc"][g]["mean"] = np.round(np.mean(events[n]["data"][g], axis=0))
        events[n]["desc"][g]["sd"] = np.round(np.std(events[n]["data"][g], axis=0))
        descmat[g, 0:6] = events[n]["desc"][g]["mean"]
    events[n]["desc"]["total"] = {}
    events[n]["desc"]["total"]["mean"] = np.mean(descmat, axis=0)
    events[n]["desc"]["total"]["sd"] = np.std(descmat, axis=0)



# %% RUN ANALYSES: DESCRIPTIVES, SPM{t}

# Calculate group ensemble descriptives from file
desc = {}
for n in scenario:
    desc[n] = {}
    for a in analyses:        
        desc[n][a] = {}
        for v in osimvars[a]:
            desc[n][a][v] = {}
            for s in range(len(subjtype)):    
                desc[n][a][v][s] = {}                      
                desc[n][a][v][s]["mean"] = np.mean(datamat[n][a][v][s], axis = 0)
                desc[n][a][v][s]["sd"] = np.std(datamat[n][a][v][s], axis = 0)


# Run SPM{t} and inference across all legs, analyses, variables and group pairs
bonferroni = 1      # 0=no, 1=yes
significance = [0.05, 0.03125]  # [treat variables as independent vs Bonferroni corrected for 16 comparisons]
spmt = {}
spmtinf = {}
for n in scenario:
    spmt[n] = {}
    spmtinf[n] = {}   
    for a in analyses: 
        spmt[n][a] = {}
        spmtinf[n][a] = {}
        for v in osimvars[a]:
            Y0 = datamat[n][a][v][0]
            Y1 = datamat[n][a][v][1]
            spmt[n][a][v] = spm1d.stats.ttest2(Y0, Y1, equal_var=False)
            spmtinf[n][a][v] = spmt[n][a][v].inference(alpha = significance[bonferroni], two_tailed=True, interp=True)


# Combine for output
sdp = {}
sdp["desc"] = desc
sdp["events"] = events
sdp["spmt"] = spmt
sdp["spmtinf"] = spmtinf
            
# Pickle it
#with open(os.path.join(outpath, outfilename + "-forpub.pkl"),"wb") as f: pk.dump(sdp, f)
    

# %% PLOT OUTPUT

# Plot parameters
eventlabels = ["PFO1", "PFS2", "NFO1", "NFS2", "PFO3", "PFS4"]
eventlabelalign = ["left", "right", "left", "right", "left", "right"]
eventlabeladjust = [0.01, -0.01, 0.01, -0.01, 0.01, -0.01]   
plotfont = {'fontname': 'Arial'}
bonftext = ["", "-bonferroni"]

# Figure headers
plotheads = {}
plotheads["ik"] = ["Hip flexion", "Hip adduction", "Hip rotation", "Knee flexion", "Ankle dorsiflexion", "Lumbar extension", "Lumbar bending", "Lumbar rotation"]
plotheads["id"] = [k + " moment" for k in plotheads["ik"]]


# Generate plots
for n in scenario:

    # Scenario parameters
    nsubjs = [np.size(datamat[n]["ik"]["time"][s], axis=0) for s in range(len(subjtype))]
    eventlist = 100 * np.round(events[n]["desc"]["total"]["mean"]) / 101    

    # Create plot area
    fig = plt.figure(constrained_layout=True, figsize=(18, 6))   
    fig.suptitle("Step-down-and-pivot: %s vs %s - %s limb" % (subjtypefulllabel[0].upper(), subjtypefulllabel[1].upper(), n.title()), fontsize=20)
    heights = [2, 1, 2]
    spec = fig.add_gridspec(nrows = 3, ncols = len(osimvars["ik"]), height_ratios = heights) 
    
    # Plot results
    for s in range(len(subjtype)): 
    
        # Create plots
        x = range(101)
        for col in range(len(osimvars["ik"])):        
            
            # Mean + stdev
            for r, row in enumerate([0, 2]):        
                 
                an = analyses[r]
                               
                # Mean
                m0 = desc[n][an][osimvars[an][col]][0]["mean"]
                m1 = desc[n][an][osimvars[an][col]][1]["mean"]
                
                # Upper
                u0 = m0 + desc[n][an][osimvars[an][col]][0]["sd"]
                u1 = m1 + desc[n][an][osimvars[an][col]][1]["sd"]
                
                # Lower
                l0 = m0 - desc[n][an][osimvars[an][col]][0]["sd"]
                l1 = m1 - desc[n][an][osimvars[an][col]][1]["sd"]       
                
                # Plot
                ax = fig.add_subplot(spec[row, col])
                ax.set_title(plotheads[an][col], fontsize = 12)
                if (row == 0) and ((col == 0) or (col == 4)):
                    ax.set_ylabel("Angle (deg)", fontsize = 12)
                elif (row == 2) and ((col == 0) or (col == 4)):
                    ax.set_ylabel("Moment (%BW*HT)", fontsize = 12)   
                ax.fill_between(x, l1, u1, alpha = 0.3, linewidth = 0.0, color = "blue")
                ax.fill_between(x, l0, u0, alpha = 0.3, linewidth = 0.0, color = "red")
                ax.plot(x, m1, label = subjtypefulllabel[1], linewidth = 2.0, color = "blue") 
                ax.plot(x, m0, label = subjtypefulllabel[0], linewidth = 2.0, color = "red")
                ax.set_xlim([x[0], x[-1]])
                ax.set_xlabel("% of step-down-and-pivot task", fontsize = 8)
                for v in range(1, 5): ax.axvline(x = eventlist[v], linewidth = 1.0, linestyle = ":", color = "k")
                if (row == 0 and col == 0): ax.legend(frameon = False, loc = "lower left")
            
                # Event labels
                #for at in range(6): ax.text((eventlist[at] / 100) + eventlabeladjust[at], 0.95, eventlabels[at], transform = ax.transAxes, horizontalalignment = eventlabelalign[at], fontsize = 8)
                
                # SPM significance shading
                issig = [1 if abs(z) > spmtinf[n][an][osimvars[an][col]].zstar else 0 for t, z in enumerate(spmtinf[n][an][osimvars[an][col]].z)]
                issigdiff = np.diff([0] + issig + [0])
                t0s = np.where(issigdiff == 1)
                t1s = np.where(issigdiff == -1)  # Should be the same length as t0s, I hope!
                if t0s[0].tolist():
                    for t in range(np.size(t0s[0])):
                        ax.axvspan(t0s[0][t], t1s[0][t], alpha = 0.3, color = "grey")               
                    
    # save to pdf
    plt.savefig(os.path.join(outpath, outfilename + bonftext[bonferroni] + "-" + n + "-forpub.pdf"))