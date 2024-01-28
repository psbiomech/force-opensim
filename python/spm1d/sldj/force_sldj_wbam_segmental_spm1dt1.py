# -*- coding: utf-8 -*-
"""
FORCE SINGLE LEG DROP JUMP: WBAM SEGMENTAL SPM{t}

@author: Prasanna Sritharan, August 2023
"""


import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import spm1d
import pickle as pk



# Data file
srcpath = r"C:\Users\Owner\Documents\data\FORCe\outputdatabase_sldj\csvfolder"
srcfile = "force_sldj_results_subject_descriptives_stability_normalised.csv"

# Output file
outpath = r"C:\Users\Owner\Documents\data\FORCe\outputdatabase_sldj\spm1d"
if not os.path.isdir(outpath): os.makedirs(outpath)
outfilename = "force_sldj_spm1dt_wbam_segmental_more_ctrl"


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
scenario = ["stance", "swing"]
trialcombo =[["stance_more", "stance"], ["stance_more", "stance"]]   # symptomatic stance leg and asymptompatic swing leg

# Get data into arrays
datamat = {}
eventmat = []
for sn, n in enumerate(scenario):
    df = df1[df1["data_leg_role"] == n]
    datamat[n] = {}
    for a in analyses:
        datamat[n][a] = {}
        for v in ["time"] + osimvars[a]:
            
            datamat[n][a][v] = {}
            
            # Group 1
            gp1data = df[(df["subj_type"] == subjtype[0]) & (df["trial_combo"] == trialcombo[sn][0]) & (df["analysis"] == a) & (df["variable"] == v)]
            datamat[n][a][v][0] = gp1data.loc[:, "t1":"t101"].to_numpy()
                   
            # Group 2
            gp2data = df[(df["subj_type"] == subjtype[1]) & (df["trial_combo"].str.contains(trialcombo[sn][1])) & (df["analysis"] == a) & (df["variable"] == v)]
            datamat[n][a][v][1] = gp2data.loc[:, "t1":"t101"].to_numpy()
        
            # Get max knee flexion angle, i.e. most negative knee extension
            if (n == "stance") and (v == "knee_angle") and (a == "ik"): 
                eventmat.append(np.argmax(gp1data.loc[:, "t1":"t101"].to_numpy(), axis=1))
                eventmat.append(np.argmax(gp2data.loc[:, "t1":"t101"].to_numpy(), axis=1))



# %% EVENTS

# Event time steps and descriptives
events = {}
events["data"] = {}
events["desc"] = {}
descmat = np.zeros(2)
for g in range(len(subjtype)):
    events["data"][g] = eventmat[g]
    events["desc"][g] = {}
    events["desc"][g]["mean"] = np.round(np.mean(events["data"][g], axis=0))
    events["desc"][g]["sd"] = np.round(np.std(events["data"][g], axis=0))
    descmat[g] = events["desc"][g]["mean"]
events["desc"]["total"] = {}
events["desc"]["total"]["mean"] = np.mean(descmat, axis=0)
events["desc"]["total"]["sd"] = np.std(descmat, axis=0)



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
            spmtinf[n][a][v] = spmt[n][a][v].inference(alpha = 0.05, two_tailed=True, interp=True)


# Combine for output
sldj = {}
sldj["desc"] = desc
sldj["events"] = events
sldj["spmt"] = spmt
sldj["spmtinf"] = spmtinf
            
# Pickle it
#with open(os.path.join(outpath, outfilename + ".pkl"),"wb") as f: pk.dump(sldj, f)
    

# %% PLOT OUTPUT

# Plot parameters
eventlabels = ["PFO1", "PFS2", "NFO1", "NFS2", "PFO3", "PFS4"]
eventlabelalign = ["left", "right", "left", "right", "left", "right"]
eventlabeladjust = [0.01, -0.01, 0.01, -0.01, 0.01, -0.01]   

# Figure headers
plotheads = {}
plotheads["ik"] = ["Hip flexion", "Hip adduction", "Hip rotation", "Knee flexion", "Ankle dorsiflexion", "Lumbar extension", "Lumbar bending", "Lumbar rotation"]
plotheads["id"] = [k + " moment" for k in plotheads["ik"]]


# Generate plots
eventlist = events["desc"]["total"]["mean"]
for n in scenario:

    # Scenario parameters
    nsubjs = [np.size(datamat[n]["ik"]["time"][s], axis=0) for s in range(len(subjtype))]

    # Create plot area
    fig = plt.figure(constrained_layout=True, figsize=(24, 10))   
    fig.suptitle("Single-leg drop jump. Stance limb: %s (n%d) vs %s (n%d). %s limb data." % (subjtypefulllabel[0].upper(), np.shape(gp1data)[0], subjtypefulllabel[1].upper(), np.shape(gp2data)[0], n.title()), fontsize=20)
    heights = [2, 1, 0.5, 2, 1]
    spec = fig.add_gridspec(nrows = 5, ncols = len(osimvars["ik"]), height_ratios = heights) 
    
    # Plot results
    for s in range(len(subjtype)): 
    
        # Create plots
        x = range(101)
        for col in range(len(osimvars["ik"])):        
            
            # Mean + stdev
            for r, row in enumerate([0, 3]):        
                 
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
                if (row == 0) and (col == 0):
                    ax.set_ylabel("Angle (deg)", fontsize = 12)
                elif (row == 3) and (col == 0):
                    ax.set_ylabel("Moment (%BW*HT)", fontsize = 12)   
                ax.fill_between(x, l1, u1, alpha = 0.3, linewidth = 0.0, color = "blue")
                ax.fill_between(x, l0, u0, alpha = 0.3, linewidth = 0.0, color = "red")
                ax.plot(x, m1, label = subjtypefulllabel[1], linewidth = 2.0, color = "blue") 
                ax.plot(x, m0, label = subjtypefulllabel[0], linewidth = 2.0, color = "red")
                ax.set_xlim([x[0], x[-1]])
                ax.axvline(x = eventlist, linewidth = 1.0, linestyle = ":", color = "k")
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
                
            # SPM inference
            for r, row in enumerate([1, 4]):  
                
                an = analyses[r]
                
                # plot
                ax = fig.add_subplot(spec[row, col])
                ax.set_xlabel("% of stance", fontsize = 12)
                if col == 0: ax.set_ylabel("SPM{t}", fontsize = 10) 
                ax.axvline(x = eventlist, linewidth = 1.0, linestyle = ":", color = "k")
                spmtinf[n][an][osimvars[an][col]].plot(plot_ylabel = False)
                    
                    
    # save to pdf
    plt.savefig(os.path.join(outpath, outfilename + "_" + n + ".pdf"))