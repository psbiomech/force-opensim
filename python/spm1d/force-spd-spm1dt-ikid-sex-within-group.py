# -*- coding: utf-8 -*-
"""
Run SPMt 1-D analyses: FORCE step-down-pivot IK ID by sex within group

@author: Prasanna Sritharan, June 2022
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import spm1d
import pickle as pk



# data file
srcpath = r"C:\Users\Owner\Documents\data\FORCe\outputdatabase\csvfolder\r-output"
srcfile = "force_sdp_results_descriptives_by_subject.csv"

# output file
outpath = r"C:\Users\Owner\Documents\data\FORCe\outputdatabase\spm1d"
if not os.path.isdir(outpath): os.makedirs(outpath)
outpkl = "force-sdp-spm1dt-ikid-sex-within-group.pkl"
outfigprefix = "force-sdp-spm1dt-ikid-"


# %% PREPARE DATA

# load XLS data into a dataframe
df0 = pd.read_csv(os.path.join(srcpath, srcfile))

# remove stdev rows
df = df0[df0["statistic"]=="mean"]

# desired variables
analyses = ["ik", "id"]
osimvars = {}
osimvars["ik"] = ["hip_flexion", "hip_adduction", "hip_rotation", "knee_angle", "ankle_angle", "lumbar_extension", "lumbar_bending", "lumbar_rotation"]
osimvars["id"] = [k + "_moment" for k in osimvars["ik"]]
leg = ["pivot", "nonpivot"]
groups = ["sym", "asym", "dom", "ndom"]
sex = [1, 2]
sexstr = ["M", "F"]


# split into smaller data frames for analysis
datamat = {}
for lg in leg:
    datamat[lg] = {}
    legdata = df[df["data_leg_role"] == lg]
    for an in analyses:        
        datamat[lg][an] = {}
        andata = legdata[legdata["analysis"] == an]
        for va in osimvars[an]:
            datamat[lg][an][va] = {}
            vadata = andata[andata["variable"] == va]
            for gp in groups:            
                datamat[lg][an][va][gp] = {}
                gpdata = vadata[vadata["data_leg_type"] == gp]
                for sx in sexstr:
                    sxdata = gpdata[gpdata["sex"] == sx]    
                    sxdata = sxdata.loc[:, "t1":"t101"]
                    datamat[lg][an][va][gp][sx] = sxdata.to_numpy()


# event times
events = {}
events["data"] = {}
events["desc"] = {}
descmat = np.zeros((4,6))
for gn, g in enumerate(groups):
    dfreduced = df.loc[(df["analysis"]=="ik") & (df["variable"]=="time") & (df["data_leg_role"]=="pivot") & (df["data_leg_type"]==g)]    
    events["data"][g] = dfreduced.loc[:, "es1_PFO1":"es6_PFS4"]
    events["desc"][g] = {}
    events["desc"][g]["mean"] = np.round(np.mean(events["data"][g].to_numpy(), axis=0))
    events["desc"][g]["sd"] = np.round(np.std(events["data"][g].to_numpy(), axis=0))
    descmat[gn, 0:6] = events["desc"][g]["mean"]
events["desc"]["total"] = {}
events["desc"]["total"]["mean"] = np.mean(descmat, axis=0)
events["desc"]["total"]["sd"] = np.std(descmat, axis=0)


# %% RUN ANALYSES :DESCRIPTIVES, SPM{t}
#
# Compare M vs F within groups



# calculate descriptives from file
desc = {}
for lg in leg:
    desc[lg] = {}
    for an in analyses:        
        desc[lg][an] = {}
        for va in osimvars[an]:
            desc[lg][an][va] = {}
            for gp in groups:    
                desc[lg][an][va][gp] = {}                      
                for sx in sexstr:
                    desc[lg][an][va][gp][sx] = {}
                    desc[lg][an][va][gp][sx]["mean"] = np.mean(datamat[lg][an][va][gp][sx], axis = 0)
                    desc[lg][an][va][gp][sx]["sd"] = np.std(datamat[lg][an][va][gp][sx], axis = 0)



# comparisons
pairs = {}
pairs["sex"] = ["M", "F"]

# run SPM{t} and inference across all legs, analyses, variables and group pairs
spmt = {}
spmtinf = {}
for lg in leg:
    spmt[lg] = {}
    spmtinf[lg] = {}
    for an in analyses: 
        spmt[lg][an] = {}
        spmtinf[lg][an] = {}
        for va in osimvars[an]:
            spmt[lg][an][va] = {}
            spmtinf[lg][an][va] = {}
            for gp in groups:
                spmt[lg][an][va][gp] = {}
                spmtinf[lg][an][va][gp] = {}                
                for pa in pairs:
                    Y0 = datamat[lg][an][va][gp][pairs[pa][0]]
                    Y1 = datamat[lg][an][va][gp][pairs[pa][1]]
                    spmt[lg][an][va][gp][pa] = spm1d.stats.ttest2(Y0, Y1, equal_var=False)
                    spmtinf[lg][an][va][gp][pa] = spmt[lg][an][va][gp][pa].inference(alpha=0.05, two_tailed=True, interp=True)


# combine for output
sdp = {}
sdp["desc"] = desc
sdp["events"] = events
sdp["spmt"] = spmt
sdp["spmtinf"] = spmtinf
            
# save to pickle
with open(os.path.join(outpath, outpkl),"wb") as f: pk.dump(sdp, f)
    

# %% PLOT OUTPUT

# plot parameters
eventlist = 100 * np.round(events["desc"]["total"]["mean"]) / 101
eventlabels = ["PFO1", "PFS2", "NFO1", "NFS3", "PFO2", "PFS4"]
eventlabelalign = ["left", "right", "left", "right", "left", "right"]
eventlabeladjust = [0.01, -0.01, 0.01, -0.01, 0.01, -0.01]
limblabel= ["pivot", "nonpivot"]
pairlabels = {}
pairlabels["sex"] = ["male", "female"]
pairnsubjs = {}
pairnsubjs["sex"] = [0, 1]
filelabels = ["male-female"]

# plots
for gp in groups:
    
    # number of M vs F per group
    nsubjs = [np.size(datamat["pivot"]["ik"]["ankle_angle"][gp][sx], axis=0) for sx in sexstr]
    
    for p, pa in enumerate(pairs):
        for ln, lmb in enumerate(leg):
        
            # plot setup
            fig = plt.figure(constrained_layout=True, figsize=(50, 15))   
            fig.suptitle("Step-down-pivot - %s (n%d vs n%d) - %s %s limb" % (filelabels[p].upper().replace("-", " vs "), nsubjs[pairnsubjs[pa][0]], nsubjs[pairnsubjs[pa][1]], gp, limblabel[ln]), fontsize=20)
            heights = [2, 1, 0.5, 2, 1]
            spec = fig.add_gridspec(nrows=5, ncols=len(osimvars["ik"]), height_ratios=heights)
            
            # create plots
            x = range(101)
            for col in range(len(osimvars["ik"])):        
                
                # mean + sd
                for r, row in enumerate([0, 3]):        
                    
                    an = analyses[r]
                    
                    # mean
                    m0 = desc[lmb][an][osimvars[an][col]][gp][pairs[pa][0]]["mean"]
                    m1 = desc[lmb][an][osimvars[an][col]][gp][pairs[pa][1]]["mean"]
                    
                    # upper
                    u0 = m0 + desc[lmb][an][osimvars[an][col]][gp][pairs[pa][0]]["sd"]
                    u1 = m1 + desc[lmb][an][osimvars[an][col]][gp][pairs[pa][1]]["sd"]
                    
                    # lower
                    l0 = m0 - desc[lmb][an][osimvars[an][col]][gp][pairs[pa][0]]["sd"]
                    l1 = m1 - desc[lmb][an][osimvars[an][col]][gp][pairs[pa][1]]["sd"]       
                    
                    # plot
                    ax = fig.add_subplot(spec[row, col], title=osimvars[an][col].replace("_", "-"))
                    if (row == 0) and (col == 0):
                        ax.set_ylabel("Angle (deg)")
                    elif (row == 3) and (col == 0):
                        ax.set_ylabel("Moment (%BW*HT)")   
                    ax.fill_between(x, l0, u0, alpha=0.4)
                    ax.fill_between(x, l1, u1, alpha=0.4)
                    ax.plot(x, m0, label=pairlabels[pa][0], linewidth=2.0)
                    ax.plot(x, m1, label=pairlabels[pa][1], linewidth=2.0) 
                    ax.set_xlim([x[0],x[-1]])
                    for v in range(1, 5): ax.axvline(x=eventlist[v], linewidth=1.0, linestyle=":", color="k")
                    if (row == 0 and col == 0): ax.legend(frameon=False, loc="lower left")
                
                    # event labels
                    for at in range(6): ax.text((eventlist[at] / 100) + eventlabeladjust[at], 0.95, eventlabels[at], transform=ax.transAxes, horizontalalignment=eventlabelalign[at], fontsize=8)
                
                # SPM inference
                for r, row in enumerate([1, 4]):  
                    
                    an = analyses[r]
                    
                    # plot
                    ax = fig.add_subplot(spec[row, col], xlabel="% task")
                    if col == 0: ax.set_ylabel("SPM{t}") 
                    for v in range(1, 5): ax.axvline(x=eventlist[v], linewidth=1.0, linestyle=":", color="k")
                    spmtinf[lmb][an][osimvars[an][col]][gp][pa].plot(plot_ylabel=False)
                    
                    
            # save to pdf
            plt.savefig(os.path.join(outpath, outfigprefix + filelabels[p] + "-" + gp.lower() + "-" + limblabel[ln] + ".pdf")) 
