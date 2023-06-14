# -*- coding: utf-8 -*-
"""
Run SPMt 1-D regression: FORCE step-down-pivot IK ID by SHOMRI

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
srcfile = "force_sdp_normalised_descriptives_subject.csv"

# output file
outpath = r"C:\Users\Owner\Documents\data\FORCe\outputdatabase\spm1d"
if not os.path.isdir(outpath): os.makedirs(outpath)
outpkl = "force-sdp-spm1dtR-ikid-shomri-group.pkl"
outfigprefix = "force-sdp-spm1dtR-ikid-"


# %% PREPARE DATA

# load XLS data into a dataframe
df0 = pd.read_csv(os.path.join(srcpath, srcfile))

# drop all rows with no SHOMRI scores (NA)
df1 = df0.dropna(axis=0, how="any")

# remove stdev rows
df = df1[df1["statistic"]=="mean"]

# desired variables
analyses = ["ik", "id"]
osimvars = {}
osimvars["ik"] = ["hip_flexion", "hip_adduction", "hip_rotation", "knee_angle", "ankle_angle", "lumbar_extension", "lumbar_bending", "lumbar_rotation"]
osimvars["id"] = [k + "_moment" for k in osimvars["ik"]]
leg = ["ipsi", "contra"]
groups = ["SYM", "ASYM", "DOM", "NDOM"]


# split into smaller data frames for analysis
datamat = {}
for lg in leg:
    datamat[lg] = {}
    legdata = df[df["data_limb"]==lg.upper()]
    for an in analyses:        
        datamat[lg][an] = {}
        andata = legdata[legdata["analysis"]==an]
        for va in osimvars[an]:
            datamat[lg][an][va] = {}
            vadata = andata[andata["variable"]==va]
            for gp in groups:            
                gpdata = vadata[vadata[lg + "_limb"]==gp]
                gpdata = gpdata.drop(gpdata.columns[range(0, 16)], axis = 1)
                datamat[lg][an][va][gp] = gpdata.to_numpy()


# get the SHOMRI scores
shomri = {}
for lg in leg:
    shomri[lg] = {}
    legdata = df[df["data_limb"]==lg.upper()]
    for an in analyses:        
        shomri[lg][an] = {}
        andata = legdata[legdata["analysis"]==an]
        for va in osimvars[an]:
            shomri[lg][an][va] = {}
            vadata = andata[andata["variable"]==va]
            for gp in groups:            
                gpdata = vadata[vadata[lg + "_limb"]==gp]
                gpdata = gpdata[["shomri_" + lg]]
                shomri[lg][an][va][gp] = gpdata.to_numpy()


# # event times
events = {}
events["data"] = {}
events["desc"] = {}
descmat = np.zeros((4,6))
for gn, g in enumerate(groups):
    dfreduced = df.loc[(df["analysis"]=="ik") & (df["variable"]=="time") & (df["data_limb"]=="IPSI") & (df["ipsi_limb"]==g)]    
    events["data"][g] = dfreduced[["E" + str(e + 1) for e in range(6)]]
    events["desc"][g] = {}
    events["desc"][g]["mean"] = np.round(np.mean(events["data"][g].to_numpy(), axis=0))
    events["desc"][g]["sd"] = np.round(np.std(events["data"][g].to_numpy(), axis=0))
    descmat[gn, 0:6] = events["desc"][g]["mean"]
events["desc"]["total"] = {}
events["desc"]["total"]["mean"] = np.mean(descmat, axis=0)
events["desc"]["total"]["sd"] = np.std(descmat, axis=0)


# %% RUN ANALYSES :DESCRIPTIVES, SPM{t} regression


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
                desc[lg][an][va][gp]["mean"] = np.mean(datamat[lg][an][va][gp], axis = 0)
                desc[lg][an][va][gp]["sd"] = np.std(datamat[lg][an][va][gp], axis = 0)


# run SPM{t} regression and inference across all legs, analyses, variables
# and group pairs
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
                Y = datamat[lg][an][va][gp]
                x = shomri[lg][an][va][gp].flatten()
                spmt[lg][an][va][gp] = spm1d.stats.regress(Y, x)
                spmtinf[lg][an][va][gp] = spmt[lg][an][va][gp].inference(alpha=0.05)


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
nsubjs = [np.size(datamat["ipsi"]["ik"]["ankle_angle"][g], axis=0) for g in groups]
eventlist = 100 * np.round(events["desc"]["total"]["mean"]) / 101
eventlabels = ["IFO1", "IFS2", "CFO1", "CFS3", "IFO2", "IFS4"]
eventlabelalign = ["left", "right", "left", "right", "left", "right"]
eventlabeladjust = [0.01, -0.01, 0.01, -0.01, 0.01, -0.01]
limblabel= ["pivot", "non-pivot"]



# plots
for g, gp in enumerate(groups):
    for ln, lmb in enumerate(leg):
    
        # plot setup
        fig = plt.figure(constrained_layout=True, figsize=(50, 15))   
        fig.suptitle("Step-down-pivot - %s vs SHOMRI (n%d) - %s limb" % (gp, nsubjs[g], limblabel[ln].title()), fontsize=20)
        heights = [2, 1, 0.5, 2, 1]
        spec = fig.add_gridspec(nrows=5, ncols=len(osimvars["ik"]), height_ratios=heights)
        
        # create plots
        x = range(101)
        for col in range(len(osimvars["ik"])):        
            
            # mean + sd
            for r, row in enumerate([0, 3]):        
                
                an = analyses[r]
                
                # mean
                m0 = desc[lmb][an][osimvars[an][col]][gp]["mean"]
                
                # upper
                u0 = m0 + desc[lmb][an][osimvars[an][col]][gp]["sd"]
                
                # lower
                l0 = m0 - desc[lmb][an][osimvars[an][col]][gp]["sd"]    
                
                # plot
                ax = fig.add_subplot(spec[row, col], title=osimvars[an][col].replace("_", "-"))
                if (row == 0) and (col == 0):
                    ax.set_ylabel("Angle (deg)")
                elif (row == 3) and (col == 0):
                    ax.set_ylabel("Moment (%BW*HT)")   
                ax.fill_between(x, l0, u0, alpha=0.4)
                ax.plot(x, m0, linewidth=2.0)
                ax.set_xlim([x[0],x[-1]])
                for v in range(1, 5): ax.axvline(x=eventlist[v], linewidth=1.0, linestyle=":", color="k")
            
                # event labels
                for at in range(6): ax.text((eventlist[at] / 100) + eventlabeladjust[at], 0.95, eventlabels[at], transform=ax.transAxes, horizontalalignment=eventlabelalign[at], fontsize=8)
            
            # SPM inference
            for r, row in enumerate([1, 4]):  
                
                an = analyses[r]
                
                # plot
                ax = fig.add_subplot(spec[row, col], xlabel="% task")
                if col == 0: ax.set_ylabel("SPM{t}") 
                for v in range(1, 5): ax.axvline(x=eventlist[v], linewidth=1.0, linestyle=":", color="k")
                spmtinf[lmb][an][osimvars[an][col]][gp].plot(plot_ylabel=False)
                
                
        # save to pdf
        plt.savefig(os.path.join(outpath, outfigprefix + gp.lower() + "-" + limblabel[ln] + ".pdf")) 
