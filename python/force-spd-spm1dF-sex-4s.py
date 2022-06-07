# -*- coding: utf-8 -*-
"""
Run SPMt 1-D analyses: FORCE step-down-pivot by sex across all groups

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
outpkl = "force-sdp-spm1dF-sex-4s.pkl"
outfigprefix = "force-sdp-spm1dF-"


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
leg = ["ipsi", "contra"]
groups = ["SYM", "DOM"]  # ignore ASYM and NDOM
sex = [1, 2]
sexstr = ["M", "F"]

# split into smaller data frames for analysis
datamat = {}
for lg in leg:
    datamat[lg] = {}
    legdata = df[df["data_limb"] == lg.upper()]
    for an in analyses:        
        datamat[lg][an] = {}
        andata = legdata[legdata["analysis"] == an]
        for va in osimvars[an]:
            datamat[lg][an][va] = {}
            vadata = andata[andata["variable"] == va]
            for gp in groups:            
                datamat[lg][an][va][gp] = {}
                gpdata = vadata[vadata[lg + "_limb"] == gp]
                for sx in sex:
                    sxdata = gpdata[gpdata["sex"] == sx]    
                    sxdata = sxdata.drop(sxdata.columns[range(0, 14)], axis = 1)
                    datamat[lg][an][va][gp][sexstr[sx - 1]] = sxdata.to_numpy()


# # event times
events = {}
events["data"] = {}
events["desc"] = {}
descmat = np.zeros((2,6))
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


# %% RUN ANALYSES :DESCRIPTIVES, SPM{F} (ONE WAY ANOVA)

# comparisons
pairs = {}
pairs["sex"] = ["M", "F"]

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
                for sxs in sexstr:
                    desc[lg][an][va][gp][sxs] = {}
                    desc[lg][an][va][gp][sxs]["mean"] = np.mean(datamat[lg][an][va][gp][sxs], axis = 0)
                    desc[lg][an][va][gp][sxs]["sd"] = np.std(datamat[lg][an][va][gp][sxs], axis = 0)


# run SPM{F} (one way ANOVA) and inference
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
            Y = []
            Ystr = []
            for gp in groups:
                spmt[lg][an][va][gp] = {}
                spmtinf[lg][an][va][gp] = {}   
                for pa in pairs:
                    Y.append(datamat[lg][an][va][gp][pairs[pa][0]])
                    Y.append(datamat[lg][an][va][gp][pairs[pa][1]])
                    Ystr.append(gp + "_" + pairs[pa][0])
                    Ystr.append(gp + "_" + pairs[pa][1])
            spmt[lg][an][va] = spm1d.stats.anova1(tuple(Y), equal_var=False)
            spmtinf[lg][an][va] = spmt[lg][an][va].inference(alpha=0.05, interp=True)


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
eventlabels = ["IFO1", "IFS2", "CFO1", "CFS3", "IFO2", "IFS4"]
eventlabelalign = ["left", "right", "left", "right", "left", "right"]
eventlabeladjust = [0.01, -0.01, 0.01, -0.01, 0.01, -0.01]
limblabel= ["pivot", "non-pivot"]
grouplabels = {}
grouplabels["sex"] = ['symptomatic male', 'symptomatic female', 'control male', 'control female']
pairnsubjs = {}
pairnsubjs["sex"] = [0, 1]
filelabels = ["male-female"]

# number of M vs F per group
nsubjs = [np.size(datamat["ipsi"]["ik"]["ankle_angle"][gp][sxs], axis=0) for gp in groups for sxs in sexstr]



# plots
for ln, lmb in enumerate(leg):
    
    # plot setup
    fig = plt.figure(constrained_layout=True, figsize=(50, 15))   
    fig.suptitle("Step-down-pivot - sym male, sym female, ctrl male, ctrl female (%d, %d, %d, %d)), %s limb" % (nsubjs[0], nsubjs[1], nsubjs[2], nsubjs[3], limblabel[ln]), fontsize=20)
    heights = [2, 1, 0.5, 2, 1]
    spec = fig.add_gridspec(nrows=5, ncols=len(osimvars["ik"]), height_ratios=heights)       
        
    # create plots
    x = range(101)
    for col in range(len(osimvars["ik"])):        
        
        # mean + sd
        for r, row in enumerate([0, 3]):        
            
            an = analyses[r]
            
            # collate waveforms
            m = []
            u = []
            l = []
            for gp in groups:
                for pa in pairs:
                    for e in range(len(pairs[pa])):
                        
                        # mean
                        m0 = desc[lmb][an][osimvars[an][col]][gp][pairs[pa][e]]["mean"]
                        m.append(m0)
                        
                        # upper
                        u0 = m0 + desc[lmb][an][osimvars[an][col]][gp][pairs[pa][0]]["sd"]
                        u.append(u0)
                        
                        # lower
                        l0 = m0 - desc[lmb][an][osimvars[an][col]][gp][pairs[pa][0]]["sd"]
                        l.append(l0)      
            
            # axis
            ax = fig.add_subplot(spec[row, col], title=osimvars[an][col].replace("_", "-"))
            if (row == 0) and (col == 0):
                ax.set_ylabel("Angle (deg)")
            elif (row == 3) and (col == 0):
                ax.set_ylabel("Moment (%BW*HT)")   
                
            # plot waveforms, set xlims
            for ms in range(len(m)):
                ax.fill_between(x, l[ms], u[ms], alpha=0.4)
                ax.plot(x, m[ms], label=grouplabels["sex"][ms], linewidth=2.0)
            ax.set_xlim([x[0],x[-1]])
            
            # event lines
            for v in range(1, 5): ax.axvline(x=eventlist[v], linewidth=1.0, linestyle=":", color="k")
        
            # event labels and legend
            if (row == 0 and col == 0): ax.legend(frameon=False, loc="lower left")
            for at in range(6): ax.text((eventlist[at] / 100) + eventlabeladjust[at], 0.95, eventlabels[at], transform=ax.transAxes, horizontalalignment=eventlabelalign[at], fontsize=8)
        
        # SPM inference
        for r, row in enumerate([1, 4]):  
            
            an = analyses[r]
            
            # plot
            ax = fig.add_subplot(spec[row, col], xlabel="% task")
            if col == 0: ax.set_ylabel("SPM{F}") 
            for v in range(1, 5): ax.axvline(x=eventlist[v], linewidth=1.0, linestyle=":", color="k")
            spmtinf[lmb][an][osimvars[an][col]].plot(plot_ylabel=False)
            
            
    # save to pdf
    plt.savefig(os.path.join(outpath, outfigprefix + filelabels[0] + "-all-" + limblabel[ln] + ".pdf")) 
