# -*- coding: utf-8 -*-
"""
Run SPMt 1-D analyses: FORCE step-down-pivot SO by sex within groups

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
outpkl = "force-sdp-spm1dt-so-group.pkl"
outfigprefix = "force-sdp-spm1dt-so-"


# %% PREPARE DATA

# load XLS data into a dataframe
df0 = pd.read_csv(os.path.join(srcpath, srcfile))

# remove stdev rows
df = df0[df0["statistic"]=="mean"]

# desired variables
analyses = ["so"]
osimvars = {}
osimvars["so"] = ["addbrev", "addlong", "addmagDist", "addmagIsch", "addmagMid", "addmagProx", "bflh", "bfsh", "edl", "ehl", "elbow_flex", "fdl", "fhl", "gaslat", "gasmed", "glmax1", "glmax2", "glmax3", "glmed1", "glmed2", "glmed3", "glmin1", "glmin2", "glmin3", "grac", "iliacus", "lumbar_bend", "lumbar_ext", "lumbar_rot", "perbrev", "perlong", "piri", "pro_sup", "psoas", "recfem", "sart", "semimem", "semiten", "shoulder_add", "shoulder_flex", "shoulder_rot", "soleus", "tfl", "tibant", "tibpost", "time", "vasint", "vaslat", "vasmed", "wrist_dev", "wrist_flex"]
leg = ["ipsi", "contra"]
groups = ["SYM", "ASYM", "DOM", "NDOM"]


# split into smaller data frames for analysis
datamat0 = {}
for lg in leg:
    datamat0[lg] = {}
    legdata = df[df["data_limb"]==lg.upper()]
    for an in analyses:        
        datamat0[lg][an] = {}
        andata = legdata[legdata["analysis"]==an]
        for va in osimvars[an]:
            datamat0[lg][an][va] = {}
            vadata = andata[andata["variable"]==va]           
            for gp in groups:            
                gpdata = vadata[vadata[lg + "_limb"]==gp]
                gpdata = gpdata.drop(gpdata.columns[range(0, 16)], axis = 1)
                datamat0[lg][an][va][gp] = gpdata.to_numpy()



# muscle groups
muscgrps = [("add", ["addbrev", "addlong", "addmagDist", "addmagIsch", "addmagMid", "addmagProx"]),
            ("bflh", ["bflh"]),
            ("bfsh", ["bfsh"]),
            ("gaslat", ["gaslat"]),
            ("gasmed", ["gasmed"]),
            ("glmax", ["glmax1", "glmax2", "glmax3"]),
            ("glmed", ["glmed1", "glmed2", "glmed3"]),
            ("glmin", ["glmin1", "glmin2", "glmin3"]),
            ("grac", ["grac"]),
            ("iliacus", ["iliacus"]),
            ("per", ["perbrev", "perlong"]),
            ("piri", ["piri"]),
            ("psoas", ["psoas"]),
            ("recfem", ["recfem"]),
            ("sart", ["sart"]),
            ("semimem", ["semimem"]),
            ("semiten", ["semiten"]),
            ("soleus", ["soleus"]),
            ("tfl", ["tfl"]),
            ("tibant", ["tibant"]),
            ("tibpost", ["tibpost"]),
            ("vas", ["vasint", "vaslat", "vasmed"])]



# sum muscles into muscle groups
mgmat = {}
for lg in leg:
    mgmat[lg] = {}
    for gp in groups:
        mgmat[lg][gp] = {}
        for mg in muscgrps:
            mgsum = np.zeros(np.shape(datamat0[lg]["so"][mg[1][0]][gp]))
            for m in mg[1]:
                mgsum = np.add(mgsum, datamat0[lg]["so"][m][gp])
            mgmat[lg][gp][mg[0]] = mgsum



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


# %% RUN ANALYSES :DESCRIPTIVES, SPM{t}


# calculate descriptives from file
desc = {}
for lg in leg:
    desc[lg] = {}      
    for gp in groups:
        desc[lg][gp] = {}
        for mg in muscgrps:               
            desc[lg][gp][mg[0]] = {}                      
            desc[lg][gp][mg[0]]["mean"] = np.mean(mgmat[lg][gp][mg[0]], axis = 0)
            desc[lg][gp][mg[0]]["sd"] = np.std(mgmat[lg][gp][mg[0]], axis = 0)



# comparisons
pairs = {}
pairs["group"] = ["DOM", "SYM"]
pairs["limb"] = ["ASYM", "SYM"]

# run SPM{t} and inference across all legs, analyses, variables and group pairs
spmt = {}
spmtinf = {}
for lg in leg:
    spmt[lg] = {}
    spmtinf[lg] = {}
    for mg in muscgrps:     
        spmt[lg][mg[0]] = {}
        spmtinf[lg][mg[0]] = {}
        for pa in pairs:
            Y0 = mgmat[lg][pairs[pa][0]][mg[0]]
            Y1 = mgmat[lg][pairs[pa][1]][mg[0]]
            spmt[lg][mg[0]][pa] = spm1d.stats.ttest2(Y0, Y1, equal_var=False)
            spmtinf[lg][mg[0]][pa] = spmt[lg][mg[0]][pa].inference(alpha=0.05, two_tailed=True, interp=True)



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
mgplots = [["add", "bflh", "bfsh", "gaslat", "gasmed", "glmax", "glmed", "glmin", "grac", "iliacus", "per"], 
           ["piri", "psoas", "recfem", "sart", "semimem", "semiten", "soleus", "tfl", "tibant", "tibpost", "vas"]]
nsubjs = [np.size(datamat0["ipsi"]["so"]["soleus"][g], axis=0) for g in groups]
eventlist = 100 * np.round(events["desc"]["total"]["mean"]) / 101
eventlabels = ["IFO1", "IFS2", "CFO1", "CFS3", "IFO2", "IFS4"]
eventlabelalign = ["left", "right", "left", "right", "left", "right"]
eventlabeladjust = [0.01, -0.01, 0.01, -0.01, 0.01, -0.01]
limblabel= ["pivot", "non-pivot"]
pairlabels = {}
pairlabels["group"] = ["control", "symptomatic"]
pairlabels["limb"] = ["asymptomatic", "symptomatic"]
pairnsubjs = {}
pairnsubjs["group"] = [2, 0]
pairnsubjs["limb"] = [1, 0]
filelabels = ["ctrl-sym", "asym-sym"]

# plots
for p, pa in enumerate(pairs):
    for ln, lmb in enumerate(leg):
    
        # plot setup
        fig = plt.figure(constrained_layout=True, figsize=(65, 15))   
        fig.suptitle("Step-down-pivot - %s (n%d vs n%d) - %s limb" % (filelabels[p].upper().replace("-", " vs "), nsubjs[pairnsubjs[pa][0]], nsubjs[pairnsubjs[pa][1]], limblabel[ln].title()), fontsize=20)
        heights = [2, 1, 0.5, 2, 1]
        spec = fig.add_gridspec(nrows=5, ncols=len(mgplots[0]), height_ratios=heights)
        
        # create plots
        x = range(101)
        for col in range(len(mgplots[0])):        
            
            # mean + sd
            for r, row in enumerate([0, 3]):        

                # mean
                m0 = desc[lmb][pairs[pa][0]][mgplots[r][col]]["mean"]
                m1 = desc[lmb][pairs[pa][1]][mgplots[r][col]]["mean"]
                
                # upper
                u0 = m0 + desc[lmb][pairs[pa][0]][mgplots[r][col]]["sd"]
                u1 = m1 + desc[lmb][pairs[pa][1]][mgplots[r][col]]["sd"]
                
                # lower
                l0 = m0 - desc[lmb][pairs[pa][0]][mgplots[r][col]]["sd"]
                l1 = m1 - desc[lmb][pairs[pa][1]][mgplots[r][col]]["sd"]      
                
                # plot
                ax = fig.add_subplot(spec[row, col], title=mgplots[r][col].replace("_", "-"))
                if (row == 0) and (col == 0):
                    ax.set_ylabel("Muscle force (BW)")
                elif (row == 3) and (col == 0):
                    ax.set_ylabel("Muscle force (BW)")   
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
                
                # plot
                ax = fig.add_subplot(spec[row, col], xlabel="% task")
                if col == 0: ax.set_ylabel("SPM{t}") 
                for v in range(1, 5): ax.axvline(x=eventlist[v], linewidth=1.0, linestyle=":", color="k")
                spmtinf[lmb][mgplots[r][col]][pa].plot(plot_ylabel=False)
                
                
        # save to pdf
        plt.savefig(os.path.join(outpath, outfigprefix + filelabels[p] + "-" + limblabel[ln] + ".pdf")) 
