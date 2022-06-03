# -*- coding: utf-8 -*-
"""
Run SPM1D analyses: FORCE step-down-pivot

@author: Prasanna Sritharan, June 2022
"""

import os
import numpy as np
import pandas as pd
import spm1d



# data file
srcpath = r"C:\Users\Owner\Documents\projects\force-moco\rstats\r-data"
srcfile = "force_sdp_results_normalised_sdp.csv"
    
# load XLS data into a dataframe
df = pd.read_csv(os.path.join(srcpath, srcfile))

# dict containing data blocks
datamat = {}

# desired variables
analyses = ["ik", "id"]
osimvars = {}
osimvars["ik"] = ["hip_flexion", "hip_adduction", "hip_rotation", "knee_angle", "ankle_angle", "lumbar_extension", "lumbar_bending", "lumbar_rotation"]
osimvars["id"] = ["hip_flexion_moment", "hip_adduction_moment", "hip_rotation_moment", "knee_angle_moment", "ankle_angle_moment", "lumbar_extension_moment", "lumbar_bending_moment", "lumbar_rotation_moment"]


# split into smaller data frames, first by limb data (IPSI, CONTRA)
leg = ["ipsi", "contra"]
groups = ["SYM", "ASYM", "DOM", "NDOM"]
for lg in leg:
    datamat[lg] = {}
    legdata = df[df["data_limb"]==lg.upper()]
    
    # split into analysis
    for an in analyses:        
        datamat[lg][an] = {}
        andata = legdata[legdata["analysis"]==an]
        
        # split into variables
        for va in osimvars[an]:
            datamat[lg][an][va] = {}
            vadata = andata[andata["variable"]==va]
            
            # split into limb type (DOM, NDOM, SYM, ASYM)
            for gp in groups:            
                gpdata = vadata[vadata[lg + "_limb"]==gp]
                gpdata = gpdata.drop(gpdata.columns[range(0, 21)], axis = 1)
                datamat[lg][an][va][gp] = gpdata.to_numpy()


# run a test
Y0 = datamat["contra"]["ik"]["hip_rotation"]["DOM"]
Y1 = datamat["contra"]["ik"]["hip_rotation"]["SYM"]
t  = spm1d.stats.ttest2(Y0, Y1, equal_var=False)

# analysis
ti = t.inference(alpha=0.05, two_tailed=True, interp=True)

# plot
print(ti)
ti.plot()


                
            
            
        