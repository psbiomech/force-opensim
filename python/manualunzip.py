# -*- coding: utf-8 -*-
"""
Manually unzip FORCe symptomatics

@author: Prasanna Sritharan
"""

import os
import glob
import zipfile


# source
srcdir = r"C:\Users\Owner\Documents\data\FORCe\inputdatabase\ForceMaster_LTUSymptomatics\3.2. OpenSim_C3D_Good"

# file prefix
fprefix = "FAILT"



# subject folder list
folderlist = glob.glob(os.path.join(srcdir, fprefix + "*"), recursive = True)
subjlist = [os.path.split(f)[1] for f in folderlist]

# parse folders
failedfiles = []
for n, subj in enumerate(subjlist):

    print("%s" % subj)    

    # zip files
    ziplist = glob.glob(os.path.join(folderlist[n], "*.zip"), recursive = True)
    
    # extract zip files into source folder
    for fzipfile in ziplist:
        
        # zip file name
        fname = os.path.split(fzipfile)[1]
        
        # unzip
        try:
            print("---> %s" % fname)
            with zipfile.ZipFile(os.path.join(fzipfile), "r") as zipf:
                zipf.extractall(folderlist[n])
        except:
            print("---> %s *** FAILED ***" % fname)
            failedfiles.append(fname)