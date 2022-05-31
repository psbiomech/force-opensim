# -*- coding: utf-8 -*-
"""
Manually set subject mass if zero

@author: Prasanna Sritharan
"""


import os
import pickle as pk
import pandas as pd


# paths
rootpath = r"C:\Users\Owner\Documents\data\FORCe"
infolder = "inputdatabase"
outfolder = "outputdatabase"

# load meta database
metadatafile = "FORCE_SDP.pkl"  
dbfilepath = os.path.join(rootpath, outfolder, metadatafile)
with open(dbfilepath, "rb") as fid:
    forcedb = pk.load(fid)

# load participant data
xlsfile = "FORCE-ParticipantData-All.xlsx"
pdata = pd.read_excel(os.path.join(rootpath, infolder, xlsfile))

# check files and add mass if required
for subj in forcedb.keys():
    subjfolder = forcedb[subj]["outpath"]
    for grp in forcedb[subj]["trials"].keys():
        for trial in forcedb[subj]["trials"][grp].keys():
                
            # load the trialkey and the osimkey
            tkeypath = os.path.join(rootpath, outfolder, subj, trial, trial + "_trialkey.pkl")
            okeypath = os.path.join(rootpath, outfolder, subj, trial, trial + "_osimkey.pkl")
            with open(tkeypath, "rb") as kfid0:
                trialkey = pk.load(kfid0)
            with open(okeypath, "rb") as ofid0:
                osimkey = pk.load(ofid0)                   
                
            # check mass and update if required
            if trialkey.mass == 0.0:
                
                # find mass in participant data
                subjmass = pdata.loc[pdata["id"] == subj]["mass"].to_list()[0]
                
                # update mass in trialkey and osimkey
                print("%s: %.3f ---> %.3f" % (trial, trialkey.mass, subjmass))
                trialkey.mass = subjmass
                osimkey.mass = subjmass
                
                # save the updated trialkey and osimkey
                with open(tkeypath, "wb") as kfid1:
                    pk.dump(trialkey, kfid1)
                with open(okeypath, "wb") as ofid1:
                    pk.dump(osimkey, ofid1)
                    
                
                
                    
