# -*- coding: utf-8 -*-
"""
Sandbox

@author: Prasanna Sritharan
"""


# %% LOAD USER SETTINGS

import usersettings as uset

user = uset.FORCESettings_SDP()



# %% LOAD META DATA

import pickle as pk
import os

dbfilepath = os.path.join(user.rootpath, user.outfolder, user.metadatafile)
with open(dbfilepath, "rb") as fid:
    forcedb = pk.load(fid)




# %% LOAD AN OSIMRESULTSKEY

import pickle as pk
import os

# file path and name prefix
fprefix = "FAILTCRT05_SDP05"
fpath = r"C:\Users\Owner\Documents\data\FORCe\outputdatabase\FAILTCRT05\FAILTCRT05_SDP05"

# OsimResultsKey
pkfile = os.path.join(fpath, fprefix + "_opensim_results.pkl")
with open(pkfile, "rb") as fid: 
    osimkey1 = pk.load(fid)
    
    

# %% CALCULATE JOINT ANGULAR IMPULSE

import analyses as an

impl = an.calculate_joint_angular_impulse(osimkey1, user)