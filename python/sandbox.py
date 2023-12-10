# -*- coding: utf-8 -*-
"""
Sandbox

@author: Prasanna Sritharan
"""


# %% LOAD USER SETTINGS

import usersettings as uset

user = uset.FORCESettings_SDP()



# %% LOAD META DATA

import usersettings as uset
import pickle as pk
import os

user = uset.FORCESettings_SLDJ()

dbfilepath = os.path.join(user.rootpath, user.outfolder, user.metadatafile)
with open(dbfilepath, "rb") as fid:
    forcedb = pk.load(fid)


# %% LOAD A TRIALKEY

import pickle as pk
import os

# file path and name prefix
fprefix = "FAILT08_SLDJ03"
fpath = r"C:\Users\Owner\Documents\data\FORCe\outputdatabase_sldj\FAILT08\FAILT08_SLDJ03"

# OsimResultsKey
pkfile = os.path.join(fpath, fprefix + "_trialkey.pkl")
with open(pkfile, "rb") as fid: 
    trialkey1 = pk.load(fid)
    
    

# %% LOAD AN OSIMKEY

import pickle as pk
import os

# file path and name prefix
fprefix = "FAILT02_SLDJ08"
fpath = r"C:\Users\Owner\Documents\data\FORCe\outputdatabase_sldj\FAILT02\FAILT02_SLDJ08"

# OsimResultsKey
pkfile = os.path.join(fpath, fprefix + "_osimkey.pkl")
with open(pkfile, "rb") as fid: 
    osimkey1 = pk.load(fid)
    
    
    
# %% LOAD A C3DKEY

import pickle as pk
import os

# file path and name prefix
fprefix = "FAILT01_SLDJ01"
fpath = r"C:\Users\Owner\Documents\data\FORCe\outputdatabase_sldj\FAILT01\FAILT01_SLDJ01"

# OsimResultsKey
pkfile = os.path.join(fpath, fprefix + "_c3dkey.pkl")
with open(pkfile, "rb") as fid: 
    c3dkey1 = pk.load(fid)
    
    


# %% LOAD AN OSIMRESULTSKEY

import pickle as pk
import os

# file path and name prefix
fprefix = "FAILT08_SLDJ03"
fpath = r"C:\Users\Owner\Documents\data\FORCe\outputdatabase_sldj\FAILT08\FAILT08_SLDJ03"

# OsimResultsKey
pkfile = os.path.join(fpath, fprefix + "_opensim_results.pkl")
with open(pkfile, "rb") as fid: 
    osimresultskey1 = pk.load(fid)
    

# %% LOAD A STABILITY KEY

import pickle as pk
import os

# file path and name prefix
fprefix = "FAILT24_SLDJ06"
fpath = r"C:\Users\Owner\Documents\data\FORCe\outputdatabase_sldj\FAILT24\FAILT24_SLDJ06"

# OsimResultsKey
pkfile = os.path.join(fpath, fprefix + "_stability.pkl")
with open(pkfile, "rb") as fid: 
    stabilitykey1 = pk.load(fid)
    


 # %% LOAD A WBAM KEY

import pickle as pk
import os

# file path and name prefix
fprefix = "FAILTCRT01_SLDJ02"
fpath = r"C:\Users\Owner\Documents\data\FORCe\outputdatabase_sldj\FAILTCRT01\FAILTCRT01_SLDJ02"

# OsimResultsKey
pkfile = os.path.join(fpath, fprefix + "_wbam.pkl")
with open(pkfile, "rb") as fid: 
    wbamkey1 = pk.load(fid)
     
     
     

# %% CALCULATE JOINT ANGULAR IMPULSE

import analyses as an

impl = an.calculate_joint_angular_impulse(osimkey1, user)