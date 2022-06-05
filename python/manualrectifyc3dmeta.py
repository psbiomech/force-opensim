# -*- coding: utf-8 -*-
"""
Manually correct event label errors

@author: Prasanna Sritharan, June 2022
"""

import os
import pickle as pk



# %% FAILTCRT01_SDP07

pkfile = r"C:\Users\Owner\Documents\data\FORCe\outputdatabase\FAILTCRT01\FAILTCRT01_SDP07\FAILTCRT01_SDP07_osimkey.pkl"
with open(pkfile, "rb") as fid:
    osimkey = pk.load(fid)

osimkey.events["labels"] = ["LFO", "LFS", "RFO", "RFS", "LFO", "LFS"]

with open(pkfile,"wb") as fid: 
  pk.dump(osimkey, fid)