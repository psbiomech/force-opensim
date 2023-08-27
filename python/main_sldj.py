# -*- coding: utf-8 -*-
"""
Process and run LASEM FORCE project data through OpenSim

@author: Prasanna Sritharan
"""


import datetime


print("\n\n\n")
print("----------------------------------------")
print("FORCE: DATA PROCESSING & OPENSIM")
print("----------------------------------------")

# start time stamp
ts0 = datetime.datetime.now();
print("Start: %s" % ts0)

print("----------------------------------------")
print("\n")


# %% SET THE LAB

import labsetup as labs

print("Loading lab info...", end="")
lasem = labs.LabKeyLasemForceSLDJ()
print("Done.\n")


# %% USER SETTINGS

import usersettings as uset

print("Loading user settings... ", end="")
user = uset.FORCESettings_SLDJ()
print("Done.\n")


# %% GET META DATABASE (BUILD NEW OR LOAD EXISTING)...


# OPTION 1: Build new database

# import builddatabase as bd

# print("Building new output database... ", end="")
# forcedb, failedfilesBD = bd.build_database(user, "sldj")
# print("Done.\n")


# OPTION 2: Load existing database

import pickle as pk
import os

print("Loading existing output database... ", end="")
dbfilepath = os.path.join(user.rootpath, user.outfolder, user.metadatafile)
with open(dbfilepath, "rb") as fid:
    forcedb = pk.load(fid)
print("Done.\n")


# %% EXTRACT C3D AND CREATE OPENSIM DATA FILES

# import c3dextract as c3dex

# print("Extracting C3D data, creating OpenSim files...\n")
# failed_files = c3dex.c3d_batch_process(user, forcedb, lasem, -1, restart = ("FAILT70", "FAILT70"))
# print("\nC3D data extract done.\n")


# %% RUN OPENSIM PIPELINE

# import opensimpipeline as osp

# print("Running OpenSim model scaling: SCALE...\n")
# osp.opensim_pipeline(forcedb, user, ["scale"])
# print("\nOpenSim model scaling (SCALE) completed.\n")

# print("Running OpenSim analyses: IK, ID...\n")
# failedfiles_ikid = osp.opensim_pipeline(forcedb, user, ["ik", "id"], restart = "FAILTCRT24")
# print("\nOpenSim analyses (IK, ID) completed.\n")

# print("Running OpenSim analyses: BK...\n")
# failedfiles_bk = osp.opensim_pipeline(forcedb, user, ["bk"], restart = "FAILT100")
# print("\nOpenSim analyses (BK) completed.\n")

# print("Running OpenSim analyses: SO...\n")
# osp.opensim_pipeline(forcedb, user, ["so"])
# print("\nOpenSim analyses (SO) completed.\n")

# print("Running OpenSim analyses: RRA, CMC...\n")
# osp.opensim_pipeline(forcedb, user, ["rra",  "cmc"])
# print("\nOpenSim analyses (RRA, CMC) completed.\n")

# print("Running OpenSim analyses: JR...\n")
# osp.opensim_pipeline(forcedb, user, ["jr"])
# print("\nOpenSim analyses (JR) completed.\n")

#****** FOR TESTING ONLY ******
# osp.opensim_pipeline(forcedb, user, ["bk"], "FAILT99")
#******************************


# %% COLLATE AND EXPORT RESULTS

# import opensimresults as osr

# print("Converting OpenSim results to Pickle...\n")
# failedfiles = osr.opensim_results_batch_process(forcedb, ["ik", "id", "bk"], user, 101)
# print("\nOpenSim results converted to Pickle.\n")

# print("Exporting OpenSim results to CSV...\n")
# osr.export_opensim_results(forcedb, user, ["ik", "id"])
# print("CSV export complete.\n")

# print("Exporting normalised OpenSim results to CSV...\n")
# osr.export_opensim_results(forcedb, user, ["ik", "id"], True)
# print("Normalised CSV export complete.\n")

# print("Exporting normalised OpenSim results to CSV (subject descriptives)...\n")
# osr.export_opensim_results_subject_mean(forcedb, user, ["ik", "id"])
# print("Normalised CSV export complete.\n")

# print("Exporting normalised OpenSim results to CSV (subject descriptives)...\n")
# osr.export_opensim_results_subject_mean(forcedb, user, ["ik", "id"], True)
# print("Normalised CSV export complete.\n")


# %% ANALYSES

# import analyses as an

# print("Running post-hoc analyses...\n")
# an.analyses_batch_process(forcedb, user)
# print("Analyses complete.\n")

# print("Exporting analysis results...\n")
# an.export_joint_angular_impulse(forcedb, user)
# print("Analyses results export complete.\n")



# %% STABILITY

import stability

# print("Running stability analyses (MoS)...\n")
# stability.batch_process_stability(user, forcedb, treadmill_speed = 0.0)
# print("Stability analyses (MoS) complete.\n")

print("Exporting stability results (MoS, WBAM) to CSV...\n")
stability.export_margin_of_stability(forcedb, user, 101, False)
print("Stability CSV export complete.\n")

print("Exporting normalised stability results (MoS) to CSV...\n")
stability.export_margin_of_stability_subject_mean(forcedb, user, 101, False)
print("Normalised stability CSV export complete.\n")




# %% END

print("\n")
print("----------------------------------------")

# end time stamp
ts1 = datetime.datetime.now();
print("End: %s" % ts1)
datetime.datetime.now()

print("----------------------------------------")
print("\n")


