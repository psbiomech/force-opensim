# -*- coding: utf-8 -*-
"""
Build output database for OpenSim projects

@author: Prasanna Sritharan
"""

import glob
import os
import re
import shutil
import pickle as pk



'''
-----------------------------------
------------- CLASSES -------------
-----------------------------------
'''

# NO CLASSES



'''
-----------------------------------
------------ FUNCTIONS ------------
-----------------------------------
'''


'''
build_database(user, task):
    Build OpenSim output database from user, return metadata dict.
'''
def build_database(user, task):
 
    # meta dict
    meta = {}
 
    # input folders
    for infolder in user.infolder:
 
        # parse database, get subject list
        inpath = os.path.join(user.rootpath, infolder, user.subjprefix + "*")
        folderlist = glob.glob(inpath, recursive = True)
        subjlist = [os.path.split(f)[1] for f in folderlist]
        
        # build metadata dict
        outpath = os.path.join(user.rootpath, user.outfolder)
        fnpatobj = re.compile(user.fnpat) 
        failedfiles = []
        for n, subj in enumerate(subjlist):        
            
            # basic info
            meta[subj] = {}
            meta[subj]["subj"] = subj
            meta[subj]["project"] = user.project
            meta[subj]["outpath"] = os.path.join(outpath, subj)
            
            # trial subfolders
            meta[subj]["trials"] = {}
            for group in user.trialgroupfolders:
                
                # parse subfolders
                meta[subj]["trials"][group] = {}
                groupinpath = os.path.join(folderlist[n], group, "*.c3d")
                groupfolderlist = glob.glob(groupinpath)
                triallist = [os.path.splitext(os.path.split(f)[1])[0] for f in groupfolderlist]
                
                # skip groups that don't have a least one static trial and at
                # least one required dynamic file
                hasstatic = any([user.staticprefix.casefold() in t.casefold() for t in triallist])
                hasdynamic = any([any([c.casefold() in t.casefold() for c in user.trialprefixes["sdp"]]) for t in triallist])
                if not(hasstatic) or not(hasdynamic): continue
                
                # trials (for selected task only)
                for m, trial in enumerate(triallist):            
                    
                    trial = trial.upper()
                    
                    # file name tokens, skip those that don't match regex
                    trialtoks = fnpatobj.fullmatch(trial)
                    if trialtoks is None: continue
                    
                    # build meta data dict
                    trialprefix = trialtoks.group(user.tasktoknum)                    
                    if (trialprefix.casefold() == user.staticprefix.casefold()) or (trialprefix.casefold() in [t.casefold() for t in user.trialprefixes[task.casefold()]]):                                   
                        meta[subj]["trials"][group][trial] = {}                    
                        meta[subj]["trials"][group][trial]["trial"] = trial
                        meta[subj]["trials"][group][trial]["c3dfile"] = trial + ".c3d"
                        meta[subj]["trials"][group][trial]["osim"] = subj.upper() + ".osim"
                        meta[subj]["trials"][group][trial]["inpath"] = os.path.split(groupfolderlist[m])[0]
                        meta[subj]["trials"][group][trial]["outpath"] = os.path.join(outpath, subj, group, trial)
                        meta[subj]["trials"][group][trial]["task"] = task
                        meta[subj]["trials"][group][trial]["condition"] = trialprefix.casefold()
                        meta[subj]["trials"][group][trial]["isstatic"] = False
                        meta[subj]["trials"][group][trial]["usedstatic"] = False
                        if trialprefix.casefold() == user.staticprefix.casefold():
                            meta[subj]["trials"][group][trial]["task"] = "static"
                            meta[subj]["trials"][group][trial]["condition"] = "static"
                            meta[subj]["trials"][group][trial]["isstatic"] = True
                            if trial.casefold().endswith(user.staticused.casefold()):
                                meta[subj]["trials"][group][trial]["usedstatic"] = True
                    
                                
                # determine which file to use as static trial in OpenSim, in most
                # cases use the static file set in the user settings, if not found,
                # then use the first available, or manually override later
                finaltriallist = meta[subj]["trials"][group].keys()           
                for trial in finaltriallist:             
                    if trial.casefold().endswith(user.staticused.casefold()):
                        meta[subj]["trials"][group][trial]["usedstatic"] = True                            
                        break
                    elif user.staticprefix.casefold() in trial.casefold():
                        meta[subj]["trials"][group][trial]["usedstatic"] = True
                        break
                    
        # create subdfolders if required, copy C3D files into output database
        if not os.path.exists(outpath): os.makedirs(outpath)
        for subj in meta:
            for group in meta[subj]["trials"]:                
                for trial in meta[subj]["trials"][group]:
                    trialoutpath = meta[subj]["trials"][group][trial]["outpath"]
                    if not os.path.exists(meta[subj]["outpath"]): os.makedirs(meta[subj]["outpath"])
                    if not os.path.exists(os.path.join(meta[subj]["outpath"], group)): os.makedirs(os.path.join(meta[subj]["outpath"], group))
                    if not os.path.exists(trialoutpath): os.makedirs(trialoutpath)
                    shutil.copy(os.path.join(meta[subj]["trials"][group][trial]["inpath"], meta[subj]["trials"][group][trial]["c3dfile"]), trialoutpath)
                
        # clean up, remove empty subject meta dict keys
        for subj in subjlist:
            if not any([bool(meta[subj]["trials"][g]) for g in meta[subj]["trials"]]):
                meta.pop(subj)
                        
    # save the metadata dict
    with open(os.path.join(outpath, user.project + ".pkl"),"wb") as fid: pk.dump(meta, fid)
    
    return meta, failedfiles



'''
-----------------------------------
-------- MANUAL OVERRIDES ---------
-----------------------------------
'''


'''
manual_usedstatic(meta, user, subj, group, actual_static_trial, reset_first):
    Some subjects may not have a static trial with the suffix designated in
    user.staticused, e.g. user.staticused = "STATIC01", but subject SUBJ9 only
    has one static SUBJ01_STATIC02, then this should be used for OpenSim. This 
    function sets the usedstatic flag for this trial in the meta dict.
'''
def manual_usedstatic(meta, user, subj, group, actual_static_trial, reset_first):
    
    # reset all static trials
    if reset_first:
        triallist = meta[subj]["trials"][group].keys()
        for trial in triallist:
            if user.staticprefix.casefold() in trial.casefold():
                meta[subj]["trials"][group][trial]["usedstatic"] = False
            
        
    # set the required static trial
    meta[subj]["trials"][group][actual_static_trial]["usedstatic"] = True
    

