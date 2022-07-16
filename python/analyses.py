# -*- coding: utf-8 -*-
"""
OpenSim post-hoc analyses

@author: Prasanna Sritharan
"""



import os
import numpy as np
import pickle as pk
import pandas as pd



'''
-----------------------------------
------------- CLASSES -------------
-----------------------------------
'''

# NO CLASSES




'''
-----------------------------------
------- FUNCTIONS: ANALYSES -------
-----------------------------------
'''


'''
analyses_batch_process(meta, user):
    Batch process post-hoc analyses.
'''
def analyses_batch_process(meta, user):
    
    # extract OpenSim data
    failedfiles = []
    for subj in meta:
    
        print("\n")
        print("%s" % "*" * 30)
        print("SUBJECT: %s" % subj)
        print("%s" % "*" * 30)
        
        for group in meta[subj]["trials"]:
            
            print("Group: %s" % group)
            print("%s" % "=" * 30)                      
            
            # process dynamic trials only
            for trial in  meta[subj]["trials"][group]:                
                
                # ignore static trials
                isstatic = meta[subj]["trials"][group][trial]["isstatic"]
                if isstatic: continue

                try:

                    # load results key
                    pklpath = meta[subj]["trials"][group][trial]["outpath"]
                    with open(os.path.join(pklpath, trial + "_opensim_results.pkl"), "rb") as fid0:
                        osimresultskey = pk.load(fid0)
                            
                    # create analyses dict key if it does not exist
                    if "analyses" not in osimresultskey.results:
                        osimresultskey.results["analyses"] = {}
                                    
                    print("Dynamic trial: %s" % trial)
                    
                    
                    
                    # ******************************
                    # ANALSYES
    
                    # joint angular impulse, incl. pelvis
                    print("---> Joint angular impulse")
                    osimresultskey = calculate_joint_angular_impulse(osimresultskey, user)                
    
    
                    #                
                    # ******************************
     
                    # pickle updated results key
                    with open(os.path.join(pklpath, trial + "_opensim_results.pkl"), "wb") as fid1:
                        pk.dump(osimresultskey, fid1)   

                except:
                    print("*** FAILED ***")
                    failedfiles.append(trial)                    
                        
                          
    print("\n")                
    
    return failedfiles    



'''
calculate_joint_angular_impulse(osimresultskey, user):
    Calculate the joint angular impulse for each time window specified by
    events, report positive and negative separately.
'''
def calculate_joint_angular_impulse(osimresultskey, user):
    
    # get event times
    events = osimresultskey.events["time"]
    
    # get joint moments
    Tdata = osimresultskey.results["split"][user.idcode]
    
    # calculate joint angular impulse on each leg
    impl = {}
    for leg in user.leg:
        
        # data
        time = Tdata[leg]["data"][:, 0]
        T = Tdata[leg]["data"][:, 1:]
        headers = Tdata[leg]["headers"][1:]
        
        # net
        impl[leg] = {}
        impl[leg]["net"] = {}
        impl[leg]["net"]["headers"] = headers
        impl[leg]["net"]["net"] = np.trapz(T, time, axis = 0)
        
        # net positive
        Tpos = T.copy()
        Tpos[Tpos < 0] = 0
        impl[leg]["net"]["pos"] = np.trapz(Tpos, time, axis = 0)
        
        # net negative
        Tneg = T.copy()
        Tneg[Tneg > 0] = 0
        impl[leg]["net"]["neg"] = np.trapz(Tneg, time, axis = 0)  
        
        # events
        impl[leg]["windows"] = {}
        impl[leg]["windows"]["headers"] = headers
        for e in range(0, len(events) - 1):
            
            # indices for time window
            idx0 = np.where(time >= events[e])[0][0]
            idx1 = np.where(time <= events[e + 1])[0][-1]
            
            # event window data
            timewin = time[idx0:idx1]
            Twin = T[idx0:idx1, :]
            
            # net window
            wlabel = "w" + str(e + 1)
            impl[leg]["windows"][wlabel] = {}
            impl[leg]["windows"][wlabel]["net"] = np.trapz(Twin, timewin, axis = 0)
            
            # net window positive
            Twinpos = Tpos[idx0:idx1, :]
            impl[leg]["windows"][wlabel]["pos"] = np.trapz(Twinpos, timewin, axis = 0)
            
            # net window negative
            Twinneg = Tneg[idx0:idx1, :]
            impl[leg]["windows"][wlabel]["neg"] = np.trapz(Twinneg, timewin, axis = 0)
           
    
    # append to results key
    osimresultskey.results["analyses"]["joint_angular_impulse"] = impl
    
    return osimresultskey
    
  
  
  
'''
-----------------------------------
-------- FUNCTIONS: EXPORT --------
-----------------------------------
'''


'''
export_joint_angular_impulse(meta, user):
    Write joint angular impulses to Excel.
'''
def export_joint_angular_impulse(meta, user):
        
    # empty output list of lists
    # (create the output table as a list of lists, then convert to dataframe
    # as iteratively appending new dataframe rows is computationally expensive)
    csvdata = []   
    
    # extract OpenSim data
    print("Collating data into lists...\n")
    failedfiles = []
    for subj in meta:
    
        print("\n")
        print("%s" % "*" * 30)
        print("SUBJECT: %s" % subj)
        print("%s" % "*" * 30)
        
        for group in meta[subj]["trials"]:
            
            print("Group: %s" % group)
            print("%s" % "=" * 30)                      
            
            # process dynamic trials only
            for trial in  meta[subj]["trials"][group]:                
                
                # ignore static trials
                isstatic = meta[subj]["trials"][group][trial]["isstatic"]
                if isstatic: continue

                try:

                    # load results key
                    pklpath = meta[subj]["trials"][group][trial]["outpath"]
                    with open(os.path.join(pklpath, trial + "_opensim_results.pkl"), "rb") as fid0:
                        osimresultskey = pk.load(fid0)
                        
                        
                    # movement task and condition
                    movement = osimresultskey.task
                    
                    # windows
                    windows = ["w" + str(e) for e in range(1, len(osimresultskey.events["labels"]))]
                    
                    # pivot leg
                    pivot_leg = osimresultskey.events["labels"][0][0].lower()
                    
                     
                    # foot
                    for f, foot in enumerate(["r","l"]):
                        
                        # pivot leg or non-pivot leg data
                        if foot == pivot_leg:
                            legtask = "pivot"
                        else:
                            legtask = "nonpivot"
                        
                        # joint angular impulse full results for leg
                        full_results = osimresultskey.results["analyses"]["joint_angular_impulse"][foot]
                        
                        # get data
                        implabels = ["net", "windows"]
                        for implabel in implabels:
                            
                            # net impulse
                            if implabel == "net":
                                winlabel = "net"
                                for v, varlabel in enumerate(full_results["net"]["headers"]):
                                    
                                    # skip pelvis and knee beta force
                                    if ("pelvis" in varlabel) or ("knee_angle_beta_force" in varlabel): continue                                    
                                    
                                    # values
                                    net = full_results["net"]["net"][v]
                                    pos = full_results["net"]["pos"][v]
                                    neg = full_results["net"]["neg"][v]

                                    # csv row
                                    csvrow = [subj, group, trial, movement, foot, legtask, implabel, winlabel, varlabel, net, pos, neg]
                                    csvdata.append(csvrow)
                                    
                            elif implabel == "windows":
                                for winlabel in windows:
                                    for v, varlabel in enumerate(full_results["windows"]["headers"]):
                                        
                                        # skip pelvis
                                        if "pelvis" in varlabel: continue
                                        
                                        # values
                                        net = full_results["windows"][winlabel]["net"][v]
                                        pos = full_results["windows"][winlabel]["pos"][v]
                                        neg = full_results["windows"][winlabel]["neg"][v]
    
                                        # csv row
                                        csvrow = [subj, group, trial, movement, foot, legtask, implabel, winlabel, varlabel, net, pos, neg]
                                        csvdata.append(csvrow)                                    
                                    
                except:
                    print("Dynamic trial: %s *** FAILED ***" % trial)
                    failedfiles.append(trial)
                else:
                    print("Dynamic trial: %s" % trial)        
    
    
    # create dataframe
    # appearing in the output
    print("\nCreating dataframe...")
    headers = ["subject", "group", "trial", "movement", "foot", "leg_task", "type", "window", "variable", "net", "positive", "negative"]
    csvdf = pd.DataFrame(csvdata, columns = headers)

    # write data to file with headers
    print("\nWriting to CSV text file...")
    csvfile = user.project.lower() + "_joint_angular_impulse_all.csv"
    fpath = os.path.join(user.rootpath, user.outfolder, user.csvfolder)
    if not os.path.exists(fpath): os.makedirs(fpath)
    csvdf.to_csv(os.path.join(fpath,csvfile), index = False)
    
    print("\n")                
    
    return None       



    
