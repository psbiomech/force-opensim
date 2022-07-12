# -*- coding: utf-8 -*-
"""
OpenSim post-hoc analyses

@author: Prasanna Sritharan
"""



import os
import numpy as np
import pickle as pk






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
        
        # net
        impl[leg] = {}
        impl[leg]["net"] = {}
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



