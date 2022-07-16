# -*- coding: utf-8 -*-
"""
Load and format OpenSim results

@author: Prasanna Sritharan
"""

import os
import pandas as pd
import numpy as np
import pickle as pk
from scipy.interpolate import interp1d


'''
-----------------------------------
------------- CLASSES -------------
-----------------------------------
'''



'''
OsimResultsKey:
    Data storage class containing all OpenSim output data, raw and normalised
    to both BW and %BW*HT. Data is resample to the desired number of samples.
'''
class OsimResultsKey():
    def __init__(self, osimkey, analyses, user, nsamp):
        self.subject = osimkey.subject
        self.trial = osimkey.trial
        self.age = osimkey.age
        self.mass = osimkey.mass
        self.model = osimkey.model
        self.lab = osimkey.lab
        self.task = osimkey.task
        self.condition = osimkey.condition
        self.events = osimkey.events
        self.outpath = osimkey.outpath
        self.__get_results_raw(osimkey, analyses, nsamp)
        self.__get_results_split(osimkey, analyses, user, nsamp)
        return None
        
    def __get_results_raw(self, osimkey, analyses, nsamp):
        
        # initialise dict
        results = {}
       
        # file suffix and extension
        filext = {}
        filext["ik"] = "_ik.mot"
        filext["id"] = "_id.sto"
        filext["so"] = "_so_force.sto"
        filext["rra"] = []
        filext["cmc"] = []
        filext["jr"] = "_jr_ReactionLoads.sto"
        
        # header rows
        # note: may differ from actual number of header rows as pandas skips
        # some blank initial rows
        headnum = {}
        headnum["ik"] = 8
        headnum["id"] = 6
        headnum["so"] = 10
        headnum["rra"] = []
        headnum["cmc"] = []
        headnum["jr"] = 9
       
        # get OpenSim data
        for ans in analyses:
            
            # skip scale
            if ans.casefold() == "scale": continue
            
            # load data 
            datafile = os.path.join(osimkey.outpath, ans, osimkey.trial + filext[ans])
            datadf = pd.read_csv(datafile, sep="\t", header=headnum[ans])
            headers = datadf.columns.tolist()
            data = datadf.to_numpy()
            
            # resample data
            datanew = resample1d(data, nsamp)
            
            # store in dict
            results[ans] = {}
            results[ans]["data"] = datanew
            results[ans]["headers"] = headers
        
        self.results = {}
        self.results["raw"] = results        
            
        return None
    
    def __get_results_split(self, osimkey, analyses, user, nsamp):
        
        # initialise dict
        results = {}
        
        # results processing parameters
        flip = user.results_flip
        columns = user.results_columns
        headers = user.results_headers
        
        # get OpenSim data
        for ans in analyses:
            
            # skip scale
            if ans.casefold() == "scale": continue
        
            # initialise dict
            results[ans] = {}        
            
            # split by feet
            for f, foot in enumerate(["r", "l"]):
                                
                # copy raw data
                data0 = None
                data0 = self.results["raw"][ans]["data"].copy()                
                
                
                # ###################################
                # ADDITIONAL PROCESSING BASED ON TASK
                #
                # Includes left leg trial flipping, as the way this is done can
                # be trial-dependent
                
                # match task and find time window for foot
                
                # static trial
                if self.task.casefold() == "static":
                    print("Static trial. Nothing to be done.")

                
                # run stride cycle on ipsilateral leg, stance on contralateral
                elif self.task.casefold() == "run_stridecycle":
                    
                    # flip columns for left leg trials
                    if f == 2:
                        data0[:, flip[ans]] = np.multiply(data0[:, flip[ans]], -1)
                    
                    # time window depends on leg task
                    if self.events["leg_task"][f] == "run_stridecycle":
                        e0 = self.events["labels"].index(foot.upper() + "FS")
                        e1 = e0 + 4
                        t0 = self.events["time"][e0]
                        t1 = self.events["time"][e1]
                    else:
                        e0 = self.events["labels"].index(foot.upper() + "FS")
                        e1 = self.events["labels"].index(foot.upper() + "FO")
                        t0 = self.events["time"][e0]
                        t1 = self.events["time"][e1]                        
                
                
                # run stance on both legs
                elif self.task.casefold() == "run_stance":
                    
                    # flip columns for left leg trials
                    if f == 2:
                        data0[:, flip[ans]] = np.multiply(data0[:, flip[ans]], -1)
                        
                    # time window
                    e0 = self.events["labels"].index(foot.upper() + "FS")
                    e1 = self.events["labels"].index(foot.upper() + "FO")
                    t0 = self.events["time"][e0]
                    t1 = self.events["time"][e1]
                    
                
                # step down and pivot
                elif self.task.casefold() == "sdp":
                    
                    # flip columns for left-foot-first left-turning trials
                    if osimkey.events["labels"][0][0].casefold() == "l":
                        data0[:, flip[ans]] = np.multiply(data0[:, flip[ans]], -1)                   
                    
                    # time window
                    e0 = 0
                    e1 = 5
                    t0 = self.events["time"][e0]
                    t1 = self.events["time"][e1]                    

                    
                #
                # ###################################

                # trim columns
                data0 = data0[:, columns[ans][f]].copy()
                
                # trim rows (time window)
                r00 = np.where(data0[:, 0] <= t0)[0]
                if r00.size == 0:
                    r0 = 0
                else:
                    r0 = r00[-1]
                r1 = np.where(data0[:, 0] <= t1)[0][-1]
                data1 = data0[r0:r1 + 1, :]
                
                # resample data, currently uses simple 1D interpolation but
                # need to find a package that emulates Matlab's resample()
                data = resample1d(data1, nsamp)
                            
                # store in dict
                results[ans][foot] = {}
                results[ans][foot]["data"] = data
                results[ans][foot]["headers"] = headers[ans]
        
        self.results["split"] = results        
            
        return None        
        
        

'''
-----------------------------------
------------ FUNCTIONS ------------
-----------------------------------
'''


'''
opensim_results_batch_process(meta, analyses, nsamp):
    Batch process OpenSim results text files to OsimResultsKeys.
'''
def opensim_results_batch_process(meta, analyses, user, nsamp):
    
    # extract OpenSim data
    osimkey = {}
    failedfiles = []
    for subj in meta:
    
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
                            
                    # load the trial OsimKey
                    c3dpath = meta[subj]["trials"][group][trial]["outpath"]
                    pkfile = os.path.join(c3dpath, trial + "_osimkey.pkl")
                    with open(pkfile, "rb") as fid:
                        osimkey = pk.load(fid)
                        
                    # get the OpenSim results
                    osimresultskey = OsimResultsKey(osimkey, analyses, user, nsamp)
                    
                    # save OsimResultsKey to file
                    with open(os.path.join(c3dpath, trial + "_opensim_results.pkl"), "wb") as f:
                        pk.dump(osimresultskey, f)
                                    
                except:
                    print("Dynamic trial: %s *** FAILED ***" % trial)
                    failedfiles.append(trial)                    
                else:
                    print("Dynamic trial: %s" % trial)
                          
    print("\n")                
    
    return failedfiles
    
    

'''
export_opensim_results(meta, user, analyses):
    Collate OpenSim results into dataframes and export to text for Rstats.
'''
def export_opensim_results(meta, user, analyses):
    
    # empty output list of lists
    # (create the output table as a list of lists, then convert to dataframe
    # as iteratively appending new dataframe rows is computationally expensive)
    csvdata = []
        
    # extract OpenSim data
    print("Collating data into lists...\n")
    failedfiles = []
    for subj in meta:
    
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
                
                    # load the trial OsimResultsKey
                    c3dpath = meta[subj]["trials"][group][trial]["outpath"]
                    pkfile = os.path.join(c3dpath,trial + "_opensim_results.pkl")
                    with open(pkfile,"rb") as fid:
                        osimresultskey = pk.load(fid)
                        
                    # trial task and condition
                    task = osimresultskey.task

                    # pivot leg
                    pivot_leg = osimresultskey.events["labels"][0][0].lower()
                    
                    # generic event labels:
                    #   PFO1: pivot limb foot off FP1
                    #   PFS2: pivot limb foot strike FP2
                    #   NFO1: non-pivot limb foot off FP1
                    #   etc...
                    events_gen_labels = ["PFO1", "PFS2", "NFO1", "NFS3", "PFO2", "PFS4"]
                    
                    # event timing
                    events_times = osimresultskey.events["time"]
                    events_steps = np.round(user.samples * (osimresultskey.events["time"] - osimresultskey.events["time"][0]) / (osimresultskey.events["time"][5] - osimresultskey.events["time"][0]))
                    
                    # foot
                    for f, foot in enumerate(["r","l"]):

                        # pivot leg or non-pivot leg data
                        if foot == pivot_leg:
                            data_leg_role = "pivot"
                        else:
                            data_leg_role = "nonpivot"
                        
                        # analysis
                        for ans in analyses:
                            
                            # ignore scaling
                            if ans.casefold() == "scale": continue
                        
                            # data array
                            data = osimresultskey.results["split"][ans][foot]["data"]
                            varheader = osimresultskey.results["split"][ans][foot]["headers"]
                        
                            # variable
                            for v, variable in enumerate(varheader):
    
                                # data for the variable (includes time)
                                drow = data[:, v]
    
                                # create new line of data
                                csvrow = [subj, trial, task, foot, data_leg_role] + events_times.tolist() + events_steps.tolist() + [ans, variable] + drow.tolist()
                                csvdata.append(csvrow)
                
                except:
                    print("Dynamic trial: %s *** FAILED ***" % trial)
                    failedfiles.append(trial)
                else:
                    print("Dynamic trial: %s" % trial)

    # create dataframe
    print("\nCreating dataframe...")
    headers = ["subject", "trial", "task", "data_leg", "data_leg_role"] + ["et" + str(e + 1) + "_" + ev for e, ev in enumerate(events_gen_labels)] + ["es" + str(e + 1) + "_" + ev for e, ev in enumerate(events_gen_labels)] + ["analysis", "variable"] + ["t" + str(n) for n in range(1,102)]
    csvdf = pd.DataFrame(csvdata, columns = headers)

    # write data to file with headers
    print("\nWriting to CSV text file...")
    csvfile = user.csvfileprefix + ".csv"
    fpath = os.path.join(user.rootpath, user.outfolder, user.csvfolder)
    if not os.path.exists(fpath): os.makedirs(fpath)
    csvdf.to_csv(os.path.join(fpath,csvfile), index = False)
    
    
    print("\n")
   
    return failedfiles



'''
resample1d(data, nsamp):
    Simple resampling by 1-D interpolation (rows = samples, cols = variable).
    Data can be a 1-D or multiple variables in a 2D array-like object.
'''
def resample1d(data, nsamp):

    # convert list to array
    if isinstance(data, list):
        data = np.array([data]).transpose()
        ny = 1        

    # data dimensions
    nx = data.shape[0]
    ny = data.shape[1]

    # old sample points
    x = np.linspace(0, nx - 1, nx)
        
    # new sample points
    xnew = np.linspace(0, nx - 1, nsamp)
    
    # parse columns
    datanew = np.zeros([nsamp, ny])
    for col in range(0, ny):
        
        # old data points
        y = data[:, col]
        
        # convert to cubic spline function
        fy = interp1d(x, y, kind = "cubic", fill_value = "extrapolate")
    
        # new data points
        ynew = fy(xnew)
    
        # store column
        datanew[:, col] = ynew
        
    
    return datanew