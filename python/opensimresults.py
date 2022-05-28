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
    def __init__(self, osimkey, analyses, nsamp):
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
        self.__get_results_split(analyses, nsamp)
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
        filext["jr"] = []
        
        # header rows
        # note: may differ from actual number of header rows as pandas skips
        # some blank initial rows
        headnum = {}
        headnum["ik"] = 8
        headnum["id"] = 6
        headnum["so"] = 10
        headnum["rra"] = []
        headnum["cmc"] = []
        headnum["jr"] = []
       
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
    
    def __get_results_split(self, analyses, nsamp):
        
        # initialise dict
        results = {}
        
        # left leg flip columns (incl. time): R, L
        flip = {}
        flip["ik"] = [3, 4, 7, 25, 26]
        flip["id"] = [3, 4, 7, 15, 16]
        flip["so"] = []
        flip["rra"] = []
        flip["cmc"] = []
        flip["jr"] = []  
        
        # foot columns (incl. time): R, L
        columns = {}
        columns["ik"] = [[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32], 
                         [0, 1, 2, 3, 4, 5, 6, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 33, 34, 35, 36, 37, 38, 39]]
        columns["id"] = [[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 13, 14, 15, 16, 17, 20, 21, 22, 26, 28, 30, 32, 34, 36, 37], 
                         [0, 1, 2, 3, 4, 5, 6, 10, 11, 12, 13, 14, 15, 18, 19, 23, 24, 25, 27, 29, 31, 33, 35, 38, 39]]
        columns["so"] = [[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90],
                         [0, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 91, 92, 93, 94, 95, 96, 97]]
        columns["rra"] = []
        columns["cmc"] = []
        columns["jr"] = []   
        
        # headers
        headers = {}
        headers["ik"] = ["time", "pelvis_tilt", "pelvis_list", "pelvis_rotation", "pelvis_tx", "pelvis_ty", "pelvis_tz", "hip_flexion", "hip_adduction", "hip_rotation", "knee_angle", "knee_angle_beta", "ankle_angle", "subtalar_angle", "mtp_angle", "lumbar_extension", "lumbar_bending", "lumbar_rotation", "arm_flex", "arm_add", "arm_rot", "elbow_flex", "pro_sup", "wrist_flex", "wrist_dev"]
        headers["id"] = ["time", "pelvis_tilt_moment", "pelvis_list_moment", "pelvis_rotation_moment", "pelvis_tx_force", "pelvis_ty_force", "pelvis_tz_force", "hip_flexion_moment", "hip_adduction_moment", "hip_rotation_moment", "lumbar_extension_moment", "lumbar_bending_moment", "lumbar_rotation_moment", "knee_angle_moment", "knee_angle_beta_force", "arm_flex_moment", "arm_add_moment", "arm_rot_moment", "ankle_angle_moment", "elbow_flex_moment", "subtalar_angle_moment", "pro_sup_moment", "mtp_angle_moment", "wrist_flex_moment", "wrist_dev_moment"]
        headers["so"] = ["time", "addbrev", "addlong", "addmagDist", "addmagIsch", "addmagMid", "addmagProx", "bflh", "bfsh", "edl", "ehl", "fdl", "fhl", "gaslat", "gasmed", "glmax1", "glmax2", "glmax3", "glmed1", "glmed2", "glmed3", "glmin1", "glmin2", "glmin3", "grac", "iliacus", "perbrev", "perlong", "piri", "psoas", "recfem", "sart", "semimem", "semiten", "soleus", "tfl", "tibant", "tibpost", "vasint", "vaslat", "vasmed", "lumbar_ext", "lumbar_bend", "lumbar_rot", "shoulder_flex", "shoulder_add", "shoulder_rot", "elbow_flex", "pro_sup", "wrist_flex", "wrist_dev"]
        headers["rra"] = []
        headers["cmc"] = []
        headers["jr"] = []       
        
        # get OpenSim data
        for ans in analyses:
            
            # skip scale
            if ans.casefold() == "scale": continue
        
            # initialise dict
            results[ans] = {}        
            
            # split by feet
            for f, foot in enumerate(["r", "l"]):
                                
                # copy raw data
                data0 = self.results["raw"][ans]["data"]
                
                # flip columns for left leg
                if f == 2:
                    data0[:, flip[ans]] = -1 * data0[:, flip[ans]]
                    
                # trim columns
                data0 = data0[:, columns[ans][f]]
                
                # ###################################
                # PROCESS EVENTS BASED ON TASK
                
                # match task and find time window for foot
                
                # static trial
                if self.task.casefold() == "static":
                    print("Static trial. Nothing to be done.")

                
                # run stride cycle on ipsilateral leg, stance on contralateral
                if self.task.casefold() == "run_stridecycle":
                    
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
                    e0 = self.events["labels"].index(foot.upper() + "FS")
                    e1 = self.events["labels"].index(foot.upper() + "FO")
                    t0 = self.events["time"][e0]
                    t1 = self.events["time"][e1]
                    
                    
                #
                # ###################################
                
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
def opensim_results_batch_process(meta, analyses, nsamp):
    
    # extract OpenSim data
    osimkey = {}
    for subj in meta:
    
        print("%s" % "*" * 30)
        print("SUBJECT: %s" % subj)
        print("%s" % "*" * 30)
        print("\n")
        
        for group in meta[subj]["trials"]:
            
            print("Group: %s" % group)
            print("%s" % "=" * 30)                      
            
            # process dynamic trials only
            for trial in  meta[subj]["trials"][group]:                
                
                # ignore static trials
                isstatic = meta[subj]["trials"][group][trial]["isstatic"]
                if isstatic: continue
            
                print("Dynamic trial: %s" % trial)
            
                # load the trial OsimKey
                c3dpath = meta[subj]["trials"][group][trial]["outpath"]
                pkfile = os.path.join(c3dpath,trial + "_osimkey.pkl")
                with open(pkfile,"rb") as fid:
                    osimkey = pk.load(fid)
                    
                # get the OpenSim results
                osimresultskey = OsimResultsKey(osimkey, analyses, nsamp)
                
                # save OsimResultsKey to file
                with open(os.path.join(c3dpath, trial + "_opensim_results.pkl"),"wb") as f:
                    pk.dump(osimresultskey, f)
                    
    
    print("\n")                
    
    return None
    
    

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
    for subj in meta:
    
        print("%s" % "*" * 30)
        print("SUBJECT: %s" % subj)
        print("%s" % "*" * 30)
        print("\n")
                
        for group in meta[subj]["trials"]:
            
            print("Group: %s" % group)
            print("%s" % "=" * 30)                      
            
            # process dynamic trials only
            for trial in  meta[subj]["trials"][group]:                
                
                # ignore static trials
                isstatic = meta[subj]["trials"][group][trial]["isstatic"]
                if isstatic: continue
            
                print("Dynamic trial: %s" % trial)
            
                # load the trial OsimResultsKey
                c3dpath = meta[subj]["trials"][group][trial]["outpath"]
                pkfile = os.path.join(c3dpath,trial + "_opensim_results.pkl")
                with open(pkfile,"rb") as fid:
                    osimresultskey = pk.load(fid)
                    
                # trial task and condition
                task = osimresultskey.task
                condition = osimresultskey.condition
                
                # foot
                for f, foot in enumerate(["r","l"]):
                    
                    # leg task
                    leg_task = osimresultskey.events["leg_task"][f]
                    
                    # analysis
                    for ans in analyses:
                        
                        # ignore scaling
                        if ans.casefold() == "scale": continue
                    
                        # data array
                        data = osimresultskey.results["split"][ans][foot]["data"]
                        varheader = osimresultskey.results["split"][ans][foot]["headers"]
                    
                        # variable
                        for v, variable in enumerate(varheader):
                            
                            # ignore time
                            if v == 0: continue

                            # data for the variable
                            drow = data[:, v]

                            # create new line of data
                            csvrow = [subj, group, trial, task, condition, foot, leg_task, ans, variable] + drow.tolist()
                            csvdata.append(csvrow)


    # create empty dataframe
    print("\nCreating dataframe...")
    headers = ["subject", "group", "trial", "task", "condition", "foot", "leg_task", "analysis", "variable"] + ["t" + str(n) for n in range(1,102)]
    csvdf = pd.DataFrame(csvdata, columns = headers)

    # write data to file with headers
    print("\nWriting to CSV text file...")
    csvfile = user.csvfileprefix + task + ".csv"
    fpath = os.path.join(user.rootpath, user.outfolder, user.csvfolder)
    if not os.path.exists(fpath): os.makedirs(fpath)
    csvdf.to_csv(os.path.join(fpath,csvfile), index = False)
    
    
    print("\n")
   
    return csvdf



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