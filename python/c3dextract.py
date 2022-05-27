# -*- coding: utf-8 -*-
"""
LASEM C3D file data extract

@author: Prasanna Sritharan
"""

import scipy.signal as signal
import scipy.interpolate as interp
import statistics as stats
import numpy as np
import pyc3dserver as c3d
import pickle as pk
import opensimpipeline as osp
import os
import re



'''
-----------------------------------
------------- CLASSES -------------
-----------------------------------
'''

'''
C3DKey:
    C3D data storage class containing all C3D file data.
'''
class C3DKey():
    def __init__(self, sname, tname, fmeta, fforces, fmarkers):
        self.subject_name = sname
        self.trial_name = tname
        self.meta = fmeta
        self.forces = fforces
        self.markers = fmarkers
        return None



'''
TrialKey:
    Trial data storage class containing all C3D data relevant to trial in the
    laboratory and force plate frames.
'''
class TrialKey():
    def __init__(self, lab, task, condition, c3dkey, xdir, static_fp_channel, mass, staticprefix):       
        self.subject_name = str(c3dkey.subject_name)
        self.trial_name = str(c3dkey.trial_name)
        self.lab_name = lab.lab_name
        self.task = task
        self.condition = condition
        self.mass = mass
        self.__set_events(c3dkey, task, static_fp_channel)
        self.__set_markers(lab, c3dkey, xdir)        
        self.__set_force_plates(lab, c3dkey, xdir, staticprefix) 
        self.__set_forces(lab, c3dkey, staticprefix)
        return None

    def __set_events(self, c3dkey, task, static_fp_channel):

        # initialise dict        
        events = {}
              
        # process events, or if no events, add general events
        if c3dkey.meta["EVENT"]["USED"] == 0:            
            events["labels"] = ["GEN", "GEN"]
            events["time"] = [c3dkey.markers["TIME"][0], c3dkey.markers["TIME"][-1]]            
        
        else:
            
            # build events list and time
            foot = [f[0].upper() for f in c3dkey.meta["EVENT"]["CONTEXTS"]]
            elabels = [foot[i] + "F" + f.split()[1][0] for i, f in enumerate(c3dkey.meta["EVENT"]["LABELS"])]
            etime = c3dkey.meta["EVENT"]["TIMES"][:,1]
            
            # sort the events list and time as sometimes the C3D stores events
            # and times out of order in its meta data
            sortidxs = np.argsort(etime)
            events["labels"] = [elabels[e] for e in sortidxs]
            events["time"] = [etime[e] for e in sortidxs]
                    
        # relative time, normalise to first frame
        events["time0"] = events["time"] - ((c3dkey.meta["TRIAL"]["ACTUAL_START_FIELD"][0] - 1) / c3dkey.meta["TRIAL"]["CAMERA_RATE"])
            
      
        
        # ###################################
        # PROCESS EVENTS BASED ON TASK
        
        # match task and perform required computations
        # (unfortunately we cannot use match-case before Python 3.10)
 
        
        # static trials
        if task.casefold() == "static":
            
            # calculate subject mass, return end frames for 25%-75% window
            mass, fidx0, fidx1 = calculate_subject_mass(c3dkey, static_fp_channel)
            self.mass = mass
            
            # time window
            events["window_time0"] = events["time0"][fidx0:fidx1 + 1]
            events["window_labels"] = events["labels"][fidx0:fidx1 + 1]
            
            # no force plate sequence
            events["fp_sequence"] = [[0, 0]]
            
            # leg task is static (R, L)
            events["leg_task"] = ["static", "static"]

            # last event index (0-based) for OpenSim analyses that require
            # kinetics (e.g., ID, SO, RRA and CMC)
            events["opensim_last_event_idx"] = -1
            
        # run full stride cycle
        elif task.casefold().startswith("run_stridecycle"):
            
            # some trials will have 2 full stride cycles, some will only have
            # one, so in the latter case the contralateral leg will report
            # stance only
            
            # calculate the time window of interest (assume 7 events means a 
            # full stride cycle is available on each leg, less then 7 events
            # assume stance on contralateral leg)
            fsidx0 = np.where(np.char.find(events["labels"],"FS")>=0)[0][0]
            if len(events["labels"]) < 7:
                fsidx1 = np.where(np.char.find(events["labels"],"FS")>=0)[0][2]
            else:
                fsidx1 = np.where(np.char.find(events["labels"],"FS")>=0)[0][3]                        
            events["window_time0"] = events["time0"][fsidx0:fsidx1 + 1]
            events["window_labels"] = events["labels"][fsidx0:fsidx1 + 1]
            
            # list the individual intervals
            events["window_intervals0"] = np.array([[t0, t1] for t0, t1 in zip(events["window_time0"][0:-1], events["window_time0"][1:])])
            
            # find force plate sequence for each interval (row) defined in the
            # array window_intervals0
            # note: sequences defined as per GaitExtract using a 2D array:
            #   rows: event intervals
            #   col1: right foot
            #   col2: left foot
            events["fp_sequence"] = [[0, 0]]
            if len(events["labels"]) < 7:
                if events["window_labels"][0][0] == "R":
                    events["fp_sequence"] = np.array([[4, 0], [0, 0], [0, 3], [0, 0]])
                else:
                    events["fp_sequence"] = np.array([[0, 4], [0, 0], [3, 0], [0, 0]])
            else:
                if events["window_labels"][0][0] == "R":
                    events["fp_sequence"] = np.array([[4, 0], [0, 0], [0, 3], [0, 0], [0, 0], [0, 0]])
                else:
                    events["fp_sequence"] = np.array([[0, 4], [0, 0], [3, 0], [0, 0], [0, 0], [0, 0]])                

            # leg task is stride cycle on ipsilateral, stance on contralateral
            if len(events["labels"]) < 7:
                if events["window_labels"][0][0] == "R":
                    events["leg_task"] = ["run_stridecycle", "run_stance"]
                else:
                    events["leg_task"] = ["run_stance", "run_stridecycle"]
            else:
                events["leg_task"] = ["run_stridecycle", "run_stridecycle"]              
            
            # last event index (0-based) for OpenSim analyses that require
            # kinetics (e.g., ID, SO, RRA and CMC)
            events["opensim_last_event_idx"] = 4
            

        # run stance phase only, both legs
        elif task.casefold().startswith("run_stance"):
                                        
            # assume first pair of stance phases (IFS ---> CFO) is the required 
            # window, thus will have one stance phase per leg.
            
            # calculate the time window of interest (assume first FS is start
            # of the trial time window, second FO is end of window)
            fsidx0 = np.where(np.char.find(events["labels"],"FS")>=0)[0][0]
            foidx1 = np.where(np.char.find(events["labels"],"FO")>=0)[0][1]
            events["window_time0"] = events["time0"][fsidx0:foidx1 + 1]
            events["window_labels"] = events["labels"][fsidx0:foidx1 + 1]
            
            # list the individual intervals
            events["window_intervals0"] = np.array([[t0, t1] for t0, t1 in zip(events["window_time0"][0:-1], events["window_time0"][1:])])
            
            # find force plate sequence for each interval (row) defined in the
            # array window_intervals0
            # note: sequences defined as per GaitExtract using a 2D array:
            #   rows: event intervals
            #   col1: right foot
            #   col2: left foot
            events["fp_sequence"] = [[0, 0]]
            if events["window_labels"][0][0] == "R":
                events["fp_sequence"] = np.array([[4, 0], [0, 0], [0, 3]])
            else:
                events["fp_sequence"] = np.array([[0, 4], [0, 0], [3, 0]])
                
            # leg task is same for both legs  (R, L)
            events["leg_task"] = ["run_stance", "run_stance"]
            
            # last event index (0-based) for OpenSim analyses that require
            # kinetics (e.g., ID, SO, RRA and CMC)
            events["opensim_last_event_idx"] = 3
            
        #
        # ###################################
        

        self.events = events
        
        return None
        
    def __set_markers(self, lab, c3dkey, xdir):
        
        # initialise dict
        markers = {}
        markers["labels"] = c3dkey.markers["LABELS"]
        
        # marker rate
        markers["rate"] = c3dkey.markers["RATE"]    
        orig_units = c3dkey.markers["UNITS"]
        
        # marker scale factor, convert to m
        markers["units"] = "m"
        if orig_units == "mm":
            markers["scale"] = 0.001
        elif orig_units == "m":
            markers["scale"] = 1
        else:
            markers["scale"] = 1

        # coordinate transforms
        markers["transforms"] = {}
        markers["transforms"]["lab_to_opensim"] = lab.transform_mat_lab_to_opensim[[1, -1, 2, -2, 3, -3].index(xdir)]

        # offset.offset_marker marker
        markers["offset_marker"] = lab.offset_marker
        
        # marker time and frames
        markers["time"] = c3dkey.markers["TIME"]
        markers["time0"] = c3dkey.markers["TIME"] - c3dkey.markers["TIME"][0]    
        markers["frames"] = c3dkey.markers["FRAME"]
        markers["frames0"] = c3dkey.markers["FRAME"] - c3dkey.markers["FRAME"][0]    

        # marker data, convert to m
        markers["data"] = {}
        for mkr in c3dkey.markers["DATA"]["POS"].keys():
            markers["data"][mkr] = c3dkey.markers["DATA"]["POS"][mkr] * markers["scale"]
        
        self.markers = markers
        
        return None
     
    def __set_force_plates(self, lab, c3dkey, xdir, staticprefix):
        
        # initialise dict
        force_plates = {}
        
        # determine used force plate depending on if static or dynamic trial
        if staticprefix.casefold() in c3dkey.trial_name.casefold():
            fpused = lab.fp_used_static
        else:
            fpused = lab.fp_used
                        
        # get force plate info for only used force plates
        force_plates["fp_used"] = [] 
        force_plates["fp_used_str"] = []        
        for f in fpused:
            
            # dict field
            dict_name = lab.fp_dict_name_prefix + str(f)
            force_plates["fp_used"].append(f)
            force_plates["fp_used_str"].append(dict_name)          
            
            # coordinate transforms
            force_plates[dict_name] = {}  
            force_plates[dict_name]["transforms"] = {}
            force_plates[dict_name]["transforms"]["lab_to_opensim"] = lab.transform_mat_lab_to_opensim[[1, -1, 2, -2, 3, -3].index(xdir)]
            force_plates[dict_name]["transforms"]["fp_to_lab"] = lab.transform_mat_fp_to_lab 
        
            # offsets
            offset_scale = self.markers["scale"]
            force_plates[dict_name]["offsets"] = {}
            force_plates[dict_name]["offsets"]["fp_centre_to_fp_origin_fp"] = -1*c3dkey.meta["FORCE_PLATFORM"]["ORIGIN"][f-1] * offset_scale
            force_plates[dict_name]["offsets"]["lab_to_fp_centre_lab"] = find_fp_centre_from_lab_origin(c3dkey.meta["FORCE_PLATFORM"]["CORNERS"][f-1]) * offset_scale
            force_plates[dict_name]["offsets"]["lab_to_fp_origin_lab"] = change_coordinates(force_plates[dict_name]["offsets"]["fp_centre_to_fp_origin_fp"], force_plates[dict_name]["transforms"]["fp_to_lab"], force_plates[dict_name]["offsets"]["lab_to_fp_centre_lab"])
        
        self.force_plates = force_plates
        
        return None
          
    def __set_forces(self, lab, c3dkey, staticprefix):
        
        # initialise dict
        forces = {}        
              
        # force plate time and frames
        forces["time"] = c3dkey.forces["TIME"]
        forces["time0"] = c3dkey.forces["TIME"] - c3dkey.forces["TIME"][0]    
        forces["frames"] = c3dkey.forces["FRAME"]
        forces["frames0"] = c3dkey.forces["FRAME"] - c3dkey.forces["FRAME"][0]    
        forces["rate"] = c3dkey.forces["RATE"]

        # determine used force plate depending on if static or dynamic trial
        if staticprefix.casefold() in c3dkey.trial_name.casefold():
            fpused = lab.fp_used_static
        else:
            fpused = lab.fp_used
    
        # get forces for only used force plates
        ns = len(forces["time"])
        for f in fpused:
            
            # dict field
            dict_name = lab.fp_dict_name_prefix + str(f)
            forces[dict_name] = {}
            
            # create force plate channel names (lazy formatting)
            channels = [s + str(f) for s in lab.fp_channel_prefixes]
            channels = [i + j for i, j in zip(channels, lab.fp_channel_suffixes)]
     
            # force plate labels: find first index, assume the others are in
            # the correct order
            indx0 = np.where(c3dkey.forces["LABELS"]==channels[0])[0][0]
            indx = list(range(indx0,indx0+lab.n_fp_channels))
            forces[dict_name]["labels"] = np.ndarray.tolist(np.array(c3dkey.forces["LABELS"])[indx])
     
            # force rate and units
            forces[dict_name]["rate"] = c3dkey.forces["RATE"]    
            old_units = np.ndarray.tolist(np.array(c3dkey.forces["UNITS"])[indx])
            
            # force scale factor, convert to N and Nm
            forces[dict_name]["scale"] = []
            forces[dict_name]["units"] = []
            for x in range(lab.n_fp_channels):
                if old_units[x] == "Nmm":
                    forces[dict_name]["scale"].append(0.001)
                    forces[dict_name]["units"].append("Nm")
                elif old_units[x] == "Nm":
                    forces[dict_name]["scale"].append(1)
                    forces[dict_name]["units"].append("Nm")
                elif old_units[x] == "N":
                    forces[dict_name]["scale"].append(1)
                    forces[dict_name]["units"].append("N")
                else:
                    forces[dict_name]["scale"].append(1)
                    forces[dict_name]["units"].append("N")
        
            # force plate data, convert to N and Nm
            forces[dict_name]["data"] = {}
            F = np.zeros([ns,3])
            M = np.zeros([ns,3])
            for c in range(3):
                F[:,c] = c3dkey.forces["DATA"][channels[c]] * forces[dict_name]["scale"][c]
                M[:,c] = c3dkey.forces["DATA"][channels[c+3]] * forces[dict_name]["scale"][c+3]                
            forces[dict_name]["data"]["F"] = F
            forces[dict_name]["data"]["M"] = M
                
            # calculate centre-of-pressure
            vc2o = self.force_plates[dict_name]["offsets"]["fp_centre_to_fp_origin_fp"]                        
            cop = calculate_centre_of_pressure_fp(ns, vc2o, F, M)
            forces[dict_name]["data"]["cop"] = cop
                
            # calculate free moments
            T = calculate_vertical_free_moment(ns, vc2o, F, M, cop)
            forces[dict_name]["data"]["T"] = T
            
            # convert to lab (Vicon) coordinates
            originvec0 = [0, 0, 0]
            rotmat = self.force_plates[dict_name]["transforms"]["fp_to_lab"] 
            originvec = self.force_plates[dict_name]["offsets"]["lab_to_fp_origin_lab"]                    
            F_lab = np.zeros([ns,3])
            M_lab = np.zeros([ns,3])
            cop_lab = np.zeros([ns,3])
            T_lab = np.zeros([ns,3])
            for n in range(ns):
                F_lab[n,:] = change_coordinates(F[n,:], rotmat, originvec0)
                M_lab[n,:] = change_coordinates(M[n,:], rotmat, originvec0)
                cop_lab[n,:] = change_coordinates(cop[n,:], rotmat, originvec)
                T_lab[n,:] = change_coordinates(T[n,:], rotmat, originvec0)
            forces[dict_name]["data"]["F_lab"] = F_lab
            forces[dict_name]["data"]["M_lab"] = M_lab
            forces[dict_name]["data"]["cop_lab"] = cop_lab
            forces[dict_name]["data"]["T_lab"] = T_lab
            
        self.forces = forces
        
        return None




'''
OpenSimKey:
    OpenSim data storage class for all data, model names, and other data for 
    processing through OpenSim.
'''
class OpenSimKey():
    def __init__(self, trialkey, ref_model, c3dpath, fp_filter_butter_order, fp_filter_cutoff, fp_filter_threshold, fp_smooth_cop_offset, fp_smooth_window, marker_filter_butter_order, marker_filter_cutoff):
        self.subject = trialkey.subject_name
        self.trial = trialkey.trial_name
        self.mass = trialkey.mass
        self.age = 0.0
        self.lab = trialkey.lab_name
        self.model = trialkey.subject_name + ".osim"
        self.task = trialkey.task
        self.condition = trialkey.condition
        self.outpath = c3dpath
        self.__set_events(trialkey)
        self.__set_markers(trialkey, marker_filter_butter_order, marker_filter_cutoff) 
        self.__set_forces(trialkey, fp_filter_butter_order, fp_filter_cutoff, fp_filter_threshold, fp_smooth_cop_offset, fp_smooth_window)      
        return None
    
    def __set_events(self, trialkey):
        
        # initialise dict
        events = {}
        
        # add events for trial
        events["time"] = trialkey.events["window_time0"]
        events["labels"] = trialkey.events["window_labels"]
        events["leg_task"] = trialkey.events["leg_task"]
        events["opensim_last_event_idx"] = trialkey.events["opensim_last_event_idx"]
    
        self.events = events
    
        return None
    
    def __set_markers(self, trialkey, marker_filter_butter_order, marker_filter_cutoff):
        
        # initialise dict
        markers = {}
        
        # time and frame vectors
        markers["time"] = trialkey.markers["time0"]
        markers["frames"] = np.arange(1,len(trialkey.markers["time0"]) + 1)             
        markers["rate"] = trialkey.markers["rate"]
        markers["units"] = trialkey.markers["units"]
        
        # get markers
        data = {}
        for mkr in trialkey.markers["data"].keys():
            
            # current marker data
            data_lab = trialkey.markers["data"][mkr]
            
            # convert data from lab coordinates to OpenSim coordinates
            ns = len(trialkey.markers["time0"])
            originvec0 = [0, 0, 0]
            rotmat = trialkey.markers["transforms"]["lab_to_opensim"]
            data[mkr] = np.zeros([ns,3])
            for n in range(ns): data[mkr][n,:] = change_coordinates(data_lab[n,:], rotmat, originvec0)
                  
        # marker offset
        markers["offset"] = [0., 0., 0.]
        if trialkey.markers["offset_marker"]: markers["offset"] = [data[trialkey.markers["offset_marker"]][0, 0], 0., data[trialkey.markers["offset_marker"]][0, 2]]
                       
        # offset X and Z trajectories (OpenSim coordinates system) of markers
        data_offset = data
        markers0 = {}
        for mkr in data.keys():
            for n in range(ns): data_offset[mkr][n,:] = data[mkr][n,:] - markers["offset"]
            markers0[mkr] = data_offset[mkr]
        
        # filter marker data
        for mkr in data.keys():
            markers[mkr] = filter_timeseries(markers0[mkr], markers["rate"], marker_filter_butter_order, marker_filter_cutoff)
                
        self.markers = markers
        
        return None
                   
    def __set_forces(self, trialkey, fp_filter_butter_order, fp_filter_cutoff, fp_filter_threshold, fp_smooth_cop_offset, fp_smooth_window):
 
        # initiliase temporary output arrays
        data = {}
        leg = ["right","left"]
        ns = len(trialkey.forces["time0"])
        for h, g in enumerate(leg):
            data[g] = {}
            data[g]["F"] = np.zeros([ns,3])
            data[g]["cop"] = np.zeros([ns,3])
            data[g]["T"] = np.zeros([ns,3])
    
        # initialise dict
        forces = {}        
      
        # time and frame vectors
        forces["time"] = trialkey.forces["time0"]
        forces["frames"] = np.arange(1,len(trialkey.forces["time0"]) + 1)
        forces["rate"] = trialkey.forces["rate"]

        # get forces for only used force plates (dynamic trials only)
        if trialkey.task != "static":
            forces["data"] = {}
            for i, fp in enumerate(trialkey.force_plates["fp_used"]):
                
                # dict field
                dict_name = trialkey.force_plates["fp_used_str"][i]
    
                # force data in lab (Vicon) coordinates
                F_lab = trialkey.forces[dict_name]["data"]["F_lab"]
                cop_lab = trialkey.forces[dict_name]["data"]["cop_lab"]
                T_lab = trialkey.forces[dict_name]["data"]["T_lab"]
    
                # convert data from lab coordinates to OpenSim coordinates
                originvec0 = [0, 0, 0]
                rotmat = trialkey.force_plates[dict_name]["transforms"]["lab_to_opensim"]
                F = np.zeros([ns,3])
                cop = np.zeros([ns,3])
                T = np.zeros([ns,3])
                for n in range(ns):
                    F[n,:] = change_coordinates(F_lab[n,:], rotmat, originvec0)
                    cop[n,:] = change_coordinates(cop_lab[n,:], rotmat, originvec0)
                    T[n,:] = change_coordinates(T_lab[n,:], rotmat, originvec0)            
                                                
                # offset the CoP using offset marker
                offset = [0., 0., 0.]
                if trialkey.markers["offset_marker"]: offset = self.markers["offset"]
                for n in range(ns): cop[n,:] = cop[n,:] - offset

                # # filter and floor the force plate data then smooth the
                # # foot-strike and foot-off edges
                # F1, T1, cop1 = filter_and_floor_fp(F, T, cop, 1, forces["rate"], filter_butter_order, filter_cutoff, filter_threshold)                
                # F2, T2, cop2 = smooth_transitions(F1, T1, cop1, 1, filter_threshold, smooth_cop_offset, smooth_window)

                # smooth the foot-strike and foot-off edges, then filter and 
                # floor the force plate data
                F1, T1, cop1 = smooth_transitions(F, T, cop, 1, fp_filter_threshold, fp_smooth_cop_offset, fp_smooth_window)
                F2, T2, cop2 = filter_and_floor_fp(F1, T1, cop1, 1, forces["rate"], fp_filter_butter_order, fp_filter_cutoff, fp_filter_threshold)
                    
                # for each force plate, add force plate data for any active
                # intervals to the output array for the relevant foot
                for h, g in enumerate(leg):
                    for n, m in enumerate(trialkey.events["fp_sequence"][:,h]):
                        if m == fp:                                     
                            idx0 = np.where(forces["time"] >= trialkey.events["window_intervals0"][n,0])[0][0]
                            idx1 = np.where(forces["time"] <= trialkey.events["window_intervals0"][n,1])[0][-1]
                            idx0s = idx0 - fp_smooth_window
                            idx1s = idx1 + fp_smooth_window - 1
                            data[g]["F"][idx0s:idx1s + 1,:] = F2[idx0s:idx1s + 1,:]
                            data[g]["cop"][idx0s:idx1s + 1,:] = cop2[idx0s:idx1s + 1,:]
                            data[g]["T"][idx0s:idx1s + 1,:] = T2[idx0s:idx1s + 1,:]
                
        # store force data
        forces["data"] = data
        
        self.forces = forces
        
        return None




'''
-----------------------------------
----- FUNCTIONS: C3D EXTRACT ------
-----------------------------------
'''



'''
c3d_batch_process(user, meta, lab, xdir, usermass):
    Batch processing for C3D data extract, and OpenSim input file write,
    obtains mass from used static trial in each group if mass = -1.
'''
def c3d_batch_process(user, meta, lab, xdir, usermass):

    # extract C3D data for OpenSim
    print("\n")
    for subj in meta:
        
        print("\n")
        print("%s" % "*" * 30)
        print("SUBJECT: %s" % subj)
        print("%s" % "*" * 30)
        print("\n")
        
        for group in meta[subj]["trials"]:
            
            print("Group: %s" % group)
            print("%s" % "=" * 30)
            
            
            # ###################################
            # PROCESS STATIC TRIAL
            
            mass = 0.0
            osimkey = {}
            for trial in meta[subj]["trials"][group]:                

                #****** FOR TESTING ONLY ******
                trialre = re.compile("FAILTCRT01_STATIC\d+")
                if not trialre.match(trial):
                    print("%s ---> SKIP" % trial)
                    continue
                #******************************
                
                # ignore dynamic trials
                isstatic = meta[subj]["trials"][group][trial]["isstatic"]
                usedstatic = meta[subj]["trials"][group][trial]["usedstatic"]
                if not(isstatic): continue
            
                print("\nStatic trial: %s" % trial)
                print("%s" % "-" * 30)
            
                # process C3D file and generate OsimKey for trial
                c3dfile = meta[subj]["trials"][group][trial]["c3dfile"]
                c3dpath = meta[subj]["trials"][group][trial]["outpath"]
                task = meta[subj]["trials"][group][trial]["task"]
                condition = meta[subj]["trials"][group][trial]["condition"]
                osimkey = c3d_extract(trial, c3dfile, c3dpath, lab, task, condition, xdir, user.refmodelfile, user.staticfpchannel, mass, user.fp_filter_butter_order, user.fp_filter_cutoff, user.fp_filter_threshold, user.fp_smooth_cop_fixed_offset, user.fp_smooth_window, user.marker_filter_butter_order, user.marker_filter_cutoff, user.staticprefix)                           
                
                # get the mass from the used static trial
                if usedstatic: mass = osimkey.mass
            
            #
            # ###################################            
            
            
            # override the mass from used static trial if user supplied
            if usermass != -1: mass = usermass
            
            
            # ###################################
            # PROCESS DYNAMIC TRIAL
            
            # process dynamic C3D files
            for trial in  meta[subj]["trials"][group]:                

                #****** FOR TESTING ONLY ******                
                trialre = re.compile("FAILTCRT01_SDP01")
                if not trialre.match(trial):
                    print("%s ---> SKIP" % trial)
                    continue
                #******************************
                
                # ignore static trials
                isstatic = meta[subj]["trials"][group][trial]["isstatic"]
                if isstatic: continue
            
                print("\nDynamic trial: %s" % trial)
                print("%s" % "-" * 30)
            
                # process C3D file and generate OsimKey for trial
                c3dfile = meta[subj]["trials"][group][trial]["c3dfile"]
                c3dpath = meta[subj]["trials"][group][trial]["outpath"]
                task = meta[subj]["trials"][group][trial]["task"]
                condition = meta[subj]["trials"][group][trial]["condition"]
                c3d_extract(trial, c3dfile, c3dpath, lab, task, condition, xdir, user.refmodelfile, user.staticfpchannel, mass, user.fp_filter_butter_order, user.fp_filter_cutoff, user.fp_filter_threshold, user.fp_smooth_cop_fixed_offset, user.fp_smooth_window, user.marker_filter_butter_order, user.marker_filter_cutoff, user.staticprefix)                           
                     
            #
            # ###################################                    

    return None



'''
c3d_extract(trial, c3dpath, c3dpath, lab, condition, xdir, ref_model, 
            static_fp_channel, mass, fp_filter_butter_order, fp_filter_cutoff, 
            fp_filter_threshold, fp_smooth_cop_offset, fp_smooth_window):
    Extract the motion data from the C3D file to arrays, and returns a dict
    containing all the relevant file metadata, force data and marker data.
'''
def c3d_extract(trial, c3dfile, c3dpath, lab, task, condition, xdir, ref_model, static_fp_channel, mass, fp_filter_butter_order, fp_filter_cutoff, fp_filter_threshold, fp_smooth_cop_offset, fp_smooth_window, marker_filter_butter_order, marker_filter_cutoff, staticprefix):
    
    # load C3D file
    itf = c3d.c3dserver()
    c3d.open_c3d(itf, c3dpath + "/" + c3dfile)       
    
    # get all file metadata, and all force plate and video C3D data
    fmeta = c3d.get_dict_groups(itf)
    fforces = c3d.get_dict_forces(itf, frame=True, time=True)
    fmarkers = c3d.get_dict_markers(itf, frame=True, time=True)
    
    # subject and trial name
    sname= fmeta["SUBJECTS"]["NAMES"][0]
    if not sname: sname = "NoName"
    tname = trial
    
    # C3D key with all data from C3D file
    c3dkey = C3DKey(sname, tname, fmeta, fforces, fmarkers)

    # trial data only from C3D key
    trialkey = TrialKey(lab, task, condition, c3dkey, xdir, static_fp_channel, mass, staticprefix)
    
    # opensim input data
    osimkey = OpenSimKey(trialkey, ref_model, c3dpath, fp_filter_butter_order, fp_filter_cutoff, fp_filter_threshold, fp_smooth_cop_offset, fp_smooth_window, marker_filter_butter_order, marker_filter_cutoff)
    
    # save key files
    with open(os.path.join(c3dpath, trial + "_c3dkey.pkl"),"wb") as f: pk.dump(c3dkey, f)
    with open(os.path.join(c3dpath, trial + "_trialkey.pkl"),"wb") as g: pk.dump(trialkey, g)
    with open(os.path.join(c3dpath, trial + "_osimkey.pkl"),"wb") as h: pk.dump(osimkey, h)

    # write OpenSim input data files
    osp.write_ground_forces_mot_file(osimkey)
    osp.write_marker_trajctory_trc_file(osimkey) 
    
    return osimkey
    




'''
-----------------------------------
--------- FUNCTIONS: MISC ---------
-----------------------------------
'''


'''
change_coordinates(oldvec, rotmat, originvec):
    Perform coordinate transformation and change of origin
'''
def change_coordinates(oldvec, rotmat, originvec):
    return originvec + np.matmul(rotmat, oldvec)



'''
find_fp_centre_from_lab_origin(corners):
    Find force plate centre in lab (Vicon) coordinates
'''
def find_fp_centre_from_lab_origin(corners):
    
    # lab origin ---> corner 1
    lo_c1 = corners[0]
    
    # lab origin ---> corner 3
    lo_c3 = corners[2]
    
    # corner 1 ---> corner 3
    c1_c3 = lo_c3 - lo_c1
    
    # corner 1 ---> centre
    c1_ct = c1_c3 / 2
    
    # lab origin ---> centre
    return lo_c1 + c1_ct



'''
calculate_centre_of_pressure_fp(n, vc2o, F, M):
    Calculate centre of pressure in force plate coordinates using equations
    from Tim Dorn's GaitExtract toolbox. Note: the GaitExtract documentation
    contains an error in the equations on page 14. The correct equations are
    given in the Matlab function getKinetics.m
        
        Given:
        
            (cop_x, cop_y, cop_z): vector from origin to centre of pressure
            (a, b, c): vector from geometric centre to origin
            (F_x, F_y, F_z): force plate force
            (M_x, M_y, M_z): force plate moment
                
        Equations:
            
            cop_x = -1 * ( (M_y + c * F_x) / F_z ) + a
            cop_y = ( (M_x - c * F_y) / F_z ) + b
            cop_z = 0 (default)        
'''
def calculate_centre_of_pressure_fp(ns, vc2o, F, M):
    cop = np.zeros([ns,3])
    for n in range(ns):
        cop[n,0] = -1 * ((M[n,1] + vc2o[2] * F[n,0]) / F[n,2]) + vc2o[0]
        cop[n,1] = ((M[n,0] - vc2o[2] * F[n,1]) / F[n,2]) + vc2o[1]
        if np.isnan(cop[n,0]) or np.isnan(cop[n,1]):
            cop[n,0] = 0
            cop[n,1] = 0
    return cop



'''
calculate_vertical_free_moment(ns, vc2o, F, M, cop):
    Calculate vertical free moment in force plate coordinates using equations
    from Tim Dorn's GaitExtract toolbox.
        
        Given:
        
            (T_x, T_y, T_z): vector of free moment
            (cop_x, cop_y, cop_z): vector from origin to centre of pressure
            (a, b, c): vector from geometric centre to origin
            (F_x, F_y, F_z): force plate force
            (M_x, M_y, M_z): force plate moment
                
        Equations:
            
            T_x = 0 (default)
            T_y = 0 (default)
            T_z = M_z - ((cop_x - a) * F_y) + ((cop_y - b) * F_x)       
'''
def calculate_vertical_free_moment(ns, vc2o, F, M, cop):
    T = np.zeros([ns,3])
    for n in range(ns):
        T[n,2] = M[n,2] - ((cop[n,0] - vc2o[0]) * F[n,1]) + ((cop[n,1] - vc2o[1]) * F[n,0])
    return T



'''
calculate_subject_mass(c3dkey, static_fp_channel):
    Calculate subject mass from a static trial using vertical force plate data 
    from the middle part of the trial.
'''
def calculate_subject_mass(c3dkey, static_fp_channel):
    
    # raw force plate data
    data = c3dkey.forces["DATA"][static_fp_channel]
    frames = list(range(0,len(data)))
    
    # trim data to a window from 25% to 75% of the datastream length
    # (arbitrary window in middle to avoid movement at start and end of trial)
    idx0 = round(np.percentile(frames,25))
    idx1 = round(np.percentile(frames,75))
    data = data[idx0:idx1 + 1]
    
    # calculate average mass
    mass = abs(stats.mean(data)) / 9.81
    
    return mass, idx0, idx1
    

    
'''
filter_and_floor_fp(F, T, CoP, vert_idx, sample_rate, butter_order, cutoff, 
                    threshold):
    Filter the force plate data, and floor data below threshold.
'''    
def filter_and_floor_fp(F, T, cop, vert_col_idx, sample_rate, butter_order, cutoff, threshold):
    
    # apply filter
    F1 = filter_timeseries(F, sample_rate, butter_order, cutoff)
    T1 = filter_timeseries(T, sample_rate, butter_order, cutoff)
    
    # apply threshold
    Fy = F1[:, vert_col_idx]
    idxs = np.where(Fy < threshold)
    F1[idxs, :] = 0.0
    T1[idxs, :] = 0.0
    cop[idxs, :] = 0.0
    
    return F1, T1, cop
    


'''
filter_timeseries(data_raw, sample_rate, butter_order, cutoff):
    Filter timeseries data. Raw data can be a list, or an array with rows
    representing time steps and columns as variables. Set cutoff < 0 if
    filtering not required.
''' 
def filter_timeseries(data_raw, sample_rate, butter_order, cutoff):
    
    # return if cutoff < 0 (i.e. filtering not required)
    if cutoff < 0: return data_raw
    
    # filter design
    Wn = sample_rate / 2
    normalised_cutoff = cutoff / Wn
    b, a = signal.butter(butter_order, normalised_cutoff, "lowpass")
    
    # apply filter
    data_filtered = signal.filtfilt(b, a, data_raw, axis = 0)

    return data_filtered
    

  
'''
smooth_transitions(F, T, CoP, vert_col_idx, threshold, cop_fixed_offset, 
                   window):
    Smooth the transitions at foot-strike and foot-off. Affix the transitioning
    forces at the heel (i.e. first valid GRF) to avoid drift of the CoP from
    the force plate centre to the foot, which can create large moments about
    the lower limb joints.
'''
def smooth_transitions(F, T, cop, vert_col_idx, threshold, cop_fixed_offset, window):
    
    # ignore if threshold is a negative number
    if threshold < 0: return F, T, cop
    
    # find threshold indices
    Fy = F[:, vert_col_idx].copy()
    idxs1 = np.where(Fy >= threshold)
    Fy[idxs1] = 1.0
    Fy = np.insert(Fy, 0, 0)
    Fy_diff = np.diff(Fy)
    idxup = np.where(Fy_diff == 1)
    idxdn = np.where(Fy_diff == -1)
    
    # rectify drift on foot strike
    for xi in idxup[0]:
        
        # x window
        x = int(xi)
        x0window = np.array([0, window - 1])
        x1window = range(0, window)
        
        # F window
        f0window = np.row_stack([[0.0, 0.0, 0.0], F[x, :]])
        fspline = interp.interp1d(x0window, f0window, kind = "linear", axis = 0) 
        f1window = fspline(x1window)
        F[x - window:x, :] = f1window
        
        # T window
        t0window = np.row_stack([[0.0, 0.0, 0.0], T[x, :]])
        tspline = interp.interp1d(x0window, t0window, kind = "linear", axis = 0) 
        t1window = tspline(x1window)
        T[x - window:x, :] = t1window      
        
        # CoP window, affix to a point reasonably ahead
        cop[x - window:x, :] = cop[x, :]

    # rectify drift on foot off
    for xi in idxdn[0]:
        
        # x window
        x = int(xi)
        x0window = np.array([0, window - 1])
        x1window = range(0, window)
        
        # F window
        f0window = np.row_stack([F[x - 1, :], [0.0, 0.0, 0.0]])
        fspline = interp.interp1d(x0window, f0window, kind = "linear", axis = 0) 
        f1window = fspline(x1window)
        F[x - 1:x + window - 1, :] = f1window
        
        # T window
        t0window = np.row_stack([T[x - 1, :], [0.0, 0.0, 0.0]])
        tspline = interp.interp1d(x0window, t0window, kind = "linear", axis = 0) 
        t1window = tspline(x1window)
        T[x - 1:x + window - 1, :] = t1window      
        
        # CoP window
        cop[x - 1:x + window - 1, :] = cop[x - 1, :]

        
    return F, T, cop
    
    