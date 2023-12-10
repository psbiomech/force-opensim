# -*- coding: utf-8 -*-
"""
Stability analyses: 
    Margin of stability
    Probability of Instability
    Whole Body Angular Momentum

@author: Prasanna Sritharan
"""

import opensim as osim
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import numpy as np
import shapely
import scipy.constants as consts
import scipy.interpolate as interp
import scipy.integrate as integ
from scipy.spatial.transform import Rotation
#import scipy.io as scio
import pickle as pk
import pandas as pd
import os


# Use non-interactive backend to avoid memory accumulation
plt.switch_backend("Agg")




'''
-----------------------------------
------------- CLASSES -------------
-----------------------------------
'''

# NO CLASSES



'''
-----------------------------------
--------- BATCH FUNCTIONS ---------
-----------------------------------
'''


'''
batch_process_stability(user, meta, artifperturb, treadmill_speed, restart):
    Batch process margin of stability, and whole body angular momentum.
'''
def batch_process_stability(user, meta, artifperturb = False, treadmill_speed = 0.0, restart = -1):

    
    failedfiles = []
    startflag = 0
    for subj in meta:
                
        # Skip to restart participant, process until last restart participant.
        # Python uses lazy evaluation so combined expressions are efficient.
        if restart != -1:
            if startflag == 1:
                if (type(restart) == tuple) and (subj == restart[1]):
                    startflag = 0            
            elif startflag == 0:
                if (type(restart) == str) and (subj == restart):
                    startflag = 1
                elif (type(restart) == tuple) and (subj == restart[0]):
                    if restart[0] != restart[1]:
                        startflag = 1
                else:
                    continue    


    # Subjects
    for subj in meta:
        
        print("-----------------------------------")
        print("SUBJECT: %s" % subj)
        #print("-----------------------------------")  
        
        # Groups
        for group in meta[subj]["trials"]:
        
            # Trials
            for trial in meta[subj]["trials"][group]:
                
                #****** FOR TESTING ONLY ******
                # if not (trial == "FAILTCRT13_SLDJ05"): continue
                #******************************
                
                # ignore static trials
                isstatic = meta[subj]["trials"][group][trial]["isstatic"]
                if isstatic: continue
                
                print("-----------------------------------") 
                print("TRIAL: %s" % trial)
                print("-----------------------------------") 
                
                try:
               
                    # Load the OpenSim results (OsimResultsKey)
                    trialpath = meta[subj]["trials"][group][trial]["outpath"]
                    with open(os.path.join(trialpath, trial + "_opensim_results.pkl"), "rb") as fid0:
                        datakey = pk.load(fid0)   
                                
                    # Load OpenSim input data (OsimKey)
                    with open(os.path.join(trialpath, trial + "_osimkey.pkl"), "rb") as fid1:
                        osimkey = pk.load(fid1)            
                    
                    
                    
                    # ******************************
                    # MARGIN OF STABILITY
                    
                    # Calculate simple margin of stability
                    print("---> Calculating margin of stability...")
                    stable = margin_of_stability(user, datakey, osimkey, treadmill_speed)
                    
                    # Visualise
                    # print("---> Generating visualisations of margin of stability...")
                    # plot_margin_of_stability(datakey, stable)
                    # visualise_stability_timehistory(datakey, stable, ylim = [0, 1.5])
                    # animate_stability_timehistory(datakey, stable)
                    
                    
                    
                    # ******************************
                    # WHOLE BODY ANGULAR MOMENTUM
                    
                    # Calculate whole body angular momentum
                    print("---> Calculating whole body angular momentum...")
                    wbam = whole_body_angular_momentum(user, datakey)
                    
                    # Visualise
                    # print("---> Generating visualisations of whole body angular momentum...")
                    # plot_whole_body_angular_momentum(datakey, wbam)
                    
                    
                except:
                    print("*** FAILED ***")
                    failedfiles.append(trial)     
                    #raise
                                    
    return failedfiles    




'''
-----------------------------------
------- STABILITY ANALYSES --------
-----------------------------------
'''


'''
whole_body_angular_momentum(user, datakey):
    Calculate whole body angular momentum as the sum of segmental angular 
    momentum with respect to the whole body centre-of-mass, from mechanical
    dynamics: L = sum(r x mv + Iw)
'''
def whole_body_angular_momentum(user, datakey):

    # Get the model
    osimfile = os.path.join(datakey.outpath, datakey.model)
    model = osim.Model(osimfile)
    
    # Get the bodies and inertia
    # Assumption: moment of inertia are provided as principal moments of inertia
    # therefore no inertia matrix products.
    bodyset = model.getBodySet()
    nbods = bodyset.getSize()
    bodies = {}
    for b in range(nbods):
        bname = bodyset.get(b).getName()
        bodies[bname] = {}
        bodies[bname]["m"] = bodyset.get(b).getMass()
        bodies[bname]["I"] = np.array([bodyset.get(b).getInertia().getMoments().get(i) for i in range(3)])

    # Get the absolute segmental position and orientation
    positions = {}
    for b in range(nbods):
        bname = bodyset.get(b).getName()
        positions[bname] = {}
        ridxs = [h for h, header in enumerate(datakey.results["raw"]["bk"]["headers"][0]) if bname in header]
        positions[bname]["r"] = datakey.results["raw"]["bk"]["data"][:, ridxs[0:3], 0]
        positions[bname]["theta"] = np.radians(datakey.results["raw"]["bk"]["data"][:, ridxs[3:6], 0])            
        
    # Get the segmental absolute linear and angular velocities
    # Convert angular velocities to rad/s
    velocities = {}
    for b in range(nbods):
        bname = bodyset.get(b).getName()
        velocities[bname] = {}
        vidxs = [h for h, header in enumerate(datakey.results["raw"]["bk"]["headers"][1]) if bname in header]
        velocities[bname]["v"] = datakey.results["raw"]["bk"]["data"][:, vidxs[0:3], 1]
        velocities[bname]["w"] = np.radians(datakey.results["raw"]["bk"]["data"][:, vidxs[3:6], 1])
    
    # Get the absolute centre-of-mass position and velocity
    com = {}
    cidxs0 = [h for h, header in enumerate(datakey.results["raw"]["bk"]["headers"][0]) if "center_of_mass" in header]
    cidxs1 = [h for h, header in enumerate(datakey.results["raw"]["bk"]["headers"][1]) if "center_of_mass" in header]
    com["r"] = datakey.results["raw"]["bk"]["data"][:, cidxs0, 0]    
    com["v"] = datakey.results["raw"]["bk"]["data"][:, cidxs1, 1]
    
    # Get the time vector from BodyKinematics
    timevec = datakey.results["raw"]["bk"]["data"][:, 0, 0]
        
    # Calculate segmental angular momentum in inertial coordinates
    nsamps = np.shape(com["r"])[0]
    L_seg = {}
    for b in range(nbods):

        # At each time step
        bname = bodyset.get(b).getName()
        L_seg[bname] = np.zeros([nsamps, 3])
        for n in range(nsamps):

            # Angular velocity relative to segment centre-of-mass in inertial
            # coordinates
            w = np.reshape(velocities[bname]["w"][n, :], [3, 1])
            
            # Linear kinematics relative to whole body centre-of-mass in 
            # inertial coordinates
            r = positions[bname]["r"][n, :] - com["r"][n, :]
            v = velocities[bname]["v"][n, :] - com["v"][n, :]               

            # Inertia relative to segment centre of mass in body coordinates
            m = bodies[bname]["m"]
            I = np.diag(bodies[bname]["I"])
            
            # Body orientation
            theta3 = positions[bname]["theta"][n, :]
            
            # Build the change of basis rotation matrix to convert inertia to
            # inertial coordinates
            rotmat = Rotation.from_rotvec(theta3, degrees=True).as_matrix()
            
            # Change basis of inertia tensor to inertial coordinates given
            # angular momentum in inertial coordinates:
            #   I_iner = R * I * R_inv = R * I * R_trans
            # Note: this is change of basis only, inertia is still about the
            # segmental centre-of-mass.
            I_iner = np.matmul(rotmat, np.matmul(I, np.linalg.inv(rotmat)))       
            
            # Instantaneous segmental angular momentum relative to centre-of mass:
            #   L_seg_t = r x mv + Iw
            L_seg_t = np.cross(r, m * v) + np.matmul(I_iner, w).T
            
            # Store
            L_seg[bname][n, :] = L_seg_t


    # Calculate whole body angular momentum
    L = np.zeros([nsamps, 3])
    for b in range(nbods): 
        L = L + L_seg[bodyset.get(b).getName()]
    
    # Range of whole body angular momentum
    L_max = np.amax(L, axis = 0)
    L_min = np.amin(L, axis = 0)
    L_range = np.absolute(L_max - L_min)

    # Integrated whole body angular momentum
    L_int = integ.simpson(L, timevec, axis = 0)
    
    # Average 3D centre of mass velocity (for normalisation)
    CoM_v_norm = np.linalg.norm(com["v"], axis = 1)
    CoM_v_mean = np.mean(CoM_v_norm)
    
    
    # Store in dict
    wbam = {}
    wbam["bodies"] = bodies
    wbam["positions"] = positions
    wbam["velocities"] = velocities
    wbam["com"] = com
    wbam["L_seg"] = L_seg
    wbam["L"] = L
    wbam["L_range"] = L_range
    wbam["L_int"] = L_int
    wbam["task"] = datakey.task
    wbam["timevec"] = timevec
    wbam["CoM_v_mean"] = CoM_v_mean

    # Pickle it
    trialname = datakey.trial
    with open(os.path.join(datakey.outpath, trialname + "_wbam.pkl"), "wb") as f:
        pk.dump(wbam, f)

    return wbam



'''
perturbed_margin_of_stability(user, datakey, pertubation, stepsize, treadmill_speed):
    
    *** TBD: DO NOT USE. NEEDS TO BE REWORKED FOR FORCE/TRAIL DATA ***
    
    Calculate the 2D margin of stability by artificially perturbing the centre
    of mass (CoM) velocity (default = 1 m/s) at regular angular increments 
    (default = 16 steps). For treadmill trials, add the treadmill speed to the 
    anteroposterior (X) centre of mass velocity (default = 0.0 m/s).
'''
def perturbed_margin_of_stability(user, datakey, perturbation = 1.0, asteps = 16, treadmill_speed = 0.0):
    
    # Number of steps should be a positive number
    if asteps <= 0: asteps = 16
    
    # Add an extra step to complete the loop
    asteps = asteps + 1
          
    # Centre of mass position (OpenSim global)
    comdata = datakey.data_osim_results[user.bkcode]["pos"]
    comxidx = comdata["headers"].index("center_of_mass_X")
    CoM_r = comdata["data"][:, comxidx:comxidx + 3]
    
    # Centre of mass velocity (OpenSim global)
    CoM_v = datakey.data_osim_results[user.bkcode]["vel"]["data"][:, comxidx:comxidx + 3].copy()   # Requires copy() or will only create new reference, not a new variable!
        
    # Adjust fore-aft velocity (for treadmill walking, default = 0.0 m/s)
    CoM_v[:, 0] = CoM_v[:, 0] + treadmill_speed
    
    # Add the 2D perturbation velocity vector counter-clockwise at regular
    # angular increments from the positive anteroposterior (X) axis
    perturb_angles = np.linspace(0, 2 * np.pi, asteps) #[0:-1]    # drop last element as it repeats the first
    CoM_v_perturb = np.zeros([asteps, np.size(CoM_v, axis = 0), 3])    
    for t in range(np.size(CoM_v, axis = 0)):
        for a, ang in enumerate(perturb_angles):
            CoM_v_perturb[a, t, 0] = np.cos(ang) * perturbation
            CoM_v_perturb[a, t, 2] = -1 * np.sin(ang) * perturbation
               
    # Instantaneous centre of mass position (m) and pendulum natural 
    # frequency (rad/s)
    CoM_l = np.linalg.norm(CoM_r, axis = 1)
    w0 = np.sqrt(consts.g / CoM_l)

    # Extrapolated centre of mass without artificial perturbation
    XCoM_r = CoM_r + np.divide(CoM_v, np.tile(np.reshape(w0, [-1, 1]), [1, 3]))
    
    # Extrapolated centre of mass with artifical perturbation
    XCoM_r_perturb = np.zeros([asteps, np.size(CoM_v, axis = 0), 3])
    for t in range(np.size(CoM_v, axis = 0)):
        for a in range(len(perturb_angles)):
            XCoM_r_perturb[a, t, :] = CoM_r[t, :] + ((CoM_v[t, :] + CoM_v_perturb[a, t, :]) / w0[t])
            
    # Create CoM and XCoM as 2D Points for artificial perturbation
    CoM = {}
    CoM["CoM_r"] = CoM_r[:, [0, 2]]
    CoM["CoM_v"] = CoM_v[:, [0, 2]]
    CoM["CoM_l"] = CoM_l
    CoM["w0"] = w0    
    CoM["CoM_v_perturb"] = CoM_v_perturb[:, :, [0, 2]]
    CoM["XCoM_r_perturb"] = XCoM_r_perturb[:, :, [0, 2]]
    CoM["XCoM_pt_perturb"] = []
    CoM["CoM_pt"] = []
    CoM["XCoM_pt"] = []
    for t in range(np.size(CoM_r, axis = 0)):
        CoM["CoM_pt"].append(shapely.Point(CoM_r[t, [0, 2]]))
        CoM["XCoM_pt"].append(shapely.Point(XCoM_r[t, [0, 2]]))
        pts_perturb = []
        for a in range(len(perturb_angles)):
            pts_perturb.append(shapely.Point(XCoM_r_perturb[a, t, [0, 2]]))
        CoM["XCoM_pt_perturb"].append(pts_perturb)
    
    
    # Base of support (OpenSim global) from foot markers and GRF
    # Assume marker data already converted to model coordinate system
    markers = datakey.data_c3dextract_osim["markers"]
    forces = datakey.data_c3dextract_osim["forces"]
    timevec = datakey.data_osim_results[user.bkcode]["pos"]["data"][:, datakey.data_osim_results[user.bkcode]["pos"]["headers"].index("time")]
    BoS = construct_base_of_support(markers, forces, timevec)  
    
    # Two-dimensional margin of stability (MoS, b) in transverse (x-z) plane
    # with artificial perturbation:
    # Minimum distance from XCoM to the edge of base of stability (BoS),
    # equivalent to shortest Cartesian distance from point to the edge of the 
    # convex hull formed by the points representing the BoS. Shapely
    # distance() calculates this automatically. Set b as negative if XCoM is 
    # outside the convex hull, positive if inside.
    MoS = {}
    MoS["b_abs_perturb"] = []
    MoS["isstable_perturb"] = []
    MoS["b_perturb"] = []
    for t in range(len(BoS["baseshape"])):
        b_abs_perturb = []
        isstable_perturb = []
        for a in range(len(perturb_angles)):
            b_abs_perturb.append(BoS["baseshape"][t].exterior.distance(CoM["XCoM_pt_perturb"][t][a]))
            isstable_perturb.append(BoS["baseshape"][t].contains(CoM["XCoM_pt_perturb"][t][a]))       
        MoS["b_abs_perturb"].append(b_abs_perturb)
        MoS["isstable_perturb"].append(isstable_perturb)   
        MoS["b_perturb"].append([MoS["b_abs_perturb"][t][bidx] if stab else -1 * MoS["b_abs_perturb"][t][bidx] for bidx, stab in enumerate(MoS["isstable_perturb"][t])])
        
    # Margins of stability along coordinates (b_x, b_z)
    MoS["b_x_perturb"] = []
    MoS["b_z_perturb"] = []
    for t in range(len(BoS["baseshape"])):
        b_x_a = []
        b_z_a = []
        for a in range(len(perturb_angles)):
        
            # Instantaneous parameters
            xcom = CoM["XCoM_pt_perturb"][t][a].coords[0]
            comv = CoM["CoM_v_perturb"][a, t, :]
            bounds = BoS["baseshape"][t].bounds   # tuple: (min x, min z, max x, max z)
            
            # Margin along coordinates depends on direction of velocity
            
            # X direction
            if (comv[0] >= 0):
                b_x = bounds[2] - xcom[0]
            else:
                b_x = xcom[0] - bounds[0]
            b_x_a.append(b_x)
               
            # Z direction
            if (comv[1] >= 0):
                b_z = bounds[3] - xcom[1]
            else:
                b_z = xcom[1] - bounds[1]        
            b_z_a.append(b_z)
            
        MoS["b_x_perturb"].append(b_x_a)
        MoS["b_z_perturb"].append(b_z_a)
            
    # Construct dict with all stability data
    perturbdict = {}
    perturbdict["MoS"] = MoS
    perturbdict["BoS"] = BoS
    perturbdict["CoM"] = CoM
    perturbdict["angles"] = perturb_angles
    perturbdict["perturbation"] = perturbation
    perturbdict["treadmill"] = treadmill_speed
        
    # Pickle it
    trialname = datakey.trial["name"]
    with open(os.path.join(datakey.trial["path"], trialname + "_stability_perturbed.pkl"), "wb") as f:
        pk.dump(perturbdict, f)     
    
    return perturbdict



'''
margin_of_stability(user, datakey, osimkey, treadmill_speed):
    Calculate 2D margin of stability (b) using the geometric method described
    by Hof et al. (2005). Use Shapely to do the vector geometry. For treadmill
    trials, add the treadmill speed to the anteroposterior centre of mass 
    velocity (default = 0.0 m/s).
'''
def margin_of_stability(user, datakey, osimkey, treadmill_speed = 0.0):
    
    # Centre of mass position (OpenSim global)   
    comdata = datakey.results["raw"][user.bkcode]["data"][:, :, 0]
    comxidx = datakey.results["raw"][user.bkcode]["headers"][0].index("center_of_mass_X")
    CoM_r = comdata[:, comxidx:comxidx + 3]
    
    # Centre of mass velocity (OpenSim global)
    CoM_v = datakey.results["raw"][user.bkcode]["data"][:, comxidx:comxidx + 3, 1].copy()  # Requires copy() or will only create new reference, not a new variable!
    
    # Adjust fore-aft velocity (for treadmill walking, default = 0.0 m/s)
    CoM_v[:, 0] = CoM_v[:, 0] + treadmill_speed
    
    # Instantaneous centre of mass position (m) and pendulum natural 
    # frequency (rad/s)
    CoM_l = np.linalg.norm(CoM_r, axis = 1)
    w0 = np.sqrt(consts.g / CoM_l)
    
    # Extrapolated centre of mass
    XCoM_r = CoM_r + np.divide(CoM_v, np.tile(np.reshape(w0, [-1, 1]), [1, 3]))
    
    # Create CoM and XCoM as 2D Points
    CoM = {}
    CoM["CoM_r"] = CoM_r[:, [0, 2]]
    CoM["CoM_v"] = CoM_v[:, [0, 2]]
    CoM["CoM_l"] = CoM_l
    CoM["w0"] = w0
    CoM["XCoM_r"] = XCoM_r[:, [0, 2]]
    CoM["CoM_pt"] = []
    CoM["XCoM_pt"] = []
    for t in range(np.size(CoM_r, axis = 0)):
        CoM["CoM_pt"].append(shapely.Point(CoM_r[t, [0, 2]]))
        CoM["XCoM_pt"].append(shapely.Point(XCoM_r[t, [0, 2]]))
    
    
    # Base of support (OpenSim global) from foot markers and GRF
    # Assume marker data already converted to model coordinate system
    markers = osimkey.markers
    forces = osimkey.forces    
    timevec = datakey.results["raw"][user.bkcode]["data"][:, datakey.results["raw"][user.bkcode]["headers"][0].index("time"), 0]
    BoS = construct_base_of_support(markers, forces, timevec)  
    
    # Two-dimensional margin of stability (MoS, b) in transverse (x-z) plane: 
    # Minimum distance from XCoM to the edge of base of stability (BoS),
    # equivalent to shortest Cartesian distance from point to the edge of the 
    # convex hull formed by the points representing the BoS. Shapely
    # distance() calculates this automatically. Set b as negative if XCoM is 
    # outside the convex hull, positive if inside.
    MoS = {}
    MoS["b_abs"] = []
    MoS["isstable"] = []
    for t in range(len(BoS["baseshape"])):
        
        # No base of support means empty Shapely GeometryColleciton created by
        # default. Report nominal values (0.0, None) for these cases.
        if type(BoS["baseshape"][t]) == shapely.geometry.collection.GeometryCollection:
            MoS["b_abs"].append(0.0)
            MoS["isstable"].append(None)
            
        # Valid base of support    
        else:
            MoS["b_abs"].append(BoS["baseshape"][t].exterior.distance(CoM["XCoM_pt"][t]))
            MoS["isstable"].append(BoS["baseshape"][t].contains(CoM["XCoM_pt"][t]))        
    
    # Append sign (pos/neg) to margin of stability: pos = stable, neg = unstable
    MoS["b"] = [MoS["b_abs"][bidx] if stab else -1 * MoS["b_abs"][bidx] for bidx, stab in enumerate(MoS["isstable"])]
    
    # Margins of stability along coordinates (b_x, b_z)
    MoS["b_x"] = []
    MoS["b_z"] = []
    for t in range(len(BoS["baseshape"])):
        
        # If no boundary of stability (flight), just return nominal values (0.0)
        if type(BoS["baseshape"][t]) == shapely.geometry.collection.GeometryCollection:
            MoS["b_x"].append(0.0)
            MoS["b_z"].append(0.0)
            continue
        
        
        # Instantaneous parameters
        xcom = CoM["XCoM_pt"][t].coords[0]
        comv = CoM["CoM_v"][t]
        bounds = BoS["baseshape"][t].bounds   # tuple: (min x, min z, max x, max z)
        
        # Margin along coordinates depends on direction of velocity
        
        # X direction
        if (comv[0] >= 0):
            b_x = bounds[2] - xcom[0]
        else:
            b_x = xcom[0] - bounds[0]
        MoS["b_x"].append(b_x)
           
        # Z direction
        if (comv[1] >= 0):
            b_z = bounds[3] - xcom[1]
        else:
            b_z = xcom[1] - bounds[1]        
        MoS["b_z"].append(b_z)        
            
    # Construct dict with all stability data
    stabledict = {}
    stabledict["MoS"] = MoS
    stabledict["BoS"] = BoS
    stabledict["CoM"] = CoM
    
    # Trial info
    stabledict["task"] = osimkey.task
    stabledict["timevec"] = timevec
        
    # Pickle it
    trialname = datakey.trial
    with open(os.path.join(datakey.outpath, trialname + "_stability.pkl"), "wb") as f:
        pk.dump(stabledict, f)     
    
    return stabledict



'''
construct_base_of_support(markers, forces, timevec):
    Construct the base of support Polygons using foot markers and the phase of
    support during the trial. The Polygon is constructed at each time step from
    the convex hull of the 2D point cloud defined by the relevant foot markers.
    Resample to the OpenSim results timesteps.
'''
def construct_base_of_support(markers, forces, timevec):
    
    # Markers defining base of support in counterclockwise order
    # Temporary: assume all foot markers always define base of support. This is
    # actually only true during foot-flat. 
    bosmarkers = [[], \
                  ["RHEEL", "RLMAL", "RMFL", "RP5MT", "RTOE", "RP1MT"], \
                  ["LHEEL", "LLMAL", "LMFL", "LP5MT", "LTOE", "LP1MT"], \
                  ["RHEEL", "RLMAL", "RMFL", "RP5MT", "RTOE", "LTOE", "LP5MT", "LMFL", "LLMAL", "LHEEL"]]
    
    # Get GRF
    # Note: OpenSim GRF files are set up as FP1: left, FP2: right, even though 
    # the lab setup is FP1: right, FP2: left.
    grftime = forces["time"]
    grfy_right = forces["data"]["right"]["F"][:, 1]
    grfy_left = forces["data"]["left"]["F"][:, 1]
        
    # Determine single or double support from GRF:
    #   0: no support
    #   1: single support, right leg
    #   2: single support, left leg
    #   3: double support
    issingle_fp_right = np.array([(fp > 0) * 1 for fp in grfy_right]).astype(int)
    issingle_fp_left = np.array([(fp > 0) * 2 for fp in grfy_left]).astype(int)
    issupport = issingle_fp_right + issingle_fp_left
    
    # Interpolation function: GRF support (zero-order)
    supportfun = interp.interp1d(grftime, issupport, kind = "zero", fill_value = "extrapolate", axis = 0)
    
    # Interpolation functions: 2D marker position (cubic spline)
    mkrfun = {}
    for m in markers:
        if m in ["frames", "offset", "rate", "time", "units"]:
            continue
        else:
            mkrfun[m] = interp.interp1d(markers["time"], markers[m][:, [0, 2]], kind = "cubic", fill_value = "extrapolate", axis = 0)
        
    
    # Dynamic base of support
    bos = {}
    bos["coordinates"] = []
    bos["baseshape"] = []
    bos["support"] = []
    bos["markers"] = []
    for t, time in enumerate(timevec):

        # Determine instantaneous GRF support type
        supportval = supportfun(time)
        bos["support"].append(supportval)
        
        # Get the instantaneous marker list
        instantmkrs = bosmarkers[int(supportval)]
        bos["markers"].append(instantmkrs)
        
        # Get the instantaneous marker coordinates
        mkrcoords = [mkrfun[m](time) for m in instantmkrs]
        bos["coordinates"].append(mkrcoords)
        
        # Build base of support Polygon, but return the convex hull
        baseshape = shapely.MultiPoint(mkrcoords).convex_hull
        bos["baseshape"].append(baseshape)
        
    return bos



'''
resample1d(data, nsamp):
    Simple resampling by 1-D interpolation (rows = samples, cols = variable).
    Data can be a 1-D or multiple variables in a 2-D array-like object.
'''
def resample1d(data, nsamp):

    # Convert list to 1D array
    if isinstance(data, list):
        data = np.array([data]).transpose()
        ny = 1        

    # Data dimensions
    nx = data.shape[0]
    ny = data.shape[1]

    # Old sample points
    x = np.linspace(0, nx - 1, nx)
        
    # New sample points
    xnew = np.linspace(0, nx - 1, nsamp)
    
    # Parse columns
    datanew = np.zeros([nsamp, ny])
    for col in range(0, ny):
        
        # Old data points
        y = data[:, col]
        
        # Convert to cubic spline function
        fy = interp.interp1d(x, y, kind = "cubic", fill_value = "extrapolate")
    
        # New data points
        ynew = fy(xnew)
    
        # Store column
        datanew[:, col] = ynew
        
    
    return datanew    



'''
-----------------------------------
-------- EXPORT FUNCTIONS ---------
-----------------------------------
'''


'''
export_stability_metrics(meta, user, nsamp, normalise):
    Write stablity time histories to CSV.
    NOTE: Currently hardcoded for SLDJ
'''
def export_stability_metrics(meta, user, nsamp, normalise = False):
   
    # empty output list of lists
    # (create the output table as a list of lists, then convert to dataframe
    # as iteratively appending new dataframe rows is computationally expensive)
    csvdata = []
        
    # extract OpenSim data
    print("Collating data into lists...\n")
    failedfiles = []
    for sn, subj in enumerate(meta):
    
        print("%s" % "*" * 30)
        print("SUBJECT: %s" % subj)
        print("%s" % "*" * 30)

        # subject index (may be different to the subject code)
        subjidx = sn

        # subject type
        if subj.startswith("FAILTCRT"):
            subj_type = "ctrl"
            subj_type_code = 1
        else:
            subj_type = "sym"
            subj_type_code = 0
                
        for group in meta[subj]["trials"]:
            
            print("Group: %s" % group)
            print("%s" % "=" * 30)                      
            
            # process dynamic trials only
            for trial in  meta[subj]["trials"][group]:                
                
                # ignore static trials
                isstatic = meta[subj]["trials"][group][trial]["isstatic"]
                if isstatic: continue
            
                try:
                
                    # load the trial StabilityKey
                    c3dpath = meta[subj]["trials"][group][trial]["outpath"]
                    pkfile = os.path.join(c3dpath,trial + "_stability.pkl")
                    with open(pkfile,"rb") as fid:
                        stabilitykey = pk.load(fid)

                    # load the trial WBAMKey
                    c3dpath = meta[subj]["trials"][group][trial]["outpath"]
                    pkfile = os.path.join(c3dpath,trial + "_wbam.pkl")
                    with open(pkfile,"rb") as fid:
                        wbamkey = pk.load(fid)
                        
                    # load the trial OsimResultsKey
                    c3dpath = meta[subj]["trials"][group][trial]["outpath"]
                    pkfile = os.path.join(c3dpath,trial + "_opensim_results.pkl")
                    with open(pkfile,"rb") as fid:
                        osimresultskey = pk.load(fid)
                                            
                    # participant data
                    age = osimresultskey.age
                    mass = osimresultskey.mass
                    height = osimresultskey.height
                    sex = osimresultskey.sex
                    dom_foot = osimresultskey.dom_foot
                    aff_side = osimresultskey.aff_side
                    shomri_r = osimresultskey.shomri[0]
                    shomri_l = osimresultskey.shomri[1]
                    
                    # for bilateral symptomatics, affected side based on shomri
                    if aff_side == 1:
                        more_aff_side = "r"
                    elif aff_side == 2:
                        more_aff_side = "l"
                    elif aff_side == 0 or aff_side == 3:
                        if np.isnan(shomri_r) or np.isnan(shomri_l):
                            more_aff_side = "r"   # default until data available
                        if shomri_r > shomri_l:
                            more_aff_side = "r"
                        else:
                            more_aff_side = "l"
                        
                        
                    # trial task type
                    task = osimresultskey.task

                    # pivot leg
                    #pivot_leg = osimresultskey.events["labels"][0][0].lower()
                    
                    # generic event labels:
                    #   FS: Foot-strike
                    #   F0: foot-off
                    #events_gen_labels = ["FS", "FO"]
                    
                    # event timing: relative time and time steps
                    #events_times = osimresultskey.events["time"] - osimresultskey.events["time"][0]
                    #events_steps = np.round(user.samples * (osimresultskey.events["time"] - osimresultskey.events["time"][0]) / (osimresultskey.events["time"][5] - osimresultskey.events["time"][0]))
                    
                    # Centre of mass velocity
                    comvidx = osimresultskey.results["raw"][user.bkcode]["headers"][1].index("center_of_mass_X")
                    comv = osimresultskey.results["raw"][user.bkcode]["data"][:, comvidx:comvidx + 3, 1]
                    
                    
                    # foot
                    for f, foot in enumerate(["r","l"]):
                        
                        
                        # NOTE: CURRENTLY HARDCODED FOR SLDJ
                        
                        # skip contralateral leg data
                        if osimresultskey.events["labels"][0][0] != foot.upper(): continue                        
                        
                        # trial leg: more affected or less affected
                        if foot == more_aff_side:
                            trial_leg = "less"
                        else:
                            trial_leg = "more"                        
                        
                                                
                        # Get the Euclidean (2D) margin of stability and check                       
                        # if it needs to be trimmed if events are not accurately
                        # labelled. Return the first and last index.
                        b = np.array(stabilitykey["MoS"]["b"])                       
                        bidx0 = 0
                        bidx1 = len(b) - 1
                        while ((b.dtype is np.dtype("bool")) and (b[bidx0] is None)) or ((b.dtype is not np.dtype("bool")) and (b[bidx0] == 0.0)): bidx0 = bidx0 + 1
                        while ((b.dtype is np.dtype("bool")) and (b[bidx1] is None)) or ((b.dtype is not np.dtype("bool")) and (b[bidx1] == 0.0)): bidx1 = bidx1 - 1
                        
                        
                        # margin of stability components
                        for ans in ["b", "b_abs", "b_x", "b_z", "isstable"]:
                        
                            # Timeseries data, convert to 1D array and trim
                            drow = np.array(stabilitykey["MoS"][ans])
                            drow = drow[bidx0:bidx1 + 1]
                                                                                                                
                            # Normalisation factors
                            # TBD: Need to determine an appropriate normalisation
                            # factor, e.g., leg length or height. Use height for
                            # the moment.
                            normfactor = 1.0
                            if normalise:
                                normfactor = 1.0 / height
                                            
                            # Normalise if required
                            if ans != "isstable":
                                drow = drow * normfactor
                            
                            # Resample if required
                            if len(drow) != nsamp:
                                drow = resample1d(np.reshape(drow, (len(drow), 1)), nsamp)                            
                            
                            # create new line of data
                            csvrow = [subj, subjidx, trial, subj_type, subj_type_code, task, foot, age, mass, height, sex, dom_foot, aff_side, shomri_r, shomri_l, more_aff_side, trial_leg, -1, -1, ans] + drow.flatten().tolist()
                            csvdata.append(csvrow)



                        # whole body angular momentum components
                        wbam_int = wbamkey["L_int"]
                        wbam_range = wbamkey["L_range"]
                        for ans in ["L"]:
                        
                            # Timeseries data
                            dmat = wbamkey[ans]
                            dmat = dmat[bidx0:bidx1 + 1, :]
                                                                                    
                            # Normalisation factors
                            normfactor = 1.0
                            if normalise:
                                normfactor = 1.0 / (mass * height * wbamkey["CoM_v_mean"])
                                            
                            # Normalise if required
                            dmat = dmat * normfactor
                            wbam_int = wbam_int * normfactor
                            wbam_range = wbam_range * normfactor

                            # Resample if required
                            if dmat.shape[0] != nsamp:
                                dmat = resample1d(dmat, nsamp)
                                
                            # create new line of data
                            for d, dim in enumerate(["X", "Y", "Z"]):
                                csvrow = [subj, subjidx, trial, subj_type, subj_type_code, task, foot, age, mass, height, sex, dom_foot, aff_side, shomri_r, shomri_l, more_aff_side, trial_leg] + [wbam_int[d], wbam_range[d]] + [ans + "_" + dim] + dmat[:, d].flatten().tolist()
                                csvdata.append(csvrow)
                                 
                
                except:
                    print("Dynamic trial: %s *** FAILED ***" % trial)
                    failedfiles.append(trial)
                    #raise
                else:
                    print("Dynamic trial: %s" % trial)

    # create dataframe
    print("\nCreating dataframe...")
    headers = ["subject", "subj_idx", "trial", "subj_type", "subj_type_code", "task", "data_leg", "age", "mass", "height", "sex", "dom_foot", "aff_side", "shomri_r", "shomri_l", "more_aff_leg", "leg_type", "var_integral", "var_range", "variable"] + ["t" + str(n) for n in range(1,102)]
    csvdf = pd.DataFrame(csvdata, columns = headers)

    # write data to file with headers
    print("\nWriting to CSV text file...")
    normalisestr = ["", "_normalised"]
    csvfile = user.csvfileprefix + "_stability" + normalisestr[int(normalise)] + ".csv"
    fpath = os.path.join(user.rootpath, user.outfolder, user.csvfolder)
    if not os.path.exists(fpath): os.makedirs(fpath)
    csvdf.to_csv(os.path.join(fpath,csvfile), index = False)
    
    
    print("\n")
   
    return failedfiles



'''
export_stability_metrics_subject_mean(meta, user, nsamp, normalise):
    Write stablity time histories to CSV.
    NOTE: Currently hardcoded for SLDJ
'''
def export_stability_metrics_subject_mean(meta, user, nsamp, normalise = False):
   
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

        # subject type
        if subj.startswith("FAILTCRT"):
            subj_type = "ctrl"
        else:
            subj_type = "sym"
                
        for group in meta[subj]["trials"]:
            
            print("Group: %s" % group)
            print("%s" % "=" * 30)                      
            
            # process dynamic trials only
            for trial in  meta[subj]["trials"][group]:                
                
                # ignore static trials
                isstatic = meta[subj]["trials"][group][trial]["isstatic"]
                if isstatic: continue
            
                try:
                
                    # load the trial StabilityKey
                    c3dpath = meta[subj]["trials"][group][trial]["outpath"]
                    pkfile = os.path.join(c3dpath,trial + "_stability.pkl")
                    with open(pkfile,"rb") as fid:
                        stabilitykey = pk.load(fid)
                        
                    # load the trial WBAMKey
                    c3dpath = meta[subj]["trials"][group][trial]["outpath"]
                    pkfile = os.path.join(c3dpath,trial + "_wbam.pkl")
                    with open(pkfile,"rb") as fid:
                        wbamkey = pk.load(fid)                        

                    # load the trial OsimResultsKey
                    c3dpath = meta[subj]["trials"][group][trial]["outpath"]
                    pkfile = os.path.join(c3dpath,trial + "_opensim_results.pkl")
                    with open(pkfile,"rb") as fid:
                        osimresultskey = pk.load(fid)
                                            
                    # participant data
                    age = osimresultskey.age
                    mass = osimresultskey.mass
                    height = osimresultskey.height
                    sex = osimresultskey.sex
                    dom_foot = osimresultskey.dom_foot
                    aff_side = osimresultskey.aff_side
                    shomri_r = osimresultskey.shomri[0]
                    shomri_l = osimresultskey.shomri[1]
                    
                    # for bilateral symptomatics, affected side based on shomri
                    if aff_side == 1:
                        more_aff_side = "r"
                    elif aff_side == 2:
                        more_aff_side = "l"
                    elif aff_side == 0 or aff_side == 3:
                        if np.isnan(shomri_r) or np.isnan(shomri_l):
                            more_aff_side = "r"   # default until data available
                        if shomri_r > shomri_l:
                            more_aff_side = "r"
                        else:
                            more_aff_side = "l"
                        
                        
                    # trial task type
                    task = osimresultskey.task
                   
                    # foot
                    for f, foot in enumerate(["r","l"]):
                        
                        
                        # NOTE: CURRENTLY HARDCODED FOR SLDJ
                        
                        # skip contralateral leg data
                        if osimresultskey.events["labels"][0][0] != foot.upper(): continue                        
                        
                        # trial leg: more affected or less affected
                        if foot == more_aff_side:
                            trial_leg = "less"
                        else:
                            trial_leg = "more"                        
                        
                                                
                        # Get the Euclidean (2D) margin of stability and check                       
                        # if it needs to be trimmed if events are not accurately
                        # labelled. Return the first and last index.
                        b = np.array(stabilitykey["MoS"]["b"])                       
                        bidx0 = 0
                        bidx1 = len(b) - 1
                        while ((b.dtype is np.dtype("bool")) and (b[bidx0] is None)) or ((b.dtype is not np.dtype("bool")) and (b[bidx0] == 0.0)): bidx0 = bidx0 + 1
                        while ((b.dtype is np.dtype("bool")) and (b[bidx1] is None)) or ((b.dtype is not np.dtype("bool")) and (b[bidx1] == 0.0)): bidx1 = bidx1 - 1
                        
                        
                        # margin of stability components
                        for ans in ["b", "b_abs", "b_x", "b_z", "isstable"]:
                        
                            # Data, convert to 1D array and trim
                            drow = np.array(stabilitykey["MoS"][ans])
                            drow = drow[bidx0:bidx1 + 1]
                                                            
                            # Resample if required
                            if len(drow) != nsamp:
                                drow = resample1d(np.reshape(drow, (len(drow), 1)), nsamp)
                                                                                                                
                            # Normalisation factors
                            # TBD: Need to determine an appropriate normalisation
                            # factor, e.g., leg length or height. Use height for
                            # the moment.
                            normfactor = 1.0
                            if normalise:
                                normfactor = 1.0 / height
                                            
                            # Normalise if required
                            if ans != "isstable":
                                drow = drow * normfactor
                            
                            # Create new line of data
                            csvrow = [subj, trial, subj_type, task, foot, age, mass, height, sex, dom_foot, aff_side, shomri_r, shomri_l, more_aff_side, trial_leg, ans] + drow.flatten().tolist()
                            csvdata.append(csvrow)
                            
                            
                        # Whole body angular momentum components
                        for ans in ["L"]:
                        
                            # Data
                            dmat = wbamkey[ans]
                            dmat = dmat[bidx0:bidx1 + 1, :]
                            
                            # Flip signs for frontal and transverse plane
                            # components for left leg trials
                            if foot == "l":
                                dmat[:, 0:2] = -1 * dmat[:, 0:2]

                                
                            # Resample if required
                            if dmat.shape[0] != nsamp:
                                dmat = resample1d(dmat, nsamp)
                                                                                                                
                            # Normalisation factors
                            normfactor = 1.0
                            if normalise:
                                normfactor = 1.0 / (mass * height * wbamkey["CoM_v_mean"])
                                            
                            # Normalise if required
                            dmat = dmat * normfactor
                            
                            # Create new line of data
                            for d, dim in enumerate(["X", "Y", "Z"]):
                                csvrow = [subj, trial, subj_type, task, foot, age, mass, height, sex, dom_foot, aff_side, shomri_r, shomri_l, more_aff_side, trial_leg] + [ans + "_" + dim] + dmat[:, d].flatten().tolist()
                                csvdata.append(csvrow)
                                
                
                except:
                    print("Dynamic trial: %s *** FAILED ***" % trial)
                    failedfiles.append(trial)
                    #raise
                else:
                    print("Dynamic trial: %s" % trial)

    # create timeseries dataframe
    print("\nCreating dataframes...")
    headers = ["subject", "trial", "subj_type", "task", "data_leg", "age", "mass", "height", "sex", "dom_foot", "aff_side", "shomri_r", "shomri_l", "more_aff_leg", "leg_type", "variable"] + ["t" + str(n) for n in range(1,102)]    
    csvdf = pd.DataFrame(csvdata, columns = headers)     

    # group
    csvdf.drop("trial", axis = 1, inplace = True)   # std() doesn't like the trial column
    csvdf_grouped = csvdf.groupby(["subject", "subj_type", "task", "data_leg", "age", "mass", "height", "sex", "dom_foot", "aff_side", "shomri_r", "shomri_l", "more_aff_leg", "leg_type", "variable"])

    # descriptives
    csvdf_grouped_mean = csvdf_grouped.mean().reset_index()
    csvdf_grouped_sd = csvdf_grouped.std().reset_index()
    
    # rearrange dataframes (much easier in dplyr with relocate()!)
    csvdf_grouped_mean["statistic"] = "mean"
    dfmean = csvdf_grouped_mean.pop("statistic")
    csvdf_grouped_mean.insert(csvdf_grouped_mean.columns.get_loc("variable") + 1, dfmean.name, dfmean)      
    csvdf_grouped_sd["statistic"] = "sd"
    dfsd = csvdf_grouped_sd.pop("statistic")
    csvdf_grouped_sd.insert(csvdf_grouped_sd.columns.get_loc("variable") + 1, dfsd.name, dfsd)    
    
    # interleave mean and sd rows
    csvdf_grouped_mean["sortidx"] = range(0, len(csvdf_grouped_mean))
    csvdf_grouped_sd["sortidx"] = range(0, len(csvdf_grouped_sd))
    csvdf_descriptives = pd.concat([csvdf_grouped_mean, csvdf_grouped_sd]).sort_values(["sortidx", "statistic"])
    csvdf_descriptives.drop(columns = "sortidx", inplace = True)

    # write data to file with headers
    print("\nWriting to CSV text file...")
    normalisestr = ["", "_normalised"]
    csvfile = user.csvdescfileprefix + "_stability" + normalisestr[int(normalise)] + ".csv"
    fpath = os.path.join(user.rootpath, user.outfolder, user.csvfolder)
    if not os.path.exists(fpath): os.makedirs(fpath)
    csvdf_descriptives.to_csv(os.path.join(fpath,csvfile), index = False)
    
    print("\n")
   
    return failedfiles




'''
export_wbam_discrete_subject_mean(meta, user, normalise):
    Write stablity time histories to CSV.
    NOTE: Currently hardcoded for SLDJ
'''
def export_wbam_discrete_subject_mean(meta, user, normalise = False):
   
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

        # subject type
        if subj.startswith("FAILTCRT"):
            subj_type = "ctrl"
        else:
            subj_type = "sym"
                
        for group in meta[subj]["trials"]:
            
            print("Group: %s" % group)
            print("%s" % "=" * 30)                      
            
            # process dynamic trials only
            for trial in  meta[subj]["trials"][group]:                
                
                # ignore static trials
                isstatic = meta[subj]["trials"][group][trial]["isstatic"]
                if isstatic: continue
            
                try:
                        
                    # load the trial WBAMKey
                    c3dpath = meta[subj]["trials"][group][trial]["outpath"]
                    pkfile = os.path.join(c3dpath,trial + "_wbam.pkl")
                    with open(pkfile,"rb") as fid:
                        wbamkey = pk.load(fid)                        

                    # load the trial OsimResultsKey
                    c3dpath = meta[subj]["trials"][group][trial]["outpath"]
                    pkfile = os.path.join(c3dpath,trial + "_opensim_results.pkl")
                    with open(pkfile,"rb") as fid:
                        osimresultskey = pk.load(fid)
                                            
                    # participant data
                    age = osimresultskey.age
                    mass = osimresultskey.mass
                    height = osimresultskey.height
                    sex = osimresultskey.sex
                    dom_foot = osimresultskey.dom_foot
                    aff_side = osimresultskey.aff_side
                    shomri_r = osimresultskey.shomri[0]
                    shomri_l = osimresultskey.shomri[1]
                    
                    # for bilateral symptomatics, affected side based on shomri
                    if aff_side == 1:
                        more_aff_side = "r"
                    elif aff_side == 2:
                        more_aff_side = "l"
                    elif aff_side == 0 or aff_side == 3:
                        if np.isnan(shomri_r) or np.isnan(shomri_l):
                            more_aff_side = "r"   # default until data available
                        if shomri_r > shomri_l:
                            more_aff_side = "r"
                        else:
                            more_aff_side = "l"
                        
                        
                    # trial task type
                    task = osimresultskey.task

                    # foot
                    for f, foot in enumerate(["r","l"]):
                        
                        
                        # NOTE: CURRENTLY HARDCODED FOR SLDJ
                        
                        # skip contralateral leg data
                        if osimresultskey.events["labels"][0][0] != foot.upper(): continue                        
                        
                        # trial leg: more affected or less affected
                        if foot == more_aff_side:
                            trial_leg = "less"
                        else:
                            trial_leg = "more"                        
                        
    
                        # Whole body angular momentum components
                        for ans in ["L_int", "L_range"]:
                        
                            # Data
                            dmat = wbamkey[ans]
                            
                            # Flip signs for frontal and transverse plane
                            # components for left leg trials
                            if (foot == "l") and (ans == "L_int"):
                                dmat[0:2] = -1 * dmat[0:2]
                                                                                                                
                            # Normalisation factors
                            normfactor = 1.0
                            if normalise:
                                normfactor = 1.0 / (mass * height * wbamkey["CoM_v_mean"])
                                            
                            # Normalise if required
                            dmat = dmat * normfactor
                            
                            # Create new line of data
                            for d, dim in enumerate(["X", "Y", "Z"]):
                                csvrow = [subj, trial, subj_type, task, foot, age, mass, height, sex, dom_foot, aff_side, shomri_r, shomri_l, more_aff_side, trial_leg] + [ans + "_" + dim] + [dmat[d]]
                                csvdata.append(csvrow)
                                
                
                except:
                    print("Dynamic trial: %s *** FAILED ***" % trial)
                    failedfiles.append(trial)
                    #raise
                else:
                    print("Dynamic trial: %s" % trial)

    # create timeseries dataframe
    print("\nCreating dataframes...")
    headers = ["subject", "trial", "subj_type", "task", "data_leg", "age", "mass", "height", "sex", "dom_foot", "aff_side", "shomri_r", "shomri_l", "more_aff_leg", "leg_type", "variable", "value"]
    csvdf = pd.DataFrame(csvdata, columns = headers)     

    # group
    csvdf.drop("trial", axis = 1, inplace = True)   # std() doesn't like the trial column
    csvdf_grouped = csvdf.groupby(["subject", "subj_type", "task", "data_leg", "age", "mass", "height", "sex", "dom_foot", "aff_side", "shomri_r", "shomri_l", "more_aff_leg", "leg_type", "variable"])

    # descriptives
    csvdf_grouped_mean = csvdf_grouped.mean().reset_index()
    csvdf_grouped_sd = csvdf_grouped.std().reset_index()
    
    # rearrange dataframes (much easier in dplyr with relocate()!)
    csvdf_grouped_mean["statistic"] = "mean"
    dfmean = csvdf_grouped_mean.pop("statistic")
    csvdf_grouped_mean.insert(csvdf_grouped_mean.columns.get_loc("variable") + 1, dfmean.name, dfmean)      
    csvdf_grouped_sd["statistic"] = "sd"
    dfsd = csvdf_grouped_sd.pop("statistic")
    csvdf_grouped_sd.insert(csvdf_grouped_sd.columns.get_loc("variable") + 1, dfsd.name, dfsd)    
    
    # interleave mean and sd rows
    csvdf_grouped_mean["sortidx"] = range(0, len(csvdf_grouped_mean))
    csvdf_grouped_sd["sortidx"] = range(0, len(csvdf_grouped_sd))
    csvdf_descriptives = pd.concat([csvdf_grouped_mean, csvdf_grouped_sd]).sort_values(["sortidx", "statistic"])
    csvdf_descriptives.drop(columns = "sortidx", inplace = True)

    # write data to file with headers
    print("\nWriting to CSV text file...")
    normalisestr = ["", "_normalised"]
    csvfile = user.csvdescfileprefix + "_wbam_discrete" + normalisestr[int(normalise)] + ".csv"
    fpath = os.path.join(user.rootpath, user.outfolder, user.csvfolder)
    if not os.path.exists(fpath): os.makedirs(fpath)
    csvdf_descriptives.to_csv(os.path.join(fpath,csvfile), index = False)
    
    print("\n")
   
    return failedfiles







'''
-----------------------------------
---------- VISUALISATION ----------
-----------------------------------
'''



'''
plot_whole_body_angular_momentum(datakey, wbam, showplot):
    Plot the whole body angular momentum.
'''
def plot_whole_body_angular_momentum(datakey, wbam, showplot = False, colours = ["black", "blue", "red"], plotlimits = [(-1, 1), (-1, 1), (-0.15, 0.15)]):
    
    subject = datakey.subject
    trial = datakey.trial
    
    figpath = os.path.join(datakey.outpath, "viz", "WBAM", "timeseries")  
    if not os.path.exists(figpath): os.makedirs(figpath)   
            
    # Global parameters
    plotlabels = ["Frontal plane (OpenSim X axis)", "Transverse plane (OpenSim Y axis)", "Sagittal plane (OpenSim Z axis)"]
    filelabels= ["Fore-aft", "Vertical", "Mediolateral"]
    #plotlimits = [(-1, 1), (-1, 1), (-0.15, 0.15)] 
    
    # Plot
    for p in range(3):
        plt.figure()
        plt.plot(wbam["L"][:, p], label = plotlabels[p], color = colours[p])
        plt.title("WBAM: " + plotlabels[p])
        plt.xlabel(xlabel = "Time step")
        plt.ylabel("WBAM (kgm$^2$/s)")

        plt.savefig(os.path.join(figpath, subject + "_" + trial + "_WBAM_" + filelabels[p] + ".png"))   
        
        if not showplot: plt.close()
    
    return None


'''
plot_margin_of_stability(datakey, stable, showplot, colours, plotlimits):
    Plot the margin of stability (b): 2D, X and Z.
'''
def plot_margin_of_stability(datakey, stable, showplot = False, colours = ["black", "blue", "red"], plotlimits = [(-1, 1), (-1, 1), (-0.2, 0.2)]):
    
    subject = datakey.subject
    trial = datakey.trial
    
    figpath = os.path.join(datakey.outpath, "viz", "MoS", "timeseries")  
    if not os.path.exists(figpath): os.makedirs(figpath)   
            
    # Global parameters
    moslabels = ["b", "b_x", "b_z"]
    plotlabels = ["Euclidean: b", "Fore-aft (OpenSim X): bx", "Mediolateral (OpenSim Z): bz"]
    filelabels = ["Euclidean", "Fore-aft", "Mediolateral"]
    
    # Plot
    for p in range(3):
        plt.figure()
        plt.plot(stable["MoS"][moslabels[p]], label = plotlabels[p], color = colours[p])
        plt.title("MoS: " + plotlabels[p])
        plt.xlabel(xlabel = "Time step")
        plt.ylabel("MoS (m)")
        plt.ylim(plotlimits[p])

        plt.savefig(os.path.join(figpath, subject + "_" + trial + "_MoS_" + filelabels[p] + ".png")) 
        
        if not showplot: plt.close()
        
    return None
    
    

'''
visualise_stability_timehistory(datakey, stable, showplot, colours, xlim, ylim):
    Visualise the base of support (BoS), centre of mass (CoM) and extrapolated 
    centre of mass (XCoM) at each time step. Export as individual frames.
'''
def visualise_stability_timehistory(datakey, stable, showplot = False, colours = ["blue", "red", "black"], xlim = [-1, 1], ylim = [-1, 1]):

    # Output folder
    figpath = os.path.join(datakey.outpath, "viz", "MoS", "slides")  
    if not os.path.exists(figpath): os.makedirs(figpath)

    # Plot time history of base of support, centre of mass, and extrapolated
    # centre of mass, save to file
    for t in range(len(stable["BoS"]["baseshape"])):
        
        fig = plt.figure()
        fig.suptitle(datakey.trial + " Step: " + str(t))
        
        # Plot
        baseshape = stable["BoS"]["baseshape"][t]
        com = stable["CoM"]["CoM_pt"][t]
        xcom = stable["CoM"]["XCoM_pt"][t]
        plot_bos_com_xcom(fig, baseshape, com, xcom, colours, xlim, ylim)     
        #plt.show()
        
        # Save
        fig.savefig(os.path.join(figpath, datakey.trial + "_" + str(t) + ".png"))
        
        if not showplot: plt.close()

    return None
    


'''
animate_stability_timehistory(datakey, stable, showplot):
    Visualise the base of support (BoS), centre of mass (CoM) and extrapolated 
    centre of mass (XCoM) at each time step. Export as MP4 animation.
'''    
def animate_stability_timehistory(datakey, stable, showplot = False, xlim = [-0.5, 0.5], ylim = [-0.5, 1.0]):
          
    fig = plt.figure()
    ax = fig.gca()
    
    # Initialise animation
    def init_func():
        ax.clear()
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)
        ax.set_aspect("equal", adjustable="box")
        
    # Animator function
    def animate(i):
    
        ax.clear()
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)
        ax.set_xlabel("Mediolateral (m)")
        ax.set_ylabel("Fore-aft (m)")
        ax.set_aspect("equal", adjustable="box")
        ax.set_title(datakey.trial + " Step: " + str(i))
        
        # Plot base of support shape if it exists (must be type Polygon)
        if type(stable["BoS"]["baseshape"][i]) == shapely.geometry.polygon.Polygon:
            shapex, shapey = stable["BoS"]["baseshape"][i].exterior.xy
            ax.plot(shapey, shapex, color = "green")
        
        # Centre of mass
        comx, comy = stable["CoM"]["CoM_pt"][i].xy
        ax.plot(comy, comx, marker = ".", markersize = 20, color = "blue")
    
        # Extrapolated centre of mass
        xcomx, xcomy = stable["CoM"]["XCoM_pt"][i].xy
        ax.plot(xcomy, xcomx, marker = ".", markersize = 20, color = "red") 
        
        # Connect with line
        ax.plot([comy, xcomy], [comx, xcomx], linestyle = "-", color = "black")

    # Animate
    anim = FuncAnimation(fig, animate, init_func=init_func, frames = len(stable["BoS"]["baseshape"]))
    
    # Save a MP4    
    figpath = os.path.join(datakey.outpath, "viz", "MoS", "animation")  
    if not os.path.exists(figpath): os.makedirs(figpath) 
    anim.save(os.path.join(figpath, datakey.trial + ".mp4"))
    
    if not showplot: plt.close()
    
    return None
    


'''
plot_stability_artificial_perturbation(datakey, perturb):
    Visualise the variation in margin of stability if baseline motion is
    artificially perturbed in all directions. Export as individual frames.
'''   
def plot_stability_artificial_perturbation(datakey, perturb):
        
    # Output folder
    figpath = os.path.join(datakey.trial["path"], "viz", "figures", "peturb")  
    if not os.path.exists(figpath): os.makedirs(figpath)    
    
    for t in range(len(perturb["BoS"]["baseshape"])):
    
        fig = plt.figure(constrained_layout=True, figsize=(12, 6))
        fig.suptitle(datakey.trial["name"] + "    Perturbation: " + str(perturb["perturbation"]) + "m/s    Step: " + str(t))
        spec = fig.add_gridspec(nrows = 1, ncols = 5, width_ratios =[0.1, 1, 0.1, 1, 0.1]) 
        ax1 = fig.add_subplot(spec[0, 1], projection = "polar")
        ax2 = fig.add_subplot(spec[0, 3])
        
        # Polar plot       
        ax1.set_theta_offset(np.pi/2)
        ax1.plot(np.linspace(0, 2 * np.pi, 100), np.zeros(100), color="black", linestyle = '--')  # Reference circle: zero margin of stability, i.e. edge of base of support    
        angles = perturb["angles"]
        mos = np.array(perturb["MoS"]["b_perturb"][t])
        ax1.plot(angles, mos, markersize = 20)
        ax1.set_rlim(2, -2)  # Reversed axis
        ax1.set_title("Margin of stability (m)")
            
        # Visualisation of baseline
        baseshape = perturb["BoS"]["baseshape"][t]
        com = perturb["CoM"]["CoM_pt"][t]
        xcom0 = perturb["CoM"]["XCoM_pt"][t]
        plot_bos_com_xcom(fig, baseshape, com, xcom0)
        ax2.set_title("Visualisation")
        ax2.legend()
        
        # Visualisation of perturbed extrapolated centre of mass
        for a in range(np.size(perturb["angles"])):
            xcom = perturb["CoM"]["XCoM_pt_perturb"][t][a]
            plot_bos_com_xcom(fig, baseshape, com, xcom, ("blue", "orange", "grey"))
            
        plt.plot()
        
        # Save
        fig.savefig(os.path.join(figpath, datakey.trial["name"] + "_perturbed_" + str(t) + ".png"))   
    
    return None



'''
stability_artificial_perturbation(datakey, perturb):
    Visualise the variation in margin of stability if baseline motion is
    artificially perturbed in all directions. Export as MP4.
'''    
def animate_stability_artificial_perturbation(datakey, perturb):
          
    fig = plt.figure()
    fig = plt.figure(constrained_layout=True, figsize=(12, 6))
    spec = fig.add_gridspec(nrows = 1, ncols = 5, width_ratios =[0.1, 1, 0.1, 1, 0.1])
    ax1 = fig.add_subplot(spec[0, 1])
    ax2 = fig.add_subplot(spec[0, 3])
    
    # Initialise animation
    def init_func():
        ax1.clear() 
        ax2.clear()
        
    # Animator function
    def animate(i):
            
        print(str(i))
        
        fig.suptitle(datakey.trial["name"] + "    Perturbation: " + str(perturb["perturbation"]) + "m/s    Step: " + str(i))
        ax1 = fig.add_subplot(spec[0, 1], projection = "polar")    
        ax2 = fig.add_subplot(spec[0, 3])
 
        # Polar plot
        ax1.set_theta_offset(np.pi/2)
        ax1.plot(np.linspace(0, 2 * np.pi, 100), np.zeros(100), color="black", linestyle = '--')  # Reference circle: zero margin of stability, i.e. edge of base of support    
        angles = perturb["angles"]
        mos = np.array(perturb["MoS"]["b_perturb"][i])
        ax1.plot(angles, mos, markersize = 20)
        ax1.set_rlim(2, -2)  # Reversed axis
        ax1.set_title("Margin of stability (m)")
            
        # Visualisation of baseline        
        baseshape = perturb["BoS"]["baseshape"][i]
        com = perturb["CoM"]["CoM_pt"][i]
        xcom0 = perturb["CoM"]["XCoM_pt"][i]
        plot_bos_com_xcom(fig, baseshape, com, xcom0)
        ax2.set_title("Visualisation")
        ax2.legend()
        
        # Visualisation of perturbed extrapolated centre of mass
        for a in range(np.size(perturb["angles"])):
            xcom = perturb["CoM"]["XCoM_pt_perturb"][i][a]
            plot_bos_com_xcom(fig, baseshape, com, xcom, ("blue", "orange", "grey"))

    # Animate
    anim = FuncAnimation(fig, animate, init_func = init_func, frames = len(perturb["BoS"]["baseshape"]))
    
    # Save a MP4    
    figpath = os.path.join(datakey.trial["path"], "viz", "animation")  
    if not os.path.exists(figpath): os.makedirs(figpath) 
    anim.save(os.path.join(figpath, datakey.trial["name"] + "_perturbed.mp4"))
    
    return None



'''
plot_bos_com_xcom(fig, baseshape, com, xcom, colours, xlim, ylim):
    Plot instantaneous base of support (BoS), centre of mass (CoM) and 
    extrapolated centre of mass (XCoM). Coordinates flipped (z-x) for easier 
    visualisation.
'''
def plot_bos_com_xcom(fig, baseshape, com, xcom, colours = ("blue", "red", "black"), xlim = [-1, 1], ylim = [-1, 1]):
    
    #fig = plt.figure()
    
    # Plot base of support shape if it exists (must be type Polygon)
    if type(baseshape) == shapely.geometry.polygon.Polygon:
        shapex, shapey = baseshape.exterior.xy
        plt.plot(shapey, shapex, color = "green")
    
    # Centre of mass
    comx, comy = com.xy
    plt.plot(comy, comx, marker = ".", markersize = 20, label = "CoM", color = colours[0])

    # Extrapolated centre of mass
    xcomx, xcomy = xcom.xy
    plt.plot(xcomy, xcomx, marker = ".", markersize = 20, label = "XCoM", color = colours[1])    
    
    # Connect with line
    plt.plot([comy, xcomy], [comx, xcomx], linestyle = "-", color = colours[2])
    
    plt.xlim(xlim)
    plt.ylim(ylim)
    plt.xlabel("Mediolateral (m)")
    plt.ylabel("Fore-aft (m)")
    
    ax = plt.gca()
    ax.set_aspect("equal", adjustable="box")
    
    return fig
    
    

    
    