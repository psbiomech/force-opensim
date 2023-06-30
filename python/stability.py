# -*- coding: utf-8 -*-
"""
Stability analyses: 
    Margin of stability
    Probability of Instability
    Whole Body Angular Momentum

@author: Prasanna Sritharan
"""

import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import numpy as np
import shapely
import scipy.constants as consts
import scipy.interpolate as interp
import pickle as pk
import os


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
batch_process_stability(user, meta, perturb, treadmill_speed, asteps, perturbation):
    Batch process margin of stability with or without artificial perturbation.
'''
def batch_process_stability(user, meta, artifperturb = False, treadmill_speed = 0.0, asteps = 16, perturbation = 1.0):

    # Subjects
    for subj in meta:
        
        print("-----------------------------------")
        print("SUBJECT: %s" % subj)
        print("-----------------------------------")  
        
        # Trials
        for trial in meta[subj]["trials"]:
            
            print("%s" % trial)
           
            # Load the data
            trialname = meta[subj]["trials"][trial]["name"]
            trialpath = meta[subj]["trials"][trialname]["path"]
            with open(os.path.join(trialpath, trialname + "_trial_data.pkl"), "rb") as fid:
                datakey = pk.load(fid)   
                
            # Calculate simple margin of stability
            stable = margin_of_stability(user, datakey, treadmill_speed)
            
            # Calculate artificially perturbed margin of stability
            if artifperturb:
                perturb = perturbed_margin_of_stability(user, datakey, perturbation, asteps, treadmill_speed)
                
            # Visualise
            plot_margin_of_stability(datakey, stable)
             
    return None    



'''
perturbed_margin_of_stability(user, datakey, pertubation, stepsize, treadmill_speed):
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
margin_of_stability(user, datakey, treadmill_speed):
    Calculate 2D margin of stability (b) using the geometric method described
    by Hof et al. (2005). Use Shapely to do the vector geometry. For treadmill
    trials, add the treadmill speed to the anteroposterior centre of mass 
    velocity (default = 0.0 m/s).
'''
def margin_of_stability(user, datakey, treadmill_speed = 0.0):
    
    # Centre of mass position (OpenSim global)
    comdata = datakey.data_osim_results[user.bkcode]["pos"]
    comxidx = comdata["headers"].index("center_of_mass_X")
    CoM_r = comdata["data"][:, comxidx:comxidx + 3]
    
    # Centre of mass velocity (OpenSim global)
    CoM_v = datakey.data_osim_results[user.bkcode]["vel"]["data"][:, comxidx:comxidx + 3].copy()  # Requires copy() or will only create new reference, not a new variable!
    
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
    markers = datakey.data_c3dextract_osim["markers"]
    forces = datakey.data_c3dextract_osim["forces"]
    timevec = datakey.data_osim_results[user.bkcode]["pos"]["data"][:, datakey.data_osim_results[user.bkcode]["pos"]["headers"].index("time")]
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
        MoS["b_abs"].append(BoS["baseshape"][t].exterior.distance(CoM["XCoM_pt"][t]))
        MoS["isstable"].append(BoS["baseshape"][t].contains(CoM["XCoM_pt"][t]))        
    MoS["b"] = [MoS["b_abs"][bidx] if stab else -1 * MoS["b_abs"][bidx] for bidx, stab in enumerate(MoS["isstable"])]
    
    # Margins of stability along coordinates (b_x, b_z)
    MoS["b_x"] = []
    MoS["b_z"] = []
    for t in range(len(BoS["baseshape"])):
        
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
        
    # Pickle it
    trialname = datakey.trial["name"]
    with open(os.path.join(datakey.trial["path"], trialname + "_stability.pkl"), "wb") as f:
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
                  ["RHEE", "RLM", "RMT5", "RMT2", "RMM"], \
                  ["LHEE", "LLM", "LMT5", "LMT2", "LMM"], \
                  ["RHEE", "RLM", "RMT5", "RMT2", "LMT2", "LMT5", "LLM", "LHEE"]]
    
    # Get GRF
    # Note: OpenSim GRF files are set up as FP1: left, FP2: right, even though 
    # the lab setup is FP1: right, FP2: left.
    grftime = forces["time0"]
    grfy_right = forces["data"]["FP2"]["F"][:, 1]
    grfy_left = forces["data"]["FP1"]["F"][:, 1]
        
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
    for m in markers["data"]:
        mkrfun[m] = interp.interp1d(markers["time0"], markers["data"][m][:, [0, 2]], kind = "cubic", fill_value = "extrapolate", axis = 0)
    
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
        baseshape = shapely. MultiPoint(mkrcoords).convex_hull
        bos["baseshape"].append(baseshape)
        
    return bos



'''
plot_margin_of_stability(stable):
    Plot the margin of stability (b): 2D, X and Z.
'''
def plot_margin_of_stability(datakey, stable):
    
    subject = datakey.subject
    trial = datakey.trial["name"]
    
    figpath = os.path.join(datakey.trial["path"], "viz", "timeseries")  
    if not os.path.exists(figpath): os.makedirs(figpath)   
            
    # Global parameters
    moslabels = ["b", "b_x", "b_z"]
    moscolours = ["black", "blue", "red"]
    plotlabels = ["2D", "X", "Z"]
    plotlimits = [(-1, 1), (-1, 1), (-0.15, 0.15)] 
    
    # Plot
    figlabels = ["b", "bx", "bz"]
    for p in range(3):
        plt.figure()
        plt.plot(stable["MoS"][moslabels[p]], label = plotlabels[p], color = moscolours[p])
        plt.title("MoS: " + plotlabels[p])
        plt.xlabel(xlabel = "Time step")
        plt.ylabel("MoS (m)")
        plt.ylim(plotlimits[p])

        plt.savefig(os.path.join(figpath, subject + "_" + trial + "_" + figlabels[p] + ".png"))    
        
    return None
    
    

'''
visualise_stability_timehistory(datakey, stable):
    Visualise the base of support (BoS), centre of mass (CoM) and extrapolated 
    centre of mass (XCoM) at each time step. Export as individual frames.
'''
def visualise_stability_timehistory(datakey, stable):

    # Output folder
    figpath = os.path.join(datakey.trial["path"], "viz", "figures", "stable")  
    if not os.path.exists(figpath): os.makedirs(figpath)

    # Plot time history of base of support, centre of mass, and extrapolated
    # centre of mass, save to file
    for t in range(len(stable["BoS"]["baseshape"])):
        
        fig = plt.figure()
        fig.suptitle(datakey.trial["name"] + " Step: " + str(t))
        
        # Plot
        baseshape = stable["BoS"]["baseshape"][t]
        com = stable["CoM"]["CoM_pt"][t]
        xcom = stable["CoM"]["XCoM_pt"][t]
        plot_bos_com_xcom(fig, baseshape, com, xcom)     
        plt.show()
        
        # Save
        fig.savefig(os.path.join(figpath, datakey.trial["name"] + "_" + str(t) + ".png"))

    return None
    


'''
animate_stability_timehistory(datakey, stable):
    Visualise the base of support (BoS), centre of mass (CoM) and extrapolated 
    centre of mass (XCoM) at each time step. Export as MP4 animation.
'''    
def animate_stability_timehistory(datakey, stable):
          
    fig = plt.figure()
    ax = fig.gca()
    
    # Initialise animation
    def init_func():
        ax.clear()
        ax.set_xlim(-0.5, 0.5)
        ax.set_ylim(-0.5, 1.0)
        ax.set_aspect("equal", adjustable="box")
        
    # Animator function
    def animate(i):
    
        ax.clear()
        ax.set_xlim(-0.5, 0.5)
        ax.set_ylim(-0.5, 1.0)
        ax.set_xlabel("Mediolateral (m)")
        ax.set_ylabel("Fore-aft (m)")
        ax.set_aspect("equal", adjustable="box")
        ax.set_title(datakey.trial["name"] + " Step: " + str(i))
        
        # Base shape
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
    figpath = os.path.join(datakey.trial["path"], "viz", "animation")  
    if not os.path.exists(figpath): os.makedirs(figpath) 
    anim.save(os.path.join(figpath, datakey.trial["name"] + ".mp4"))
    
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
plot_bos_com_xcom(fig, baseshape, com, xcom, colours):
    Plot instantaneous base of support (BoS), centre of mass (CoM) and 
    extrapolated centre of mass (XCoM). Coordinates flipped (z-x) for easier 
    visualisation.
'''
def plot_bos_com_xcom(fig, baseshape, com, xcom, colours = ("blue", "red", "black")):
    
    #fig = plt.figure()
    
    # Base shape
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
    
    plt.xlim([-0.75, 0.75])
    plt.ylim([-0.5, 1.0])
    plt.xlabel("Mediolateral (m)")
    plt.ylabel("Fore-aft (m)")
    
    ax = plt.gca()
    ax.set_aspect("equal", adjustable="box")
    
    return fig
    
    
'''
resample1d(data, nsamp):
    Simple resampling by 1-D interpolation (rows = samples, cols = variable).
    Data can be a 1-D or multiple variables in a 2D array-like object.
'''
def resample1d(data, nsamp):

    # Convert list to array
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