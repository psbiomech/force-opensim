# -*- coding: utf-8 -*-
"""
MOCO: WALKING 2D EXAMPLE - STEP ASYMMETRY

@author: Prasanna Sritharan, June 2022


This is a Python implementation step asymmetry functions originally written 
by Russell T. Johnson for Matlab. 

See Moco example: example2DWalkingStepAsymmetry.m
"""


import opensim as osim
import numpy as np


"""
computeStepAsymmetry():
    
Calculate the values of the step length and step time asymmetry from the
results of the simulation.
"""
def computeStepAsymmetry(model_file, threshold, solution_file, grf_file):
    
    # Load model
    model = osim.Model(model_file)
    
    # Load predicted GRF
    grf = np.genfromtxt(grf_file, skip_header=5, delimiter="\t")
    
    # GRF time vector and vertical components
    tvec = grf[:, 0]


    # **********************************
    # STEP TIME ASYMMETRY
    
    # Find index of heel strike on each leg
    hs_idxR = findHeelStrikeIndex(grf[:, 2], threshold)
    hs_idxL = findHeelStrikeIndex(grf[:, 8], threshold)
   
    # Compute step time on each leg
    hs_timeR = tvec[hs_idxR]
    hs_timeL = tvec[hs_idxL]
    if hs_timeR < hs_timeL:
        step_timeL = hs_timeL - hs_timeR
        step_timeR = tvec[-1] - step_timeL
    else:
        step_timeR = hs_timeR - hs_timeL
        step_timeL = tvec[-1] - step_timeR  
    
    # Calculate step time asymmetry (%)
    step_time_asym = 100 * (step_timeR - step_timeL) \
                                    / (step_timeR + step_timeL)
    
    
    # **********************************
    # STEP LENGTH ASYMMETRY
    
    # Get the states for each limb at the instant of heel strike on that limb
    states_traj = osim.StatesTrajectory().createFromStatesTable(model, 
                        osim.TimeSeriesTable(solution_file), False, True, True)
    statesR = states_traj.get(hs_idxR)
    statesL = states_traj.get(hs_idxL)
    
    # Calculate the step length
    step_lengthR = calculateStepLength(model, statesR)
    step_lengthL = calculateStepLength(model, statesL)
    
    # Calculate step length asymmetry (%)
    step_length_asym = 100 * (step_lengthR - step_lengthL) \
                                    / (step_lengthR + step_lengthL)
    
    
    return step_time_asym, step_length_asym



"""
findHeelStrikeIndex():
    
Find heel strike index by determining foot contact on-off instances. If no
heel strike is found, then assume first index is heel strike. This 
implementation differs from that of Russell T. Johnson's Matlab version, but
follows the same prinicples.
"""
def findHeelStrikeIndex(grfy, threshold):
    
    # Find windows representing ground contact
    is_contact = (grfy > threshold).astype(int)

    # Calculate transition points, i.e. heel strike (1) and foot off (-1)
    contact_diff = np.diff(np.insert(is_contact, 0, 1))
    
    # Extract heel strike and foot off indices. If no heel strike found, 
    # assume first index is heel strike.
    idxs = np.where(contact_diff == 1)[0]
    if idxs.size == 0:
        idx = 0
    else:
        idx = idxs[0]
        
    return int(idx)
    
    

"""
calculateStepLength():
    
Find step length by configuring the model at heel strike, then compute distance
between contact spheres along the fore-aft coordinate.
"""    
def calculateStepLength(model, state):
    
    # Configure the model at heel strike
    model.initSystem()
    model.realizePosition(state)
    
    # Get the heel contact spheres
    contact_r = model.getContactGeometrySet().get("heel_r")
    contact_l = model.getContactGeometrySet().get("heel_l")
    
    # Find the positions of the contact spheres in the global frame
    pos_r = contact_r.getFrame().getPositionInGround(state)
    pos_l = contact_l.getFrame().getPositionInGround(state)
    
    # Step length is the difference between global position of the left and
    # right along the fore-aft coordinate (x)
    step_length = abs(pos_r.get(0) - pos_l.get(0))
    
    return step_length
    
