# -*- coding: utf-8 -*-
"""
MOCO: WALKING 2D TRACKING EXAMPLE

@author: Prasanna Sritharan, June 2022

This is a Python implementation of an example optimal control
problem (2-D walking) orginally created in C++ by Antoine Falisse
(see: example2DWalking.cpp). This Python script is itself adapted from a
Matlab implementation by Brian Umberger.
"""

from math import pi
import opensim as osim



# %% DEFINE THE OPTIMAL CONTROL PROBLEM


# Create tracking problem
track = osim.MocoTrack()
track.setName("gait_tracking")

# Construct a ModelProcessor and set it on the tool
modelProcessor = osim.ModelProcessor("2D_gait.osim")
track.setModel(modelProcessor)


# Get the coordinates into a TableProcessor
tableProcessor = osim.TableProcessor("referenceCoordinates.sto")
tableProcessor.append(osim.TabOpLowPassFilter(6.0))

# Set the coordinates as the reference states
track.setStatesReference(tableProcessor)
track.set_allow_unused_references(True)

# Only coordinates are provided so derive to get speeds and track these too
track.set_track_reference_position_derivatives(True)

# use coordinates and speeds for initial guess
track.set_apply_tracked_states_to_guess(True)

# Define the global state tracking weight
track.set_states_global_tracking_weight(1.0)


# Set initial time, final time and mesh interval
track.set_initial_time(0.0)
track.set_final_time(0.47008941)

# Call initialize() to receive a pre-configured MocoStudy object based on the 
# settings above. Use this to customize the problem beyond the MocoTrack
# interface. 
study = track.initialize()

# Get the MocoProblem from the MocoStudy to perform more customisation, this is
# a writable reference
problem = study.updProblem()



# %% GOALS


# Symmetry: This goal allows us to simulate only one step with left-right
# symmetry that we can then double to create a full gait cycle.
symmetryGoal = osim.MocoPeriodicityGoal("symmetryGoal")
problem.addGoal(symmetryGoal)


# Enforce periodic symmetry
model = modelProcessor.process()
model.initSystem()
state_names = [model.getStateVariableNames().getitem(sn) for sn in range(model.getNumStateVariables())]
for sn in state_names:
    
    # Symmetric coordinate values and speeds (except for pelvis_tx value):
    # Here we constrain final coordinate values of one leg to match the initial
    # value of the other leg.
    if sn.startswith("/jointset"):
        if "_r" in sn:
            sn1 = sn.replace("_r", "_l")
        elif "_l" in sn:
            sn1 = sn.replace("_l", "_r")
        elif "pelvis_tx" not in sn:
            sn1 = sn
        symmetryGoal.addStatePair(osim.MocoPeriodicityGoalPair(sn, sn1))

    # Symmetric muscle activations:
    # Here, we constrain final muscle activation values of one leg to match the
    # initial activation values of the other leg.
    if sn.startswith("/activation"):
        if "_r" in sn:
            sn1 = sn.replace("_r", "_l")
        elif "_l" in sn:
            sn1 = sn.replace("_l", "_r")
        pair = osim.MocoPeriodicityGoalPair(sn, sn1)
        symmetryGoal.addStatePair(pair)        

# The pelvis_tx speed has periodic symmetry
symmetryGoal.addStatePair(osim.MocoPeriodicityGoalPair("/jointset/groundPelvis/pelvis_tx/speed"))                

# Lumbar control has periodic symmetry. The other controls don't need 
# symmetry enforces as they are all muscle excitations. Their behaviour will be
# contstrained by the periodicity imposed on their respective activations.
symmetryGoal.addControlPair(osim.MocoPeriodicityGoalPair('/lumbarAct'))


# Get a reference to the MocoControlGoal that is added to a MocoTrack problem
# by default and adjust the weight
effort = osim.MocoControlGoal().safeDownCast(problem.updGoal("control_effort"))
effort.setWeight(10.0)


# Add a contact tracking goal:
# Track the left and right vertical and fore-aft GRFs
grf_tracking_weight = 1
contact_r = ["contactHeel_r", "contactFront_r"]
contact_l = ["contactHeel_l", "contactFront_l"]
if grf_tracking_weight != 0:
    
    # Create a contact tracking goal
    contactTracking = osim.MocoContactTrackingGoal("contact", grf_tracking_weight)
    contactTracking.setExternalLoadsFile("referenceGRF.sto")
    
    # Add contact groups. The sum of all the contact forces in each group
    # should track the force data from a single ExternalForce
    contactTracking.addContactGroup(contact_r, "Right_GRF")
    contactTracking.addContactGroup(contact_l, "Left_GRF")
    
    # 2D walking problem so consider force errors in the sagittal plane only
    contactTracking.setProjection("plane")
    contactTracking.setProjectionVector(osim.Vec3(0, 0, 1))
    
    # add the goal to the MocoProblem
    problem.addGoal(contactTracking)
 
    
# %% BOUNDS
  
    
# Coordinate bounds as dict
coord_bounds = {}
coord_bounds["/jointset/groundPelvis/pelvis_tilt/value"] = [-20*pi/180, -10*pi/180]
coord_bounds["/jointset/groundPelvis/pelvis_tilt/value"] = [0.0, 1.0]
coord_bounds["/jointset/groundPelvis/pelvis_ty/value"] = [0.75, 1.25]
coord_bounds["/jointset/hip_r/hip_flexion_r/value"] = [-10*pi/180, 60*pi/180]
coord_bounds["/jointset/hip_l/hip_flexion_l/value"] = [-10*pi/180, 60*pi/180]
coord_bounds["/jointset/knee_l/knee_angle_r/value"] = [-50*pi/180, 0]
coord_bounds["/jointset/knee_l/knee_angle_l/value"] = [-50*pi/180, 0]
coord_bounds["/jointset/ankle_l/ankle_angle_r/value"] = [-15*pi/180, 25*pi/180]
coord_bounds["/jointset/ankle_l/ankle_angle_l/value"] = [-15*pi/180, 25*pi/180]
coord_bounds["/jointset/lumbar/lumbar/value"] = [0, 20*pi/180]

# Set bounds
for bnd in coord_bounds:
    problem.setStateInfo(bnd, coord_bounds[bnd]);
    

# %% SOLVE

# Solve the problem for a single step
solution = study.solve()

# Create a full stride cycle from the single step solution, see API 
# documentation for use of this utility function
full_stride_cycle_solution = osim.createPeriodicTrajectory(solution)
full_stride_cycle_solution.write("walk_2D_stride_cycle_predicted_solution.sto")

# Also extract the predicted ground forces, see API documentation for use of
# this utility function
contact_forces_table = osim.createExternalLoadsTableForGait(model, full_stride_cycle_solution, contact_r, contact_l)
osim.STOFileAdapter().write(contact_forces_table, "walk_2D_stride_cycle_predicted_ground_forces.sto")





