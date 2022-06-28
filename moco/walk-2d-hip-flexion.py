# -*- coding: utf-8 -*-
"""
WALKING 2D PREDICTION: MINIMISE HIP FORCE

@author: Prasanna Sritharan
"""


from math import pi
import opensim as osim
import os
import stepasymmetry as stas



# %% DEFINE OPTIMAL CONTROL PROBLEM



# Create MocoStudy
study = osim.MocoStudy()
study.setName("walk_2d_hip force")


# Get the model
model_processor = osim.ModelProcessor("2D_gait.osim")
model_processor.append(osim.ModOpTendonComplianceDynamicsModeDGF("implicit"))

# Get the MocoProblem and set the model on it
problem = study.updProblem()
problem.setModelProcessor(model_processor)



# %% SET GOALS


# Periodicity:
# All states are periodic except pevlix_tx value, lumbar control is periodic
periodicity_goal = osim.MocoPeriodicityGoal("periodicity")
model = model_processor.process()
model.initSystem()
state_names = [model.getStateVariableNames().getitem(sn) for sn in range(model.getNumStateVariables())]
for sn in state_names:
    if "pelvis_tx/value" not in sn:
        periodicity_goal.addStatePair(osim.MocoPeriodicityGoalPair(sn))
periodicity_goal.addStateControlPair(osim.MocoPeriodicityGoalPair("/lumbarAct"))
problem.addGoal(periodicity_goal)

# Average gait speed:
# Arbitrary "normal" walking speed
speed_goal = osim.MocoAverageSpeedGoal("average_speed")
speed_goal.set_desired_average_speed(1.0)
problem.addGoal(speed_goal)

# Hip joint angles:
hip_goal = osim.MocoSumSquaredStateGoal()
hip_goal.setPattern("/jointset/hip_r/hip_flexion_r/value")
problem.addGoal(hip_goal)

# Control effort:
# Minimise overall control effort
effort_goal = osim.MocoControlGoal("effort", 1.0)
effort_goal.setExponent(3)
effort_goal.setDivideByDisplacement(True)
problem.addGoal(effort_goal)



# %% BOUNDS


# Set time bounds
problem.setTimeBounds(0.0, 0.94)

# Coordinate bounds as dict
coord_bounds = {}
coord_bounds["/jointset/groundPelvis/pelvis_tilt/value"] = [-20*pi/180, 20*pi/180]
coord_bounds["/jointset/groundPelvis/pelvis_tx/value"] = [0, 2]
coord_bounds["/jointset/groundPelvis/pelvis_ty/value"] = [0.75, 1.25]
coord_bounds["/jointset/hip_r/hip_flexion_r/value"] = [-10*pi/180, 60*pi/180]
coord_bounds["/jointset/hip_l/hip_flexion_l/value"] = [-10*pi/180, 60*pi/180]
coord_bounds["/jointset/knee_r/knee_angle_r/value"] = [-50*pi/180, 0]
coord_bounds["/jointset/knee_l/knee_angle_l/value"] = [-50*pi/180, 0]
coord_bounds["/jointset/ankle_r/ankle_angle_r/value"] = [-15*pi/180, 25*pi/180]
coord_bounds["/jointset/ankle_l/ankle_angle_l/value"] = [-15*pi/180, 25*pi/180]
coord_bounds["/jointset/lumbar/lumbar/value"] = [0, 20*pi/180]

# Set coordinate bounds
for bnd in coord_bounds:
    problem.setStateInfo(bnd, coord_bounds[bnd])




# %% SOLVE


# Configure the solver
solver = study.initCasADiSolver()
solver.set_num_mesh_intervals(100)
solver.set_verbosity(2)
solver.set_optim_convergence_tolerance(1e-4)
solver.set_optim_constraint_tolerance(1e-4)
solver.set_optim_max_iterations(2000)

# Set the initial guess to be the symmetric two-steps tracking solution
# from example2DWalking.py. Run this first, or proceed without a guess.
guess_file = "walk_2D_two_steps_tracking_solution.sto"
if os.path.isfile(guess_file):
    solver.setGuessFile(guess_file)
        
# Print the Study to file
study.printToXML("walking_2d_hip_force.omoco")

# Solve
solution = study.solve()
solution.write("walking_2d_hip_force_solution.sto")


# Write predicted GRF to file
contact_r = ["contactHeel_r", "contactFront_r"]
contact_l = ["contactHeel_l", "contactFront_l"]
contact_forces_table = osim.createExternalLoadsTableForGait(model, solution, contact_r, contact_l)
osim.STOFileAdapter().write(contact_forces_table, "walking_2d_hip_force_ground_forces.sto")

# Calculate step time asymmetry
step_time_asym, step_length_asym = stas.computeStepAsymmetry(model, 25, "walk_2d_hip_force_solution.sto", "walk_2d_hip_force_ground_forces.sto")
print("\nStep time asymmetry: \t%f\n" % step_time_asym)
print("\nStep length asymmetry:\t%f\n" % step_length_asym)