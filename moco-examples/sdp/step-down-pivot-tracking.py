# -*- coding: utf-8 -*-
"""
Created on Fri Jun 17 06:22:23 2022

@author: Owner
"""

import os
import opensim as osim


# source files
srcdir = r"C:\Users\Owner\Documents\data\FORCe\outputdatabase\FAILTCRT01\FAILTCRT01_SDP01"
ikfile = r"ik\FAILTCRT01_SDP01_ik.mot"
extfile = "FAILTCRT01_SDP01_ExternalLoads.xml"
modelfile = "FAILTCRT01.osim"


# Create and name an instance of the MocoTrack tool.
track = osim.MocoTrack()
track.setName("sdp_tracking")

# Construct a ModelProcessor and set it on the tool. The default
# muscles in the model are replaced with optimization-friendly
# DeGrooteFregly2016Muscles, and adjustments are made to the default muscle
# parameters.
modelProcessor = osim.ModelProcessor(os.path.join(srcdir, modelfile))
modelProcessor.append(osim.ModOpAddExternalLoads(os.path.join(srcdir, extfile)))
modelProcessor.append(osim.ModOpIgnoreTendonCompliance())
modelProcessor.append(osim.ModOpReplaceMusclesWithDeGrooteFregly2016())
# Only valid for DeGrooteFregly2016Muscles.
modelProcessor.append(osim.ModOpIgnorePassiveFiberForcesDGF())
# Only valid for DeGrooteFregly2016Muscles.
modelProcessor.append(osim.ModOpScaleActiveFiberForceCurveWidthDGF(1.5))
track.setModel(modelProcessor)

# Construct a TableProcessor of the coordinate data and pass it to the 
# tracking tool. TableProcessors can be used in the same way as
# ModelProcessors by appending TableOperators to modify the base table.
# A TableProcessor with no operators, as we have here, simply returns the
# base table.
track.setStatesReference(osim.TableProcessor(os.path.join(srcdir, ikfile)))
track.set_states_global_tracking_weight(10)

# This setting allows extra data columns contained in the states
# reference that don't correspond to model coordinates.
track.set_allow_unused_references(True)

# Since there is only coordinate position data in the states references,
# this setting is enabled to fill in the missing coordinate speed data using
# the derivative of splined position data.
track.set_track_reference_position_derivatives(True)

# Initial time, final time, and mesh interval.
track.set_initial_time(1.05)
track.set_final_time(2.93)
track.set_mesh_interval(0.08)

# Instead of calling solve(), call initialize() to receive a pre-configured
# MocoStudy object based on the settings above. Use this to customize the
# problem beyond the MocoTrack interface.
study = track.initialize()

# Get a reference to the MocoControlCost that is added to every MocoTrack
# problem by default.
problem = study.updProblem()
effort = osim.MocoControlGoal.safeDownCast(problem.updGoal("control_effort"))

# Put a large weight on the pelvis CoordinateActuators, which act as the
# residual, or 'hand-of-god', forces which we would like to keep as small
# as possible.
model = modelProcessor.process()
model.initSystem()
forceSet = model.getForceSet()
for i in range(forceSet.getSize()):
    forcePath = forceSet.get(i).getAbsolutePathString()
    if 'pelvis' in str(forcePath):
        effort.setWeightForControl(forcePath, 10)

# Solve and visualize.
solution = study.solve()
study.visualize(solution)