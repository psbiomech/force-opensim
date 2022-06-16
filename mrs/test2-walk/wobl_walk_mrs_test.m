%% FORCE STEP-DOWN PIVOT: MRS TEST

% clear variables and command window
clear all; clc; close all;

%% Input information
% Add here the paths of IK, ID , US and EMG data trials you want to work with
Misc.IKfile  = {fullfile(pwd,'PFOA1_COND1_TRIAL6_ik.mot')};
Misc.IDfile  = {fullfile(pwd,'PFOA1_COND1_TRIAL6_id.sto')};
Misc.model_path   = {fullfile(pwd,'PFOA1.osim')};
Misc.OutPath      = fullfile(pwd,'test-out');                    % folder to store results
time = [0.05 1.01]; 

%% Settings

% name of the results file
Misc.OutName ='PFOA1_COND1_TRIAL2_';

% select degrees of freedom
Misc.DofNames_Input={'ankle_angle_r','knee_angle_r','hip_flexion_r','hip_adduction_r','hip_rotation_r'};    % select the DOFs you want to include in the optimization

% Plotter Bool: Boolean to select if you want to plot lots of output information of intermediate steps in the script
Misc.PlotBool = 1;
% MRS Bool: Select if you want to run the generic muscle redundancy solver
Misc.MRSBool = 1;
% Validation Bool: Select if you want to run the muscle redundancy solver with the optimized parameters
Misc.ValidationBool = 0; 	% TO DO: we should report results of EMG driven simulation as well

%% Run muscle tendon estimator:
[Results,DatStore,Misc] = solveMuscleRedundancy(time,Misc);

