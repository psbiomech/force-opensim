 function user = getUserScriptSettings()

%GETUSERSCRIPTSETTINGS User settings for processing
%   Prasanna Sritharan, Februaru 2019



%% FOLDER PATHS


% folder paths: Lenovo laptop
user.CODEROOT = 'C:\Users\Owner\Documents\data\FORCe\'; 
user.OUTPATH1 = 'C:\Users\Owner\Documents\data\FORCe\pcaDatabase\by-leg\';   % location of output database
user.OUTPATH2 = 'C:\Users\Owner\Documents\data\FORCe\pcaDatabase\both-legs\';   % location of output database
user.SRCPATH = 'C:\Users\Owner\Documents\data\FORCe\outputDatabase\csvfolder';   % location of source data

% additional folder paths: three-group analysis
user.GROUP3OUTPATH = fullfile(user.OUTPATH2, 'three-group\');


%% GENERAL
% ------------------------------

% general parameters
user.LIMBS = {'pivot', 'nonpivot'};
user.GROUPS = {'sym','ctrl'};
user.SUBJPREFIX = {'FAILT', 'FAILTCRT'};
user.TRIALCOMBO = {{'pivot_more', 'pivot'}, {'pivot_less', 'pivot'}};
user.FOOT = {'r','l'};
user.RESAMP = 101;
user.GRAVITY = 9.81;    % m/s2
user.DATACOLS = 32:132;

% general parameters: three-group analysis
user.GROUPS3 = {'moresym', 'lesssym', 'ctrl'};
user.SUBJPREFIX3 = {'FAILT', 'FAILT', 'FAILTCRT'};
user.TRIALCOMBO3 = {'pivot_more', 'pivot_less', 'pivot'};

% user.STATICPREFIX = 'STATIC';  % trial code prefix for static trials (eg. for trial code WALK01, the trialprefix is 'WALK')
% user.DYNAMICPREFIX = 'HOPP';  % trial code prefix for dynamic trials (eg. for trial code WALK01, the trialprefix is 'WALK')
% user.SEPARATOR = '_';   % file name separator between subject code and trial code (eg. for file FAILT01_WALK01.c3d, the separator is '_')
% user.TRIALTYPE = {'static','dynamic'};  % trial type: static calibration trial, or dynamic gait trial


%% FEATURE SELECTION PARAMETERS: BY LEGS

% Rajagopal model
user.feature.rajagopal.ik.label = 'angle';
user.feature.rajagopal.ik.headers = {'hip_flexion','hip_adduction','hip_rotation','knee_angle','ankle_angle','lumbar_extension','lumbar_bending','lumbar_rotation'};
user.feature.rajagopal.ik.idx = [8:10 11 13 16:18];

user.feature.rajagopal.id.label = 'moment';
user.feature.rajagopal.id.headers = {'hip_flexion','hip_adduction','hip_rotation_moment','knee_angle_moment','ankle_angle_moment','lumbar_extension_moment','lumbar_bending_moment','lumbar_rotation_moment'};
user.feature.rajagopal.id.idx = [8:10 14 19 11:13];

% user.feature.gait2392.so.pcacols = [];
% user.feature.gait2392.so.flipcols = [];  % column indices *before* trimming to only PCA cols
% user.feature.gait2392.so.mergecols = {{21:23,2:4,9:11,29,30:32,33:34,35},{64:66,45:47,51:54,72,73:75,76:77,78}};    % column indices *before* trimming to only PCA cols
% user.feature.gait2392.so.label = 'force';
% user.feature.gait2392.so.headers = {'gmax','gmed','hams','rf','vas','gas','sol'};
% 
% user.feature.gait2392.rra.pcacols = {[8:14 2:4 26:28],[17:23 2:4 26:28]};
% user.feature.gait2392.rra.flipcols = [3 4 7 27 28];  % column indices *before* trimming to only PCA cols
% user.feature.gait2392.rra.mergecols = [];       % column indices *before* trimming to only PCA cols
% user.feature.gait2392.rra.label = 'angle';
% user.feature.gait2392.rra.headers = {'hip_flex','hip_add','hip_rot','knee_angle','knee_rot','knee_add','ankle_angle','pelvis_tilt','pelvis_list','pelvis_rot','lumbar_ext','lumbar_bend','lumbar_rot'};
% 
% user.feature.gait2392.cmc.pcacols = [];
% user.feature.gait2392.cmc.flipcols = [];  % column indices *before* trimming to only PCA cols
% user.feature.gait2392.cmc.mergecols = {{21:23,2:4,9:11,29,30:32,33:34,35},{64:66,45:47,51:54,72,73:75,76:77,78}};    % column indices *before* trimming to only PCA cols
% user.feature.gait2392.cmc.label = 'force';
% user.feature.gait2392.cmc.headers = {'gmax','gmed','hams','rf','vas','gas','sol'};


%% FEATURE SELECTION PARAMETERS: BOTH LEGS

% Rajagopal model
user.feature2.rajagopal.ik.label = 'angle';
user.feature2.rajagopal.ik.headers = {'hip_adduction','hip_flexion','hip_rotation','knee_angle','ankle_angle','lumbar_bending','lumbar_extension','lumbar_rotation', 'hip_adduction_np','hip_flexion_np','hip_rotation_np','knee_angle_np','ankle_angle_np'};
user.feature2.rajagopal.ik.idx = [8:10 11 13 16:18 33:35 36 38];

user.feature2.rajagopal.id.label = 'moment';
user.feature2.rajagopal.id.headers = {'hip_adduction','hip_flexion','hip_rotation','knee_angle','ankle_angle','lumbar_bending','lumbar_extension','lumbar_rotation', 'hip_adduction_np','hip_flexion_np','hip_rotation_np','knee_angle_np','ankle_angle_np'};
user.feature2.rajagopal.id.idx = [8:10 14 19 11:13 33:35 39 44];

% user.feature.gait2392.so.pcacols = [];
% user.feature.gait2392.so.flipcols = [];  % column indices *before* trimming to only PCA cols
% user.feature.gait2392.so.mergecols = {{21:23,2:4,9:11,29,30:32,33:34,35},{64:66,45:47,51:54,72,73:75,76:77,78}};    % column indices *before* trimming to only PCA cols
% user.feature.gait2392.so.label = 'force';
% user.feature.gait2392.so.headers = {'gmax','gmed','hams','rf','vas','gas','sol'};
% 
% user.feature.gait2392.rra.pcacols = {[8:14 2:4 26:28],[17:23 2:4 26:28]};
% user.feature.gait2392.rra.flipcols = [3 4 7 27 28];  % column indices *before* trimming to only PCA cols
% user.feature.gait2392.rra.mergecols = [];       % column indices *before* trimming to only PCA cols
% user.feature.gait2392.rra.label = 'angle';
% user.feature.gait2392.rra.headers = {'hip_flex','hip_add','hip_rot','knee_angle','knee_rot','knee_add','ankle_angle','pelvis_tilt','pelvis_list','pelvis_rot','lumbar_ext','lumbar_bend','lumbar_rot'};
% 
% user.feature.gait2392.cmc.pcacols = [];
% user.feature.gait2392.cmc.flipcols = [];  % column indices *before* trimming to only PCA cols
% user.feature.gait2392.cmc.mergecols = {{21:23,2:4,9:11,29,30:32,33:34,35},{64:66,45:47,51:54,72,73:75,76:77,78}};    % column indices *before* trimming to only PCA cols
% user.feature.gait2392.cmc.label = 'force';
% user.feature.gait2392.cmc.headers = {'gmax','gmed','hams','rf','vas','gas','sol'};




end
