function [pcadata, pcaweights, pcaout, pcainfo] = performWeightedPCA(refmodel)

%PERFORMWEIGHTEDPCA Undertake PCA analysis for FORCe SDP
%   Prasanna Sritharan, February 2022
%
% Based on PCA scripts by Prasanna Sritharan for ACLR hopping 
% (published AnnBiomedEng 2022)
%
% Procedure:
%   1. Pool all data into 3D matrices
%   2. Calculate weightings
%   3. Apply weighted PCA

% User settings
addpath('..');
user = getUserScriptSettings();
csvpath = user.SRCPATH;
outpath = user.OUTPATH2;
groups = user.GROUPS;
subjprefix = user.SUBJPREFIX;
trialcombo = user.TRIALCOMBO;
modelparams = user.feature2.(refmodel);

fprintf('Perform Weighted PCA on raw data waveforms.\n');
fprintf('------------------------------------------------\n'); 

% Load FORCe SDP data into tables
fprintf('Loading raw data...\n');
limbdata = readtable(fullfile(csvpath, 'force_sdp_results_all_trials_ikid_bothlegs.csv'));

% Data tables
pcadata.ik = [];
pcadata.id = [];
pcainfo = struct;
pcaout = struct;
qs = [];
x = 1;  % total observations
ntrials = [0 0];    % number of observations per group
subjlist = {{}, {}};    % list of subjects in each group


fprintf('Collating data into matrices...\n');

% Groups
for g=1:2

    % Subjects
    allsubjects = unique(limbdata.subject);
    subjects = allsubjects(~cellfun(@isempty, regexp(allsubjects, [subjprefix{g} '\d+'])));
    for s=1:length(subjects)

        % Trials
        q = 0;
        trials = unique(limbdata.trial(strcmpi(limbdata.subject, subjects{s}) & contains(limbdata.trial_combo, trialcombo{g})));
        for t=1:length(trials)

            % Increment trial counter
            q = q + 1;
            ntrials(g) = ntrials(g) + 1;

            % IK data
            ikrows0 = limbdata(strcmpi(limbdata.subject, subjects{s}) & strcmpi(limbdata.trial, trials{t}) & strcmpi(limbdata.analysis, 'ik'), :);
            ikrows = ikrows0(~(strcmpi(ikrows0.data_leg_role, 'nonpivot') & contains(ikrows0.variable, 'lumbar')), :);   % remove non-pivot lumbar variables
            ikdata0 = ikrows{:, 30:130}';
            ikdata = ikdata0(:, modelparams.ik.idx);
            pcadata.ik(x+q-1, :, :) = ikdata;

            % ID data
            idrows0 = limbdata(strcmpi(limbdata.subject, subjects{s}) & strcmpi(limbdata.trial, trials{t}) & strcmpi(limbdata.analysis, 'id'), :);
            idrows = idrows0(~(strcmpi(idrows0.data_leg_role, 'nonpivot') & contains(idrows0.variable, 'lumbar')), :);   % remove non-pivot lumbar variables
            iddata0 = idrows{:, 30:130}';
            iddata = iddata0(:, modelparams.id.idx);
            pcadata.id(x+q-1, :, :) = iddata;   
                
            % Add subject names to list, append limb suffix
            subjlist{g} = [subjlist{g}, subjects{s}];
                
        end

        % For each row of data, list number of trials per subject
        % (i.e. each row is 1 trial; if a subject has q valid trials in
        % total, then for each row associated with that subject, record q)
        if (q > 0)
            qs(x:x+q-1) = q;
            x = x + q;
        end            

    end

end

% Store database info
pcainfo.observations.total = x-1;
pcainfo.observations.(groups{1}) = ntrials(1);
pcainfo.observations.(groups{2}) = ntrials(2);
pcainfo.variables = size(pcadata.ik, 2);
pcainfo.subjects.all.(groups{1}) = subjlist{1};
pcainfo.subjects.all.(groups{2}) = subjlist{2};
pcainfo.subjects.unique.(groups{1}) = unique(subjlist{1});
pcainfo.subjects.unique.(groups{2}) = unique(subjlist{2});
for d={'ik','id'}
    pcainfo.(d{1}).label = modelparams.(d{1}).label;
    pcainfo.(d{1}).varnames = modelparams.(d{1}).headers;
end
pcainfo.(['is' groups{1}]) = [true(pcainfo.observations.(groups{1}), 1); false(pcainfo.observations.(groups{2}), 1)];


% Calculate weighting vector
fprintf('\nBuilding weighting vector...\n');
Tis = 1./qs';
pcaweights = Tis/sum(Tis);

% Undertake weighted PCA
fprintf('Performing weighted PCA...\n');
pcaout.ik = [];
pcaout.id = [];
for d={'ik','id'}
    for c=1:size(pcadata.(d{1}), 3)
        [coeff, score, ~, ~, explained] = pca(squeeze(pcadata.(d{1})(:,:,c)), 'Weights', pcaweights);
        pcaout.(d{1}).coeff(:, :, c) = coeff;
        pcaout.(d{1}).score(:, :, c) = score;
        pcaout.(d{1}).explained(:, c) = explained;
    end
end

% Save PCA inputs and outputs
fprintf('\nSaving PCA inputs and outputs...\n');
if ~exist(outpath,'dir'), mkdir(outpath); end
save([outpath 'pca.mat'],'pcadata','pcaweights','pcaout','pcainfo');

fprintf('------------------------------------------------\n');

    
end


