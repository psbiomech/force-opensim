function [pcadata,pcaweights,pcaout,pcainfo] = performWeightedPCA(refmodel)

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
outpath = user.OUTPATH;
limbs = user.LIMBS;
groups = user.GROUPS;
subjprefix = user.SUBJPREFIX;
modelparams = user.feature.(refmodel);


fprintf('------------------------------------------------\n'); 

% Load FORCe SDP data into tables
fprintf('Loading data into tables:\n');
fprintf('---> Pivot leg\n');
sdp.pivot = readtable(fullfile(csvpath, 'force_sdp_results_all_trials_ikid_pivot.csv'));
fprintf('---> Non-pivot leg\n');
sdp.nonpivot = readtable(fullfile(csvpath, 'force_sdp_results_all_trials_ikid_nonpivot.csv'));

% Data tables
pcadata = struct;
pcainfo = struct;
pcaweights = struct;
pcaout = struct;

% Limb
for b=1:2
    
    fprintf('\nWeighted PCA: %s LIMB\n', upper(limbs{b}))

    % Initialise output data tables
    pcadata.(limbs{b}).ik = [];
    pcadata.(limbs{b}).id = [];

    % Data from relevant table
    limbdata = sdp.(limbs{b});

    % Groups
    qs = [];
    x = 1;  % total observations
    ntrials = [0 0];    % number of observations per group
    subjlist = {{}, {}};
    for g=1:2

        % Subjects
        allsubjects = unique(sdp.(limbs{b}).subject);
        subjects = allsubjects(~cellfun(@isempty, regexp(allsubjects, [subjprefix{g} '\d+'])));
        for s=1:length(subjects)
    
            % Trials
            q = 0;
            trials = unique(limbdata.trial(strcmpi(limbdata.subject, subjects{s})));
            for t=1:length(trials)

                % Increment trial counter
                q = q + 1;
                ntrials(g) = ntrials(g) + 1;

                % IK data
                ikrows = limbdata(strcmpi(limbdata.subject, subjects{s}) & strcmpi(limbdata.trial, trials{t}) & strcmpi(limbdata.analysis, 'ik'), :);
                ikdata0 = ikrows{:, 20:120}';
                ikdata = ikdata0(:, modelparams.ik.idx);
                pcadata.(limbs{b}).ik(x+q-1, :, :) = ikdata;

                % ID data
                idrows = limbdata(strcmpi(limbdata.subject, subjects{s}) & strcmpi(limbdata.trial, trials{t}) & strcmpi(limbdata.analysis, 'id'), :);
                iddata0 = idrows{:, 20:120}';
                iddata = iddata0(:, modelparams.id.idx);
                pcadata.(limbs{b}).id(x+q-1, :, :) = iddata;   

                % Add subject names to list
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
    pcainfo.(limbs{b}).observations.(groups{1}) = ntrials(1);
    pcainfo.(limbs{b}).observations.(groups{2}) = ntrials(2);
    pcainfo.(limbs{b}).variables = size(pcadata.(limbs{b}).ik, 2);
    pcainfo.(limbs{b}).subjects.all.(groups{1}) = subjlist{1};
    pcainfo.(limbs{b}).subjects.all.(groups{2}) = subjlist{2};
    pcainfo.(limbs{b}).subjects.unique.(groups{1}) = unique(subjlist{1});
    pcainfo.(limbs{b}).subjects.unique.(groups{2}) = unique(subjlist{2});
    for d={'ik','id'}
        pcainfo.(limbs{b}).(d{1}).label = modelparams.(d{1}).label;
        pcainfo.(limbs{b}).(d{1}).varnames = modelparams.(d{1}).headers;
    end
    pcainfo.(['is' groups{1}]) = [true(pcainfo.(limbs{b}).observations.(groups{1}), 1); false(pcainfo.(limbs{b}).observations.(groups{2}), 1)];


    % Calculate weighting vector
    fprintf('---> Building weighting vector...\n');
    Tis = 1./qs';
    pcaweights.(limbs{b}) = Tis/sum(Tis);

    % Undertake weighted PCA
    fprintf('---> Performing weighted PCA...\n');
    pcaout.(limbs{b}).ik = [];
    pcaout.(limbs{b}).id = [];
    for d={'ik','id'}
        for c=1:size(pcadata.(limbs{b}).(d{1}),3)
            [coeff,score,~,~,explained] = pca(squeeze(pcadata.(limbs{b}).(d{1})(:,:,c)), 'Weights', pcaweights.(limbs{b}));
            pcaout.(limbs{b}).(d{1}).coeff(:,:,c) = coeff;
            pcaout.(limbs{b}).(d{1}).score(:,:,c) = score;
            pcaout.(limbs{b}).(d{1}).explained(:,c) = explained;
        end
    end

end

% Save PCA inputs and outputs
fprintf('\nSaving PCA inputs and outputs...\n');
if ~exist(outpath,'dir'), mkdir(outpath); end
save([outpath 'pca.mat'],'pcadata','pcaweights','pcaout','pcainfo');

fprintf('------------------------------------------------\n');

    
end



