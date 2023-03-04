function aclr = collateDataMAT(aclr,fsource,qsource,rra_string)


%COLLATEDATA Trim, resample and save all data to MAT file
%   Prasanna Sritharan, March 2020
%
% Parameters:
%   fsource: string indicating muscle forces source data: 'so','cmc'
%       (default='so')
%   qsource: string indicating joint angle source data: 'ik','rra'
%       (default='so')
%   rra_string: file name suffix for RRA model name and kinematics file,
%       must be provided if use_rra_model=true otherwise empty string is
%       assumed, e.g. if the model name after RRA is
%       'SUBJECT_RRA_Model.osim, then rra_string='_RRA' (optional,
%       default='')

% check parameters
if nargin==3
    rra_string='';
end


% user settings
addpath('..');
user = getUserScriptSettings();
outpath = user.OUTPATH;
gravity = user.GRAVITY;
groups = user.GROUPS;
numsamples = user.RESAMP;
dynamicprefix = user.DYNAMICPREFIX;
ikfolder = user.opensim.ik.CODE;
idfolder = user.opensim.id.CODE;
sofolder = user.opensim.so.CODE;
rrafolder = user.opensim.rra.CODE;
cmcfolder = user.opensim.cmc.CODE;
ikhead = user.opensim.ik.NUMHEADERLINES;
idhead = user.opensim.id.NUMHEADERLINES;
sohead = user.opensim.so.NUMHEADERLINES;
rrahead = user.opensim.rra.NUMHEADERLINES;
cmchead = user.opensim.cmc.NUMHEADERLINES;
ikparams = user.feature.gait2392.ik;
idparams = user.feature.gait2392.id;
soparams = user.feature.gait2392.so;
rraparams = user.feature.gait2392.rra;
cmcparams = user.feature.gait2392.cmc;


% groups
fprintf('Trimming, flipping, resampling and collating data into MAT files...\n');
fprintf('------------------------------------------------\n'); 
n = 1;
for g=1:length(groups)

    % subjects
    subjs = fieldnames(aclr.(groups{g}));
    for s=1:length(subjs)

        % subject parameters
        mass = aclr.(groups{g}).(subjs{s}).mass;
        height = aclr.(groups{g}).(subjs{s}).height;
        
        % normalisation factors
        bw = mass*gravity;
        bwht = (bw*height)/100;
        
        % trials
        trials = fieldnames(aclr.(groups{g}).(subjs{s}));
        for t=1:length(trials)        
            
            % check if struct field is really a trial
            if ~startsWith(trials{t},dynamicprefix), continue; end         

            fprintf('%s %s\t%s',n,groups{g},subjs{s},trials{t});              
            
            % trial folder location and file name prefix
            dest = aclr.(groups{g}).(subjs{s}).(trials{t}).fullpath; 
            fprefix = aclr.(groups{g}).(subjs{s}).(trials{t}).fileprefix;
                       
            % time window for landing phase
            twindow = aclr.(groups{g}).(subjs{s}).(trials{t}).timevec.landing;

            % trial leg
            leg = aclr.(groups{g}).(subjs{s}).(trials{t}).leg;            
            
            % create ignore flag
            aclr.(groups{g}).(subjs{s}).(trials{t}).ignore = false;                

            try
 
                n = n + 1;
                fprintf(' ---> Row ID: %d\n',n);
                                
                % ********************
                % OpenSim data
               
                % muscle forces
                % note: struct name is 'so' even if cmc is used for the
                % forces, need to change this to a generic name (TBD)
                if strcmpi(fsource,'so')
                    sopath0 = fullfile(dest,sofolder,[fprefix '_StaticOptimization_force.sto']);
                    sohead0 = sohead;
                    soparams0 = soparams;
                elseif strcmpi(fsource,'cmc')
                    sopath0 = fullfile(dest,cmcfolder,'RRA_Actuation_force.sto');
                    sohead0 = cmchead;
                    soparams0 = cmcparams;
                end   
                soopts = detectImportOptions(sopath0,'FileType','text','VariableNamesLine',sohead0);
                sodata = readmatrix(sopath0,'FileType','text','NumHeaderLines',sohead0);
                sodata = resampleData(trimData(sodata,twindow(1),twindow(2),1),numsamples);
                so.data = sodata;
                so.normalised = [sodata(:,1) sodata(:,2:end)/bw];
                so.pca = so.normalised;
                so.pca = cell2mat(cellfun(@(g)sum(so.pca(:,g),2),soparams0.mergecols{leg},'UniformOutput',false));
                so.headers = soopts.VariableNames;
                so.pcaheaders = soparams0.headers;                
                
                % joint angles
                % note: struct name is 'ik' even if rra is used for the
                % angles, need to change this to a generic name (TBD)                
                if strcmpi(qsource,'ik')
                    ikpath0 = fullfile(dest,ikfolder,[fprefix '_ik.mot']);
                    ikhead0 = ikhead;
                    ikparams0 = ikparams;
                elseif strcmpi(qsource,'rra')
                    ikpath0 = fullfile(dest,rrafolder,[fprefix rra_string '_Kinematics_q.sto']);
                    ikhead0 = rrahead;
                    ikparams0 = rraparams;
                end
                ikopts = detectImportOptions(ikpath0,'FileType','text','VariableNamesLine',ikhead0);
                ikdata = readmatrix(ikpath0,'FileType','text','NumHeaderLines',ikhead0);
                ikdata = resampleData(trimData(ikdata,twindow(1),twindow(2),1),numsamples);
                ik.data = ikdata;
                ik.normalised = ikdata;   
                ik.pca = ik.normalised;
                if (leg==2), ik.pca(:,ikparams0.flipcols) = -ik.pca(:,ikparams0.flipcols); end
                ik.pca = ik.pca(:,ikparams0.pcacols{leg});
                ik.headers = ikopts.VariableNames;  
                ik.pcaheaders = ikparams0.headers;

                % joint moments
                idpath = fullfile(dest,idfolder,[fprefix '_id.sto']);
                idopts = detectImportOptions(idpath,'FileType','text','VariableNamesLine',idhead);
                iddata = readmatrix(idpath,'FileType','text','NumHeaderLines',idhead);
                iddata = resampleData(trimData(iddata,twindow(1),twindow(2),1),numsamples);
                id.data = iddata;
                id.normalised = [iddata(:,1) iddata(:,2:end)/bwht];
                id.pca = id.normalised;
                if (leg==2), id.pca(:,idparams.flipcols) = -id.pca(:,idparams.flipcols); end
                id.pca = id.pca(:,idparams.pcacols{leg});   
                id.headers = idopts.VariableNames;
                id.pcaheaders = idparams.headers;

                
                % ********************
                % Subject and trial data            

                % get trial info
                trial = aclr.(groups{g}).(subjs{s}).(trials{t});

                % get subj info, delete trial info
                subj = aclr.(groups{g}).(subjs{s});
                subjfields = fieldnames(subj);
                for j=1:length(subjfields)
                    if startsWith(subjfields(j),dynamicprefix)
                        subj = rmfield(subj,subjfields{j});
                    end
                end


                % save data as MAT file  
                save(fullfile(dest,[fprefix '.mat']),'subj','trial','ik','id','so');                                      
                                
            catch
                fprintf('---> ERROR. Ignoring trial.\n');
                aclr.(groups{g}).(subjs{s}).(trials{t}).ignore = true;                
                fid = fopen('aclr_collate_ignore_trials.txt','a');
                fprintf(fid,'%s %s %s %s\n',aclr.(groups{g}).(subjs{s}).(trials{t}).id,groups{g},subjs{s},trials{t});
                fclose(fid);
                n = n - 1; % decrease counter
                
            end
            
        end

    end
    
end

% save metastruct
save([outpath '\aclr.mat'],'aclr');

fprintf('------------------------------------------------\n');

    
end



