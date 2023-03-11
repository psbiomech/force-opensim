function [pcsexplained, totalvalidpcs, totalvariance, paselected, paquantl] = performParallelAnalysis(pcadata, pcaweights, pcaout, pcainfo)

%PERFORMPARALLELANALYSIS Undertake Horn's Parallel Analysis for FORCe SDP
%   Prasanna Sritharan, February 2022
%
% Based on PCA scripts by Prasanna Sritharan for ACLR hopping 
% (published AnnBiomedEng 2022)


% User settings
addpath('..');
user = getUserScriptSettings();
outpath = user.OUTPATH2;


fprintf('Horn''s Parallel Analysis on PCA scores.\n');
fprintf('------------------------------------------------\n'); 


% Load the 95th percentiles of the eigenvalues of the normally-distributed
% random variables if it exists, otherwise create one and save to file
eigfile = fullfile(outpath, 'randvareig.mat');
if exist(eigfile,'file')
    fprintf('Random variable file exists. Loading...\n');
    aux = load(eigfile);
    randvareigpcts = aux.randvareigpcts;
else
    fprintf('Random variable file does not exist. Generating random variables...\n');
    obs = pcainfo.observations.total;
    vars = pcainfo.variables;
    [randvareigmeans, randvareigpcts] = generateRandomSet(obs,vars);
    save(eigfile, 'randvareigmeans', 'randvareigpcts');
end

% Data tables
pcsexplained = {};
totalvariance = {};
totalvalidpcs = 0;
paselected = struct;
paquantl = struct;

% Perform Parallel Analysis per variable
fprintf('Performing Parallel Analysis...\n');

% Analysis
inc = 1;
for d={'ik','id'}
    
    % Retain eigenvalues of R, the weighted correlation matrix of data X,
    % that exceed the equivalent eigenvalue from the random variable set at
    % the 95th percentile
    dataset = pcadata.(d{1}); 
    label = pcainfo.(d{1}).label;
    for v=1:size(dataset, 3)
        

        % Find scores
        data = squeeze(dataset(:, :, v));        
        latent = sort(eig(weightedcorrs(data,pcaweights)), 'descend');
        numvalidscores = find(latent>randvareigpcts, 1, 'last');
                
        % Record table of results for Parallel Analysis
        % (how much each selected PC explains the associated waveform)
        explained = pcaout.(d{1}).explained(:,v);
        varname = pcainfo.(d{1}).varnames{v};
        pcsexplained{inc,1} = [label '_' varname];
        for bb=1:numvalidscores
            pcsexplained{inc,1+bb} = explained(bb)./100;
        end
        totalvariance{inc} = sum(cell2mat(pcsexplained(inc, 2:end)));
        totalvalidpcs = totalvalidpcs + numvalidscores;
        inc = inc + 1;
        
        % Record the eigenvector matrix (coefficients), and also the scores
        % for the selected PCs for each variable
        paselected.(d{1}).(varname).coeff = squeeze(pcaout.(d{1}).coeff(:, :, v));
        paselected.(d{1}).(varname).score = squeeze(pcaout.(d{1}).score(:, 1:numvalidscores, v));        
                
        % Get quantiles and indices of scores associated with top and
        % bottom quantiles
        paquantl.(d{1}).(varname).quantl = quantile(paselected.(d{1}).(varname).score, [.025 .25 .50 .75 .975]);
        paquantl.(d{1}).(varname).bottom.idx = paselected.(d{1}).(varname).score <= paquantl.(d{1}).(varname).quantl(2, :);
        paquantl.(d{1}).(varname).top.idx = paselected.(d{1}).(varname).score >= paquantl.(d{1}).(varname).quantl(4, :);        
        
        % Get the data associated with the PC quantiles
        for bb=1:numvalidscores
            paquantl.(d{1}).(varname).bottom.mean(:, bb) = mean(pcadata.(d{1})(paquantl.(d{1}).(varname).bottom.idx(:,bb), :, v), 1);
            paquantl.(d{1}).(varname).bottom.std(:, bb) = std(pcadata.(d{1})(paquantl.(d{1}).(varname).bottom.idx(:,bb), :, v), [], 1);        
            paquantl.(d{1}).(varname).top.mean(:, bb) = mean(pcadata.(d{1})(paquantl.(d{1}).(varname).top.idx(:,bb), :, v), 1);
            paquantl.(d{1}).(varname).top.std(:, bb) = std(pcadata.(d{1})(paquantl.(d{1}).(varname).top.idx(:,bb), :, v), [], 1);
            paquantl.(d{1}).(varname).normpc(:, :, bb) = weightedcorrs([paselected.(d{1}).(varname).score(:,bb) squeeze(pcadata.(d{1})(:, :, v))], pcaweights);
            paquantl.(d{1}).(varname).normpc2(:, bb) = squeeze(paquantl.(d{1}).(varname).normpc(2:end, 1, bb)).^2;
        end
        
    end   

end


% Save results
fprintf('\nSaving Parallel Analysis inputs and outputs...\n');
if ~exist(outpath,'dir'), mkdir(outpath); end
save(fullfile(outpath,'parallelanalysis.mat'),'pcsexplained','totalvalidpcs','totalvariance','paselected','paquantl');

fprintf('------------------------------------------------\n'); 

end

