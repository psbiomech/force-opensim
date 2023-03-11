function ttable = performTTest(final, pcainfo)


%PERFORMTTEST Perform T-test on final feature set
%   Prasanna Sritharan, March 2022
%
% Based on PCA scripts by Prasanna Sritharan for ACLR hopping 
% (published AnnBiomedEng 2022)


% User settings
addpath('..');
user = getUserScriptSettings();
outpath = user.OUTPATH2;
groups = user.GROUPS;


fprintf('T-tests, effect size and descriptives for final PC set.\n');
fprintf('------------------------------------------------'); 

% Data tables
ttable = struct;

% Run tests and build output table
npcs = length(final.labels);
ttable.table = cell(npcs, 8);
for n=1:npcs
    
    fprintf('%s\n',final.labels{n});
    
    % Descriptives
    fprintf('Calculating descriptives...\n');
    ttable.table{n,1} = final.labels{n};
    ttable.table{n,2} = mean(final.data(pcainfo.(['is' groups{1}]), n));
    ttable.table{n,3} = std(final.data(pcainfo.(['is' groups{1}]), n));
    ttable.table{n,4} = mean(final.data(~pcainfo.(['is' groups{1}]), n));
    ttable.table{n,5} = std(final.data(~pcainfo.(['is' groups{1}]), n));    
    
    % T-test (t and p)
    fprintf('Performing two-sample t-tests...\n');
    [~, pval, ~, stats] = ttest2(final.data(pcainfo.(['is' groups{1}]), n), final.data(~pcainfo.(['is' groups{1}]), n));
    ttable.table{n,6} = stats.tstat;
    ttable.table{n,7} = pval;
    
    % Effect size (Hedges g)
    fprintf('---> Performing effect size using Hedges g...\n');
    mesout = mes(final.data(pcainfo.(['is' groups{1}]), n), final.data(~pcainfo.(['is' groups{1}]), n), 'hedgesg');
    ttable.table{n,8} = mesout.hedgesg;
    
    % Headers
    ttable.headers = {'feature', 'mean1', 'sd1', 'mean2', 'sd2', 't', 'p', 'g'};

end
    


% Save results
fprintf('\nSaving T-tests inputs and outputs...\n');
if ~exist(outpath,'dir'), mkdir(outpath); end
save(fullfile(outpath,'ttable.mat'),'ttable');

fprintf('------------------------------------------------\n');

end

