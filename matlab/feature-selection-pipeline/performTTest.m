function ttable = performTTest(final,pcainfo)


%PERFORMTTEST Perform T-test on final feature set
%   Prasanna Sritharan, Mario Andrés Muñoz, August 2020
%
% Based on original scripts written by Mario Andrés Muñoz, 2013
%
% Tests performed at 99.9% significance (alpha=0.001). Also calculate
% Hedges g effect size and descriptives.


% user settings
addpath('..');
user = getUserScriptSettings();
outpath = user.FEATPATH;


fprintf('T-tests, effect size and descriptives for final PC set.\n');
fprintf('------------------------------------------------\n'); 

% run tests and build output table
npcs = length(final.labels);
ttable = cell(npcs,8);
for n=1:npcs
    
    fprintf('%s\n',final.labels{n});
    
    % descriptives
    fprintf('\tCalculating descriptives...\n');
    ttable{n,1} = final.labels{n};
    ttable{n,2} = mean(final.data(pcainfo.isACLR,n));
    ttable{n,3} = std(final.data(pcainfo.isACLR,n));
    ttable{n,4} = mean(final.data(~pcainfo.isACLR,n));
    ttable{n,5} = std(final.data(~pcainfo.isACLR,n));    
    
    % t-test (t and p)
    fprintf('\tPerforming two-sample t-tests...\n');
    [~,pval,~,stats] = ttest2(final.data(pcainfo.isACLR,n),final.data(~pcainfo.isACLR,n));
    ttable{n,6} = stats.tstat;
    ttable{n,7} = pval;
    
    % effect size (Hedges g)
    fprintf('\tPerforming effect size using Hedges g...\n');
    mesout = mes(final.data(pcainfo.isACLR,n),final.data(~pcainfo.isACLR,n),'hedgesg');
    ttable{n,8} = mesout.hedgesg;
    
end


% save results
save(fullfile(outpath,'ttable.mat'),'ttable');

fprintf('------------------------------------------------\n');

end

