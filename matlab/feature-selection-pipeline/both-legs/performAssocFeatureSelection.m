function [acorrtable, associdx] = performAssocFeatureSelection(final, training, pcainfo)

%PERFORMASSOCFEATURESSELECTION Find PCs associated with primary features
%   Prasanna Sritharan, Mario Andrés Muñoz, August 2020
%
% Based on PCA scripts by Prasanna Sritharan for ACLR hopping 
% (published AnnBiomedEng 2022)


% user settings
addpath('..');
user = getUserScriptSettings();
outpath = user.OUTPATH2;
groups = user.GROUPS;

fprintf('Associated features analysis.\n');
fprintf('------------------------------------------------\n'); 

% determine correlation between retained PCs and all PCs after Parallel
% Analysis
corrassoc = corr(training.data, training.data(:, final.trainidx));

% determine associated features
x = 1;
acorrtable = {};
for f=1:length(final.trainidx)

    % correlation coefficients
    corrassoc(final.trainidx(f), f) = NaN;   % ignore self correlations
    iscorr = abs(corrassoc(:, f))>0.5;
    associdx = find(iscorr)';
    
    fprintf('%s\n', final.labels{f});
    
    % if no associated featues, then record empty row
    if isempty(associdx)
        fprintf('----> None\n');
        acorrtable{x,1} = final.labels{f};
        acorrtable{x,2} = 'None';
        x = x + 1;
        continue;
    end
    
    % build output table
    for j=associdx
        
        fprintf('----> %s\n',training.labels{j});
        
        % descriptives
        fprintf('\t\tCalculating descriptives...\n');
        acorrtable{x,1} = final.labels{f};
        acorrtable{x,2} = training.labels{j};
        acorrtable{x,3} = mean(training.data(pcainfo.(['is' groups{1}]),j));
        acorrtable{x,4} = std(training.data(pcainfo.(['is' groups{1}]),j));
        acorrtable{x,5} = mean(training.data(~pcainfo.(['is' groups{1}]),j));
        acorrtable{x,6} = std(training.data(~pcainfo.(['is' groups{1}]),j));
        
        % correlation coefficient (rho)
        acorrtable{x,7} = corrassoc(j,f);
        
        % t-test (t and p)
        fprintf('\t\tPerforming two-sample t-tests...\n');
        [~, pval, ~, stats] = ttest2(training.data(pcainfo.(['is' groups{1}]),j), training.data(~pcainfo.(['is' groups{1}]), j));
        acorrtable{x,8} = stats.tstat;
        acorrtable{x,9} = pval;
        
        % effect size (Hedges g)
        fprintf('\t\tPerforming effect size using Hedges g...\n');        
        mesvals = mes(training.data(pcainfo.(['is' groups{1}]),j), training.data(~pcainfo.(['is' groups{1}]),j), 'hedgesg');
        acorrtable{x,10} = mesvals.hedgesg;
        
        % increment counter
        x = x + 1;
        
    end
end


% save results
save(fullfile(outpath,'assocfeatures.mat'), 'acorrtable', 'associdx');


fprintf('------------------------------------------------\n'); 

end

