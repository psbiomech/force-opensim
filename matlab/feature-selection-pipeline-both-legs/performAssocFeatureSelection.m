function [acorrtable,associdx] = performAssocFeatureSelection(final,training,pcainfo)

%PERFORMASSOCFEATURESSELECTION Find PCs associated with primary features
%   Prasanna Sritharan, Mario Andrés Muñoz, August 2020
%
% Based on original scripts written by Mario Andrés Muñoz, 2013


% user settings
addpath('..');
user = getUserScriptSettings();
outpath = user.FEATPATH;

fprintf('Associated features analysis.\n');
fprintf('------------------------------------------------\n'); 

% determine correlation between retained PCs and all PCs after Parallel
% Analysis
corrassoc = corr(training.data,training.data(:,final.trainidx));

% determine associated features
x = 1;
acorrtable = {};
for f=1:length(final.trainidx)

    % correlation coefficients
    corrassoc(final.trainidx(f),f) = NaN;   % ignore self correlations
    iscorr = abs(corrassoc(:,f))>0.5;
    associdx = find(iscorr)';
    
    fprintf('%s\n',final.labels{f});
    
    % if no associated featues, then record empty row
    if isempty(associdx)
        fprintf('\t----> None\n');
        acorrtable{x,1} = final.labels{f};
        acorrtable{x,2} = 'None';
        x = x + 1;
        continue;
    end
    
    % build output table
    for j=associdx
        
        fprintf('\t----> %s\n',training.labels{j});
        
        % descriptives
        fprintf('\t\t\tCalculating descriptives...\n');
        acorrtable{x,1} = final.labels{f};
        acorrtable{x,2} = training.labels{j};
        acorrtable{x,3} = mean(training.data(pcainfo.isACLR,j));
        acorrtable{x,4} = std(training.data(pcainfo.isACLR,j));
        acorrtable{x,5} = mean(training.data(~pcainfo.isACLR,j));
        acorrtable{x,6} = std(training.data(~pcainfo.isACLR,j));
        
        % correlation coefficient (rho)
        acorrtable{x,7} = corrassoc(j,f);
        
        % t-test (t and p)
        fprintf('\t\t\tPerforming two-sample t-tests...\n');
        [~,pval,~,stats] = ttest2(training.data(pcainfo.isACLR,j),training.data(~pcainfo.isACLR,j));
        acorrtable{x,8} = stats.tstat;
        acorrtable{x,9} = pval;
        
        % effect size (Hedges g)
        fprintf('\t\t\tPerforming effect size using Hedges g...\n');        
        mesvals = mes(training.data(pcainfo.isACLR,j),training.data(~pcainfo.isACLR,j),'hedgesg');
        acorrtable{x,10} = mesvals.hedgesg;
        
        % increment counter
        x = x + 1;
        
    end
end


% save results
save(fullfile(outpath,'assocfeatures.mat'),'acorrtable','associdx');


fprintf('------------------------------------------------\n'); 

end

