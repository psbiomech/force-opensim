function tbls = generateOutputTables(pcsexplained,totalvariance,ttable,acorrtable)


%GENERATEOUTPUTTABLES Generate output tables
%   Prasanna Sritharan, Mario Andrés Muñoz, August 2020
%
% Based on original scripts written by Mario Andrés Muñoz, 2013


% user settings
addpath('..');
user = getUserScriptSettings();
outpath = user.FEATPATH;


fprintf('Generating output tables.\n');
fprintf('------------------------------------------------\n');

% output folder
tabledir = fullfile(outpath,'tables');
if ~exist(tabledir,'dir'), mkdir(tabledir); end

% Table 1: variance explained by PCs
fprintf('Creating Table 1: Variance explained by PCs...\n');
pcsexplained_percent_1decpl = pcsexplained;
pcsexplained_percent_1decpl(:,2:end) = cellfun(@(x) round(x*100,1),pcsexplained_percent_1decpl(:,2:end),'UniformOutput',false);
totalvariance_percent_1decpl = cellfun(@(x) round(x*100,1),totalvariance,'UniformOutput',false);
headers1 = ['variable',cellfun(@(x) ['PC' num2str(x)],num2cell(1:size(pcsexplained_percent_1decpl,2)-1),'UniformOutput',false),'total'];
tbls.table1 = [headers1; [pcsexplained_percent_1decpl,totalvariance_percent_1decpl']];
writecell(tbls.table1,fullfile(tabledir,'table1_pcsexplained.csv'));

% Table 2: t-test, effect size and descriptives
fprintf('Creating Table 2: T-test, effect size and descriptives for final selected PCs...\n');
headers2 = {'feature','mean_aclr','sd_aclr','mean_ctrl','sd_ctrl','t','P','g'};
tbls.table2 = [headers2; ttable];
writecell(tbls.table2,fullfile(tabledir,'table2_ttestout.csv'));

% Table 3: qualitative interpretation of PCs
% N/A, descriptive table
fprintf('Creating Table 3: Qualitative interpretation of PCs...\n');
tbls.table3 = 'N/A, qualitative interpretation of PCs, descriptive table';

% Table 4: associated features
fprintf('Creating Table 4: Associated features...\n');
headers4 = {'main_feature','associated_feature','mean_aclr','sd_aclr','mean_ctrl','sd_ctrl','rho','t','P','g'};
tbls.table4 = [headers4; acorrtable];
writecell(tbls.table4,fullfile(tabledir,'table4_assocfeatures.csv'));

fprintf('------------------------------------------------\n');

end

