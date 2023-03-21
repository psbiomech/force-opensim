%% FEATURE SELECTION PIPELINE - FORCE SDP
% Prasanna Sritharan, March 2022
%
% Note: back up the results of any previous run of this script as the
% Sequential Feature Selection analysis involves random selection of
% principal components, and therefore results from the new run can
% differ substantially from those obtained on a previous run.

% User settings
addpath('..');
user = getUserScriptSettings();
outpath = user.OUTPATH1;


%% FEATURE SELECTION PIPELINE

% Weighted PCA
[pcadata, pcaweights, pcaout, pcainfo] = performWeightedPCA('rajagopal');

% Parallel Analysis
[pcsexplained, totalvalidpcs, totalvariance, paselected, paquantl] = performParallelAnalysis(pcadata, pcaweights, pcaout, pcainfo);

% Feature selection
[final, weissind, training] = performFeatureSelection(paselected, pcainfo);

% t-tests, effect size and descriptives
ttable = performTTest(final, pcainfo);

% correlation of PCs with original waveforms
[wavecorr, dataidx] = performWaveformCorr(final, pcainfo, pcadata, pcaweights);
wcorr.wavecorr = wavecorr;
wcorr.dataidx = dataidx;

% associated features
[acorrtable,associdx] = performAssocFeatureSelection(final,training,pcainfo);


%% TABLES AND FIGURES

% generate output tables
tbls = generateOutputTables(pcsexplained, totalvariance, ttable, acorrtable);

% generate output figures
generateOutputFigures(paselected, paquantl, final);


%% SAVE OUTPUTS

% pcaw.pcadata = pcadata;
% pcaw.pcaweights = pcaweights;
% pcaw.pcaout = pcaout;
% pcaw.pcainfo = pcainfo;
% 
% pa.pcsexplained = pcsexplained;
% pa.totalvalidpcs = totalvalidpcs;
% pa.totalvariance = totalvariance;
% pa.paselected = paselected;
% pa.paquantl = paquantl;
% 
% fs.final = final;
% fs.weissind = weissind;
% fs.training = training;
% 
% tt.ttable = ttable;
% 
% acorr.acorrtable = wavecorr;
% acorr.associdx = dataidx;

% % save all structs
% save(fullfile(outpath,'fsp_all_outputs.mat'),'pcaw','pa','fs','tt','wcorr','acorr','tbls');
