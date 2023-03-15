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
outpath = user.OUTPATH2;



%% FEATURE SELECTION PIPELINE

% Weighted PCA
[pcadata, pcaweights, pcaout, pcainfo] = performWeightedPCA('rajagopal');
pcaw.pcadata = pcadata;
pcaw.pcaweights = pcaweights;
pcaw.pcaout = pcaout;
pcaw.pcainfo = pcainfo;

% Parallel Analysis
[pcsexplained, totalvalidpcs, totalvariance, paselected, paquantl] = performParallelAnalysis(pcadata, pcaweights, pcaout, pcainfo);
pa.pcsexplained = pcsexplained;
pa.totalvalidpcs = totalvalidpcs;
pa.totalvariance = totalvariance;
pa.paselected = paselected;
pa.paquantl = paquantl;

% Feature selection
[final, weissind, training] = performFeatureSelection(paselected, pcainfo);
fs.final = final;
fs.weissind = weissind;
fs.training = training;

% t-tests, effect size and descriptives
ttable = performTTest(final, pcainfo);
tt.ttable = ttable;

% correlation of PCs with original waveforms
[wavecorr, dataidx] = performWaveformCorr(final, pcainfo, pcadata, pcaweights);
wcorr.wavecorr = wavecorr;
wcorr.dataidx = dataidx;

% associated features
[acorrtable, associdx] = performAssocFeatureSelection(final, training, pcainfo);
acorr.acorrtable = wavecorr;
acorr.associdx = dataidx;


%% TABLES AND FIGURES

% % generate output tables
% tbls = generateOutputTables(pcsexplained,totalvariance,ttable,acorrtable);
% 
% % generate output figures
% generateOutputFigures(paselected,paquantl,final);


%% SAVE OUTPUTS

% % save all structs
% save(fullfile(outpath,'fsp_outputs_all.mat'),'pcaw','pa','fs','tt','wcorr','acorr','tbls');
