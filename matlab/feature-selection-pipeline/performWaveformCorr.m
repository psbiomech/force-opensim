function [wavecorr,dataidx] = performWaveformCorr(final,pcainfo,pcadata,pcaweights)


%PERFORMPCWAVEFORMCORR Correlation of PC scores against original waveforms
%   Prasanna Sritharan, Mario Andrés Muñoz, August 2020
%
% Based on original scripts written by Mario Andrés Muñoz, 2013


% user settings
addpath('..');
user = getUserScriptSettings();
outpath = user.FEATPATH;

fprintf('Correlations of PC scores against original waveform data.\n');
fprintf('------------------------------------------------\n'); 

% perform correlation against waveforms
wavcorr = struct;
dataidx = cell(length(final.labels),2);
labels = {'so','ik','id'};
alias = {'force','angle','moment'};
pcregexp = '(force|moment|angle)_(\w+)_\w+';
for s=1:length(final.labels)
    
    fprintf('%s\n',final.labels{s});
    
    % get indices of original waveform data
    tokens = regexpi(final.labels{s},pcregexp,'tokens');
    dataidx{s,1} = labels{strcmpi(alias,tokens{1}{1})};
    dataidx{s,2} = find(strcmpi(pcainfo.(dataidx{s,1}).varnames,tokens{1}{2}),1);

    % calculate waveform correlation coefficients
    fprintf('\tCalculating correlation coefficients...\n');
    wavecorr.coeffs(:,:,s) = weightedcorrs([final.data(:,s) squeeze(pcadata.(dataidx{s,1})(:,:,dataidx{s,2}))],pcaweights);
    wavecorr.norms(:,s) = squeeze(wavecorr.coeffs(2:end,1,s)).^2; 
    
end


% save results
save(fullfile(outpath,'waveformcorr.mat'),'wavecorr','dataidx');


fprintf('------------------------------------------------\n'); 

end

