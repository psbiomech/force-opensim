function [wavecorr, dataidx] = performWaveformCorr(final, pcainfo, pcadata, pcaweights)


%PERFORMPCWAVEFORMCORR Correlation of PC scores against original waveforms
%   Prasanna Sritharan, March 2022
%
% Based on PCA scripts by Prasanna Sritharan for ACLR hopping 
% (published AnnBiomedEng 2022)


% User settings
addpath('..');
user = getUserScriptSettings();
outpath = user.OUTPATH2;


fprintf('Correlations of PC scores against original waveform data.\n');
fprintf('------------------------------------------------\n'); 

% Data tables
wavecorr = struct;
dataidx = cell(length(final.labels), 2);

% Perform correlation against waveforms
labels = {'ik','id'};
alias = {'angle','moment'};
pcregexp = '(moment|angle)_(\w+)_PC\d+$';
for s=1:length(final.labels)
    
    fprintf('%s\n', final.labels{s});
    
    % Get indices of original waveform data
    tokens = regexpi(final.labels{s}, pcregexp, 'tokens');
    dataidx{s, 1} = labels{strcmpi(alias, tokens{1}{1})};
    dataidx{s, 2} = find(strcmpi(pcainfo.(dataidx{s, 1}).varnames, tokens{1}{2}), 1);

    % Calculate waveform correlation coefficients
    fprintf('Calculating correlation coefficients...\n');
    wavecorr.coeffs(:, :, s) = weightedcorrs([final.data(:,s) squeeze(pcadata.(dataidx{s, 1})(:, :, dataidx{s, 2}))], pcaweights);
    wavecorr.norms(:, s) = squeeze(wavecorr.coeffs(2:end,1,s)).^2; 
    
end

% Save results
fprintf('\nSaving Waveform Correlation inputs and outputs...\n');
if ~exist(outpath,'dir'), mkdir(outpath); end
save(fullfile(outpath,'waveformcorr.mat'),'wavecorr','dataidx');


fprintf('------------------------------------------------\n'); 

end

