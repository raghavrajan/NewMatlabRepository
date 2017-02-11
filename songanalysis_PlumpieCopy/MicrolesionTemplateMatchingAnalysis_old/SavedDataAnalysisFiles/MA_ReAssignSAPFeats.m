function [Feats, FeatureNames, FeatureLabels] = MA_ReAssignSAPFeats(FeatValues, SAPFeatureFile)

% Features to be used : Duration, LogAmplitude, Entropy, MeanFrequency,
% PitchGoodness, FrequencyModulation, EntropyVariance

FeatureNames{1} = 'Duration';
FeatureNames{2} = 'LogAmplitude';
FeatureNames{3} = 'Entropy';
FeatureNames{4} = 'MeanFrequency';
FeatureNames{5} = 'PitchGoodness';
FeatureNames{6} = 'FrequencyModulation';
FeatureNames{7} = 'EntropyVariance';
% FeatureNames{8} = 'FundamentalFrequency';

FeatureLabels{1} = 'Duration (sec)';
FeatureLabels{2} = 'Log Amplitude (dB)';
FeatureLabels{3} = 'Entropy';
FeatureLabels{4} = 'Mean Frequency (Hz)';
FeatureLabels{5} = 'Pitch Goodness';
FeatureLabels{6} = 'Frequency Modulation';
FeatureLabels{7} = 'Entropy Variance';
% FeatureLabels{8} = 'Fundamental Frequency (Hz)';

Feats = [];
for i = 1:length(FeatureNames),
    Feats(:,i) = eval(['FeatValues.', FeatureNames{i}]);
end

save(SAPFeatureFile, 'Feats', 'FeatureNames', 'FeatureLabels');