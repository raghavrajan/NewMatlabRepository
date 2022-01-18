function [] = Harini_UndirHFMicChanges(IndividualBirds, BirdNames, BirdOption, FeatureToPlot_ColNo, PlotType, DifferentMicrophones)

% ======= Feature variation across days ===================================
% Part of the scripts used to make plots for Harini's data
% This function plots the chosen feature for all motif syllables. The data 
% is sorted by day and by microphone type.
% =========================================================================

Colours = 'rgbcmk';
Symbols = '+o<sd';

for i = 1:length(IndividualBirds),
    [IndividualBirds(i).RawSyllLogAmplitudeMeanValue, IndividualBirds(i).AdjustedSyllLogAmplitudeMeanValue] = Harini_AdjustLogAmplitude(IndividualBirds(i), BirdNames{i});
end

disp('Finished');