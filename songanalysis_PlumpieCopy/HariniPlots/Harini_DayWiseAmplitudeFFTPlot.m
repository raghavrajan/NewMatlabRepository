function [] = Harini_DayWiseAmplitudeFFTPlot(IndividualBirds, BirdNames, BirdOption, FeatureToPlot_ColNo, PlotType, DifferentMicrophones)

% ======= Feature variation across days ===================================
% Part of the scripts used to make plots for Harini's data
% This function plots the chosen feature for all motif syllables. The data 
% is sorted by day and by microphone type.
% =========================================================================

Colours = 'rgbcmk';
Symbols = '+o<sd';

[BaselineSubtractedMeanFFTLogAmpValues, RawMeanFFTLogAmpValues] = Harini_CalcMeanMaxFFTLogAmplitude(IndividualBirds, BirdNames);

for i = 1:length(BaselineSubtractedMeanFFTLogAmpValues),
    IndividualBirds(i).BaselineSubtractedMeanFFTLogAmpValues = BaselineSubtractedMeanFFTLogAmpValues(i).MeanFFTLogAmpValues;
    IndividualBirds(i).BaselineSubtractedMaxFFTLogAmpValues = BaselineSubtractedMeanFFTLogAmpValues(i).MaxFFTLogAmpValues;
end

for i = 1:length(RawMeanFFTLogAmpValues),
    IndividualBirds(i).RawMeanFFTLogAmpValues = RawMeanFFTLogAmpValues(i).MeanFFTLogAmpValues;
    IndividualBirds(i).RawMaxFFTLogAmpValues = RawMeanFFTLogAmpValues(i).MaxFFTLogAmpValues;
end

% Now to plot for each bird and each motif syllable within that bird - this
% is already sorted in the order of days, so it can just be plotted in
% order

for i = 1:length(BirdNames),
    MicrophoneTypes = unique(IndividualBirds(i).AllMicrophoneIndices);
    for MicrophoneType = MicrophoneTypes(:)',
        for j = 1:length(IndividualBirds(i).SortedBirdParameters(1).MotifLabels),
            Matches = find((char(IndividualBirds(i).AllSyllableData(:,1)) == IndividualBirds(i).SortedBirdParameters(1).MotifLabels(j)) & (IndividualBirds(i).AllMicrophoneIndices(:) == MicrophoneType));
            UniqueRecordingDays = unique(IndividualBirds(i).AllRecordingDayIndices(Matches));
            figure;
            hold on;
            for RecordingDay = UniqueRecordingDays(:)',
                Matches = find((char(IndividualBirds(i).AllSyllableData(:,1)) == IndividualBirds(i).SortedBirdParameters(1).MotifLabels(j)) & (IndividualBirds(i).AllRecordingDayIndices(:) == RecordingDay) & (IndividualBirds(i).AllMicrophoneIndices(:) == MicrophoneType));
                UniqueConditions = unique(IndividualBirds(i).AllConditionIndices(Matches));
                for k = UniqueConditions(:)',
                % for k = 6,
                    Matches = find((char(IndividualBirds(i).AllSyllableData(:,1)) == IndividualBirds(i).SortedBirdParameters(1).MotifLabels(j)) & (IndividualBirds(i).AllRecordingDayIndices(:) == RecordingDay) & (IndividualBirds(i).AllConditionIndices(:) == k) & (IndividualBirds(i).AllMicrophoneIndices(:) == MicrophoneType));
                    if (~isempty(Matches))
                        errorbar(RecordingDay + 0.15*(k-1), nanmean(IndividualBirds(i).RawMaxFFTLogAmpValues(Matches)), nanstd(IndividualBirds(i).RawMaxFFTLogAmpValues(Matches))/sqrt(length(find(~isnan(IndividualBirds(i).RawMaxFFTLogAmpValues(Matches))))), [Colours(k), 's-']); 
                    end
                end
            end
            title([BirdNames{i}, ': Syll ', IndividualBirds(i).SortedBirdParameters(1).MotifLabels(j)]);
        end
    end
end
