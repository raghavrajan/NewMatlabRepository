function [] = Harini_TimeOfDayPlot(IndividualBirds, BirdNames, BirdOption, FeatureToPlot_ColNo, PlotType, DifferentMicrophones)

% ======= Feature variation across different times in the day =============
% Part of the scripts used to make plots for Harini's data
% This function plots the chosen feature for all motif syllables. The data 
% is sorted by time of day and by microphone type
% =========================================================================

Colours = 'rgbcmk';
Symbols = '+o<sd';


for i = 1:length(BirdNames),
    if ~((strcmp(BirdOption, BirdNames{i})) || strcmp(BirdOption, 'All'))
        continue;
    end
    for MicrophoneType = 1:length(DifferentMicrophones),
    % Sorting plots by microphone type
        MicrophoneIndex = find(cellfun(@length, strfind(IndividualBirds(i).Microphone, DifferentMicrophones{MicrophoneType})));
        if (isempty(MicrophoneIndex))
            continue;
        end

        TempBirdParameters = IndividualBirds(i).SortedBirdParameters;
        % Keep only parts that correspond to the microphone type being
        % analyzed
        TempBirdParameters = TempBirdParameters(find([TempBirdParameters.MicrophoneIndices] == MicrophoneIndex));
        for j = 1:length(TempBirdParameters(1).MotifLabels),
            figure(j);
            subplot(length(IndividualBirds(i).Microphone),1, MicrophoneIndex);
            hold on;
            for k = 1:length(TempBirdParameters),
                SyllableIndices = find((char(TempBirdParameters(k).SyllableData(:,1)) == TempBirdParameters(1).MotifLabels(j)) & (~isnan(TempBirdParameters(k).SAPFeatsMatrix(:,FeatureToPlot_ColNo)))); % The second condition here ensures that the log amplitude is not NaN
                plot(TempBirdParameters(k).SyllableTimes(SyllableIndices), TempBirdParameters(k).SAPFeatsMatrix(SyllableIndices,FeatureToPlot_ColNo), [Colours(TempBirdParameters(k).ConditionIndices), Symbols(TempBirdParameters(k).MicrophoneIndices)], 'MarkerSize', 6);
            end
            title([BirdNames{i}, '.SyllLabel.', TempBirdParameters(1).MotifLabels(j), '.TimeofDay.vs.', TempBirdParameters(1).SAPFeat_FieldNames{FeatureToPlot_ColNo}, '.', DifferentMicrophones{MicrophoneType}]);
            xlabel('Time of day (hours)');
            ylabel(TempBirdParameters(1).SAPFeat_FieldNames{FeatureToPlot_ColNo});
        end
    end
end