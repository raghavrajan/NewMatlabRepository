function [] = Harini_ConditionWiseFeaturePlot(IndividualBirds, BirdNames, BirdOption, FeatureToPlot_ColNo, PlotType, DifferentMicrophones)

% ======= Feature variation across conditinos =============================
% Part of the scripts used to make plots for Harini's data
% This function plots the chosen feature for all motif syllables. The data 
% is sorted by condition and by microphone type
% =========================================================================

Colours = 'rgbcmk';
Symbols = '+o<sd';

% Now to plot for each bird and each motif syllable within that bird
for i = 1:length(BirdNames),

    if ~((strcmp(BirdOption, BirdNames{i})) || strcmp(BirdOption, 'All'))
        continue;
    end

    TempBirdParameters = IndividualBirds(i).SortedBirdParameters;
    for j = 1:length(TempBirdParameters(1).MotifLabels),
        figure;
        hold on;
        XLabelString = [];
        PlotIndex = 1;
        for MicrophoneIndex = 1:length(IndividualBirds(i).Microphone),
            for ConditionIndex = 1:length(IndividualBirds(i).Conditions),
                MatchingConditions = find(([IndividualBirds(i).SortedBirdParameters.ConditionIndices] == ConditionIndex) & ([IndividualBirds(i).SortedBirdParameters.MicrophoneIndices] == MicrophoneIndex));
                for k = MatchingConditions(:)',
                    SyllableIndices = find((char(TempBirdParameters(k).SyllableData(:,1)) == TempBirdParameters(1).MotifLabels(j)) & (~isnan(TempBirdParameters(k).SAPFeatsMatrix(:,FeatureToPlot_ColNo)))); % The second condition here ensures that the log amplitude is not NaN
                    plot(ones(length(SyllableIndices))*(PlotIndex-0.2), TempBirdParameters(k).SAPFeatsMatrix(SyllableIndices,FeatureToPlot_ColNo), 'k.', 'MarkerSize', 6);
                    errorbar(PlotIndex, eval([PlotType, '(TempBirdParameters(k).SAPFeatsMatrix(SyllableIndices,FeatureToPlot_ColNo))']), std(TempBirdParameters(k).SAPFeatsMatrix(SyllableIndices,FeatureToPlot_ColNo)), 'rs', 'LineWidth', 2, 'MarkerSize', 6);
                    XLabelString{end+1} = [TempBirdParameters(k).DataLabel, '.', TempBirdParameters(k).Condition, '.', TempBirdParameters(k).Microphone];
                    PlotIndex = PlotIndex + 1;
                end
            end
        end
        set(gca, 'XTick', 1:1:length(TempBirdParameters), 'XTickLabel', XLabelString, 'XTickLabelRotation', 45);
        title([BirdNames{i}, '.SyllLabel.', TempBirdParameters(1).MotifLabels(j), '.AcrossDay.', TempBirdParameters(1).SAPFeat_FieldNames{FeatureToPlot_ColNo}]);
        ylabel(TempBirdParameters(1).SAPFeat_FieldNames{FeatureToPlot_ColNo});
    end
end