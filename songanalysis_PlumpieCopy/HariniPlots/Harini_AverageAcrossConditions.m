function [] = Harini_AverageAcrossConditions(IndividualBirds, BirdNames, BirdOption, FeatureToPlot_ColNo, PlotType, DifferentMicrophones)

% ======= Average across conditions =======================================
% Plots the average feature values - averaged over different repetitions of
% the same condition over multiple days. 
% =========================================================================

Colours = 'rgbcmk';
Symbols = '+o<sd';

% Now to plot for each bird and each motif syllable within that bird
for i = 1:length(BirdNames),
    if ~((strcmp(BirdOption, BirdNames{i})) || strcmp(BirdOption, 'All'))
        continue;
    end
    for j = 1:length(IndividualBirds(i).SortedBirdParameters(1).MotifLabels),
        figure;
        hold on;
        XLabelString = [];
        PlotIndex = 1;
        for MicrophoneIndex = 1:length(IndividualBirds(i).Microphone),
            for ConditionIndex = 1:length(IndividualBirds(i).Conditions),
                MatchingSyllables = find((char(IndividualBirds(i).AllSyllableData(:,1)) == IndividualBirds(i).SortedBirdParameters(1).MotifLabels(j)) & (IndividualBirds(i).AllConditionIndices == ConditionIndex) & (IndividualBirds(i).AllMicrophoneIndices == MicrophoneIndex) & (~isnan(IndividualBirds(i).AllSyllableFeatValues(:,FeatureToPlot_ColNo))));
                SyllFeatsToPlot = IndividualBirds(i).AllSyllableFeatValues(MatchingSyllables,FeatureToPlot_ColNo);
                plot(ones(length(MatchingSyllables),1)*(PlotIndex-0.2), SyllFeatsToPlot, 'k.', 'MarkerSize', 6);
                errorbar(PlotIndex, eval([PlotType, '(SyllFeatsToPlot)']), std(SyllFeatsToPlot), 'rs', 'LineWidth', 2, 'MarkerSize', 6);
                MeanWithinCondition{i}{MicrophoneIndex}(j,ConditionIndex) = eval([PlotType, '(SyllFeatsToPlot)']); 
                XLabelString{end+1} = [IndividualBirds(i).Conditions{ConditionIndex}, '.', IndividualBirds(i).Microphone{MicrophoneIndex}];
                PlotIndex = PlotIndex + 1;
            end
        end
        set(gca, 'XTick', 1:1:length(XLabelString), 'XTickLabel', XLabelString, 'XTickLabelRotation', 45);
        title([BirdNames{i}, '.SyllLabel.', IndividualBirds(i).SortedBirdParameters(1).MotifLabels(j), '.AverageAcrossConditions.', IndividualBirds(i).SortedBirdParameters(1).SAPFeat_FieldNames{FeatureToPlot_ColNo}]);
        ylabel(IndividualBirds(i).SortedBirdParameters(1).SAPFeat_FieldNames{FeatureToPlot_ColNo});
    end
end

% Now to plot this for all birds - plot separately for microphone and for
% backpack

if (strcmp(BirdOption, 'All'))
    DifferentMicrophones = [{'BP'} {'MM'} {'HF'}];

    for i = 1:length(DifferentMicrophones),
        figure;
        hold on;
        for j = 1:length(IndividualBirds),
            Index = find(cellfun(@length, strfind(IndividualBirds(j).Microphone, DifferentMicrophones{i})));
            if (~isempty(Index))
                plot((MeanWithinCondition{j}{Index}./repmat(MeanWithinCondition{j}{Index}(:,1), 1, length(IndividualBirds(i).Conditions)))', [Colours(mod(j,length(Colours)) + 1), 's-']);
            end
        end
        set(gca, 'XTick', 1:1:length(IndividualBirds(1).Conditions), 'XTickLabel', IndividualBirds(1).Conditions);
        ylabel(['Syllable ', BirdParameters(1).SAPFeat_FieldNames{FeatureToPlot_ColNo}, ' normalized to L1']);
        title(['Syllable ', BirdParameters(1).SAPFeat_FieldNames{FeatureToPlot_ColNo}, ' normalized to L1 for ', DifferentMicrophones{i}]);
        axis tight;
        Temp = axis;
        Temp = [0 5.5 0.95*Temp(3) 1.05*Temp(4)];
        axis(Temp);
        for j = 1:length(IndividualBirds),
            Index = find(cellfun(@length, strfind(IndividualBirds(j).Microphone, DifferentMicrophones{i})));
            if (~isempty(Index))
                text(0.1, Temp(3) + (Temp(4) - Temp(3))*(j/20), BirdNames{j}, 'Color', Colours(mod(j,length(Colours)) + 1));
            end
        end
    end
end

