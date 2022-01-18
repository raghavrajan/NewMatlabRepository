function [] = Harini_ChangeWithinSameDay(IndividualBirds, BirdNames, BirdOption, FeatureToPlot_ColNo, PlotType, DifferentMicrophones)

% ======= Feature variation for experiments done on the same day ==========
% Part of the scripts used to make plots for Harini's data
% This function plots the chosen feature for all motif syllables. The data 
% is sorted by experiments done on the same day and by microphone type.
% =========================================================================

Colours = 'rgbcmk';
Symbols = '+o<sd';

% If there are two conditions done on the same day, then join those by
% lines and plot

% Now to plot for each bird and each motif syllable within that bird
for i = 1:length(BirdNames),
    if ~((strcmp(BirdOption, BirdNames{i})) || strcmp(BirdOption, 'All'))
        continue;
    end
    TempBirdParameters = IndividualBirds(i).SortedBirdParameters;
    for j = 1:length(TempBirdParameters(1).MotifLabels),
        figure;
        hold on;
        PlotIndex = 1;
        XTickNums = [];
        XLabelString = [];
        for k = 1:length(IndividualBirds(i).RecordingDays),
            for Microphone = 1:length(IndividualBirds(i).Microphone),
                MatchingDays = find(([TempBirdParameters.RecordingDayIndex] == k) & ([TempBirdParameters.MicrophoneIndices] == Microphone));
                AmplitudeWithinDay = [];
                for Days = MatchingDays(:)',
                    SyllableIndices = find((char(TempBirdParameters(Days).SyllableData(:,1)) == TempBirdParameters(1).MotifLabels(j)) & (~isnan(TempBirdParameters(Days).SAPFeatsMatrix(:,FeatureToPlot_ColNo)))); % The second condition here ensures that the log amplitude is not NaN
                    AmplitudeWithinDay(end+1,:) = [TempBirdParameters(Days).ConditionIndices eval([PlotType, '(TempBirdParameters(Days).SAPFeatsMatrix(SyllableIndices,FeatureToPlot_ColNo))']) std(TempBirdParameters(Days).SAPFeatsMatrix(SyllableIndices,FeatureToPlot_ColNo))];
                end
                if (size(AmplitudeWithinDay,1) > 1)
                    if (length(unique(AmplitudeWithinDay(:,1))) < length(AmplitudeWithinDay(:,1)))
                        [UniqueVals, UniqueIndices] = unique(AmplitudeWithinDay(:,1));
                        AmplitudeWithinDay = AmplitudeWithinDay(UniqueIndices,:);
                    end
                    plot((PlotIndex + AmplitudeWithinDay(:,1)), (AmplitudeWithinDay(:,2) - AmplitudeWithinDay(1,2))/abs(AmplitudeWithinDay(1,2)), [Colours(Microphone), 's-']);
                    for SessionNo = 1:size(AmplitudeWithinDay,1),
                        XLabelString{end+1} = IndividualBirds(i).Conditions{(AmplitudeWithinDay(SessionNo,1))};
                        XTickNums(end+1) = PlotIndex + AmplitudeWithinDay(SessionNo, 1);
                    end
                    PlotIndex = PlotIndex + max(AmplitudeWithinDay(:,1));
                end
            end
        end    
        set(gca, 'XTick', XTickNums, 'XTickLabel', XLabelString, 'XTickLabelRotation', 45);
        title([BirdNames{i}, '.SyllLabel.', TempBirdParameters(1).MotifLabels(j), '.WithinDay.', TempBirdParameters(1).SAPFeat_FieldNames{FeatureToPlot_ColNo}, '.Comparisons']);
        ylabel([PlotType, ' ', TempBirdParameters(1).SAPFeat_FieldNames{FeatureToPlot_ColNo}, ' change relative to first session on that day']);
    end
end