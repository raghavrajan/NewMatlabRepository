function [] = Harini_DayWiseFFAutocorrKaoPlot(IndividualBirds, BirdNames, BirdOption, FeatureToPlot_ColNo, PlotType, DifferentMicrophones)

% ======= Feature variation across days ===================================
% Part of the scripts used to make plots for Harini's data
% This function plots the chosen feature for all motif syllables. The data 
% is sorted by day and by microphone type.
% =========================================================================

Colours = 'rgbcmk';
Symbols = '+o<sd';

% Now to plot for each bird and each motif syllable within that bird - this
% is already sorted in the order of days, so it can just be plotted in
% order

for MicrophoneType = 1:length(DifferentMicrophones),
    % Sorting plots by microphone type
    for i = 1:length(BirdNames),
        MicrophoneIndex = find(cellfun(@length, strfind(IndividualBirds(i).Microphone, DifferentMicrophones{MicrophoneType})));
        if (isempty(MicrophoneIndex))
            continue;
        end

        if ~((strcmp(BirdOption, BirdNames{i})) || strcmp(BirdOption, 'All'))
            continue;
        end

        TempBirdParameters = IndividualBirds(i).SortedBirdParameters;
        % Keep only parts that correspond to the microphone type being
        % analyzed
        TempBirdParameters = TempBirdParameters(find([TempBirdParameters.MicrophoneIndices] == MicrophoneIndex));

        SessionMeans = [];
        SessionCVs = [];
        for j = 1:length(TempBirdParameters(1).MotifLabels),
            AllFeatValues = [];
            figure;
            hold on;
            for k = 1:length(TempBirdParameters),
                SyllableIndices = find((char(TempBirdParameters(k).SyllableData(:,1)) == TempBirdParameters(1).MotifLabels(j)) & (~isnan(TempBirdParameters(k).FFAutocorrKao(:)))); % The second condition here ensures that the log amplitude is not NaN
                TempFFValues = TempBirdParameters(k).FFAutocorrKao(SyllableIndices);
                TempSyllDurs = TempBirdParameters(k).SAPFeatsMatrix(SyllableIndices,1);
                
                % Remove outliers based on syllable duration and based on
                % FF
%                 DurOutlierThreshold = [(prctile(TempSyllDurs,25) - 3*iqr(TempSyllDurs)) (prctile(TempSyllDurs,75) + 3*iqr(TempSyllDurs))];
%                 DurOutlierIndices = find((TempSyllDurs < DurOutlierThreshold(1)) | (TempSyllDurs > DurOutlierThreshold(2)));
%                 TempFFValues(DurOutlierIndices) = [];
%                 disp(['Removed ', num2str(length(DurOutlierIndices)), ' based on syllable duration']);
%                 
%                 FFOutlierThreshold = [(prctile(TempFFValues,25) - 3*iqr(TempFFValues)) (prctile(TempFFValues,75) + 3*iqr(TempFFValues))];
%                 FFOutlierIndices = find((TempFFValues < FFOutlierThreshold(1)) | (TempFFValues > FFOutlierThreshold(2)));
%                 TempFFValues(FFOutlierIndices) = [];
%                 disp(['Removed ', num2str(length(FFOutlierIndices)), ' based on syllable FF']);
%                 
                AllFeatValues = [AllFeatValues; [ones(length(TempFFValues),1)*k TempFFValues]];
                
                plot(ones(length(TempFFValues))*(k-0.2), TempFFValues, 'k.', 'MarkerSize', 6);
                errorbar(k, mean(TempFFValues), std(TempFFValues), 'rs', 'LineWidth', 2, 'MarkerSize', 6);
                SessionMeans{j}(k,:) = [median(TempFFValues) std(TempFFValues)];
                SessionCVs{j}(k,:) = [iqr(TempFFValues)/median(TempFFValues)];
                XLabelString{k} = [TempBirdParameters(k).DataLabel, '.', TempBirdParameters(k).Condition, '.', TempBirdParameters(k).Microphone];
            end
            set(gca, 'XTick', 1:1:length(TempBirdParameters), 'XTickLabel', XLabelString, 'XTickLabelRotation', 45);
            title([BirdNames{i}, '.SyllLabel.', TempBirdParameters(1).MotifLabels(j), '.AcrossDay.MeanAmplitude.', DifferentMicrophones{MicrophoneType}]);
            ylabel('Mean FF (Hz)');
            %[p, anovatabl, stats] = anova1(AllFeatValues(:,2), AllFeatValues(:,1), 'off');
            %figure;
            %multcompare(stats);
            %title([BirdNames{i}, '.SyllLabel.', TempBirdParameters(1).MotifLabels(j), '.AcrossDay.', TempBirdParameters(1).SAPFeat_FieldNames{FeatureToPlot_ColNo}]);
        end
        % Now to plot all syllables for a bird together
        CVFig = figure;
        hold on;
        MeanFig = figure;
        hold on;
        for j = 1:length(SessionMeans),
            figure(MeanFig);
            errorbar(SessionMeans{j}(:,1), SessionMeans{j}(:,2), [Colours(mod(j,length(Colours)) + 1), 's-']);
            figure(CVFig);
            plot(SessionCVs{j}(:,1), [Colours(mod(j,length(Colours)) + 1), 's-']);
        end
        figure(MeanFig);
        legend(mat2cell(IndividualBirds(i).SortedBirdParameters(1).MotifLabels, 1, ones(length(IndividualBirds(i).SortedBirdParameters(1).MotifLabels), 1)));
        set(gca, 'XTick', 1:1:length(TempBirdParameters), 'XTickLabel', XLabelString, 'XTickLabelRotation', 45);
        title([BirdNames{i}, '.AllSyllables.', TempBirdParameters(1).MotifLabels(j), '.AcrossDay.MeanAmplitude.', DifferentMicrophones{MicrophoneType}]);
        ylabel('Mean FF (Hz)');
        
        figure(CVFig);
        Temp = axis;
        plot(Temp(1:2), ones(1,2)*0.03, 'k--');
        plot(Temp(1:2), ones(1,2)*0.025, 'k--');
        legend(mat2cell(IndividualBirds(i).SortedBirdParameters(1).MotifLabels, 1, ones(length(IndividualBirds(i).SortedBirdParameters(1).MotifLabels), 1)));
        set(gca, 'XTick', 1:1:length(TempBirdParameters), 'XTickLabel', XLabelString, 'XTickLabelRotation', 45);
        title([BirdNames{i}, '.AllSyllables.', TempBirdParameters(1).MotifLabels(j), '.AcrossDay.MeanAmplitude.', DifferentMicrophones{MicrophoneType}]);
        ylabel('CV of FF');
    end
end
disp('Finished');