function [] = Harini_PlotBoutFeatures(IndividualBirds, BirdNames, BirdOption, FeatureToPlot_ColNo, PlotType, DifferentMicrophones)

% == PCA analysis for all bouts ===========================================
% The plan is to take different features and represent each bout by one
% value for each of these features. THis way each bout will be a vector
% with a set of values. Then we can normalize each axis and then plot the
% different bouts along with their distance to see if we see a correlation.
% We also can check how these bout feature values tie up with the video
% scoring which is currently a 3 point score - D (for directed), DUN (for
% ambiguous, but with some characteristics of directed) and UN (for
% completely undirected)
% Features that will be included
% 1) Bout length
% 2) # of INs
% 3) Mean first motif duration
% 4) Mean FF
% 5) CV FF
% 6) Mean FM
% 7) Mean Frequency
% For the last 4 features, we will currently use only the syllables that
% have harmonic portions and are being used for FF measurements.
% =========================================================================

Colours = 'rgbcmk';
ColourMatrix = [1 0 0; 0 1 0; 0 0 1; 0 1 1; 1 0 1; 0 0 0];
Symbols = '+o<sd^>*pxvh';
Conditions = [{'L0'} {'L1'} {'L2'} {'L3'} {'L4'} {'UN'}];
MinTrialNo = 5;

FeaturesToUse = [{'NumINs'} {'BoutLength'} {'FirstMotifDur'} {'AverageSyllFF'} {'CVSyllFF'} {'AverageSyllMeanFM'} {'AverageSyllMeanFreq'}];

for i = 1:length(IndividualBirds),
    for j = 1:length(FeaturesToUse),
        ColumnIndex = find(strcmp(FeaturesToUse{j}, IndividualBirds(i).BoutStatisticsColumnNames));
        FeatVals{i}(:,j) = IndividualBirds(i).BoutStatistics(:,ColumnIndex);
    end
    
    % Find nan rows
    [NanRows, NanCols] = find(isnan(FeatVals{i}));
    
    % Make all values nans in nan rows
    FeatVals{i}(unique(NanRows),:) = NaN;
    
    % Convert this into a z-score
    % FeatVals{i} = (FeatVals{i} - repmat(nanmean(FeatVals{i}), size(FeatVals{i},1),1))./repmat(nanstd(FeatVals{i}), size(FeatVals{i},1),1);
    
    % Next plot the points for different bouts based on the condition for
    % all possible combinations of features
    TotalSubPlots = (length(FeaturesToUse)*(length(FeaturesToUse) - 1))/2;
    figure;
    SubPlotIndex = 1;
    
    for j = 1:length(FeaturesToUse),
        for k = j+1:length(FeaturesToUse),
            subplot(5, ceil(TotalSubPlots/5), SubPlotIndex);
            ConditionColumnIndex = find(strcmp('Condition', IndividualBirds(i).BoutStatisticsColumnNames));
            scatter(FeatVals{i}(:,j), FeatVals{i}(:,k), 5, ColourMatrix(IndividualBirds(i).BoutStatistics(:,ConditionColumnIndex),:));
            xlabel(FeaturesToUse{j});
            ylabel(FeaturesToUse{k});
            SubPlotIndex = SubPlotIndex + 1;
        end
    end

    % Next plot the points for L0 and UN bouts based on the condition for
    % all possible combinations of features
    TotalSubPlots = (length(FeaturesToUse)*(length(FeaturesToUse) - 1))/2;
    figure;
    SubPlotIndex = 1;
    ConditionColumnIndex = find(strcmp('Condition', IndividualBirds(i).BoutStatisticsColumnNames));
    L0_UNSongs = find((IndividualBirds(i).BoutStatistics(:,ConditionColumnIndex) == 1) | (IndividualBirds(i).BoutStatistics(:,ConditionColumnIndex) == 6));
    
    for j = 1:length(FeaturesToUse),
        for k = j+1:length(FeaturesToUse),
            subplot(5, ceil(TotalSubPlots/5), SubPlotIndex);
            scatter(FeatVals{i}(L0_UNSongs,j), FeatVals{i}(L0_UNSongs,k), 5, ColourMatrix(IndividualBirds(i).BoutStatistics(L0_UNSongs,ConditionColumnIndex),:));
            xlabel(FeaturesToUse{j});
            ylabel(FeaturesToUse{k});
            SubPlotIndex = SubPlotIndex + 1;
        end
    end
    
%     % Next plot the mean and 1 * std dev as confidence ellipses for the data
%     figure;
%     hold on;
%     for j = min(IndividualBirds(i).BoutStatistics(:,end)):1:max(IndividualBirds(i).BoutStatistics(:,end)),
%         Songs = find(((IndividualBirds(i).BoutStatistics(:,end)) == j) & (~isnan(PCA_Score(:,1))));
%         PlotConfidenceEllipse(PCA_Score(Songs,1:2), Colours(j), 1)
%         plot(nanmean(PCA_Score(Songs, 1)), nanmean(PCA_Score(Songs, 2)), [Colours(j), '+'], 'MarkerSize', 10, 'LineWidth', 2);
%     end
end

disp('Finished PCA analysis on Bout features');