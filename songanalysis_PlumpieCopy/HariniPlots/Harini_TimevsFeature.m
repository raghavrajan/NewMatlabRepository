function [] = Harini_TimevsFeature(IndividualBirds, BirdNames, BirdOption, FeatureToPlot_ColNo, PlotType, DifferentMicrophones)

% ================= Time vs. features =====================================
% The idea here is to look at how features of directed song change with
% passage of time since the introduction of the female.

% =========================================================================

Colours = 'rgbcmk';
Symbols = '+o<sd^>*pxvh';
Conditions = [{'L0'} {'L1'} {'L2'} {'L3'} {'L4'} {'UN'}];
MinTrialNo = 5;

% Note: time for each bout relative to female introduction time is in hours

for i = 1:length(IndividualBirds),
    NumConditions(i) = length(IndividualBirds(i).ConditionDistances);
end

[MaxNumConditions, MaxNumConditionsIndex] = max(NumConditions);

ConditionDistances = IndividualBirds(MaxNumConditionsIndex).ConditionDistances;

% Fid = fopen('INNumber_Correlations.txt', 'w');
% fprintf(Fid, 'Bird Name\tPearsons correlation\t\t\t\t\t\tSpearmans correlation\t\t\t\t\t\n');
% fprintf(Fid, '\tAllBouts\t\tDirected bouts\t\tUndirected bouts\t\tAllBouts\t\tDirected bouts\t\tUndirected bouts\t\n');
% fprintf(Fid, '\tr\tp\tr\tp\tr\tp\tr\tp\tr\tp\tr\tp\n');
% 
% % First collect data for all 3 types of bouts listed above
for i = 1:length(IndividualBirds),
    figure;
    hold on;
    ColumnIndex = find(strcmp('BoutLength', IndividualBirds(i).BoutStatisticsColumnNames));
    ConditionColumnIndex = find(strcmp('Condition', IndividualBirds(i).BoutStatisticsColumnNames));
    
    IndividualBird_AllBouts{i} = [];
    IndividualBird_DirectedBouts{i} = [];
    IndividualBird_UnDirectedBouts{i} = [];

    PositiveTimeIndices = find(IndividualBirds(i).BoutStatisticsTimeRelativeToStart > 0);
    [r, p] = corr(IndividualBirds(i).BoutStatisticsTimeRelativeToStart(PositiveTimeIndices)', IndividualBirds(i).BoutStatistics(PositiveTimeIndices, ColumnIndex), 'rows', 'complete');
    disp([BirdNames{i}, ': All conditions : r=', num2str(r), '; p=', num2str(p)]);
    
    for j = 1:length(IndividualBirds(i).Conditions),
        % First find all the bouts that correspond to the given condition 
        ConditionIndices = find(IndividualBirds(i).BoutStatistics(:,ConditionColumnIndex) == j);
        if (length(ConditionIndices) >= MinTrialNo)
            %plot(IndividualBirds(i).BoutStatisticsTimeRelativeToStart(ConditionIndices), IndividualBirds(i).BoutStatistics(ConditionIndices,ColumnIndex)+(rand(length(ConditionIndices),1)/5 - 0.1), [Colours(j), 'o']);
            plot(IndividualBirds(i).BoutStatisticsTimeRelativeToStart(ConditionIndices), IndividualBirds(i).BoutStatistics(ConditionIndices,ColumnIndex), [Colours(j), 'o']);
            % Now to check for linear correlation between bout time and num
            % INs only for bout times that are +ve, meaning they occured
            % after the female was introduced
            PositiveTimeIndices = ConditionIndices(find(IndividualBirds(i).BoutStatisticsTimeRelativeToStart(ConditionIndices) > 0));
            [r, p] = corr(IndividualBirds(i).BoutStatisticsTimeRelativeToStart(PositiveTimeIndices)', IndividualBirds(i).BoutStatistics(PositiveTimeIndices, ColumnIndex), 'rows', 'complete');
            disp([BirdNames{i}, ': ', IndividualBirds(i).Conditions{j}, ': r=', num2str(r), '; p=', num2str(p)]);
        end
    end    
end
disp('Finished plotting');