function [] = Harini_PlotINNumber_AllPossibilities(IndividualBirds, BirdNames, BirdOption, FeatureToPlot_ColNo, PlotType, DifferentMicrophones)

% ======= IN number in all possible ways ==================================
% 1. For the different conditions irrespective of whether the songs were
% directed or undirected (based on video scoring data)
% 2. For only the DIRECTED songs
% 3. For only the UNDIRECTED songs
% 4. Have to do stats on all of the above
%   a) Kruskal-wallis to check if any of them are different
%   b) Correlation to see if distance and # of INs is related - spearman's
%   rank correlation 
% =========================================================================

Colours = 'rgbcmk';
Symbols = '+o<sd';
Conditions = [{'L0'} {'L1'} {'L2'} {'L3'} {'L4'} {'UN'}];
MinTrialNo = 5;

for i = 1:length(IndividualBirds),
    NumConditions(i) = length(IndividualBirds(i).ConditionDistances);
end

[MaxNumConditions, MaxNumConditionsIndex] = max(NumConditions);

ConditionDistances = IndividualBirds(MaxNumConditionsIndex).ConditionDistances;

% First collect data for all 3 types of bouts listed above
for i = 1:length(IndividualBirds),
    ColumnIndex = find(strcmp('NumINs', IndividualBirds(i).BoutStatisticsColumnNames));
    ConditionColumnIndex = find(strcmp('Condition', IndividualBirds(i).BoutStatisticsColumnNames));
    
    IndividualBird_AllBouts{i} = [];
    IndividualBird_DirectedBouts{i} = [];
    IndividualBird_UnDirectedBouts{i} = [];
    
    for j = 1:length(IndividualBirds(i).Conditions),
        % First find all the bouts that correspond to the given condition 
        ConditionIndices = find(IndividualBirds(i).BoutStatistics(:,ConditionColumnIndex) == j);
        if (length(ConditionIndices) >= MinTrialNo)
            % Add to all bouts
            eval(['AllBouts_', IndividualBirds(i).Conditions{j}, '(', num2str(i), ') = ', num2str(nanmean(IndividualBirds(i).BoutStatistics(ConditionIndices,ColumnIndex))), ';']);
            eval(['AllBouts_SEM_', IndividualBirds(i).Conditions{j}, '(', num2str(i), ') = ', num2str(nanstd(IndividualBirds(i).BoutStatistics(ConditionIndices,ColumnIndex))/sqrt(length(find(~isnan(IndividualBirds(i).BoutStatistics(ConditionIndices,ColumnIndex)))))), ';']);
            IndividualBird_AllBouts{i} = [IndividualBird_AllBouts{i}; [IndividualBirds(i).BoutStatistics(ConditionIndices, ColumnIndex) IndividualBirds(i).BoutStatistics(ConditionIndices, ConditionColumnIndex)]];
            
            % Find all directed only bouts
            if (isempty(strfind(IndividualBirds(i).Conditions{j}, 'UN')))
                DirSongIndices = ConditionIndices(find(strcmp('D', IndividualBirds(i).BoutCategorisation(ConditionIndices))));
            else
                DirSongIndices = ConditionIndices;
            end
            
            if (length(DirSongIndices) >= MinTrialNo)
                eval(['DirectedBouts_', IndividualBirds(i).Conditions{j}, '(', num2str(i), ') = ', num2str(nanmean(IndividualBirds(i).BoutStatistics(DirSongIndices,ColumnIndex))), ';']);
                eval(['DirectedBouts_SEM_', IndividualBirds(i).Conditions{j}, '(', num2str(i), ') = ', num2str(nanstd(IndividualBirds(i).BoutStatistics(DirSongIndices,ColumnIndex))/sqrt(length(find(~isnan(IndividualBirds(i).BoutStatistics(DirSongIndices,ColumnIndex)))))), ';']);
                IndividualBird_DirectedBouts{i} = [IndividualBird_DirectedBouts{i}; [IndividualBirds(i).BoutStatistics(DirSongIndices, ColumnIndex) IndividualBirds(i).BoutStatistics(DirSongIndices, ConditionColumnIndex)]];
            else
                eval(['DirectedBouts_', IndividualBirds(i).Conditions{j}, '(', num2str(i), ') = ', num2str(NaN), ';']);
                eval(['DirectedBouts_SEM_', IndividualBirds(i).Conditions{j}, '(', num2str(i), ') = ', num2str(NaN), ';']);
            end
            
            % Find all undirected only bouts
            if (isempty(strfind(IndividualBirds(i).Conditions{j}, 'UN')))
                UnDirSongIndices = ConditionIndices(find(strcmp('UN', IndividualBirds(i).BoutCategorisation(ConditionIndices))));
            else
                UnDirSongIndices = ConditionIndices;
            end
            
            if (length(UnDirSongIndices) >= MinTrialNo)
                eval(['UnDirectedBouts_', IndividualBirds(i).Conditions{j}, '(', num2str(i), ') = ', num2str(nanmean(IndividualBirds(i).BoutStatistics(UnDirSongIndices,ColumnIndex))), ';']);
                eval(['UnDirectedBouts_SEM_', IndividualBirds(i).Conditions{j}, '(', num2str(i), ') = ', num2str(nanstd(IndividualBirds(i).BoutStatistics(UnDirSongIndices,ColumnIndex))/sqrt(length(find(~isnan(IndividualBirds(i).BoutStatistics(UnDirSongIndices,ColumnIndex)))))), ';']);
                IndividualBird_UnDirectedBouts{i} = [IndividualBird_UnDirectedBouts{i}; [IndividualBirds(i).BoutStatistics(UnDirSongIndices, ColumnIndex) IndividualBirds(i).BoutStatistics(UnDirSongIndices, ConditionColumnIndex)]];
            else
                eval(['UnDirectedBouts_', IndividualBirds(i).Conditions{j}, '(', num2str(i), ') = ', num2str(NaN), ';']);
                eval(['UnDirectedBouts_SEM_', IndividualBirds(i).Conditions{j}, '(', num2str(i), ') = ', num2str(NaN), ';']);
            end
        else
            eval(['AllBouts_', IndividualBirds(i).Conditions{j}, '(', num2str(i), ') = ', num2str(NaN), ';']);
            eval(['AllBouts_SEM_', IndividualBirds(i).Conditions{j}, '(', num2str(i), ') = ', num2str(NaN), ';']);
        end
    end
end

% Put all data together
AllBouts = [AllBouts_L0(:) AllBouts_L1(:) AllBouts_L2(:) AllBouts_L3(:) AllBouts_L4(:) AllBouts_UN(:)];
AllBouts_SEM = [AllBouts_SEM_L0(:) AllBouts_SEM_L1(:) AllBouts_SEM_L2(:) AllBouts_SEM_L3(:) AllBouts_SEM_L4(:) AllBouts_SEM_UN(:)];
AllBouts(find(AllBouts == 0)) = NaN;

DirectedBouts = [DirectedBouts_L0(:) DirectedBouts_L1(:) DirectedBouts_L2(:) DirectedBouts_L3(:) DirectedBouts_L4(:) DirectedBouts_UN(:)];
DirectedBouts_SEM = [DirectedBouts_SEM_L0(:) DirectedBouts_SEM_L1(:) DirectedBouts_SEM_L2(:) DirectedBouts_SEM_L3(:) DirectedBouts_SEM_L4(:) DirectedBouts_SEM_UN(:)];
DirectedBouts(find(DirectedBouts == 0)) = NaN;

UnDirectedBouts = [UnDirectedBouts_L0(:) UnDirectedBouts_L1(:) UnDirectedBouts_L2(:) UnDirectedBouts_L3(:) UnDirectedBouts_L4(:) UnDirectedBouts_UN(:)];
UnDirectedBouts_SEM = [UnDirectedBouts_SEM_L0(:) UnDirectedBouts_SEM_L1(:) UnDirectedBouts_SEM_L2(:) UnDirectedBouts_SEM_L3(:) UnDirectedBouts_SEM_L4(:) UnDirectedBouts_SEM_UN(:)];
UnDirectedBouts(find(UnDirectedBouts == 0)) = NaN;


% Now the plots and the statistics

% First the distributions along with the means for all 3 types of bouts in
% 3 subplots of a figure

for i = 1:length(IndividualBirds),
    % Now plot the distributions
    figure;
    set(gcf, 'Color', 'w');
    set(gcf, 'Position', [680 166 1000 850]);
    hold on;
    
    subplot(2,2,1);
    hold on;
    PlotLegend{i} = [];
    UniqueConditions = unique(IndividualBird_AllBouts{i}(:,2));
    
    for j = UniqueConditions(:)',
        ConditionIndices = find(IndividualBird_AllBouts{i}(:,2) == j);
        plot(0:1:max(IndividualBirds(i).BoutStatistics(:, ColumnIndex)), histc(IndividualBird_AllBouts{i}(ConditionIndices, 1), 0:1:max(IndividualBirds(i).BoutStatistics(:,ColumnIndex)))/length(ConditionIndices), [Colours(mod(j, length(Colours)) + 1), 'o-'], 'LineWidth', 2);
        PlotLegend{i}{end+1} = IndividualBirds(i).Conditions{j};
    end
    
    for j = 1:length(IndividualBirds(i).Conditions),
        ConditionIndices = find(IndividualBird_AllBouts{i}(:,2) == j);
        plot(mean(IndividualBird_AllBouts{i}(ConditionIndices,1)), 1 + j*0.05, [Colours(mod(j, length(Colours)) + 1), 's'], 'MarkerSize', 8, 'LineWidth', 1.5);
        plot([(mean(IndividualBird_AllBouts{i}(ConditionIndices, 1)) - std(IndividualBird_AllBouts{i}(ConditionIndices, 1))) (mean(IndividualBird_AllBouts{i}(ConditionIndices, 1)) + std(IndividualBird_AllBouts{i}(ConditionIndices, 1)))], [(1 + j*0.05) (1 + j*0.05)], Colours(mod(j, length(Colours)) + 1), 'LineWidth', 1.5);
    end
    
    Temp = [0 (max(IndividualBirds(i).BoutStatistics(:,ColumnIndex)) + 0.5) 0 1.35];
    axis(Temp);
    legend(PlotLegend{i});
    xlabel('# of INs');
    ylabel('Fraction of bouts');
    title([BirdNames{i}, ': All bouts']);
    
    subplot(2,2,2);
    hold on;
    PlotLegend{i} = [];
    UniqueConditions = unique(IndividualBird_DirectedBouts{i}(:,2));
    for j = UniqueConditions(:)',
        ConditionIndices = find(IndividualBird_DirectedBouts{i}(:,2) == j);
        plot(0:1:max(IndividualBirds(i).BoutStatistics(:, ColumnIndex)), histc(IndividualBird_DirectedBouts{i}(ConditionIndices, 1), 0:1:max(IndividualBirds(i).BoutStatistics(:,ColumnIndex)))/length(ConditionIndices), [Colours(mod(j, length(Colours)) + 1), 'o-'], 'LineWidth', 2);
        PlotLegend{i}{end+1} = IndividualBirds(i).Conditions{j};
    end
    
    for j = 1:length(IndividualBirds(i).Conditions),
        ConditionIndices = find(IndividualBird_DirectedBouts{i}(:,2) == j);
        plot(mean(IndividualBird_DirectedBouts{i}(ConditionIndices,1)), 1 + j*0.05, [Colours(mod(j, length(Colours)) + 1), 's'], 'MarkerSize', 8, 'LineWidth', 1.5);
        plot([(mean(IndividualBird_DirectedBouts{i}(ConditionIndices, 1)) - std(IndividualBird_DirectedBouts{i}(ConditionIndices, 1))) (mean(IndividualBird_DirectedBouts{i}(ConditionIndices, 1)) + std(IndividualBird_DirectedBouts{i}(ConditionIndices, 1)))], [(1 + j*0.05) (1 + j*0.05)], Colours(mod(j, length(Colours)) + 1), 'LineWidth', 1.5);
    end
    
    Temp = [0 (max(IndividualBirds(i).BoutStatistics(:,ColumnIndex)) + 0.5) 0 1.35];
    axis(Temp);
    legend(PlotLegend{i});
    xlabel('# of INs');
    ylabel('Fraction of bouts');
    title([BirdNames{i}, ': Directed bouts']);
    
    subplot(2,2,3);
    hold on;
    PlotLegend{i} = [];
    UniqueConditions = unique(IndividualBird_UnDirectedBouts{i}(:,2));
    for j = UniqueConditions(:)',
        ConditionIndices = find(IndividualBird_UnDirectedBouts{i}(:,2) == j);
        plot(0:1:max(IndividualBirds(i).BoutStatistics(:, ColumnIndex)), histc(IndividualBird_UnDirectedBouts{i}(ConditionIndices, 1), 0:1:max(IndividualBirds(i).BoutStatistics(:,ColumnIndex)))/length(ConditionIndices), [Colours(mod(j, length(Colours)) + 1), 'o-'], 'LineWidth', 2);
        PlotLegend{i}{end+1} = IndividualBirds(i).Conditions{j};
    end
    
    for j = 1:length(IndividualBirds(i).Conditions),
        ConditionIndices = find(IndividualBird_UnDirectedBouts{i}(:,2) == j);
        plot(mean(IndividualBird_UnDirectedBouts{i}(ConditionIndices,1)), 1 + j*0.05, [Colours(mod(j, length(Colours)) + 1), 's'], 'MarkerSize', 8, 'LineWidth', 1.5);
        plot([(mean(IndividualBird_UnDirectedBouts{i}(ConditionIndices, 1)) - std(IndividualBird_UnDirectedBouts{i}(ConditionIndices, 1))) (mean(IndividualBird_UnDirectedBouts{i}(ConditionIndices, 1)) + std(IndividualBird_UnDirectedBouts{i}(ConditionIndices, 1)))], [(1 + j*0.05) (1 + j*0.05)], Colours(mod(j, length(Colours)) + 1), 'LineWidth', 1.5);
    end
    
    Temp = [0 (max(IndividualBirds(i).BoutStatistics(:,ColumnIndex)) + 0.5) 0 1.35];
    axis(Temp);
    legend(PlotLegend{i});
    xlabel('# of INs');
    ylabel('Fraction of bouts');
    title([BirdNames{i}, ': Undirected bouts']);
    
    subplot(2,2,4);
    hold on;
    errorbar(ConditionDistances, AllBouts(i,:), AllBouts_SEM(i,:), 'ko-', 'LineWidth', 2, 'MarkerSize', 7);
    errorbar(ConditionDistances, DirectedBouts(i,:), DirectedBouts_SEM(i,:), 'ro-', 'LineWidth', 2, 'MarkerSize', 7);
    errorbar(ConditionDistances, UnDirectedBouts(i,:), UnDirectedBouts_SEM(i,:), 'bo-', 'LineWidth', 2, 'MarkerSize', 7);
    axis tight;
    Temp = axis;
    Temp = [(ConditionDistances(1) - 1) (ConditionDistances(end) + 1) 0 1.02*Temp(4)];
    axis(Temp);
    legend('All bouts', 'Directed bouts', 'Undirected bouts', 'Location', 'southeast');
    ylabel('Mean # of INs');
    xlabel('Distance of female (cm)');
    title([BirdNames{i}]);
end

% Make group plots for all bouts
AllBouts = [AllBouts_L0(:) AllBouts_L1(:) AllBouts_L2(:) AllBouts_L3(:) AllBouts_L4(:) AllBouts_UN(:)];
AllBouts(find(AllBouts == 0)) = NaN;

% Normalise as % of undirected and also as delta change relative to
% undirected
for i = 1:size(AllBouts,1),
    PercentChange_AllBouts(i,:) = ((AllBouts(i,:) - AllBouts(i,end)) * 100)/AllBouts(i,end);
    DeltaChange_AllBouts(i,:) = (AllBouts(i,:) - AllBouts(i,end));
end

% First for percent change
figure;
for i = 1:size(AllBouts,2),
    AllBout_Bar(i) = bar(ConditionDistances(i), nanmean(PercentChange_AllBouts(:,i)));
    set(AllBout_Bar(i), 'FaceColor', 'none', 'EdgeColor', 'k', 'BarWidth', 10, 'LineWidth', 2);
    hold on;
    errorbar(ConditionDistances(i), nanmean(PercentChange_AllBouts(:,i)), nanstd(PercentChange_AllBouts(:,i))/sqrt(length(find(~isnan(PercentChange_AllBouts(:,i))))), 'ks-', 'MarkerSize', 8, 'LineWidth', 2);
end
plot(repmat(ConditionDistances(:)', size(AllBouts,1), 1)', PercentChange_AllBouts', 'LineWidth', 0.5);
xlabel('Distance of female cage (cm)');
ylabel('% change in mean # of INs relative to Undirected song');
title(['All bouts (n=', num2str(length(BirdNames)), ' birds']);
axis tight;
Temp = axis;
Temp = [(ConditionDistances(1) - 12) (ConditionDistances(end) + 12) 1.02*Temp(3) 1.02*Temp(4)];
axis(Temp);

% Now for delta change
figure;
for i = 1:size(AllBouts,2),
    AllBout_Bar(i) = bar(ConditionDistances(i), nanmean(DeltaChange_AllBouts(:,i)));
    set(AllBout_Bar(i), 'FaceColor', 'none', 'EdgeColor', 'k', 'BarWidth', 10, 'LineWidth', 2);
    hold on;
    errorbar(ConditionDistances(i), nanmean(DeltaChange_AllBouts(:,i)), nanstd(DeltaChange_AllBouts(:,i))/sqrt(length(find(~isnan(DeltaChange_AllBouts(:,i))))), 'ks-', 'MarkerSize', 8, 'LineWidth', 2);
end
plot(repmat(ConditionDistances(:)', size(AllBouts,1), 1)', DeltaChange_AllBouts', 'LineWidth', 0.5);
xlabel('Distance of female cage (cm)');
ylabel('Difference in mean # of INs relative to Undirected song');
title(['All bouts (n=', num2str(length(BirdNames)), ' birds']);
axis tight;
Temp = axis;
Temp = [(ConditionDistances(1) - 12) (ConditionDistances(end) + 12) 1.02*Temp(3) 1.02*Temp(4)];
axis(Temp);

disp('Finished plotting');