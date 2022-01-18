function [DirBout_Stats, UndirBout_Stats, DirUndirBout_Stats] = Harini_ProcessINNos_PlotTimesBoutTypes(IndividualBirds, BirdNames, BirdOption, FeatureToPlot_ColNo, PlotType, DifferentMicrophones)

% ====== IN Number processing =============================================
% THis is to check all aspects related to IN number in all the birds.

% =========================================================================

FinalFigureDir = '/home/raghav/StudentRelated/Harini/Manuscript/results';
MinTrialNo = 3;

for i = 1:length(IndividualBirds),
    disp(['Bird #', num2str(i)]);
    ConditionColumnIndex = find(cellfun(@length,strfind(IndividualBirds(i).BoutStatisticsColumnNames, 'Condition')));
    BoutIndexColumnIndex = find(cellfun(@length,strfind(IndividualBirds(i).BoutStatisticsColumnNames, 'BoutIndex')));
    
    INLabels = IndividualBirds(i).SortedBirdParameters(1).INLabels;
    MotifLabels = IndividualBirds(i).SortedBirdParameters(1).MotifLabels;
    CommonMotifs = IndividualBirds(i).SortedBirdParameters(1).CommonMotifs;
    
    ValidSongBouts = find((IndividualBirds(i).Bouts(:,7) == 1) & (IndividualBirds(i).Bouts(:,8) > 0) & (IndividualBirds(i).Bouts(:,9) > 1));
    disp(['% of song bouts that have enough data in front and at the back = ', num2str(length(ValidSongBouts)), '/', num2str(length(find(IndividualBirds(i).Bouts(:,7) == 1))), '(', num2str(100 * length(ValidSongBouts)/length(find(IndividualBirds(i).Bouts(:,7) == 1))), '%)']);
    for j = min(IndividualBirds(i).BoutStatistics(:,ConditionColumnIndex)):max(IndividualBirds(i).BoutStatistics(:,ConditionColumnIndex)),
        Bouts{i}{j} = IndividualBirds(i).BoutStatistics(find(IndividualBirds(i).BoutStatistics(:,ConditionColumnIndex) == j), BoutIndexColumnIndex);
        Bouts{i}{j} = intersect(Bouts{i}{j}, ValidSongBouts);
        
        DirBouts{i}{j} = intersect(IndividualBirds(i).BoutStatistics(strmatch('D', IndividualBirds(i).BoutCategorisation, 'exact'), BoutIndexColumnIndex), Bouts{i}{j});
        DirUndirBouts{i}{j} = intersect(IndividualBirds(i).BoutStatistics(strmatch('DUN', IndividualBirds(i).BoutCategorisation, 'exact'), BoutIndexColumnIndex), Bouts{i}{j});
        if (j ~= 6)
            UndirBouts{i}{j} = intersect(IndividualBirds(i).BoutStatistics(strmatch('UN', IndividualBirds(i).BoutCategorisation, 'exact'), BoutIndexColumnIndex), Bouts{i}{j});
        else
            UndirBouts{i}{j} = Bouts{i}{j}; % since these bouts are anyway only in the absence of the female.
        end
        % Now for each of these bouts, get the # of INs, the last set of INs
        % with < 500ms gap and the last set of INs independent of the gap
        % In addition, for each bout, also get the total number of complete
        % motifs, partial motifs, all motif durations, first motif duration
        % Also, get the date for that particular dataset
    
        % Get Bout times
        if (~isempty(DirBouts{i}{j}))
            DirBoutTimes{i}{j} = IndividualBirds(i).BoutTimesRelativeToStart(DirBouts{i}{j}) * 60; % in minutes
        end
        if (~isempty(DirUndirBouts{i}{j}))
            DirUndirBoutTimes{i}{j} = IndividualBirds(i).BoutTimesRelativeToStart(DirUndirBouts{i}{j}) * 60; % in minutes
        end
        
        if (~isempty(UndirBouts{i}{j}))
            UndirBoutTimes{i}{j} = IndividualBirds(i).BoutTimesRelativeToStart(UndirBouts{i}{j}) * 60; % in minutes
        end
        
        if (~isempty(DirBouts{i}{j}))
            [DirBout_Stats{i}{j}] = Harini_CalculateBoutINMotifNumbers_Duration(IndividualBirds(i), DirBouts{i}{j}, MotifLabels, INLabels, CommonMotifs);
        else
            DirBout_Stats{i}{j} = [];
        end
        if (~isempty(DirUndirBouts{i}{j}))
            [DirUndirBout_Stats{i}{j}] = Harini_CalculateBoutINMotifNumbers_Duration(IndividualBirds(i), DirUndirBouts{i}{j}, MotifLabels, INLabels, CommonMotifs);
        else
            DirUndirBout_Stats{i}{j} = [];
        end
        if (~isempty(UndirBouts{i}{j}))
            [UndirBout_Stats{i}{j}] = Harini_CalculateBoutINMotifNumbers_Duration(IndividualBirds(i), UndirBouts{i}{j}, MotifLabels, INLabels, CommonMotifs);
        else
            UndirBout_Stats{i}{j} = [];
        end
    end
end


% Now to do all the stats and plots

% First for all IN nos. - one set with full IN nos and another with only IN
% numbers that are separated by <= 500ms
Distances = [0 20 60 110 165 200];

MinTrialNo = 3;
StructString = 'TotalINNumber_500ms';
LabelString = '# of INs';
ParametricOrNot = 0;
[OutputStats(1)] = Harini_DoStats_AndPlots_LM(DirBout_Stats, DirUndirBout_Stats, UndirBout_Stats, StructString, MinTrialNo, Distances, LabelString, BirdNames, ParametricOrNot);

MinTrialNo = 3;
StructString = 'CompleteMotifNumber';
LabelString = '# of complete motifs';
ParametricOrNot = 0;
[OutputStats(2)] = Harini_DoStats_AndPlots_LM(DirBout_Stats, DirUndirBout_Stats, UndirBout_Stats, StructString, MinTrialNo, Distances, LabelString, BirdNames, ParametricOrNot);

MinTrialNo = 3;
StructString = 'FirstMotifDuration';
LabelString = 'First motif duration (msec)';
ParametricOrNot = 1;
[OutputStats(3)] = Harini_DoStats_AndPlots_LM(DirBout_Stats, DirUndirBout_Stats, UndirBout_Stats, StructString, MinTrialNo, Distances, LabelString, BirdNames, ParametricOrNot);

MinTrialNo = 10;
StructString = 'LogAmplitude';
LabelString = 'Log Amplitude (dB)';
ParametricOrNot = 1;
[Amp_OutputStats(1)] = Harini_DoStats_AndPlots_LogAmp_LM(DirBout_Stats, DirUndirBout_Stats, UndirBout_Stats, StructString, MinTrialNo, Distances, LabelString, BirdNames, ParametricOrNot);

ResponseCategories = {'INTERACTION' 'DISTANCE' 'SONG TYPE' 'DISTANCE + SONG TYPE' 'NO CHANGE'};
NaNBirds = union(find(isnan(OutputStats(1).LMFit_ModelPValue)), find(isnan(OutputStats(1).LMFit_DistanceContext_Interactions(:,2))));

for i = 1:length(OutputStats),
    for j = 1:length(ResponseCategories),
        switch (ResponseCategories{j})
            case 'INTERACTION'
                Responses{i}{j} = find((OutputStats(i).LMFit_ModelPValue' < 0.05) & (OutputStats(i).LMFit_DistanceContext_Interactions(:,2) < 0.05));
                
            case 'DISTANCE'
                Responses{i}{j} = find((OutputStats(i).LMFit_ModelPValue' < 0.05) & (OutputStats(i).LMFit_Distance(:,2) < 0.05) & (OutputStats(i).LMFit_Context(:,2) >= 0.05) & (OutputStats(i).LMFit_DistanceContext_Interactions(:,2) >= 0.05));
                
            case 'SONG TYPE'
                Responses{i}{j} = find((OutputStats(i).LMFit_ModelPValue' < 0.05) & (OutputStats(i).LMFit_Distance(:,2) >= 0.05) & (OutputStats(i).LMFit_Context(:,2) < 0.05) & (OutputStats(i).LMFit_DistanceContext_Interactions(:,2) >= 0.05));
    
            case 'DISTANCE + SONG TYPE'
                Responses{i}{j} = find((OutputStats(i).LMFit_ModelPValue' < 0.05) & (OutputStats(i).LMFit_Distance(:,2) < 0.05) & (OutputStats(i).LMFit_Context(:,2) < 0.05) & (OutputStats(i).LMFit_DistanceContext_Interactions(:,2) >= 0.05));
                
            case 'NO CHANGE'
                Responses{i}{j} = union(find((OutputStats(i).LMFit_ModelPValue' >= 0.05)), find((OutputStats(i).LMFit_Distance(:,2) >= 0.05) & (OutputStats(i).LMFit_Context(:,2) >= 0.05) & (OutputStats(i).LMFit_DistanceContext_Interactions(:,2) >= 0.05)));
        end
        Responses{i}{j} = setdiff(Responses{i}{j}, NaNBirds);
    end
end

for i = 1:length(Responses),
    NumResponses(:,i) = cellfun(@length, Responses{i});
end

TempAllBirds = [];

for i = 1:length(OutputStats),
    TempResponses = [];
    for j = 1:length(Responses{i}),
        TempResponses = [TempResponses; Responses{i}{j}(:)];
    end
    TempAllBirds(:,i) = TempResponses;
end

BirdsToPlot = [8 8 3; 2 2 1; 11 5 11]; % examples to plot for # of INs, # of motifs and FMD

AllBirds = [];
for i = 1:length(OutputStats),
    AllBirds(:,i) = setdiff(TempAllBirds(:,i), BirdsToPlot(:,i), 'rows', 'stable');
end

OutputStatsYLabels = {'# of INs' '# of motifs/bout', 'FMD (msec)'};

% Now to plot all bout level features across all birds in one panel
figure;
set(gcf, 'Color', 'w');
p = panel();
p.pack(size(AllBirds,1),length(OutputStats));

p.de.margin = 12;

for i = 1:length(OutputStats),
    Index = 0;
    for j = (AllBirds(:,i))',
        Index = Index + 1;
        p(Index,i).select();
        hold on;
        NanDistances = find(~isnan(OutputStats(i).DistanceDirMeans(j,1:5)));
        errorbar(Distances(1:5), OutputStats(i).DistanceDirMeans(j,1:5), OutputStats(i).DistanceDirSEMs(j,1:5), 'ro', 'MarkerFaceColor', 'r', 'LineWidth', 0.5, 'MarkerSize', 6);
        % plot(Distances(1:5), OutputStats(i).DistanceDirMeans(j,1:5), 'rs-', 'LineWidth', 1);
        % plot(repmat(Distances(1:5)', 1, 2)', OutputStats(i).DistanceDirCIs{j}(1:5,:)', 'r')
        plot(Distances(NanDistances), polyval([OutputStats(i).LMFit_Distance(j,1) OutputStats(i).LMFit_Intercept(j,1)], Distances(NanDistances)), 'r', 'LineWidth', 1);
        
        errorbar(Distances(1:5), OutputStats(i).DistanceUnDirMeans(j,1:5), OutputStats(i).DistanceUnDirSEMs(j,1:5), 'bo', 'MarkerFaceColor', 'b', 'LineWidth', 0.5, 'MarkerSize', 6);
        NanDistances = find(~isnan(OutputStats(i).DistanceUnDirMeans(j,1:5)));
        plot(Distances(NanDistances), polyval([OutputStats(i).LMFit_Distance(j,1)+OutputStats(i).LMFit_DistanceContext_Interactions(j,1) OutputStats(i).LMFit_Intercept(j,1)+OutputStats(i).LMFit_Context(j,1)], Distances(NanDistances)), 'b', 'LineWidth', 1);
        
        % plot(Distances(1:5), OutputStats(i).DistanceUnDirMeans(j,1:5), 'bs-', 'LineWidth', 1);
        % plot(repmat(Distances(1:5)', 1, 2)', OutputStats(i).DistanceUnDirCIs{j}(1:5,:)', 'b');
        % plot(repmat(Distances(1:5)',1,2), repmat(OutputStats(i).DistanceUnDirCIs{j}(end,:), length(Distances)-1, 1), 'k--');
        plot([Distances(1) Distances(5)], ones(1,2)*OutputStats(i).DistanceUnDirMeans(j,end) + OutputStats(i).DistanceUnDirSEMs(j,end), 'k--');
        plot([Distances(1) Distances(5)], ones(1,2)*OutputStats(i).DistanceUnDirMeans(j,end) - OutputStats(i).DistanceUnDirSEMs(j,end), 'k--');
        
        axis tight;
        Temp = axis;
        Temp = [-5 185 0.99*Temp(3) 1.05*Temp(4)];
        axis(Temp);
        
        % Write the model p-value and the adjusted r^2 on the plot
        text(0, Temp(4), ['p = ', num2str(OutputStats(i).LMFit_ModelPValue(j)), '; Adj. r^{2} = ', num2str(OutputStats(i).Models{j}.mdl.Rsquared.Adjusted)], 'FontSize', 6)
        
        
        title(['Bird #', num2str(j)], 'FontSize', 10);
        if (ceil(Temp(3)) == floor(Temp(4)))
            set(gca, 'YTick', [ceil(Temp(3)*10)/10 floor(Temp(4)*10)/10]);
        else
            set(gca, 'YTick', [ceil(Temp(3)) floor(Temp(4))]);
        end
        
%        ylabel(OutputStatsYLabels{i});
        if (i > 1)
            p(Index,i).marginleft = 20;
        end
        
        set(gca, 'XTick', [0 75 150], 'XTickLabel', [{'0' '75' '150'}]);
        if (j == AllBirds(end,i))
            xlabel('Distance from female (cm)');
        end
    end
end

p.fontsize = 12;
% p.fontname = 'Times';
%p.margintop = 10;
%p.marginleft = 15;
%p.marginright = 5;
%p.marginbottom = 25;

set(gcf, 'ReSize', 'off');
set(gcf, 'Color', 'w');
set(gcf, 'Units', 'inches');
set(gcf, 'Position', [7 0.5 8.5 11]);
set(gcf, 'PaperPositionMode', 'auto');
print(fullfile(FinalFigureDir, ['AllBirdsExceptExamples.eps']), '-depsc2', '-r600');
print(fullfile(FinalFigureDir, ['AllBirdsExceptExamples.png']), '-dpng', '-r600');

% Now to plot all bout level features across all birds in one panel
figure;
set(gcf, 'Color', 'w');
p = panel();
p.pack({7/8 1/8});
p(1).pack(size(BirdsToPlot,1),length(OutputStats));
p(2).pack(1, size(BirdsToPlot,1));

p.de.margin = 20;

for i = 1:length(OutputStats),
    Index = 0;
    for j = (BirdsToPlot(:,i))',
        Index = Index + 1;
        p(1,Index,i).select();
        hold on;
        NanDistances = find(~isnan(OutputStats(i).DistanceDirMeans(j,1:5)));
        errorbar(Distances(1:5), OutputStats(i).DistanceDirMeans(j,1:5), OutputStats(i).DistanceDirSEMs(j,1:5), 'ro', 'MarkerFaceColor', 'r', 'LineWidth', 0.5, 'MarkerSize', 6);
        % plot(Distances(1:5), OutputStats(i).DistanceDirMeans(j,1:5), 'rs-', 'LineWidth', 1);
        % plot(repmat(Distances(1:5)', 1, 2)', OutputStats(i).DistanceDirCIs{j}(1:5,:)', 'r')
        plot(Distances(NanDistances), polyval([OutputStats(i).LMFit_Distance(j,1) OutputStats(i).LMFit_Intercept(j,1)], Distances(NanDistances)), 'r', 'LineWidth', 1);
        
        errorbar(Distances(1:5), OutputStats(i).DistanceUnDirMeans(j,1:5), OutputStats(i).DistanceUnDirSEMs(j,1:5), 'bo', 'MarkerFaceColor', 'b', 'LineWidth', 0.5, 'MarkerSize', 6);
        NanDistances = find(~isnan(OutputStats(i).DistanceUnDirMeans(j,1:5)));
        plot(Distances(NanDistances), polyval([OutputStats(i).LMFit_Distance(j,1)+OutputStats(i).LMFit_DistanceContext_Interactions(j,1) OutputStats(i).LMFit_Intercept(j,1)+OutputStats(i).LMFit_Context(j,1)], Distances(NanDistances)), 'b', 'LineWidth', 1);
        
        % plot(Distances(1:5), OutputStats(i).DistanceUnDirMeans(j,1:5), 'bs-', 'LineWidth', 1);
        % plot(repmat(Distances(1:5)', 1, 2)', OutputStats(i).DistanceUnDirCIs{j}(1:5,:)', 'b');
        % plot(repmat(Distances(1:5)',1,2), repmat(OutputStats(i).DistanceUnDirCIs{j}(end,:), length(Distances)-1, 1), 'k--');
        plot([Distances(1) Distances(5)], ones(1,2)*OutputStats(i).DistanceUnDirMeans(j,end) + OutputStats(i).DistanceUnDirSEMs(j,end), 'k--');
        plot([Distances(1) Distances(5)], ones(1,2)*OutputStats(i).DistanceUnDirMeans(j,end) - OutputStats(i).DistanceUnDirSEMs(j,end), 'k--');
        
        axis tight;
        Temp = axis;
        Temp = [-5 185 0.99*Temp(3) 1.05*Temp(4)];
        axis(Temp);
        
        % Write the model p-value and the adjusted r^2 on the plot
        text(0, Temp(4), ['p = ', num2str(OutputStats(i).LMFit_ModelPValue(j)), '; Adj. r^{2} = ', num2str(OutputStats(i).Models{j}.mdl.Rsquared.Adjusted)], 'FontSize', 8)
        
        
        title(['Bird #', num2str(j)]);
        if (ceil(Temp(3)) == floor(Temp(4)))
            set(gca, 'YTick', [ceil(Temp(3)*10)/10 floor(Temp(4)*10)/10]);
        else
            set(gca, 'YTick', [ceil(Temp(3)) floor(Temp(4))]);
        end
        
        ylabel(OutputStatsYLabels{i});
        if (i > 1)
            p(1,Index,i).marginleft = 20;
        end
        
        set(gca, 'XTick', [0 75 150], 'XTickLabel', [{'0' '75' '150'}]);
        if (j == BirdsToPlot(end,i))
            xlabel('Distance from female (cm)');
        end
    end
end

for i = 1:size(BirdsToPlot,1),
    p(2,1,i).select();
    TempBar = bar(NumResponses(:,i));
    set(TempBar, 'FaceColor', 'none', 'EdgeColor', 'k', 'LineWidth', 1.5);
    set(gca, 'Box', 'off');
    if (i == 1)
        ylabel('# of birds');
    end
    set(gca, 'XTick', 1:1:length(ResponseCategories), 'XTickLabel', ResponseCategories, 'XTickLabelRotation', 45);
    axis([0.5 (length(ResponseCategories)+0.5) 0 max(NumResponses(:))+1]);
end

p.fontsize = 12;
% p.fontname = 'Times';
p.margintop = 10;
p.marginleft = 15;
p.marginright = 5;
p.marginbottom = 45;

% set(gcf, 'ReSize', 'off');
set(gcf, 'Color', 'w');
set(gcf, 'Units', 'inches');
set(gcf, 'Position', [7 0.5 8 11]);
set(gcf, 'PaperPositionMode', 'auto');
print(fullfile(FinalFigureDir, ['Fig.3.ResponseExamples_Counts.eps']), '-depsc2', '-r600');
print(fullfile(FinalFigureDir, ['Fig.3.ResponseExamples_Counts.png']), '-dpng', '-r600');

ActualResponses = Responses;

% Now to plot the slopes for each of the categories.
figure;
p = panel();
p.pack(2,1);
p(1,1).select();
hold on;
Symbols = 's^*';
for i = 1:length(OutputStats),
    for j = 1:length(ResponseCategories),
        Slopes{j} = [];
    end
end

for i = 1:length(OutputStats),
    for j = 1:length(ResponseCategories)-1,
        if (isempty(Responses{i}{j}))
            continue;
        end
        Changes = setdiff(Responses{i}{j}, BirdsToPlot(:,i));
        plot(repmat([1.2 1.8] + (j-1)*2, length(Changes), 1)', [OutputStats(i).LMFit_ZScore_Distance(Changes,1) (OutputStats(i).LMFit_ZScore_Distance(Changes,1) + OutputStats(i).LMFit_ZScore_DistanceContext_Interactions(Changes,1))]', 'k');
        plot(repmat([1.2] + (j-1)*2, length(Changes), 1)', [OutputStats(i).LMFit_ZScore_Distance(Changes,1)]', ['r', Symbols(i)], 'MarkerSize', 6);
        plot(repmat([1.8] + (j-1)*2, length(Changes), 1)', [(OutputStats(i).LMFit_ZScore_Distance(Changes,1) + OutputStats(i).LMFit_ZScore_DistanceContext_Interactions(Changes,1))]', ['b', Symbols(i)], 'MarkerSize', 6);

        Changes = intersect(Responses{i}{j}, BirdsToPlot(:,i));
        plot(repmat([1.2 1.8] + (j-1)*2, length(Changes), 1)', [OutputStats(i).LMFit_ZScore_Distance(Changes,1) (OutputStats(i).LMFit_ZScore_Distance(Changes,1) + OutputStats(i).LMFit_ZScore_DistanceContext_Interactions(Changes,1))]', 'k');
        plot(repmat([1.2] + (j-1)*2, length(Changes), 1)', [OutputStats(i).LMFit_ZScore_Distance(Changes,1)]', ['r', Symbols(i)], 'MarkerSize', 9);
        plot(repmat([1.8] + (j-1)*2, length(Changes), 1)', [(OutputStats(i).LMFit_ZScore_Distance(Changes,1) + OutputStats(i).LMFit_ZScore_DistanceContext_Interactions(Changes,1))]', ['b', Symbols(i)], 'MarkerSize', 9);
    
        Slopes{j} = [Slopes{j}; [OutputStats(i).LMFit_ZScore_Distance(Responses{i}{j},1) (OutputStats(i).LMFit_ZScore_Distance(Responses{i}{j},1) + OutputStats(i).LMFit_ZScore_DistanceContext_Interactions(Responses{i}{j},1))]];
    end
end
for i = 1:length(Slopes),
    if (isempty(Slopes{i}))
        continue;
    end
    ResponseBar(i) = bar([1 2] + (i-1)*2, mean(Slopes{i}));
    set(ResponseBar(i), 'FaceColor', 'none', 'LineWidth', 0.5);
    errorbar([1 2] + (i-1)*2, mean(Slopes{i}), std(Slopes{i})/sqrt(size(Slopes{i},1)), 'ks', 'LineWidth', 0.5, 'MarkerSize', 8);
    [H, ResponseSlope_DiffProb(i)] = ttest(abs(Slopes{i}(:,1)), abs(Slopes{i}(:,2)));
    disp(['Comparing absolute slopes for Dir and undir songs for ', ResponseCategories{i}, ': p = ', num2str(ResponseSlope_DiffProb(i)), ' paired t-test']);
end

axis tight;
Temp = axis;
Temp = [(round(Temp(1)) - 0.5) (round(Temp(2)) + 0.5) 1.05*Temp(3) 1.05*Temp(4)];
axis(Temp);
ylabel('Slope');
sigstar({[1 2]}, ResponseSlope_DiffProb(1));
set(gca, 'XTick', 1:1:Temp(2), 'XTickLabel', repmat({'D' 'UN'}, 1, length(Slopes)));

p(2,1).select();
hold on;
% For interaction and distance birds, plot slopes for L0 > UN and slopes
% for L0 < UN separately for directed songs and undirected songs
DirL0_MoreThan_UN_Slopes = [];
DirUN_MoreThan_L0_Slopes = [];

UnDirL0_MoreThan_UN_Slopes = [];
UnDirUN_MoreThan_L0_Slopes = [];

SlopesToConsider = find(cellfun(@length, strfind(ResponseCategories, 'DISTANCE')));
SlopesToConsider = [1; SlopesToConsider(:)];

for i = 1:length(OutputStats),
    for k = (SlopesToConsider(:))',
        for j = Responses{i}{k}(:)',
            if (OutputStats(i).DistanceDirMeans(j,1) >= OutputStats(i).DistanceUnDirMeans(j,end))
                plot(1, OutputStats(i).LMFit_ZScore_Distance(j,1), ['r', Symbols(i)], 'LineWidth', 1, 'MarkerSize', 6);
                DirL0_MoreThan_UN_Slopes(end+1) = OutputStats(i).LMFit_ZScore_Distance(j,1);

                plot(3, (OutputStats(i).LMFit_ZScore_Distance(j,1) + OutputStats(i).LMFit_ZScore_DistanceContext_Interactions(j,1)), ['b', Symbols(i)], 'LineWidth', 1, 'MarkerSize', 6);
                UnDirL0_MoreThan_UN_Slopes(end+1) = (OutputStats(i).LMFit_ZScore_Distance(j,1) + OutputStats(i).LMFit_ZScore_DistanceContext_Interactions(j,1));
            else
                plot(2, OutputStats(i).LMFit_ZScore_Distance(j,1), ['r', Symbols(i)], 'LineWidth', 1, 'MarkerSize', 6);
                DirUN_MoreThan_L0_Slopes(end+1) = OutputStats(i).LMFit_ZScore_Distance(j,1);

                plot(4, (OutputStats(i).LMFit_ZScore_Distance(j,1) + OutputStats(i).LMFit_ZScore_DistanceContext_Interactions(j,1)), ['b', Symbols(i)], 'LineWidth', 1, 'MarkerSize', 6);
                UnDirUN_MoreThan_L0_Slopes(end+1) = (OutputStats(i).LMFit_ZScore_Distance(j,1) + OutputStats(i).LMFit_ZScore_DistanceContext_Interactions(j,1));
            end
        end
    end
end

errorbar(1.2, mean(DirL0_MoreThan_UN_Slopes), std(DirL0_MoreThan_UN_Slopes)/sqrt(length(DirL0_MoreThan_UN_Slopes)), 'rs-', 'MarkerSize', 8, 'LineWidth', 0.5, 'MarkerFaceColor', 'r');
errorbar(1.8, mean(DirUN_MoreThan_L0_Slopes), std(DirUN_MoreThan_L0_Slopes)/sqrt(length(DirUN_MoreThan_L0_Slopes)), 'rs-', 'MarkerSize', 8, 'LineWidth', 0.5, 'MarkerFaceColor', 'r');

errorbar(3.2, mean(UnDirL0_MoreThan_UN_Slopes), std(UnDirL0_MoreThan_UN_Slopes)/sqrt(length(UnDirL0_MoreThan_UN_Slopes)), 'bs-', 'MarkerSize', 8, 'LineWidth', 0.5, 'MarkerFaceColor', 'b');
errorbar(3.8, mean(UnDirUN_MoreThan_L0_Slopes), std(UnDirUN_MoreThan_L0_Slopes)/sqrt(length(UnDirUN_MoreThan_L0_Slopes)), 'bs-', 'MarkerSize', 8, 'LineWidth', 0.5, 'MarkerFaceColor', 'b');

disp(['Dir L0 /geq UN: fraction of dir song slopes < 0 = ', num2str(length(find(DirL0_MoreThan_UN_Slopes < 0))), '/', num2str(length(DirL0_MoreThan_UN_Slopes)), '(', num2str(100*length(find(DirL0_MoreThan_UN_Slopes < 0))/length(DirL0_MoreThan_UN_Slopes)), '%)']);
disp(['Dir L0 < UN: fraction of dir song slopes < 0 = ', num2str(length(find(DirUN_MoreThan_L0_Slopes < 0))), '/', num2str(length(DirUN_MoreThan_L0_Slopes)), '(', num2str(100*length(find(DirUN_MoreThan_L0_Slopes < 0))/length(DirUN_MoreThan_L0_Slopes)), '%)']);

disp(['Undir L0 /geq UN: fraction of dir song slopes < 0 = ', num2str(length(find(UnDirL0_MoreThan_UN_Slopes < 0))), '/', num2str(length(UnDirL0_MoreThan_UN_Slopes)), '(', num2str(100*length(find(UnDirL0_MoreThan_UN_Slopes < 0))/length(UnDirL0_MoreThan_UN_Slopes)), '%)']);
disp(['Undir L0 < UN: fraction of dir song slopes < 0 = ', num2str(length(find(UnDirUN_MoreThan_L0_Slopes < 0))), '/', num2str(length(UnDirUN_MoreThan_L0_Slopes)), '(', num2str(100*length(find(UnDirUN_MoreThan_L0_Slopes < 0))/length(UnDirUN_MoreThan_L0_Slopes)), '%)']);

axis tight;
Temp = axis;
Temp = [0.75 4.25 1.05*Temp(3) 1.05*Temp(4)];
axis(Temp);
plot([0.75 4.25], [0 0], 'k--');

set(gca, 'XTick', [1 2 3 4], 'XTickLabel', {'Dir L0 \geq UN' 'Dir L0 < UN' 'Dir L0 \geq UN' 'Dir L0 < UN'}, 'XTickLabelRotation', 45);
ylabel('Slope');

p.fontsize = 12;
p.de.margin = 20;
p.marginbottom = 25;
set(gcf, 'Color', 'w');
set(gcf, 'Units', 'inches');
set(gcf, 'Position', [7 0.5 4 7]);
set(gcf, 'PaperPositionMode', 'auto');

print(fullfile(FinalFigureDir, ['ResponseSlopes.eps']), '-depsc2', '-r600');
print(fullfile(FinalFigureDir, ['ResponseSlopes.png']), '-dpng', '-r600');


% 
% 
% MinTrialNo = 3;
% StructString = 'AllMotifDuration';
% LabelString = 'All motif durations (msec)';
% ParametricOrNot = 1;
% % [TotalINNumber_L0Dir_Undir_PValue] = Harini_DoStats_AndPlots(DirBout_Stats, DirUndirBout_Stats, UndirBout_Stats, StructString, MinTrialNo, Distances, LabelString, BirdNames, ParametricOrNot);
% 
% MinTrialNo = 3;
% StructString = 'Bout_FFCV';
% LabelString = 'CV of FF';
% ParametricOrNot = 1;
% [OutputStats(4)] = Harini_DoStats_AndPlots(DirBout_Stats, DirUndirBout_Stats, UndirBout_Stats, StructString, MinTrialNo, Distances, LabelString, BirdNames, ParametricOrNot);
% 
% % For FF and amplitude, min trial no corresponds to minimum # of syllables
% % so we will change that to 15
% MinTrialNo = 15;
% 
% % StructString = 'FF';
% % LabelString = 'FF (Hz)';
% % ParametricOrNot = 1;
% % [OutputStats(4)] = Harini_DoStats_AndPlots_FF(DirBout_Stats, DirUndirBout_Stats, UndirBout_Stats, StructString, MinTrialNo, Distances, LabelString, BirdNames, ParametricOrNot);
% 
% StructString = 'LogAmplitude';
% LabelString = 'Log Amplitude (dB)';
% ParametricOrNot = 1;
% [OutputStats(5)] = Harini_DoStats_AndPlots_FF(DirBout_Stats, DirUndirBout_Stats, UndirBout_Stats, StructString, MinTrialNo, Distances, LabelString, BirdNames, ParametricOrNot);
% 
% % plot all the different correlation co-efficients and their significances
% % for each of the features
% figure;
% hold on;
% for i = 1:length(OutputStats),
%     Indices = find(OutputStats(i).Distance_DirSong_Corr(:,2) < 0.05);
%     plot(ones(size(Indices))*(i-0.15), OutputStats(i).Distance_DirSong_Corr(Indices,1), 'ro', 'MarkerFaceColor', 'r');
%     
%     Indices = find((OutputStats(i).Distance_DirSong_Corr(:,2) >= 0.05) & (cell2mat(OutputStats(i).KWP_DirSong)' < 0.05));
%     plot(ones(size(Indices))*(i-0.15), OutputStats(i).Distance_DirSong_Corr(Indices,1), 'ko', 'MarkerFaceColor', 'k');
% 
%     Indices = find((OutputStats(i).Distance_DirSong_Corr(:,2) >= 0.05) & (cell2mat(OutputStats(i).KWP_DirSong)' >= 0.05));
%     plot(ones(size(Indices))*(i-0.15), OutputStats(i).Distance_DirSong_Corr(Indices,1), 'ko');
%     
%     Indices = find(OutputStats(i).Distance_UnDirSong_Corr(:,2) < 0.05);
%     plot(ones(size(Indices))*(i+0.15), OutputStats(i).Distance_UnDirSong_Corr(Indices,1), 'bo', 'MarkerFaceColor', 'b');
%     
%     Indices = find((OutputStats(i).Distance_UnDirSong_Corr(:,2) >= 0.05) & (cell2mat(OutputStats(i).KWP_UnDirSong)' < 0.05));
%     plot(ones(size(Indices))*(i+0.15), OutputStats(i).Distance_UnDirSong_Corr(Indices,1), 'ko', 'MarkerFaceColor', 'k');
% 
%     Indices = find((OutputStats(i).Distance_UnDirSong_Corr(:,2) >= 0.05) & (cell2mat(OutputStats(i).KWP_UnDirSong)' >= 0.05));
%     plot(ones(size(Indices))*(i+0.15), OutputStats(i).Distance_UnDirSong_Corr(Indices,1), 'ko');
% end
% ylabel('Correlation co-efficient');
% set(gca, 'XTick', [1 2 3], 'XTickLabel', {'# of INs' '# of motifs/bout', 'FMD'});
% Temp = axis;
% Temp = [0.5 3.5 -1 1];
% axis(Temp);
% plot([0.5 3.5], [0 0], 'k--');
% 
% % Now to plot the correlation co-efficient with L0/UN
% % First considering only significant correlations
% Symbols = 's^*+';
% % First for directed
% figure;
% hold on;
% for i = 1:2,
%     Indices = find(OutputStats(i).Distance_DirSong_Corr(:,2) < 0.05);
%     plot(OutputStats(i).DistanceDirMeans(Indices,1)/100, OutputStats(i).Distance_DirSong_Corr(Indices,1), ['r', Symbols(i)], 'LineWidth', 1, 'MarkerFaceColor', 'r', 'MarkerSize', 8);
% end
% 
% % For FMD and FF, r-values are flipped and L0/UN is also flipped to be
% % consistent with the expectations for changes between L0 and UN across all
% % features
% 
% for i = 3:3,
%     Indices = find(OutputStats(i).Distance_DirSong_Corr(:,2) < 0.05);
%     plot(1./(OutputStats(i).DistanceDirMeans(Indices,1)/100), -OutputStats(i).Distance_DirSong_Corr(Indices,1), ['r', Symbols(i)], 'LineWidth', 1, 'MarkerFaceColor', 'r', 'MarkerSize', 8);
% end
% 
% % First for undirected
% for i = 1:2,
%     Indices = find(OutputStats(i).Distance_UnDirSong_Corr(:,2) < 0.05);
%     plot(OutputStats(i).DistanceDirMeans(Indices,1)/100, OutputStats(i).Distance_UnDirSong_Corr(Indices,1), ['b', Symbols(i)], 'LineWidth', 1, 'MarkerFaceColor', 'b', 'MarkerSize', 8);
% end
% 
% % For FMD and FF, r-values are flipped and L0/UN is also flipped to be
% % consistent with the expectations for changes between L0 and UN across all
% % features
% 
% for i = 3:3,
%     Indices = find(OutputStats(i).Distance_UnDirSong_Corr(:,2) < 0.05);
%     plot(1./(OutputStats(i).DistanceDirMeans(Indices,1)/100), -OutputStats(i).Distance_UnDirSong_Corr(Indices,1), ['b', Symbols(i)], 'LineWidth', 1, 'MarkerFaceColor', 'b', 'MarkerSize', 8);
% end
% axis tight;
% Temp = axis;
% Temp = [0.95*Temp(1) 1.05*Temp(2) -1 1];
% axis(Temp);
% plot(Temp(1:2), [0 0], 'k--');
% plot([1 1], Temp(3:4), 'k--');

MinTrialNo = 3;
StructString = 'LogAmplitude';
LabelString = 'Log Amplitude (dB)';

Harini_AllFeatureStats_Plots_WithoutCI(DirBout_Stats, DirUndirBout_Stats, UndirBout_Stats, StructString, MinTrialNo, Distances, LabelString, BirdNames);

disp('Finished calculating bout stats');