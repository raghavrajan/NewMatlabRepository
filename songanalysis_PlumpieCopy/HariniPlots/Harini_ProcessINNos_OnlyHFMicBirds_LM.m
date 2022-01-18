function [DirBout_Stats, UndirBout_Stats, DirUndirBout_Stats] = Harini_ProcessINNos_OnlyHFMicBirds_LM(IndividualBirds, BirdNames, BirdOption, FeatureToPlot_ColNo, PlotType, DifferentMicrophones)

% ====== IN Number processing =============================================
% THis is to check all aspects related to IN number in all the birds.

% =========================================================================

FinalFigureDir = '/home/raghav/StudentRelated/Harini/Manuscript/results/OutliersRemoved';
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

MinTrialNo = 10;
StructString = 'LogAmplitude';
LabelString = 'Log Amplitude (dB)';
ParametricOrNot = 1;
[OutputStats(1)] = Harini_DoStats_AndPlots_LogAmp_LM(DirBout_Stats, DirUndirBout_Stats, UndirBout_Stats, StructString, MinTrialNo, Distances, LabelString, BirdNames, ParametricOrNot);

AllSyll_LMFit_ModelPValue = [];
AllSyll_LMFit_Distance = [];
AllSyll_LMFit_Context = [];
AllSyll_LMFit_DistanceContext_Interactions = [];
AllSyll_LMFit_Intercept = [];

AllSyll_LMFit_ZScore_Distance = [];
AllSyll_LMFit_ZScore_Context = [];
AllSyll_LMFit_ZScore_DistanceContext_Interactions = [];
AllSyll_LMFit_ZScore_Intercept = [];
AllSyll_BirdSyllIndex = [];

for i = 1:length(OutputStats(1).LMFit_ModelPValue),
    AllSyll_LMFit_ModelPValue = [AllSyll_LMFit_ModelPValue; OutputStats(1).LMFit_ModelPValue{i}(:)];
    AllSyll_BirdSyllIndex = [AllSyll_BirdSyllIndex; [ones(size(OutputStats(1).LMFit_ModelPValue{i}(:)))*i (1:1:length(OutputStats(1).LMFit_ModelPValue{i}))']];
    AllSyll_LMFit_Distance = [AllSyll_LMFit_Distance; OutputStats(1).LMFit_Distance{i}];
    AllSyll_LMFit_Context = [AllSyll_LMFit_Context; OutputStats(1).LMFit_Context{i}];
    AllSyll_LMFit_DistanceContext_Interactions = [AllSyll_LMFit_DistanceContext_Interactions; OutputStats(1).LMFit_DistanceContext_Interactions{i}];
    AllSyll_LMFit_Intercept = [AllSyll_LMFit_Intercept; OutputStats(1).LMFit_Intercept{i}];
    
    AllSyll_LMFit_ZScore_Distance = [AllSyll_LMFit_ZScore_Distance; OutputStats(1).LMFit_ZScore_Distance{i}];
    AllSyll_LMFit_ZScore_Context = [AllSyll_LMFit_ZScore_Context; OutputStats(1).LMFit_ZScore_Context{i}];
    AllSyll_LMFit_ZScore_DistanceContext_Interactions = [AllSyll_LMFit_ZScore_DistanceContext_Interactions; OutputStats(1).LMFit_ZScore_DistanceContext_Interactions{i}];
    AllSyll_LMFit_ZScore_Intercept = [AllSyll_LMFit_ZScore_Intercept; OutputStats(1).LMFit_ZScore_Intercept{i}];
end


ResponseCategories = {'INTERACTION' 'DISTANCE' 'SONG TYPE' 'DISTANCE + SONG TYPE' 'NO CHANGE'};
NaNBirds = union(find(isnan(AllSyll_LMFit_ModelPValue)), find(isnan(AllSyll_LMFit_DistanceContext_Interactions(:,2))));

for j = 1:length(ResponseCategories),
    switch (ResponseCategories{j})
        case 'INTERACTION'
            Responses{j} = find((AllSyll_LMFit_ModelPValue(:) < 0.05) & (AllSyll_LMFit_DistanceContext_Interactions(:,2) < 0.05));

        case 'DISTANCE'
            Responses{j} = find((AllSyll_LMFit_ModelPValue(:) < 0.05) & (AllSyll_LMFit_Distance(:,2) < 0.05) & (AllSyll_LMFit_Context(:,2) >= 0.05) & (AllSyll_LMFit_DistanceContext_Interactions(:,2) >= 0.05));

        case 'SONG TYPE'
            Responses{j} = find((AllSyll_LMFit_ModelPValue(:) < 0.05) & (AllSyll_LMFit_Distance(:,2) >= 0.05) & (AllSyll_LMFit_Context(:,2) < 0.05) & (AllSyll_LMFit_DistanceContext_Interactions(:,2) >= 0.05));

        case 'DISTANCE + SONG TYPE'
            Responses{j} = find((AllSyll_LMFit_ModelPValue(:) < 0.05) & (AllSyll_LMFit_Distance(:,2) < 0.05) & (AllSyll_LMFit_Context(:,2) < 0.05) & (AllSyll_LMFit_DistanceContext_Interactions(:,2) >= 0.05));

        case 'NO CHANGE'
            Responses{j} = union(find((AllSyll_LMFit_ModelPValue(:) >= 0.05)), find((AllSyll_LMFit_Distance(:,2) >= 0.05) & (AllSyll_LMFit_Context(:,2) >= 0.05) & (AllSyll_LMFit_DistanceContext_Interactions(:,2) >= 0.05)));
    end
    Responses{j} = setdiff(Responses{j}, NaNBirds);
end

NumResponses = cellfun(@length, Responses);

TempAllBirds = [];

TempResponses = [];
for j = 1:length(Responses),
    TempResponses = [TempResponses; Responses{j}(:)];
end
TempAllBirds = TempResponses;

BirdsToPlot = [1 9 4 3]'; % examples to plot for # of INs, # of motifs and FMD

AllBirds = setdiff(TempAllBirds, BirdsToPlot, 'rows', 'stable');


% Now to plot all bout level features across all birds in one panel
figure;
set(gcf, 'Color', 'w');
p = panel();
p.pack(4, 4);

p.de.margin = 12;

Index = 0;
RowNo = 0;
for j = AllBirds(:)',
    Index = Index + 1;
    if (mod(Index,4) == 1)
        RowNo = RowNo + 1;
    end
    p(RowNo, Index - (RowNo - 1)*4).select();
    hold on;
    NanDistances = find(~isnan(OutputStats(1).DistanceDirMeans(j,1:5)));
    errorbar(Distances(1:5), OutputStats(1).DistanceDirMeans(j,1:5), OutputStats(1).DistanceDirSEMs(j,1:5), 'ro', 'MarkerFaceColor', 'r', 'LineWidth', 0.5, 'MarkerSize', 6);
    plot(Distances(NanDistances), polyval([AllSyll_LMFit_Distance(j,1) AllSyll_LMFit_Intercept(j,1)], Distances(NanDistances)), 'r', 'LineWidth', 1);

    errorbar(Distances(1:5), OutputStats(1).DistanceUnDirMeans(j,1:5), OutputStats(1).DistanceUnDirSEMs(j,1:5), 'bo', 'MarkerFaceColor', 'b', 'LineWidth', 0.5, 'MarkerSize', 6);
    NanDistances = find(~isnan(OutputStats(1).DistanceUnDirMeans(j,1:5)));
    plot(Distances(NanDistances), polyval([AllSyll_LMFit_Distance(j,1)+AllSyll_LMFit_DistanceContext_Interactions(j,1) AllSyll_LMFit_Intercept(j,1)+AllSyll_LMFit_Context(j,1)], Distances(NanDistances)), 'b', 'LineWidth', 1);

    plot([Distances(1) Distances(5)], ones(1,2)*OutputStats(1).DistanceUnDirMeans(j,end) + OutputStats(1).DistanceUnDirSEMs(j,end), 'k--');
    plot([Distances(1) Distances(5)], ones(1,2)*OutputStats(1).DistanceUnDirMeans(j,end) - OutputStats(1).DistanceUnDirSEMs(j,end), 'k--');

    axis tight;
    Temp = axis;
    Temp = [-5 185 0.99*Temp(3) 1.05*Temp(4)];
    axis(Temp);

    % Write the model p-value and the adjusted r^2 on the plot
    text(0, Temp(4), ['p = ', num2str(AllSyll_LMFit_ModelPValue(j)), '; Adj. r^{2} = ', num2str(OutputStats(1).Models{AllSyll_BirdSyllIndex(j,1)}{AllSyll_BirdSyllIndex(j,2)}.mdl.Rsquared.Adjusted)], 'FontSize', 6)

    title(['Bird #', num2str(AllSyll_BirdSyllIndex(j,1)), ': Syll #', num2str(AllSyll_BirdSyllIndex(j,2))], 'FontSize', 10);
    if (ceil(Temp(3)) == floor(Temp(4)))
        set(gca, 'YTick', [ceil(Temp(3)*10)/10 floor(Temp(4)*10)/10]);
    else
        set(gca, 'YTick', [ceil(Temp(3)) floor(Temp(4))]);
    end

%        ylabel(OutputStatsYLabels{i});
    if (i > 1)
        p(RowNo, Index - (RowNo - 1)*4).marginleft = 20;
    end

    set(gca, 'XTick', [0 75 150], 'XTickLabel', [{'0' '75' '150'}]);
    if (j == AllBirds(end))
        xlabel('Distance from female (cm)');
    end
end

p.fontsize = 12;
% p.fontname = 'Times';
%p.margintop = 10;
%p.marginleft = 15;
%p.marginright = 5;
%p.marginbottom = 25;
p.de.margin = 25;

%set(gcf, 'ReSize', 'off');
set(gcf, 'Color', 'w');
set(gcf, 'Units', 'inches');
set(gcf, 'Position', [7 0.5 8.5 8]);
set(gcf, 'PaperPositionMode', 'auto');
print(fullfile(FinalFigureDir, ['LogAmplitude_AllBirdsExceptExamples.eps']), '-depsc2', '-r600');
print(fullfile(FinalFigureDir, ['LogAmplitude_AllBirdsExceptExamples.png']), '-dpng', '-r600');

% Now to plot all bout level features across all birds in one panel
figure;
set(gcf, 'Color', 'w');
p = panel();
p.pack({2/5 3/5});
p(1).pack(1, length(BirdsToPlot));
p(2).pack(1, 3);

p.de.margin = 20;

Index = 0;
for j = (BirdsToPlot(:))',
    Index = Index + 1;
    p(1,1,Index).select();
    hold on;
    NanDistances = find(~isnan(OutputStats(1).DistanceDirMeans(j,1:5)));
    errorbar(Distances(1:5), OutputStats(1).DistanceDirMeans(j,1:5), OutputStats(1).DistanceDirSEMs(j,1:5), 'ro', 'MarkerFaceColor', 'r', 'LineWidth', 0.5, 'MarkerSize', 6);
    plot(Distances(NanDistances), polyval([AllSyll_LMFit_Distance(j,1) AllSyll_LMFit_Intercept(j,1)], Distances(NanDistances)), 'r', 'LineWidth', 1);

    errorbar(Distances(1:5), OutputStats(1).DistanceUnDirMeans(j,1:5), OutputStats(1).DistanceUnDirSEMs(j,1:5), 'bo', 'MarkerFaceColor', 'b', 'LineWidth', 0.5, 'MarkerSize', 6);
    NanDistances = find(~isnan(OutputStats(1).DistanceUnDirMeans(j,1:5)));
    plot(Distances(NanDistances), polyval([AllSyll_LMFit_Distance(j,1)+AllSyll_LMFit_DistanceContext_Interactions(j,1) AllSyll_LMFit_Intercept(j,1)+AllSyll_LMFit_Context(j,1)], Distances(NanDistances)), 'b', 'LineWidth', 1);

    plot([Distances(1) Distances(5)], ones(1,2)*OutputStats(1).DistanceUnDirMeans(j,end) + OutputStats(1).DistanceUnDirSEMs(j,end), 'k--');
    plot([Distances(1) Distances(5)], ones(1,2)*OutputStats(1).DistanceUnDirMeans(j,end) - OutputStats(1).DistanceUnDirSEMs(j,end), 'k--');

    axis tight;
    Temp = axis;
    Temp = [-5 185 0.99*Temp(3) 1.05*Temp(4)];
    axis(Temp);

    % Write the model p-value and the adjusted r^2 on the plot
    text(0, Temp(4)/1.01, {['p = ', num2str(AllSyll_LMFit_ModelPValue(j))]; ['Adj. r^{2} = ', num2str(OutputStats(1).Models{AllSyll_BirdSyllIndex(j,1)}{AllSyll_BirdSyllIndex(j,2)}.mdl.Rsquared.Adjusted)]}, 'FontSize', 8)

    title(['Bird #', num2str(AllSyll_BirdSyllIndex(j,1)), ': Syll #', num2str(AllSyll_BirdSyllIndex(j,2))], 'FontSize', 10);
    if (ceil(Temp(3)) == floor(Temp(4)))
        set(gca, 'YTick', [ceil(Temp(3)*10)/10 floor(Temp(4)*10)/10]);
    else
        set(gca, 'YTick', [ceil(Temp(3)) floor(Temp(4))]);
    end

    if (Index == 1)
        ylabel('Log Amplitude (dB)');
        p(1, 1, Index).marginleft = 20;
    end
    
    set(gca, 'XTick', [0 75 150], 'XTickLabel', [{'0' '75' '150'}]);
    % xlabel('Distance from female (cm)');
end

p(2,1,1).select();
TempBar = bar(NumResponses);
set(TempBar, 'FaceColor', 'none', 'EdgeColor', 'k', 'LineWidth', 1);
set(gca, 'Box', 'off');
ylabel('# of syllables');
set(gca, 'XTick', 1:1:length(ResponseCategories), 'XTickLabel', ResponseCategories, 'XTickLabelRotation', 45);
axis([0.5 (length(ResponseCategories)+0.5) 0 max(NumResponses(:))+1]);

% Now to plot the slopes for each of the categories.
p(2,1,2).select();
hold on;
Symbols = 's^*';

for j = 1:length(ResponseCategories),
    Slopes{j} = [];
end

for j = 1:length(ResponseCategories)-1,
    if (length(Responses{j}) <= 1)
        continue;
    end
    Changes = setdiff(Responses{j}, BirdsToPlot);
    plot(repmat([1.2 1.8] + (j-1)*2, length(Changes), 1)', [AllSyll_LMFit_ZScore_Distance(Changes,1) (AllSyll_LMFit_ZScore_Distance(Changes,1) + AllSyll_LMFit_ZScore_DistanceContext_Interactions(Changes,1))]', 'k');
    plot(repmat([1.2] + (j-1)*2, length(Changes), 1)', [AllSyll_LMFit_ZScore_Distance(Changes,1)]', ['r', Symbols(1)], 'MarkerSize', 6);
    plot(repmat([1.8] + (j-1)*2, length(Changes), 1)', [(AllSyll_LMFit_ZScore_Distance(Changes,1) + AllSyll_LMFit_ZScore_DistanceContext_Interactions(Changes,1))]', ['b', Symbols(1)], 'MarkerSize', 6);

    Changes = intersect(Responses{j}, BirdsToPlot);
    plot(repmat([1.2 1.8] + (j-1)*2, length(Changes), 1)', [AllSyll_LMFit_ZScore_Distance(Changes,1) (AllSyll_LMFit_ZScore_Distance(Changes,1) + AllSyll_LMFit_ZScore_DistanceContext_Interactions(Changes,1))]', 'k');
    plot(repmat([1.2] + (j-1)*2, length(Changes), 1)', [AllSyll_LMFit_ZScore_Distance(Changes,1)]', ['r', Symbols(1)], 'MarkerSize', 9);
    plot(repmat([1.8] + (j-1)*2, length(Changes), 1)', [(AllSyll_LMFit_ZScore_Distance(Changes,1) + AllSyll_LMFit_ZScore_DistanceContext_Interactions(Changes,1))]', ['b', Symbols(1)], 'MarkerSize', 9);

    Slopes{j} = [Slopes{j}; [AllSyll_LMFit_ZScore_Distance(Responses{j},1) (AllSyll_LMFit_ZScore_Distance(Responses{j},1) + AllSyll_LMFit_ZScore_DistanceContext_Interactions(Responses{j},1))]];
end

for i = 1:length(Slopes),
    if (size(Slopes{i},1) <= 1)
        continue;
    end
    ResponseBar(i) = bar([1 2] + (i-1)*2, mean(Slopes{i}));
    set(ResponseBar(i), 'FaceColor', 'none', 'LineWidth', 0.5);
    errorbar([1 2] + (i-1)*2, mean(Slopes{i}), std(Slopes{i})./sqrt(size(Slopes{i},1)), 'ks', 'LineWidth', 0.5, 'MarkerSize', 8);
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

p(2,1,3).select();
hold on;
% For interaction and distance birds, plot slopes for L0 > UN and slopes
% for L0 < UN separately for directed songs and undirected songs
% For interaction and distance birds, plot slopes for L0 > UN and slopes
% for L0 < UN separately for directed songs and undirected songs
DirL0_MoreThan_UN_Slopes = [];
DirUN_MoreThan_L0_Slopes = [];

UnDirL0_MoreThan_UN_Slopes = [];
UnDirUN_MoreThan_L0_Slopes = [];

SlopesToConsider = find(cellfun(@length, strfind(ResponseCategories, 'DISTANCE')));
SlopesToConsider = [1; SlopesToConsider(:)];

for k = (SlopesToConsider(:))',
    for j = Responses{k}(:)',
        if (isnan(OutputStats(1).DistanceDirMeans(j,1)) || isnan(OutputStats(1).DistanceUnDirMeans(j,end)))
            continue;
        end
        if (OutputStats(1).DistanceDirMeans(j,1) >= OutputStats(1).DistanceUnDirMeans(j,end))
            plot(1, AllSyll_LMFit_ZScore_Distance(j,1), ['r', Symbols(1)], 'LineWidth', 1, 'MarkerSize', 6);
            DirL0_MoreThan_UN_Slopes(end+1) = AllSyll_LMFit_ZScore_Distance(j,1);

            plot(3, (AllSyll_LMFit_ZScore_Distance(j,1) + AllSyll_LMFit_ZScore_DistanceContext_Interactions(j,1)), ['b', Symbols(1)], 'LineWidth', 1, 'MarkerSize', 6);
            UnDirL0_MoreThan_UN_Slopes(end+1) = (AllSyll_LMFit_ZScore_Distance(j,1) + AllSyll_LMFit_ZScore_DistanceContext_Interactions(j,1));
        else
            plot(2, AllSyll_LMFit_ZScore_Distance(j,1), ['r', Symbols(1)], 'LineWidth', 1, 'MarkerSize', 6);
            DirUN_MoreThan_L0_Slopes(end+1) = AllSyll_LMFit_ZScore_Distance(j,1);

            plot(4, (AllSyll_LMFit_ZScore_Distance(j,1) + AllSyll_LMFit_ZScore_DistanceContext_Interactions(j,1)), ['b', Symbols(1)], 'LineWidth', 1, 'MarkerSize', 6);
            UnDirUN_MoreThan_L0_Slopes(end+1) = (AllSyll_LMFit_ZScore_Distance(j,1) + AllSyll_LMFit_ZScore_DistanceContext_Interactions(j,1));
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

set(gca, 'XTick', [1 2 3 4], 'XTickLabel', {'L0-DIR \geq ALONE' 'LO-DIR < ALONE' 'L0-DIR \geq ALONE' 'LO-DIR < ALONE'}, 'XTickLabelRotation', 45);
ylabel('Slope');

p.fontsize = 12;
% p.fontname = 'Times';
p.margintop = 10;
p.marginleft = 15;
p.marginright = 5;
p.marginbottom = 45;

% set(gcf, 'ReSize', 'off');
set(gcf, 'Color', 'w');
set(gcf, 'Units', 'inches');
set(gcf, 'Position', [7 0.5 8.5 7]);
set(gcf, 'PaperPositionMode', 'auto');
print(fullfile(FinalFigureDir, ['LogAmplitude_Fig.3.ResponseExamples_Counts.eps']), '-depsc2', '-r600');
print(fullfile(FinalFigureDir, ['LogAmplitude_Fig.3.ResponseExamples_Counts.png']), '-dpng', '-r600');

figure;
set(gcf, 'Position', [11 388 504 372]);
p = panel();
p.pack(1,1);
p(1,1).select();
hold on;
for i = 1:size(AllSyll_BirdSyllIndex,1),
    L0DirIndices = find((OutputStats(1).GroupData{AllSyll_BirdSyllIndex(i,1)}{AllSyll_BirdSyllIndex(i,2)}(:,2) == 1) & (OutputStats(1).GroupData{AllSyll_BirdSyllIndex(i,1)}{AllSyll_BirdSyllIndex(i,2)}(:,3) == 1));
    AloneIndices = find((OutputStats(1).GroupData{AllSyll_BirdSyllIndex(i,1)}{AllSyll_BirdSyllIndex(i,2)}(:,2) == 6) & (OutputStats(1).GroupData{AllSyll_BirdSyllIndex(i,1)}{AllSyll_BirdSyllIndex(i,2)}(:,3) == 3));
    if (~isempty(L0DirIndices) && ~isempty(AloneIndices))
        [H, Prob_Diff(i)] = ttest2(OutputStats(1).GroupData{AllSyll_BirdSyllIndex(i,1)}{AllSyll_BirdSyllIndex(i,2)}(L0DirIndices,1), OutputStats(1).GroupData{AllSyll_BirdSyllIndex(i,1)}{AllSyll_BirdSyllIndex(i,2)}(AloneIndices,1));
    else
        Prob_Diff(i) = NaN;
    end
    if (OutputStats(1).DistanceDirMeans(i,1) >= OutputStats(1).DistanceUnDirMeans(i,end))
        MarkerColour = 'r';
    else
        MarkerColour = 'b';
    end
    
    if (Prob_Diff(i) < 0.05)
        MarkerFaceColorYesorNo = MarkerColour;
    else
        MarkerFaceColorYesorNo = 'none';
    end
    plot(AllSyll_BirdSyllIndex(i,1)*2+(AllSyll_BirdSyllIndex(i,2)/7), OutputStats(1).DistanceDirMeans(i,1)/OutputStats(1).DistanceUnDirMeans(i,end), 'ko', 'Color', MarkerColour, 'MarkerFaceColor', MarkerFaceColorYesorNo);
    AmplitudeRatios(i) = OutputStats(1).DistanceDirMeans(i,1)/OutputStats(1).DistanceUnDirMeans(i,end);
end
ylabel('Log Amplitude ratio (L0-DIR/ALONE)');
axis tight;
Temp = axis;
Temp = [1.75 (2*AllSyll_BirdSyllIndex(end,1) + 0.25) 0.99*Temp(3) 1.01*Temp(4)];
axis(Temp);
plot(Temp(1:2), [1 1], 'k--');
p.fontsize = 12;
p.marginleft = 20;
p.marginbottom = 25;
set(gca, 'XTick', [2:2:14], 'XTickLabel', {'Bird #1', 'Bird #2', 'Bird #3', 'Bird #4', 'Bird #5', 'Bird #6', 'Bird #7'}, 'XTickLabelRotation', 45); 
set(gcf, 'Color', 'w');
set(gcf, 'PaperPositionMode', 'auto');
print(fullfile(FinalFigureDir, ['LogAmplitude_Ratio.eps']), '-depsc2', '-r600');
print(fullfile(FinalFigureDir, ['LogAmplitude_Ratio.png']), '-dpng', '-r600');


AdjustedRsquared = [];
for i = 1:length(OutputStats(1).Models),
    for j = 1:length(OutputStats(1).Models{i}),
        AdjustedRsquared = [AdjustedRsquared; OutputStats(1).Models{i}{j}.mdl.Rsquared.Adjusted];
    end
end
SigPValues = find(cell2mat(OutputStats(1).LMFit_ModelPValue) < 0.05);
NonNanPValues = find(~isnan(cell2mat(OutputStats(1).LMFit_ModelPValue)));
disp(['# of sig p-values = ', num2str(length(SigPValues)), ' out of ', num2str(length(NonNanPValues)), ' valid (non-NaN) syllables']);
disp(['Variance explained by sig. fits: mean = ', num2str(100*mean(AdjustedRsquared(SigPValues))), ' %, median = ', num2str(100*median(AdjustedRsquared(SigPValues))), '%, min = ', num2str(100*min(AdjustedRsquared(SigPValues))), '%, max = ', num2str(100*max(AdjustedRsquared(SigPValues))), '%']);


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

% MinTrialNo = 3;
% StructString = 'LogAmplitude';
% LabelString = 'Log Amplitude (dB)';
% 
% Harini_AllFeatureStats_Plots_WithoutCI(DirBout_Stats, DirUndirBout_Stats, UndirBout_Stats, StructString, MinTrialNo, Distances, LabelString, BirdNames);

disp('Finished calculating bout stats');