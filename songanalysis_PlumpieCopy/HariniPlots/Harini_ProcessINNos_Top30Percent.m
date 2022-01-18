function [DirBout_Stats, UndirBout_Stats, DirUndirBout_Stats] = Harini_ProcessINNos_Top30Percent(IndividualBirds, BirdNames, BirdOption, FeatureToPlot_ColNo, PlotType, DifferentMicrophones)

% ====== IN Number processing =============================================
% THis is to check all aspects related to IN number in all the birds.

% =========================================================================

MinTrialNo = 3;
FractionToConsider = 0.3;

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
    
        AllBoutStats{i}{j}.TotalINNumber_500ms = [];
        AllBoutStats{i}{j}.CompleteMotifNumber = [];
        AllBoutStats{i}{j}.FirstMotifDuration = [];
        
        if (~isempty(DirBouts{i}{j}))
            [DirBout_Stats{i}{j}] = Harini_CalculateBoutINMotifNumbers_Duration(IndividualBirds(i), DirBouts{i}{j}, MotifLabels, INLabels, CommonMotifs);
            AllBoutStats{i}{j}.TotalINNumber_500ms = [AllBoutStats{i}{j}.TotalINNumber_500ms DirBout_Stats{i}{j}.TotalINNumber_500ms];
            AllBoutStats{i}{j}.CompleteMotifNumber = [AllBoutStats{i}{j}.CompleteMotifNumber DirBout_Stats{i}{j}.CompleteMotifNumber];
            AllBoutStats{i}{j}.FirstMotifDuration = [AllBoutStats{i}{j}.FirstMotifDuration DirBout_Stats{i}{j}.FirstMotifDuration];
        else
            DirBout_Stats{i}{j} = [];
        end
        if (~isempty(DirUndirBouts{i}{j}))
            [DirUndirBout_Stats{i}{j}] = Harini_CalculateBoutINMotifNumbers_Duration(IndividualBirds(i), DirUndirBouts{i}{j}, MotifLabels, INLabels, CommonMotifs);
            AllBoutStats{i}{j}.TotalINNumber_500ms = [AllBoutStats{i}{j}.TotalINNumber_500ms DirUndirBout_Stats{i}{j}.TotalINNumber_500ms];
            AllBoutStats{i}{j}.CompleteMotifNumber = [AllBoutStats{i}{j}.CompleteMotifNumber DirUndirBout_Stats{i}{j}.CompleteMotifNumber];
            AllBoutStats{i}{j}.FirstMotifDuration = [AllBoutStats{i}{j}.FirstMotifDuration DirUndirBout_Stats{i}{j}.FirstMotifDuration];
        else
            DirUndirBout_Stats{i}{j} = [];
        end
        if (~isempty(UndirBouts{i}{j}))
            [UndirBout_Stats{i}{j}] = Harini_CalculateBoutINMotifNumbers_Duration(IndividualBirds(i), UndirBouts{i}{j}, MotifLabels, INLabels, CommonMotifs);
            AllBoutStats{i}{j}.TotalINNumber_500ms = [AllBoutStats{i}{j}.TotalINNumber_500ms UndirBout_Stats{i}{j}.TotalINNumber_500ms];
            AllBoutStats{i}{j}.CompleteMotifNumber = [AllBoutStats{i}{j}.CompleteMotifNumber UndirBout_Stats{i}{j}.CompleteMotifNumber];
            AllBoutStats{i}{j}.FirstMotifDuration = [AllBoutStats{i}{j}.FirstMotifDuration UndirBout_Stats{i}{j}.FirstMotifDuration];
        else
            UndirBout_Stats{i}{j} = [];
        end
        
        if (j ~= 6)
            [SortedData, SortedDataIndices] = sort(AllBoutStats{i}{j}.TotalINNumber_500ms);
            AllBoutStats{i}{j}.TotalINNumber_500ms = AllBoutStats{i}{j}.TotalINNumber_500ms(SortedDataIndices);
            AllBoutStats{i}{j}.CompleteMotifNumber = AllBoutStats{i}{j}.CompleteMotifNumber(SortedDataIndices);
            AllBoutStats{i}{j}.FirstMotifDuration = AllBoutStats{i}{j}.FirstMotifDuration(SortedDataIndices);

            if (length(SortedDataIndices) >= MinTrialNo*3)
                DirBout_Stats{i}{j}.TotalINNumber_500ms = AllBoutStats{i}{j}.TotalINNumber_500ms(round((1-FractionToConsider)*length(SortedDataIndices)):end);
                UnDirBout_Stats{i}{j}.TotalINNumber_500ms = AllBoutStats{i}{j}.TotalINNumber_500ms(1:round(FractionToConsider*length(SortedDataIndices)));

                DirBout_Stats{i}{j}.CompleteMotifNumber = AllBoutStats{i}{j}.CompleteMotifNumber(round((1-FractionToConsider)*length(SortedDataIndices)):end);
                UnDirBout_Stats{i}{j}.CompleteMotifNumber = AllBoutStats{i}{j}.CompleteMotifNumber(1:round(FractionToConsider*length(SortedDataIndices)));

                DirBout_Stats{i}{j}.FirstMotifDuration = AllBoutStats{i}{j}.FirstMotifDuration(round((1-FractionToConsider)*length(SortedDataIndices)):end);
                UnDirBout_Stats{i}{j}.FirstMotifDuration = AllBoutStats{i}{j}.FirstMotifDuration(1:round(FractionToConsider*length(SortedDataIndices)));
            else
                DirBout_Stats{i}{j} = [];
                UnDirBout_Stats{i}{j} = [];
            end
        end
    end
end


% Now to do all the stats and plots

% First for all IN nos. - one set with full IN nos and another with only IN
% numbers that are separated by <= 500ms
Distances = [0 20 60 110 165 200];

RValuesFig = figure;
RValuevsL0UNDiffFig = figure;

MinTrialNo = 3;
StructString = 'TotalINNumber';
LabelString = '# of INs';
ParametricOrNot = 0;
% [TotalINNumber_L0Dir_Undir_PValue] = Harini_DoStats_AndPlots(DirBout_Stats, DirUndirBout_Stats, UndirBout_Stats, StructString, MinTrialNo, Distances, LabelString, BirdNames, ParametricOrNot);

MinTrialNo = 3;
StructString = 'TotalINNumber_500ms';
LabelString = '# of INs_500ms';
ParametricOrNot = 0;
[OutputStats(1)] = Harini_DoStats_AndPlots(DirBout_Stats, DirUndirBout_Stats, UndirBout_Stats, StructString, MinTrialNo, Distances, LabelString, BirdNames, ParametricOrNot);

MinTrialNo = 3;
StructString = 'BoutLength';
LabelString = 'Bout length (msec)';
ParametricOrNot = 1;
% [TotalINNumber_L0Dir_Undir_PValue] = Harini_DoStats_AndPlots(DirBout_Stats, DirUndirBout_Stats, UndirBout_Stats, StructString, MinTrialNo, Distances, LabelString, BirdNames, ParametricOrNotRValuesFig, RValuevsL0UNDiffFig);

MinTrialNo = 3;
% Next for complete motif nos.
StructString = 'CompleteMotifNumber';
LabelString = '# of complete motifs';
ParametricOrNot = 0;
[OutputStats(2)] = Harini_DoStats_AndPlots(DirBout_Stats, DirUndirBout_Stats, UndirBout_Stats, StructString, MinTrialNo, Distances, LabelString, BirdNames, ParametricOrNot);

MinTrialNo = 3;
StructString = 'PartialMotifNumber';
LabelString = '# of partial motifs';
ParametricOrNot = 0;
% [TotalINNumber_L0Dir_Undir_PValue] = Harini_DoStats_AndPlots(DirBout_Stats, DirUndirBout_Stats, UndirBout_Stats, StructString, MinTrialNo, Distances, LabelString, BirdNames, ParametricOrNotRValuesFig, RValuevsL0UNDiffFig);

MinTrialNo = 3;
StructString = 'FirstMotifDuration';
LabelString = 'First motif duration (msec)';
ParametricOrNot = 1;
[OutputStats(3)] = Harini_DoStats_AndPlots(DirBout_Stats, DirUndirBout_Stats, UndirBout_Stats, StructString, MinTrialNo, Distances, LabelString, BirdNames, ParametricOrNot);

OutputStatsYLabels = {'# of INs' '# of motifs/bout', 'FMD (msec)'};

% Now to plot all bout level features across all birds in one panel
figure;
set(gcf, 'Color', 'w');
p = panel();
NonNanBirds = find(~isnan(OutputStats(1).Distance_DirSong_Corr(:,1)));
p.pack(length(NonNanBirds),length(OutputStats));
p.de.margin = 8;

for i = 1:length(OutputStats),
    Index = 0;
    for j = NonNanBirds(:)',
        Index = Index + 1;
        p(Index,i).select();
        hold on;
        plot(Distances(1:5), OutputStats(i).DistanceDirMeans(j,1:5), 'rs-', 'LineWidth', 1);
        plot(repmat(Distances(1:5)', 1, 2)', OutputStats(i).DistanceDirCIs{j}(1:5,:)', 'r')

        plot(Distances(1:5), OutputStats(i).DistanceUnDirMeans(j,1:5), 'bs-', 'LineWidth', 1);
        plot(repmat(Distances(1:5)', 1, 2)', OutputStats(i).DistanceUnDirCIs{j}(1:5,:)', 'b');
        plot(repmat(Distances(1:5)',1,2), repmat(OutputStats(i).DistanceUnDirCIs{j}(end,:), length(Distances)-1, 1), 'k--');
        
        axis tight;
        Temp = axis;
        Temp = [-5 185 0.99*Temp(3) 1.01*Temp(4)];
        axis(Temp);
        
        if (OutputStats(i).Distance_DirSong_Corr(j,2) < 0.05)
            text(170, Temp(4)/1.01, '*', 'FontSize', 12, 'Color', 'r');
        else
            if (OutputStats(i).KWP_DirSong{j} < 0.05)
                text(170, Temp(4)/1.01, '+', 'FontSize', 12, 'Color', 'r');
            else
                text(170, Temp(4)/1.01, '=', 'FontSize', 12, 'Color', 'r');
            end
        end
        
        if (OutputStats(i).Distance_UnDirSong_Corr(j,2) < 0.05)
            text(180, Temp(4)/1.01, '*', 'FontSize', 12, 'Color', 'b');
        else
            if (OutputStats(i).KWP_UnDirSong{j} < 0.05)
                text(180, Temp(4)/1.01, '+', 'FontSize', 12, 'Color', 'b');
            else
                text(170, Temp(4)/1.01, '=', 'FontSize', 12, 'Color', 'b');    
            end
        end
        
        if (ceil(Temp(3)) == floor(Temp(4)))
            set(gca, 'YTick', [ceil(Temp(3)*10)/10 floor(Temp(4)*10)/10]);
        else
            set(gca, 'YTick', [ceil(Temp(3)) floor(Temp(4))]);
        end
        if (j == NonNanBirds(round(length(NonNanBirds)/2)))
            ylabel(OutputStatsYLabels{i});
        end
        if (i > 1)
            p(Index,i).marginleft = 20;
        end
        if (j ~= NonNanBirds(end))
            set(gca, 'XTick', []);
        else
            set(gca, 'XTick', [0 75 150], 'XTickLabel', [{'0' '75' '150'}]);
            if (i == 2)
                xlabel('Distance from female (cm)');
            end
        end
    end
end

p.fontsize = 10;
% p.fontname = 'Times';
p.margintop = 10;
p.marginleft = 20;
p.marginright = 5;

set(gcf, 'ReSize', 'off');
set(gcf, 'Color', 'w');
set(gcf, 'Units', 'inches');
set(gcf, 'Position', [7 0.5 8 11]);
set(gcf, 'PaperPositionMode', 'auto');
FinalFigureDir = pwd;
print(fullfile(FinalFigureDir, ['AllBirds.BoutLevelFeatures.eps']), '-depsc2', '-r600');
print(fullfile(FinalFigureDir, ['AllBirds.BoutLevelFeatures.png']), '-dpng', '-r600');


MinTrialNo = 3;
StructString = 'AllMotifDuration';
LabelString = 'All motif durations (msec)';
ParametricOrNot = 1;
% [TotalINNumber_L0Dir_Undir_PValue] = Harini_DoStats_AndPlots(DirBout_Stats, DirUndirBout_Stats, UndirBout_Stats, StructString, MinTrialNo, Distances, LabelString, BirdNames, ParametricOrNot);

MinTrialNo = 3;
StructString = 'Bout_FFCV';
LabelString = 'CV of FF';
ParametricOrNot = 1;
[OutputStats(4)] = Harini_DoStats_AndPlots(DirBout_Stats, DirUndirBout_Stats, UndirBout_Stats, StructString, MinTrialNo, Distances, LabelString, BirdNames, ParametricOrNot);

% For FF and amplitude, min trial no corresponds to minimum # of syllables
% so we will change that to 15
MinTrialNo = 15;

StructString = 'FF';
LabelString = 'FF (Hz)';
ParametricOrNot = 1;
[OutputStats(4)] = Harini_DoStats_AndPlots_FF(DirBout_Stats, DirUndirBout_Stats, UndirBout_Stats, StructString, MinTrialNo, Distances, LabelString, BirdNames, ParametricOrNot);

StructString = 'LogAmplitude';
LabelString = 'Log Amplitude (dB)';
ParametricOrNot = 1;
[OutputStats(5)] = Harini_DoStats_AndPlots_FF_LogAmp(DirBout_Stats, DirUndirBout_Stats, UndirBout_Stats, StructString, MinTrialNo, Distances, LabelString, BirdNames, ParametricOrNot);

% plot all the different correlation co-efficients and their significances
% for each of the features
figure;
hold on;
for i = 1:length(OutputStats),
    Indices = find(OutputStats(i).Distance_DirSong_Corr(:,2) < 0.05);
    plot(ones(size(Indices))*(i-0.15), OutputStats(i).Distance_DirSong_Corr(Indices,1), 'ro', 'MarkerFaceColor', 'r');
    
    Indices = find((OutputStats(i).Distance_DirSong_Corr(:,2) >= 0.05) & (cell2mat(OutputStats(i).KWP_DirSong)' < 0.05));
    plot(ones(size(Indices))*(i-0.15), OutputStats(i).Distance_DirSong_Corr(Indices,1), 'ko', 'MarkerFaceColor', 'k');

    Indices = find((OutputStats(i).Distance_DirSong_Corr(:,2) >= 0.05) & (cell2mat(OutputStats(i).KWP_DirSong)' >= 0.05));
    plot(ones(size(Indices))*(i-0.15), OutputStats(i).Distance_DirSong_Corr(Indices,1), 'ko');
    
    Indices = find(OutputStats(i).Distance_UnDirSong_Corr(:,2) < 0.05);
    plot(ones(size(Indices))*(i+0.15), OutputStats(i).Distance_UnDirSong_Corr(Indices,1), 'bo', 'MarkerFaceColor', 'b');
    
    Indices = find((OutputStats(i).Distance_UnDirSong_Corr(:,2) >= 0.05) & (cell2mat(OutputStats(i).KWP_UnDirSong)' < 0.05));
    plot(ones(size(Indices))*(i+0.15), OutputStats(i).Distance_UnDirSong_Corr(Indices,1), 'ko', 'MarkerFaceColor', 'k');

    Indices = find((OutputStats(i).Distance_UnDirSong_Corr(:,2) >= 0.05) & (cell2mat(OutputStats(i).KWP_UnDirSong)' >= 0.05));
    plot(ones(size(Indices))*(i+0.15), OutputStats(i).Distance_UnDirSong_Corr(Indices,1), 'ko');
end
ylabel('Correlation co-efficient');
set(gca, 'XTick', [1 2 3], 'XTickLabel', {'# of INs' '# of motifs/bout', 'FMD'});
Temp = axis;
Temp = [0.5 3.5 -1 1];
axis(Temp);
plot([0.5 3.5], [0 0], 'k--');

% Now to plot the correlation co-efficient with L0/UN
% First considering only significant correlations
Symbols = 's^*+';
% First for directed
figure;
hold on;
for i = 1:2,
    Indices = find(OutputStats(i).Distance_DirSong_Corr(:,2) < 0.05);
    plot(OutputStats(i).DistanceDirMeans(Indices,1)/100, OutputStats(i).Distance_DirSong_Corr(Indices,1), ['r', Symbols(i)], 'LineWidth', 1, 'MarkerFaceColor', 'r', 'MarkerSize', 8);
end

% For FMD and FF, r-values are flipped and L0/UN is also flipped to be
% consistent with the expectations for changes between L0 and UN across all
% features

for i = 3:3,
    Indices = find(OutputStats(i).Distance_DirSong_Corr(:,2) < 0.05);
    plot(1./(OutputStats(i).DistanceDirMeans(Indices,1)/100), -OutputStats(i).Distance_DirSong_Corr(Indices,1), ['r', Symbols(i)], 'LineWidth', 1, 'MarkerFaceColor', 'r', 'MarkerSize', 8);
end

% First for undirected
for i = 1:2,
    Indices = find(OutputStats(i).Distance_UnDirSong_Corr(:,2) < 0.05);
    plot(OutputStats(i).DistanceDirMeans(Indices,1)/100, OutputStats(i).Distance_UnDirSong_Corr(Indices,1), ['b', Symbols(i)], 'LineWidth', 1, 'MarkerFaceColor', 'b', 'MarkerSize', 8);
end

% For FMD and FF, r-values are flipped and L0/UN is also flipped to be
% consistent with the expectations for changes between L0 and UN across all
% features

for i = 3:3,
    Indices = find(OutputStats(i).Distance_UnDirSong_Corr(:,2) < 0.05);
    plot(1./(OutputStats(i).DistanceDirMeans(Indices,1)/100), -OutputStats(i).Distance_UnDirSong_Corr(Indices,1), ['b', Symbols(i)], 'LineWidth', 1, 'MarkerFaceColor', 'b', 'MarkerSize', 8);
end
axis tight;
Temp = axis;
Temp = [0.95*Temp(1) 1.05*Temp(2) -1 1];
axis(Temp);
plot(Temp(1:2), [0 0], 'k--');
plot([1 1], Temp(3:4), 'k--');

Harini_AllFeatureStats_Plots(DirBout_Stats, DirUndirBout_Stats, UndirBout_Stats, StructString, MinTrialNo, Distances, LabelString, BirdNames);

disp('Finished calculating bout stats');