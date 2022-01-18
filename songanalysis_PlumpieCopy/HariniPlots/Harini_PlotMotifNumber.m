function [] = Harini_PlotMotifNumber_AllPossibilities(IndividualBirds, BirdNames, BirdOption, FeatureToPlot_ColNo, PlotType, DifferentMicrophones)

% ======= IN number in all possible ways ==================================
% 1. For the different conditions irrespective of whether the songs were
% directed or undirected (based on video scoring data)
% 2. For only the DIRECTED songs
% 3. For only the UNDIRECTED songs
% 4. Have to do stats on all of the above
%   a) Kruskal-wallis to check if any of them are different
%   b) Correlation to see if distance and # of Motifs is related - spearman's
%   rank correlation 
% =========================================================================

Colours = 'rgbcmk';
Symbols = '+o<sd^>*pxvh';
Conditions = [{'L0'} {'L1'} {'L2'} {'L3'} {'L4'} {'UN'}];
MinTrialNo = 5;

for i = 1:length(IndividualBirds),
    NumConditions(i) = length(IndividualBirds(i).ConditionDistances);
end

[MaxNumConditions, MaxNumConditionsIndex] = max(NumConditions);

ConditionDistances = IndividualBirds(MaxNumConditionsIndex).ConditionDistances;

Fid = fopen('MotifNumber_Correlations.txt', 'w');
fprintf(Fid, 'Bird Name\tPearsons correlation\t\t\t\t\t\tSpearmans correlation\t\t\t\t\t\n');
fprintf(Fid, '\tAllBouts\t\tDirected bouts\t\tUndirected bouts\t\tAllBouts\t\tDirected bouts\t\tUndirected bouts\t\n');
fprintf(Fid, '\tr\tp\tr\tp\tr\tp\tr\tp\tr\tp\tr\tp\n');

% First collect data for all 3 types of bouts listed above
for i = 1:length(IndividualBirds),
    ColumnIndex = find(strcmp('NumMotifs', IndividualBirds(i).BoutStatisticsColumnNames));
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
    % Now to do the stats for All bouts - kruskal-wallis and then
    % multcompare
    [p, tabl, stats] = kruskalwallis(IndividualBird_AllBouts{i}(:,1), IndividualBird_AllBouts{i}(:,2), 'off');
    SignificanceMatrix{i} = ones(length(IndividualBirds(i).Conditions))*NaN;
    if (p < 0.05)
        Comparisons = multcompare(stats, 'display', 'off');
        for k = 1:size(Comparisons,1),
            if (Comparisons(k,end) < 0.05)
                if ((Comparisons(k,3) < 0) && (Comparisons(k,4) < 0)) 
                    SignificanceMatrix{i}(Comparisons(k,1), Comparisons(k,2)) = -1;
                else
                    SignificanceMatrix{i}(Comparisons(k,1), Comparisons(k,2)) = 1;
                end
            else
                SignificanceMatrix{i}(Comparisons(k,1), Comparisons(k,2)) = 0;
            end
        end
    end

    % ANd then do the correlation linear
    % First Pearson's
    Indices = find(IndividualBird_AllBouts{i}(:,2) < (length(IndividualBirds(i).Conditions) - 1));
    [r, p] = corrcoef(IndividualBird_AllBouts{i}(Indices,1), IndividualBird_AllBouts{i}(Indices,2));
    Pearsons_R(i) = r(1,2);
    Pearsons_P(i) = p(1,2);

    % Also non-parametric Spearmans
    [r, p] = corr(IndividualBird_AllBouts{i}(Indices,1), ConditionDistances(IndividualBird_AllBouts{i}(Indices,2))', 'type', 'Spearman', 'rows', 'complete');
    Spearman_R(i) = r;
    Spearman_P(i) = p;

    % Now for the directed and undirected bouts, but only for birds
    % that have enough data from atleast 3 different female positions

    Directed_Pearsons_R(i) = NaN;
    Directed_Pearsons_P(i) = NaN;
    Directed_Spearman_R(i) = NaN;
    Directed_Spearman_P(i) = NaN;

    UnDirected_Pearsons_R(i) = NaN;
    UnDirected_Pearsons_P(i) = NaN;
    UnDirected_Spearman_R(i) = NaN;
    UnDirected_Spearman_P(i) = NaN;

    if (length(unique(IndividualBird_DirectedBouts{i}(:,2))) > 3)
        Indices = find(IndividualBird_DirectedBouts{i}(:,2) < (length(IndividualBirds(i).Conditions) - 1));
        [r, p] = corrcoef(IndividualBird_DirectedBouts{i}(Indices,1), IndividualBird_DirectedBouts{i}(Indices,2));
        Directed_Pearsons_R(i) = r(1,2);
        Directed_Pearsons_P(i) = p(1,2);

        % Also non-parametric Spearmans
        [r, p] = corr(IndividualBird_DirectedBouts{i}(Indices,1), ConditionDistances(IndividualBird_DirectedBouts{i}(Indices,2))', 'type', 'Spearman', 'rows', 'complete');
        Directed_Spearman_R(i) = r;
        Directed_Spearman_P(i) = p;
    end

    if (length(unique(IndividualBird_UnDirectedBouts{i}(:,2))) > 3)
        Indices = find(IndividualBird_UnDirectedBouts{i}(:,2) < (length(IndividualBirds(i).Conditions) - 1));
        [r, p] = corrcoef(IndividualBird_UnDirectedBouts{i}(Indices,1), IndividualBird_UnDirectedBouts{i}(Indices,2));
        UnDirected_Pearsons_R(i) = r(1,2);
        UnDirected_Pearsons_P(i) = p(1,2);

        % Also non-parametric Spearmans
        [r, p] = corr(IndividualBird_UnDirectedBouts{i}(Indices,1), ConditionDistances(IndividualBird_UnDirectedBouts{i}(Indices,2))', 'type', 'Spearman', 'rows', 'complete');
        UnDirected_Spearman_R(i) = r;
        UnDirected_Spearman_P(i) = p;
    end
    
    % Now write all correlation co-efficients and p-values to file
    fprintf(Fid, '%s\t%1.3f\t%1.3f\t%1.3f\t%1.3f\t%1.3f\t%1.3f\t%1.3f\t%1.3f\t%1.3f\t%1.3f\t%1.3f\t%1.3f\n', BirdNames{i}, Pearsons_R(i), Pearsons_P(i), Directed_Pearsons_R(i), Directed_Pearsons_P(i), UnDirected_Pearsons_R(i), UnDirected_Pearsons_P(i), Spearman_R(i), Spearman_P(i), Directed_Spearman_R(i), Directed_Spearman_P(i), UnDirected_Spearman_R(i), UnDirected_Spearman_P(i));
    
end
fclose(Fid);

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
    
    subplot(2,3,1);
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
    xlabel('# of motifs');
    ylabel('Fraction of bouts');
    title([BirdNames{i}, ': All bouts']);
    
    subplot(2,3,2);
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
    xlabel('# of motifs');
    ylabel('Fraction of bouts');
    title([BirdNames{i}, ': Directed bouts']);
    
    subplot(2,3,3);
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
    xlabel('# of motifs');
    ylabel('Fraction of bouts');
    title([BirdNames{i}, ': Undirected bouts']);
    
    subplot(2,3,4);
    hold on;
    errorbar(ConditionDistances, AllBouts(i,:), AllBouts_SEM(i,:), 'ko-', 'LineWidth', 2, 'MarkerSize', 7);
    errorbar(ConditionDistances, DirectedBouts(i,:), DirectedBouts_SEM(i,:), 'ro-', 'LineWidth', 2, 'MarkerSize', 7);
    errorbar(ConditionDistances, UnDirectedBouts(i,:), UnDirectedBouts_SEM(i,:), 'bo-', 'LineWidth', 2, 'MarkerSize', 7);
    axis tight;
    Temp = axis;
    Temp = [(ConditionDistances(1) - 1) (ConditionDistances(end) + 1) 0 1.02*Temp(4)];
    axis(Temp);
    legend('All bouts', 'Directed bouts', 'Undirected bouts', 'Location', 'southeast');
    ylabel('Mean # of Motifs');
    xlabel('Distance of female (cm)');
    title([BirdNames{i}]);
    
    subplot(2,3,[5 6]);
    TempImage = imagesc(SignificanceMatrix{i});
    set(TempImage, 'AlphaData', ~isnan(SignificanceMatrix{i}));
    if (min(SignificanceMatrix{i}(:)) < 0)
        colormap([0 0 1; 0.5 0.5 0.5; 1 0 0]);
    else
        colormap([0.5 0.5 0.5; 1 0 0]);
    end
    colorbar;
    for j = 1:length(IndividualBirds(i).Conditions),
        if (isempty(strfind(IndividualBirds(i).Conditions{j}, 'UN')))
            AxisLabelString{j} = [num2str(IndividualBirds(i).ConditionDistances(j)), ' cm'];
        else
            AxisLabelString{j} = IndividualBirds(i).Conditions{j};
        end
    end
    set(gca, 'XTick', [1:1:length(IndividualBirds(i).Conditions)], 'XTickLabel', AxisLabelString);
    set(gca, 'YTick', [1:1:length(IndividualBirds(i).Conditions)], 'YTickLabel', AxisLabelString);
end


% Normalise as % of undirected and also as delta change relative to
% undirected
for i = 1:size(AllBouts,1),
    PercentChange_AllBouts(i,:) = ((AllBouts(i,:) - AllBouts(i,end)) * 100)/AllBouts(i,end);
    DeltaChange_AllBouts(i,:) = (AllBouts(i,:) - AllBouts(i,end));
end

% First for percent change
figure;
set(gcf, 'Color', 'w');
set(gcf, 'Position', [680 289 750 675]);
for i = 1:size(AllBouts,2),
    AllBout_Bar(i) = bar(ConditionDistances(i), nanmean(PercentChange_AllBouts(:,i)));
    set(AllBout_Bar(i), 'FaceColor', 'none', 'EdgeColor', 'k', 'BarWidth', 10, 'LineWidth', 2);
    hold on;
    errorbar(ConditionDistances(i), nanmean(PercentChange_AllBouts(:,i)), nanstd(PercentChange_AllBouts(:,i))/sqrt(length(find(~isnan(PercentChange_AllBouts(:,i))))), 'ks-', 'MarkerSize', 10, 'LineWidth', 2);
end
plot(ConditionDistances, nanmean(PercentChange_AllBouts), 'k', 'LineWidth', 2);

for i = 1:size(AllBouts,1),
    PercentChangeLine(i) = plot(ConditionDistances, PercentChange_AllBouts(i,:), [Symbols(mod(i,length(Symbols))+1), Colours(ceil(i/length(Symbols))), '-'], 'LineWidth', 0.5, 'MarkerSize', 5);
end
legend(PercentChangeLine, BirdNames);

clear AxisLabelString;
for j = 1:length(ConditionDistances),
    if (isempty(strfind(IndividualBirds(MaxNumConditionsIndex).Conditions{j}, 'UN')))
        AxisLabelString{j} = num2str(ConditionDistances(j));
    else
        AxisLabelString{j} = 'UN';
    end
end

set(gca, 'XTick', ConditionDistances, 'XTickLabel', AxisLabelString);
xlabel('Distance of female cage (cm)');
ylabel('% change in mean # of Motifs relative to Undirected song');
title(['All bouts: (n=', num2str(length(BirdNames)), ') birds']);
axis tight;
Temp = axis;
Temp = [(ConditionDistances(1) - 12) (ConditionDistances(end) + 12) 1.02*Temp(3) 1.2*Temp(4)];
axis(Temp);

ConditionIndexMatrix = repmat(1:1:size(PercentChange_AllBouts,2), size(PercentChange_AllBouts, 1), 1);
% Omit nan rows for doing kruskalwallis
[Nanrows, Nancols] = find(isnan(PercentChange_AllBouts));
RowsToUse = setdiff(1:1:size(PercentChange_AllBouts,1), unique(Nanrows));

DataMatrix = PercentChange_AllBouts(RowsToUse,:);
IndexMatrix = ConditionIndexMatrix(RowsToUse,:);
[p, tabl, stats] = kruskalwallis(DataMatrix(:), IndexMatrix(:), 'off');
if (p < 0.05)
    Comparisons = multcompare(stats, 'display', 'off');
    % only plotting * for anything that is significantly different from
    % undirected
    SigIndices = find((Comparisons(:,end) < 0.05) & (Comparisons(:,2) == length(ConditionDistances)));
    for j = 1:length(SigIndices),
        plot(ConditionDistances(Comparisons(SigIndices(j,1))), 1.15*(Temp(4))/1.2, 'k*', 'MarkerSize', 10, 'LineWidth', 2);
    end
end

% Now calculate correlation to the directed data alone and for the birds that have
% all distances
DataMatrix = DataMatrix(1:(length(ConditionDistances)-1),:);
IndexMatrix = IndexMatrix(1:(length(ConditionDistances) - 1), :);
[AllBird_PearsonR, AllBird_PearsonP] = corr(ConditionDistances(IndexMatrix(:))', DataMatrix(:));
[AllBird_SpearmanR, AllBirds_SpearmanP] = corr(ConditionDistances(IndexMatrix(:))', DataMatrix(:), 'type', 'Spearman', 'rows', 'complete');
text(30, Temp(4)/1.2, ['r = ', num2str(AllBird_PearsonR), '; p = ', num2str(AllBird_PearsonP)], 'FontSize', 12, 'FontWeight', 'bold');

% First for delta change
figure;
set(gcf, 'Color', 'w');
set(gcf, 'Position', [680 289 750 675]);
for i = 1:size(AllBouts,2),
    AllBout_Bar(i) = bar(ConditionDistances(i), nanmean(DeltaChange_AllBouts(:,i)));
    set(AllBout_Bar(i), 'FaceColor', 'none', 'EdgeColor', 'k', 'BarWidth', 10, 'LineWidth', 2);
    hold on;
    errorbar(ConditionDistances(i), nanmean(DeltaChange_AllBouts(:,i)), nanstd(DeltaChange_AllBouts(:,i))/sqrt(length(find(~isnan(DeltaChange_AllBouts(:,i))))), 'ks-', 'MarkerSize', 10, 'LineWidth', 2);
end
plot(ConditionDistances, nanmean(DeltaChange_AllBouts), 'k', 'LineWidth', 2);

for i = 1:size(AllBouts,1),
    DeltaChangeLine(i) = plot(ConditionDistances, DeltaChange_AllBouts(i,:), [Symbols(mod(i,length(Symbols))+1), Colours(ceil(i/length(Symbols))), '-'], 'LineWidth', 0.5, 'MarkerSize', 5);
end
legend(DeltaChangeLine, BirdNames);

clear AxisLabelString;
for j = 1:length(ConditionDistances),
    if (isempty(strfind(IndividualBirds(MaxNumConditionsIndex).Conditions{j}, 'UN')))
        AxisLabelString{j} = num2str(ConditionDistances(j));
    else
        AxisLabelString{j} = 'UN';
    end
end

set(gca, 'XTick', ConditionDistances, 'XTickLabel', AxisLabelString);
xlabel('Distance of female cage (cm)');
ylabel('% change in mean # of Motifs relative to Undirected song');
title(['All bouts: (n=', num2str(length(BirdNames)), ') birds']);
axis tight;
Temp = axis;
Temp = [(ConditionDistances(1) - 12) (ConditionDistances(end) + 12) 1.02*Temp(3) 1.2*Temp(4)];
axis(Temp);

ConditionIndexMatrix = repmat(1:1:size(DeltaChange_AllBouts,2), size(DeltaChange_AllBouts, 1), 1);
% Omit nan rows for doing kruskalwallis
[Nanrows, Nancols] = find(isnan(DeltaChange_AllBouts));
RowsToUse = setdiff(1:1:size(DeltaChange_AllBouts,1), unique(Nanrows));

DataMatrix = DeltaChange_AllBouts(RowsToUse,:);
IndexMatrix = ConditionIndexMatrix(RowsToUse,:);
[p, tabl, stats] = kruskalwallis(DataMatrix(:), IndexMatrix(:), 'off');
if (p < 0.05)
    Comparisons = multcompare(stats, 'display', 'off');
    % only plotting * for anything that is significantly different from
    % undirected
    SigIndices = find((Comparisons(:,end) < 0.05) & (Comparisons(:,2) == length(ConditionDistances)));
    for j = 1:length(SigIndices),
        plot(ConditionDistances(Comparisons(SigIndices(j,1))), 1.15*(Temp(4))/1.2, 'k*', 'MarkerSize', 10, 'LineWidth', 2);
    end
end

% Now calculate correlation to the directed data alone and for the birds that have
% all distances
DataMatrix = DataMatrix(1:(length(ConditionDistances)-1),:);
IndexMatrix = IndexMatrix(1:(length(ConditionDistances) - 1), :);
[AllBird_PearsonR, AllBird_PearsonP] = corr(ConditionDistances(IndexMatrix(:))', DataMatrix(:));
[AllBird_SpearmanR, AllBirds_SpearmanP] = corr(ConditionDistances(IndexMatrix(:))', DataMatrix(:), 'type', 'Spearman', 'rows', 'complete');
text(30, Temp(4)/1.2, ['r = ', num2str(AllBird_PearsonR), '; p = ', num2str(AllBird_PearsonP)], 'FontSize', 12, 'FontWeight', 'bold');

% Now to do the same thing for all directed only bouts

% Normalise as % of undirected and also as delta change relative to
% undirected
for i = 1:size(DirectedBouts,1),
    PercentChange_DirectedBouts(i,:) = ((DirectedBouts(i,:) - DirectedBouts(i,end)) * 100)/DirectedBouts(i,end);
    DeltaChange_DirectedBouts(i,:) = (DirectedBouts(i,:) - DirectedBouts(i,end));
end

% First for percent change
figure;
set(gcf, 'Color', 'w');
set(gcf, 'Position', [680 289 750 675]);
for i = 1:size(DirectedBouts,2),
    AllBout_Bar(i) = bar(ConditionDistances(i), nanmean(PercentChange_DirectedBouts(:,i)));
    set(AllBout_Bar(i), 'FaceColor', 'none', 'EdgeColor', 'k', 'BarWidth', 10, 'LineWidth', 2);
    hold on;
    errorbar(ConditionDistances(i), nanmean(PercentChange_DirectedBouts(:,i)), nanstd(PercentChange_DirectedBouts(:,i))/sqrt(length(find(~isnan(PercentChange_DirectedBouts(:,i))))), 'ks-', 'MarkerSize', 10, 'LineWidth', 2);
end
plot(ConditionDistances, nanmean(PercentChange_DirectedBouts), 'k', 'LineWidth', 2);

for i = 1:size(DirectedBouts,1),
    PercentChangeLine(i) = plot(ConditionDistances, PercentChange_DirectedBouts(i,:), [Symbols(mod(i,length(Symbols))+1), Colours(ceil(i/length(Symbols))), '-'], 'LineWidth', 0.5, 'MarkerSize', 5);
end
legend(PercentChangeLine, BirdNames);

clear AxisLabelString;
for j = 1:length(ConditionDistances),
    if (isempty(strfind(IndividualBirds(MaxNumConditionsIndex).Conditions{j}, 'UN')))
        AxisLabelString{j} = num2str(ConditionDistances(j));
    else
        AxisLabelString{j} = 'UN';
    end
end

set(gca, 'XTick', ConditionDistances, 'XTickLabel', AxisLabelString);
xlabel('Distance of female cage (cm)');
ylabel('% change in mean # of Motifs relative to Undirected song');
title(['Directed bouts: (n=', num2str(length(BirdNames)), ') birds']);
axis tight;
Temp = axis;
Temp = [(ConditionDistances(1) - 12) (ConditionDistances(end) + 12) 1.02*Temp(3) 1.2*Temp(4)];
axis(Temp);

ConditionIndexMatrix = repmat(1:1:size(PercentChange_DirectedBouts,2), size(PercentChange_DirectedBouts, 1), 1);
% Omit nans for doing kruskalwallis
%DataMatrix = PercentChange_DirectedBouts(RowsToUse,:);
IndexMatrix = ConditionIndexMatrix;
NotNaNIndices = find(~isnan(DataMatrix));
[p, tabl, stats] = kruskalwallis(DataMatrix(NotNaNIndices), IndexMatrix(NotNaNIndices), 'off');
if (p < 0.05)
    Comparisons = multcompare(stats, 'display', 'off');
    % only plotting * for anything that is significantly different from
    % undirected
    SigIndices = find((Comparisons(:,end) < 0.05) & (Comparisons(:,2) == length(ConditionDistances)));
    for j = 1:length(SigIndices),
        plot(ConditionDistances(Comparisons(SigIndices(j,1))), 1.15*(Temp(4))/1.2, 'k*', 'MarkerSize', 10, 'LineWidth', 2);
    end
end

% Now calculate correlation to the directed data alone and for the birds that have
% all distances
% Take only rows that have at least 3 not-nan values in the first 3 columns

[Nanrows, Nancols] = find(isnan(PercentChange_DirectedBouts(:,1:3)));
RowsToUse = setdiff(1:1:size(PercentChange_DirectedBouts,1), unique(Nanrows));

DataMatrix = DeltaChange_AllBouts(RowsToUse,1:3);
IndexMatrix = ConditionIndexMatrix(RowsToUse,1:3);

[AllBird_DirectedPearsonR, AllBird_DirectedPearsonP] = corr(ConditionDistances(IndexMatrix(:))', DataMatrix(:));
[AllBird_DirectedSpearmanR, AllBirds_DirectedSpearmanP] = corr(ConditionDistances(IndexMatrix(:))', DataMatrix(:), 'type', 'Spearman', 'rows', 'complete');
text(30, Temp(4)/1.2, ['r = ', num2str(AllBird_DirectedPearsonR), '; p = ', num2str(AllBird_DirectedPearsonP)], 'FontSize', 12, 'FontWeight', 'bold');


% Now to do the same thing for all undirected only bouts

% Normalise as % of undirected and also as delta change relative to
% undirected
for i = 1:size(UnDirectedBouts,1),
    PercentChange_UnDirectedBouts(i,:) = ((UnDirectedBouts(i,:) - UnDirectedBouts(i,end)) * 100)/UnDirectedBouts(i,end);
    DeltaChange_UnDirectedBouts(i,:) = (UnDirectedBouts(i,:) - UnDirectedBouts(i,end));
end

% First for percent change
figure;
set(gcf, 'Color', 'w');
set(gcf, 'Position', [680 289 750 675]);
for i = 1:size(UnDirectedBouts,2),
    AllBout_Bar(i) = bar(ConditionDistances(i), nanmean(PercentChange_UnDirectedBouts(:,i)));
    set(AllBout_Bar(i), 'FaceColor', 'none', 'EdgeColor', 'k', 'BarWidth', 10, 'LineWidth', 2);
    hold on;
    errorbar(ConditionDistances(i), nanmean(PercentChange_UnDirectedBouts(:,i)), nanstd(PercentChange_UnDirectedBouts(:,i))/sqrt(length(find(~isnan(PercentChange_UnDirectedBouts(:,i))))), 'ks-', 'MarkerSize', 10, 'LineWidth', 2);
end
plot(ConditionDistances, nanmean(PercentChange_UnDirectedBouts), 'k', 'LineWidth', 2);

for i = 1:size(UnDirectedBouts,1),
    PercentChangeLine(i) = plot(ConditionDistances, PercentChange_UnDirectedBouts(i,:), [Symbols(mod(i,length(Symbols))+1), Colours(ceil(i/length(Symbols))), '-'], 'LineWidth', 0.5, 'MarkerSize', 5);
end
legend(PercentChangeLine, BirdNames);

clear AxisLabelString;
for j = 1:length(ConditionDistances),
    if (isempty(strfind(IndividualBirds(MaxNumConditionsIndex).Conditions{j}, 'UN')))
        AxisLabelString{j} = num2str(ConditionDistances(j));
    else
        AxisLabelString{j} = 'UN';
    end
end

set(gca, 'XTick', ConditionDistances, 'XTickLabel', AxisLabelString);
xlabel('Distance of female cage (cm)');
ylabel('% change in mean # of Motifs relative to Undirected song');
title(['UnDirected bouts: (n=', num2str(length(BirdNames)), ') birds']);
axis tight;
Temp = axis;
Temp = [(ConditionDistances(1) - 12) (ConditionDistances(end) + 12) 1.02*Temp(3) 1.2*Temp(4)];
axis(Temp);

ConditionIndexMatrix = repmat(1:1:size(PercentChange_UnDirectedBouts,2), size(PercentChange_UnDirectedBouts, 1), 1);
% Omit nans for doing kruskalwallis
%DataMatrix = PercentChange_UnDirectedBouts(RowsToUse,:);
IndexMatrix = ConditionIndexMatrix;
NotNaNIndices = find(~isnan(DataMatrix));
[p, tabl, stats] = kruskalwallis(DataMatrix(NotNaNIndices), IndexMatrix(NotNaNIndices), 'off');
if (p < 0.05)
    Comparisons = multcompare(stats, 'display', 'off');
    % only plotting * for anything that is significantly different from
    % undirected
    SigIndices = find((Comparisons(:,end) < 0.05) & (Comparisons(:,2) == length(ConditionDistances)));
    for j = 1:length(SigIndices),
        plot(ConditionDistances(Comparisons(SigIndices(j,1))), 1.15*(Temp(4))/1.2, 'k*', 'MarkerSize', 10, 'LineWidth', 2);
    end
end

% Now calculate correlation to the directed data alone and for the birds that have
% all distances
% Take only rows that have at least 3 not-nan values in the last 3 female
% presentation columns

[Nanrows, Nancols] = find(isnan(PercentChange_UnDirectedBouts(:,(length(ConditionDistances) - 3):(length(ConditionDistances) - 1))));
RowsToUse = setdiff(1:1:size(PercentChange_UnDirectedBouts,1), unique(Nanrows));

DataMatrix = DeltaChange_AllBouts(RowsToUse,(length(ConditionDistances) - 3):(length(ConditionDistances) - 1));
IndexMatrix = ConditionIndexMatrix(RowsToUse,(length(ConditionDistances) - 3):(length(ConditionDistances) - 1));

[AllBird_DirectedPearsonR, AllBird_DirectedPearsonP] = corr(ConditionDistances(IndexMatrix(:))', DataMatrix(:));
[AllBird_DirectedSpearmanR, AllBirds_DirectedSpearmanP] = corr(ConditionDistances(IndexMatrix(:))', DataMatrix(:), 'type', 'Spearman', 'rows', 'complete');
text(30, Temp(4)/1.2, ['r = ', num2str(AllBird_DirectedPearsonR), '; p = ', num2str(AllBird_DirectedPearsonP)], 'FontSize', 12, 'FontWeight', 'bold');
disp('Finished plotting');