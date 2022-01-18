function [L0Dir_Undir_PValue] = Harini_AllFeatureStats_Plots_TopMatched(DirBout_Stats, DirUndirBout_Stats, UndirBout_Stats, StructString, MinTrialNo, Distances, LabelString, BirdNames)

FinalFigureDir = '/home/raghav/StudentRelated/Harini/Manuscript/results';

% Now to do the stats
% Now, we want to normalize all of the features across all the data for
% each bird. Then plot the normalized changes for directed songs for all
% features for a given bird in the same plot.
% Normalization is just a z-score across all data for that bird and that
% feature

FeaturesToConsider = {'TotalINNumber_500ms' 'CompleteMotifNumber' 'FirstMotifDuration'};

Colours = 'rbkcmg';

for i = 1:length(DirBout_Stats)
    disp(['Bird #', num2str(i), ': ', BirdNames{i}]);
    GroupData{i} = [];
    % GroupData has 3 columns, the first one has the data corresponding to
    % the variable that that I'm looking at like TotalINNumber, etc. 2nd
    % column has the index of the distance at which the recordings were
    % done (0 for L0, 1 for L1 ... 6 for UN) and 3rd column has 1 for Dir,
    % 2 for Dir/Undir - DUN - and 3 for Undir
    
    for j = 1:length(DirBout_Stats{i}),
        if (~isempty(DirBout_Stats{i}{j}))
            TempGroupData = [];
            for k = 1:length(FeaturesToConsider),
                if (k == length(FeaturesToConsider))
                    TempGroupData = [TempGroupData [eval(['DirBout_Stats{i}{j}.', FeaturesToConsider{k}, '(:)']) ones(size(eval(['DirBout_Stats{i}{j}.', FeaturesToConsider{k}, '(:)'])))*j ones(size(eval(['DirBout_Stats{i}{j}.', FeaturesToConsider{k}, '(:)'])))*1]];
                else
                    TempGroupData = [TempGroupData eval(['DirBout_Stats{i}{j}.', FeaturesToConsider{k}, '(:)'])];
                end
            end
            GroupData{i} = [GroupData{i}; TempGroupData];
        end
    end
    
    for j = 1:length(DirUndirBout_Stats{i}),
        if (~isempty(DirUndirBout_Stats{i}{j}))
            TempGroupData = [];
            for k = 1:length(FeaturesToConsider),
                if (k == length(FeaturesToConsider))
                    TempGroupData = [TempGroupData [eval(['DirUndirBout_Stats{i}{j}.', FeaturesToConsider{k}, '(:)']) ones(size(eval(['DirUndirBout_Stats{i}{j}.', FeaturesToConsider{k}, '(:)'])))*j ones(size(eval(['DirUndirBout_Stats{i}{j}.', FeaturesToConsider{k}, '(:)'])))*2]];
                else
                    TempGroupData = [TempGroupData eval(['DirUndirBout_Stats{i}{j}.', FeaturesToConsider{k}, '(:)'])];
                end
            end
            GroupData{i} = [GroupData{i}; TempGroupData];
        end
    end
    
    for j = 1:length(UndirBout_Stats{i}),
        if (~isempty(UndirBout_Stats{i}{j}))
            TempGroupData = [];
            for k = 1:length(FeaturesToConsider),
                if (k == length(FeaturesToConsider))
                    TempGroupData = [TempGroupData [eval(['UndirBout_Stats{i}{j}.', FeaturesToConsider{k}, '(:)']) ones(size(eval(['UndirBout_Stats{i}{j}.', FeaturesToConsider{k}, '(:)'])))*j ones(size(eval(['UndirBout_Stats{i}{j}.', FeaturesToConsider{k}, '(:)'])))*3]];
                else
                    TempGroupData = [TempGroupData eval(['UndirBout_Stats{i}{j}.', FeaturesToConsider{k}, '(:)'])];
                end
            end
            GroupData{i} = [GroupData{i}; TempGroupData];
        end
    end

    % Remove NaN values
    NaNIndices = find(isnan(GroupData{i}(:,1)));
    if (~isempty(NaNIndices))
        GroupData{i}(NaNIndices,:) = [];
    end
    
    GroupData{i}(:,1:length(FeaturesToConsider)) = zscore(GroupData{i}(:,1:length(FeaturesToConsider)));
    
    % Now plot the distributions and means for all birds across all
    % distances
end

% Now to plot the mean and std across distances for D, D/UN and UN for all
% birds together as percent change relative to undirected.
% This is done only for birds where the means are already different between
% L0 directed and undirected as measured by the earlier p-value

DistanceDirMeans = [];
DistanceDirUnDirMeans = [];
DistanceUnDirMeans = [];
Symbols = 'so^+hp><*dxv.so';

Index = 0;

DirUnDirComparisonFig = figure;
set(gcf, 'Color', 'w');
set(gcf, 'Position', [188 208 1600 800]);
for i = 1:length(BirdNames),
    Index = i;
    L0DirSongIndices = find((GroupData{i}(:,length(FeaturesToConsider)+2) == 1) & (GroupData{i}(:,length(FeaturesToConsider)+1) == 1));
    UnDirUnDirSongIndices = find((GroupData{i}(:,length(FeaturesToConsider)+2) == 3) & (GroupData{i}(:,length(FeaturesToConsider)+1) == 6));
    
    % First for directed songs
    for j = 1:6,
        DirSongIndices = find((GroupData{i}(:,length(FeaturesToConsider)+2) == 1) & (GroupData{i}(:,length(FeaturesToConsider)+1) == j));
        if (length(DirSongIndices) >= MinTrialNo)
            for k = 1:length(FeaturesToConsider),
                DistanceDirMeans{k}(Index, j) = mean(GroupData{i}(DirSongIndices,k));
                DistanceDirSEMs{k}(Index, j) = std(GroupData{i}(DirSongIndices,k))/sqrt(length(DirSongIndices));
%                DistanceDirCIs{k}{Index}(j,:) = bootci(10000,@mean, GroupData{i}(DirSongIndices,k));
            end
        else
            for k = 1:length(FeaturesToConsider),
                DistanceDirMeans{k}(Index, j) = NaN;
                DistanceDirSEMs{k}(Index,j) = NaN;
            end
        end
        
        DirUnDirSongIndices = find((GroupData{i}(:,length(FeaturesToConsider)+2) == 2) & (GroupData{i}(:,length(FeaturesToConsider)+1) == j));
        if (length(DirUnDirSongIndices) >= MinTrialNo)
            for k = 1:length(FeaturesToConsider),
                DistanceDirUnDirMeans{k}(Index, j) = mean(GroupData{i}(DirUnDirSongIndices,k));
                DistanceDirUnDirSEMs{k}(Index, j) = std(GroupData{i}(DirUnDirSongIndices,k))/sqrt(length(DirUnDirSongIndices));
            end
        else
            for k = 1:length(FeaturesToConsider),
                DistanceDirUnDirMeans{k}(Index, j) = NaN;
                DistanceDirUnDirSEMs{k}(Index, j) = NaN;
            end
        end
        
        UnDirSongIndices = find((GroupData{i}(:,length(FeaturesToConsider)+2) == 3) & (GroupData{i}(:,length(FeaturesToConsider)+1) == j));
        if (length(UnDirSongIndices) >= MinTrialNo)
            for k = 1:length(FeaturesToConsider),
                DistanceUnDirMeans{k}(Index, j) = mean(GroupData{i}(UnDirSongIndices,k));
                DistanceUnDirSEMs{k}(Index, j) = std(GroupData{i}(UnDirSongIndices,k))/sqrt(length(UnDirSongIndices));
            end
        else
            for k = 1:length(FeaturesToConsider),
                DistanceUnDirMeans{k}(Index, j) = NaN;
                DistanceUnDirSEMs{k}(Index, j) = NaN;
            end
        end
    end
    
    figure(DirUnDirComparisonFig);
    subplot(3, ceil(length(DirBout_Stats)/3),Index);
    hold on;
    for k = 1:length(FeaturesToConsider),
        plot(Distances, DistanceDirMeans{k}(Index,:), [Colours(k), 's-'], 'LineWidth', 1);
        %plot(repmat(Distances(:), 1, 2)', DistanceDirCIs{k}{Index}', Colours(k))
      
        % plot(Distances, DistanceUnDirMeans{k}(Index,:), 'bs-', 'LineWidth', 1);
        %plot(repmat(Distances(:), 1, 2)', DistanceUnDirCIs{Index}', 'b');
        %plot(repmat(Distances(:),1,2), repmat(DistanceUnDirCIs{k}{Index}(end,:), length(Distances), 1), [Colours(k), '--']);
    end
    
    xlabel('Distance from female (cm)');
    axis tight;
    Temp = axis;
    Temp = [-0.5 250 0.95*Temp(3) 1.05*Temp(4)];
    axis(Temp);
    ylabel(LabelString);
    title(BirdNames{i}, 'Color', 'r');
    % legend(FeaturesToConsider);
end

YLabelStrings = {'# of INs (zscore)', '# of motifs/bout (zscore)', 'First motif duration (zscore)'};

figure;
set(gcf, 'Color', 'w');
set(gcf, 'Position', [680 162 750 1050]);
p = panel();
p.pack(3,2);
p.de.margin = 15;

% Plot dir undir group data for 
for i = 1:3,
    clear MeanSEM;
    p(i,1).select();
    hold on;
    plot(Distances(1:5), DistanceDirMeans{i}(:,1:5)', 'ko-', 'LineWidth', 0.25, 'MarkerSize', 4);
    for j = 1:5,
        if (~isempty(find(~isnan(DistanceDirMeans{i}(:,j)))))
            MeanSEM(j,:) = [nanmean(DistanceDirMeans{i}(:,j)) nanstd(DistanceDirMeans{i}(:,j))/sqrt(length(find(~isnan(DistanceDirMeans{i}(:,j)))))];
        else
            MeanSEM(j,:) = [NaN NaN];
        end
    end
    errorbar(Distances(1:5)+1, MeanSEM(:,1), MeanSEM(:,2), 'rs-', 'LineWidth', 1.5, 'MarkerSize', 8);
    plot(Distances(end), DistanceUnDirMeans{i}(:,end), 'ko', 'MarkerSize', 4, 'LineWidth', 0.25);
    errorbar(Distances(end)-10, nanmean(DistanceUnDirMeans{i}(:,end)), nanstd(DistanceUnDirMeans{i}(:,end))/sqrt(length(find(~isnan(DistanceUnDirMeans{i}(:,end))))), 'ks-', 'LineWidth', 1.5, 'MarkerSize', 8);
    axis tight;
    Temp1 = axis;
    
    set(gca, 'XTick', [0 20 60 110 165 195], 'XTickLabel', {'0' '20' '60' '110' '165' 'UN'});
    clear MeanSEM;
    ylabel(YLabelStrings{i});
    set(gca, 'YTick', [-3:1:3]);
    
    p(i,2).select();
    hold on;
    plot(Distances(1:5), DistanceUnDirMeans{i}(:,1:5)', 'ko-', 'LineWidth', 0.25, 'MarkerSize', 4);
    for j = 1:5,
        if (~isempty(find(~isnan(DistanceUnDirMeans{i}(:,j)))))
            MeanSEM(j,:) = [nanmean(DistanceUnDirMeans{i}(:,j)) nanstd(DistanceUnDirMeans{i}(:,j))/sqrt(length(find(~isnan(DistanceUnDirMeans{i}(:,j)))))];
        else
            MeanSEM(j,:) = [NaN NaN];
        end
    end
    errorbar(Distances(1:5)+1, MeanSEM(:,1), MeanSEM(:,2), 'bs-', 'LineWidth', 1.5, 'MarkerSize', 8);
    plot(Distances(end), DistanceUnDirMeans{i}(:,end), 'ko', 'MarkerSize', 4, 'LineWidth', 0.25);
    errorbar(Distances(end)-10, nanmean(DistanceUnDirMeans{i}(:,end)), nanstd(DistanceUnDirMeans{i}(:,end))/sqrt(length(find(~isnan(DistanceUnDirMeans{i}(:,end))))), 'ks-', 'LineWidth', 1.5, 'MarkerSize', 8);
    axis tight;
    Temp2 = axis;
    Temp = [-10 Distances(end)+5 1.05*min(Temp1(3), Temp2(3)) 1.1*max(Temp1(4), Temp2(4))];
    axis(Temp);
    plot(Temp(1:2), [0 0], 'k--');
    set(gca, 'XTick', [0 20 60 110 165 195], 'XTickLabel', {'0' '20' '60' '110' '165' 'UN'});
    clear MeanSEM;
    set(gca, 'YTick', [-3:1:3]);
    
    p(i,1).select();
    axis(Temp);
    plot(Temp(1:2), [0 0], 'k--');
    
    % Now to do the correlation and put the r and p values in the corner
    [CorrR, CorrP] = corr(Distances(reshape(repmat(1:1:5, 14, 1), 14*5, 1))', reshape(DistanceDirMeans{i}(:,1:5), 14*5, 1), 'rows', 'complete', 'type', 'Pearson');
    text(130, Temp(4)/1.08, {['r = ', num2str(CorrR)]; ['p = ', num2str(CorrP)]});
    if (i == 1)
        title('DIRECTED SONGS', 'Color', 'r', 'FontSize', 12);
    end
    
    p(i,2).select();
    % Now to do the correlation and put the r and p values in the corner
    [CorrR, CorrP] = corr(Distances(reshape(repmat(1:1:5, 14, 1), 14*5, 1))', reshape(DistanceUnDirMeans{i}(:,1:5), 14*5, 1), 'rows', 'complete', 'type', 'Pearson');
    text(130, Temp(4)/1.08, {['r = ', num2str(CorrR)]; ['p = ', num2str(CorrP)]});
    if (i == 1)
        title('UNDIRECTED SONGS', 'Color', 'b', 'FontSize', 12);
    end
end

p.fontsize = 12;
p.marginleft = 20;
p.margintop = 15;
set(gcf, 'ReSize', 'off');
set(gcf, 'PaperPositionMode', 'auto');
print(fullfile(FinalFigureDir, ['Fig.4.GroupDataZScores.', LabelString, '.eps']), '-depsc2', '-r600');
print(fullfile(FinalFigureDir, ['Fig.4.GroupDataZScores.', LabelString, '.png']), '-dpng', '-r600');    



disp('Done with stats and plots');
