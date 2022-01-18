function [L0Dir_Undir_PValue] = Harini_AllBoutMotifSyllFeatureAnalysis(DirBout_Stats, DirUndirBout_Stats, UndirBout_Stats, StructString, MinTrialNo, Distances, LabelString, BirdNames)

FinalFigureDir = '/home/raghav/StudentRelated/Harini/Manuscript/results';

% We want to do the analysis per bird and for the group as a whole. The
% things that we will analyze are (1) # of INs at the beginning based on
% the 500ms criterion, (2) # of complete motifs in the bout, (3) all and
% first motif duration, (4), CV of FF and (5) amplitude of individual
% motif syllables only for HF mic birds

FeaturesToConsider = {'TotalINNumber_500ms' 'CompleteMotifNumber' 'FirstMotifDuration'};

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
    for k = 1:length(FeaturesToConsider),
        NaNIndices = find(isnan(GroupData{i}(:,k)));
        if (~isempty(NaNIndices))
            GroupData{i}(NaNIndices,:) = [];
        end
    end
    
    UnNormalizedGroupData{i} = GroupData{i};
    DirIndices = find(GroupData{i}(:,end) == 1);
    [r, p] = corrcoef(GroupData{i}(DirIndices, 1:(length(FeaturesToConsider)+1)));
    
    Corr_R_Values(i,:) = r(end,1:length(FeaturesToConsider));
    Corr_P_Values(i,:) = p(end, 1:length(FeaturesToConsider));

    for Feat = 1:length(FeaturesToConsider),
        fun = @(x)sseval(x, Distances(GroupData{i}(DirIndices,end-1)), GroupData{i}(DirIndices, Feat), 'linear');
        bestx = fminsearch(fun, rand(2,1));
        [Rsq, F, Prob] = CalculateGoodnessofLinearFit(bestx, Distances(GroupData{i}(DirIndices, 4))', GroupData{i}(DirIndices, 1));
        Fmin_Corr_Rsq_Values(i,Feat) = Rsq;
        Fmin_Corr_P_Values(i,Feat) = Prob;
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

UnNormalizedDistanceDirMeans = [];
UnNormalizedDistanceDirUnDirMeans = [];
UnNormalizedDistanceUnDirMeans = [];
Symbols = 'so^+hp><*dxv.so';

Index = 0;

for i = 1:length(BirdNames),
    Index = i;
    
    % First for directed songs
    for j = 1:6,
        DirSongIndices = find((GroupData{i}(:,length(FeaturesToConsider)+2) == 1) & (GroupData{i}(:,length(FeaturesToConsider)+1) == j));
        if (length(DirSongIndices) >= MinTrialNo)
            for k = 1:length(FeaturesToConsider),
                DistanceDirMeans{k}(Index, j) = mean(GroupData{i}(DirSongIndices,k));
                DistanceDirSEMs{k}(Index, j) = std(GroupData{i}(DirSongIndices,k))/sqrt(length(DirSongIndices));
                UnNormalizedDistanceDirMeans{k}(Index, j) = mean(UnNormalizedGroupData{i}(DirSongIndices,k));
                UnNormalizedDistanceDirSEMs{k}(Index, j) = std(UnNormalizedGroupData{i}(DirSongIndices,k))/sqrt(length(DirSongIndices));
            end
        else
            for k = 1:length(FeaturesToConsider),
                DistanceDirMeans{k}(Index, j) = NaN;
                DistanceDirSEMs{k}(Index,j) = NaN;
                UnNormalizedDistanceDirMeans{k}(Index, j) = NaN;
                UnNormalizedDistanceDirSEMs{k}(Index,j) = NaN;
            end
        end
        
        DirUnDirSongIndices = find((GroupData{i}(:,length(FeaturesToConsider)+2) == 2) & (GroupData{i}(:,length(FeaturesToConsider)+1) == j));
        if (length(DirUnDirSongIndices) >= MinTrialNo)
            for k = 1:length(FeaturesToConsider),
                DistanceDirUnDirMeans{k}(Index, j) = mean(GroupData{i}(DirUnDirSongIndices,k));
                DistanceDirUnDirSEMs{k}(Index, j) = std(GroupData{i}(DirUnDirSongIndices,k))/sqrt(length(DirUnDirSongIndices));
                UnNormalizedDistanceDirUnDirMeans{k}(Index, j) = mean(UnNormalizedGroupData{i}(DirUnDirSongIndices,k));
                UnNormalizedDistanceDirUnDirSEMs{k}(Index, j) = std(UnNormalizedGroupData{i}(DirUnDirSongIndices,k))/sqrt(length(DirUnDirSongIndices));
            end
        else
            for k = 1:length(FeaturesToConsider),
                DistanceDirUnDirMeans{k}(Index, j) = NaN;
                DistanceDirUnDirSEMs{k}(Index, j) = NaN;
                UnNormalizedDistanceDirUnDirMeans{k}(Index, j) = NaN;
                UnNormalizedDistanceDirUnDirSEMs{k}(Index, j) = NaN;
            end
        end
        
        UnDirSongIndices = find((GroupData{i}(:,length(FeaturesToConsider)+2) == 3) & (GroupData{i}(:,length(FeaturesToConsider)+1) == j));
        if (length(UnDirSongIndices) >= MinTrialNo)
            for k = 1:length(FeaturesToConsider),
                DistanceUnDirMeans{k}(Index, j) = mean(GroupData{i}(UnDirSongIndices,k));
                DistanceUnDirSEMs{k}(Index, j) = std(GroupData{i}(UnDirSongIndices,k))/sqrt(length(UnDirSongIndices));
                UnNormalizedDistanceUnDirMeans{k}(Index, j) = mean(UnNormalizedGroupData{i}(UnDirSongIndices,k));
                UnNormalizedDistanceUnDirSEMs{k}(Index, j) = std(UnNormalizedGroupData{i}(UnDirSongIndices,k))/sqrt(length(UnDirSongIndices));
            end
        else
            for k = 1:length(FeaturesToConsider),
                DistanceUnDirMeans{k}(Index, j) = NaN;
                DistanceUnDirSEMs{k}(Index, j) = NaN;
                UnNormalizedDistanceUnDirMeans{k}(Index, j) = NaN;
                UnNormalizedDistanceUnDirSEMs{k}(Index, j) = NaN;
            end
        end
    end
end

Colours = 'rbkcmg';

MinTrialNo = 10;

for i = 1:length(DirBout_Stats)
    disp(['Bird #', num2str(i), ': ', BirdNames{i}]);
    % GroupData has 3 columns, the first one has the data corresponding to
    % the variable that that I'm looking at like TotalINNumber, etc. 2nd
    % column has the index of the distance at which the recordings were
    % done (0 for L0, 1 for L1 ... 6 for UN) and 3rd column has 1 for Dir,
    % 2 for Dir/Undir - DUN - and 3 for Undir
    if (~isfield(UndirBout_Stats{i}{end}, 'FF_SyllLabel'))
        continue;
    end
    for SyllIndex = 1:length(UndirBout_Stats{i}{end}.FF_SyllLabel),
        FFGroupData{i}{SyllIndex} = [];
        for j = 1:length(DirBout_Stats{i}),
            if (~isempty(DirBout_Stats{i}{j}))
                FFGroupData{i}{SyllIndex} = [FFGroupData{i}{SyllIndex}; [DirBout_Stats{i}{j}.FF{SyllIndex}(:) ones(size(DirBout_Stats{i}{j}.FF{SyllIndex}(:)))*j ones(size(DirBout_Stats{i}{j}.FF{SyllIndex}(:)))*1]];
            end
        end
    
        for j = 1:length(DirUndirBout_Stats{i}),
            if (~isempty(DirUndirBout_Stats{i}{j}))
                FFGroupData{i}{SyllIndex} = [FFGroupData{i}{SyllIndex}; [DirUndirBout_Stats{i}{j}.FF{SyllIndex}(:) ones(size(DirUndirBout_Stats{i}{j}.FF{SyllIndex}(:)))*j ones(size(DirUndirBout_Stats{i}{j}.FF{SyllIndex}(:)))*2]];
            end
        end
        
        for j = 1:length(UndirBout_Stats{i}),
            if (~isempty(UndirBout_Stats{i}{j}))
                FFGroupData{i}{SyllIndex} = [FFGroupData{i}{SyllIndex}; [UndirBout_Stats{i}{j}.FF{SyllIndex}(:) ones(size(UndirBout_Stats{i}{j}.FF{SyllIndex}(:)))*j ones(size(UndirBout_Stats{i}{j}.FF{SyllIndex}(:)))*3]];
            end
        end

        % Remove NaN values
        NaNIndices = find(isnan(FFGroupData{i}{SyllIndex}(:,1)));
        if (~isempty(NaNIndices))
            FFGroupData{i}{SyllIndex}(NaNIndices,:) = [];
        end
    
        % FFGroupData{i}{SyllIndex}(:,1) = zscore(FFGroupData{i}{SyllIndex}(:,1));
    
    end
end

DistanceFFDirMeans = [];
DistanceFFDirUnDirMeans = [];
DistanceFFUnDirMeans = [];
DistanceFFDirCVs = [];
DistanceFFDirUnDirCVs = [];
DistanceFFUnDirCVs = [];

Index = 0;

for i = 1:length(BirdNames),
    Index = i;
    for SyllIndex = 1:length(FFGroupData{i}),
        % First for directed songs
        for j = 1:6,
            DirSongIndices = find((FFGroupData{i}{SyllIndex}(:,3) == 1) & (FFGroupData{i}{SyllIndex}(:,2) == j));
            if (length(DirSongIndices) >= MinTrialNo)
                DistanceFFDirMeans{i}(SyllIndex, j) = mean(FFGroupData{i}{SyllIndex}(DirSongIndices,1));
                DistanceFFDirCVs{i}(SyllIndex, j) = std(FFGroupData{i}{SyllIndex}(DirSongIndices,1))/mean(FFGroupData{i}{SyllIndex}(DirSongIndices,1));
                DistanceFFDirSEMs{i}(SyllIndex, j) = std(FFGroupData{i}{SyllIndex}(DirSongIndices,1))/sqrt(length(DirSongIndices));
            else
                DistanceFFDirMeans{i}(SyllIndex, j) = NaN;
                DistanceFFDirCVs{i}(SyllIndex, j) = NaN;
                DistanceFFDirSEMs{i}(SyllIndex,j) = NaN;
            end

            DirUnDirSongIndices = find((FFGroupData{i}{SyllIndex}(:,3) == 2) & (FFGroupData{i}{SyllIndex}(:,2) == j));
            if (length(DirUnDirSongIndices) >= MinTrialNo)
                DistanceFFDirUnDirMeans{i}(SyllIndex, j) = mean(FFGroupData{i}{SyllIndex}(DirUnDirSongIndices,1));
                DistanceFFDirUnDirCVs{i}(SyllIndex, j) = std(FFGroupData{i}{SyllIndex}(DirUnDirSongIndices,1))/mean(FFGroupData{i}{SyllIndex}(DirUnDirSongIndices,1));
                DistanceFFDirUnDirSEMs{i}(SyllIndex, j) = std(FFGroupData{i}{SyllIndex}(DirUnDirSongIndices,1))/sqrt(length(DirUnDirSongIndices));
            else
                DistanceFFDirUnDirMeans{i}(SyllIndex, j) = NaN;
                DistanceFFDirUnDirCVs{i}(SyllIndex, j) = NaN;
                DistanceFFDirUnDirSEMs{i}(SyllIndex,j) = NaN;
            end
            
            UnDirSongIndices = find((FFGroupData{i}{SyllIndex}(:,3) == 3) & (FFGroupData{i}{SyllIndex}(:,2) == j));
            if (length(UnDirSongIndices) >= MinTrialNo)
                DistanceFFUnDirMeans{i}(SyllIndex, j) = mean(FFGroupData{i}{SyllIndex}(UnDirSongIndices,1));
                DistanceFFUnDirCVs{i}(SyllIndex, j) = std(FFGroupData{i}{SyllIndex}(UnDirSongIndices,1))/mean(FFGroupData{i}{SyllIndex}(UnDirSongIndices,1));
                DistanceFFUnDirSEMs{i}(SyllIndex, j) = std(FFGroupData{i}{SyllIndex}(UnDirSongIndices,1))/sqrt(length(UnDirSongIndices));
            else
                DistanceFFUnDirMeans{i}(SyllIndex, j) = NaN;
                DistanceFFUnDirCVs{i}(SyllIndex, j) = NaN;
                DistanceFFUnDirSEMs{i}(SyllIndex,j) = NaN;
            end
        end
    end
    if (~isempty(FFGroupData{i}))
        if (size(DistanceFFDirCVs{i}, 1) > 1)
            AllBirdDistanceFFDirCVs(i,:) = mean(DistanceFFDirCVs{i});
            AllBirdDistanceFFUnDirCVs(i,:) = mean(DistanceFFUnDirCVs{i});
        else
            AllBirdDistanceFFDirCVs(i,:) = DistanceFFDirCVs{i};
            AllBirdDistanceFFUnDirCVs(i,:) = DistanceFFUnDirCVs{i};
        end
    else
        AllBirdDistanceFFDirCVs(i,:) = ones(1,6)*NaN;
        AllBirdDistanceFFUnDirCVs(i,:) = ones(1,6)*NaN;
    end
end

YLabelStrings = {'# of INs (zscore)', '# of motifs/bout (zscore)', 'First motif duration (zscore)'};

DistancesMatrix = repmat(Distances(1:5), size(DistanceDirMeans{1},1), 1);

figure;
set(gcf, 'Color', 'w');
set(gcf, 'Position', [680 162 440 1050]);
p = panel();
p.pack(4,1);
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
    errorbar(Distances(1:5)+7, MeanSEM(:,1), MeanSEM(:,2), 'rs', 'LineWidth', 1.5, 'MarkerSize', 8);
    % Fit a line to the means and plot the fit line
    Coeffs = polyfit(DistancesMatrix(find(~isnan(DistanceDirMeans{i}(:,1:5)))), DistanceDirMeans{i}(find(~isnan(DistanceDirMeans{i}(:,1:5)))), 1);
    plot(Distances(1:5), polyval(Coeffs, Distances(1:5)), 'r', 'LineWidth', 1.5);
    
    plot(Distances(end), DistanceUnDirMeans{i}(:,end), 'ko', 'MarkerSize', 4, 'LineWidth', 0.25);
    errorbar(Distances(end)-10, nanmean(DistanceUnDirMeans{i}(:,end)), nanstd(DistanceUnDirMeans{i}(:,end))/sqrt(length(find(~isnan(DistanceUnDirMeans{i}(:,end))))), 'ks-', 'LineWidth', 1.5, 'MarkerSize', 8);
    axis tight;
    Temp = axis;
    Temp = [-10 Distances(end)+5 1.15*Temp(3) 1.15*Temp(4)];
    axis(Temp);
    
    set(gca, 'XTick', [0 20 60 110 165 195], 'XTickLabel', {'0' '20' '60' '110' '165' 'UN'});
    clear MeanSEM;
    ylabel(YLabelStrings{i});
    set(gca, 'YTick', [-3:1:3]);
    
    plot(Temp(1:2), nanmean(DistanceUnDirMeans{i}(:,end))*ones(1,2), 'k--');
    
    % Now to do the correlation and put the r and p values in the corner
    [CorrR, CorrP] = corr(Distances(reshape(repmat(1:1:5, length(BirdNames), 1), length(BirdNames)*5, 1))', reshape(DistanceDirMeans{i}(:,1:5), length(BirdNames)*5, 1), 'rows', 'complete', 'type', 'Pearson');
    text(130, Temp(4)/1.08, {['r = ', num2str(CorrR)]; ['p = ', num2str(CorrP)]});
end

clear MeanSEM;
p(4,1).select();
hold on;
plot(Distances(1:5), AllBirdDistanceFFDirCVs(:,1:5)', 'ko-', 'LineWidth', 0.25, 'MarkerSize', 4);
for j = 1:5,
    if (~isempty(find(~isnan(AllBirdDistanceFFDirCVs(:,j)))))
        MeanSEM(j,:) = [nanmean(AllBirdDistanceFFDirCVs(:,j)) nanstd(AllBirdDistanceFFDirCVs(:,j))/sqrt(length(find(~isnan(AllBirdDistanceFFDirCVs(:,j)))))];
    else
        MeanSEM(j,:) = [NaN NaN];
    end
end
errorbar(Distances(1:5)+7, MeanSEM(:,1), MeanSEM(:,2), 'rs', 'LineWidth', 1.5, 'MarkerSize', 8);
% Fit a line to the means and plot the fit line
Coeffs = polyfit(DistancesMatrix(find(~isnan(AllBirdDistanceFFDirCVs(:,1:5)))), AllBirdDistanceFFDirCVs(find(~isnan(AllBirdDistanceFFDirCVs(:,1:5)))), 1);
plot(Distances(1:5), polyval(Coeffs, Distances(1:5)), 'r', 'LineWidth', 1.5);
    
plot(Distances(end), AllBirdDistanceFFUnDirCVs(:,end), 'ko', 'MarkerSize', 4, 'LineWidth', 0.25);
errorbar(Distances(end)-10, nanmean(AllBirdDistanceFFUnDirCVs(:,end)), nanstd(AllBirdDistanceFFUnDirCVs(:,end))/sqrt(length(find(~isnan(AllBirdDistanceFFUnDirCVs(:,end))))), 'ks-', 'LineWidth', 1.5, 'MarkerSize', 8);
axis tight;
Temp = axis;
Temp = [-10 Distances(end)+5 0 1.25*Temp(4)];
axis(Temp);
plot(Temp(1:2), ones(1,2)*nanmean(AllBirdDistanceFFUnDirCVs(:,end)), 'k--');
set(gca, 'XTick', [0 20 60 110 165 195], 'XTickLabel', {'0' '20' '60' '110' '165' 'UN'});
clear MeanSEM;
xlabel('Distance from female (cm)');
ylabel('CV of FF');

% Now to do the correlation and put the r and p values in the corner
[CorrR, CorrP] = corr(Distances(reshape(repmat(1:1:5, length(BirdNames), 1), length(BirdNames)*5, 1))', reshape(AllBirdDistanceFFDirCVs(:,1:5), length(BirdNames)*5, 1), 'rows', 'complete', 'type', 'Pearson');
text(130, Temp(4)/1.08, {['r = ', num2str(CorrR)]; ['p = ', num2str(CorrP)]});

p.fontsize = 12;
p.marginleft = 25;
p.margintop = 15;
set(gcf, 'ReSize', 'off');
set(gcf, 'PaperPositionMode', 'auto');

% Now to align the y-axis labels for Rasters and PSTs

for i = 1:3,
    p(4,1).select();
    RasterYLabelPos = get(get(gca, 'ylabel'), 'Position');
    
    p(i,1).select();
    PSTYLabelPos = get(get(gca, 'ylabel'), 'Position');
    PSTYLabelPos(1) = RasterYLabelPos(1);
    set(get(gca, 'ylabel'), 'Position', PSTYLabelPos);
end

print(fullfile(FinalFigureDir, ['GroupDataZScores_OnlyDir.', LabelString, '.eps']), '-depsc2', '-r600');
print(fullfile(FinalFigureDir, ['GroupDataZScores_OnlyDir.', LabelString, '.png']), '-dpng', '-r600');    

% Now plot a version where I plot the absolute differences between the
% values and undirected songs in the absence of the female. This should
% help to show that it does become more like undirected song

for i = 1:length(DistanceDirMeans),
    FeatureDiffs{i} = abs(DistanceDirMeans{i} - repmat(DistanceUnDirMeans{i}(:,end), 1, size(DistanceDirMeans{i},2)));
end

FFCVDiffs = (AllBirdDistanceFFDirCVs - repmat(AllBirdDistanceFFUnDirCVs(:,end), 1, size(AllBirdDistanceFFDirCVs, 2))).^2;

figure;
set(gcf, 'Color', 'w');
set(gcf, 'Position', [680 162 440 1050]);
p = panel();
p.pack(4,1);
p.de.margin = 15;

% Plot dir undir group data for 
for i = 1:3,
    clear MeanSEM;
    p(i,1).select();
    hold on;
    plot(Distances(1:5), FeatureDiffs{i}(:,1:5)', 'ko-', 'LineWidth', 0.25, 'MarkerSize', 4);
    for j = 1:5,
        if (~isempty(find(~isnan(FeatureDiffs{i}(:,j)))))
            MeanSEM(j,:) = [nanmean(FeatureDiffs{i}(:,j)) nanstd(FeatureDiffs{i}(:,j))/sqrt(length(find(~isnan(FeatureDiffs{i}(:,j)))))];
        else
            MeanSEM(j,:) = [NaN NaN];
        end
    end
    errorbar(Distances(1:5)+1, MeanSEM(:,1), MeanSEM(:,2), 'rs', 'LineWidth', 1.5, 'MarkerSize', 8);
    % Fit a line to the means and plot the fit line
    Coeffs = polyfit(DistancesMatrix(find(~isnan(FeatureDiffs{i}(:,1:5)))), FeatureDiffs{i}(find(~isnan(FeatureDiffs{i}(:,1:5)))), 1);
    plot(Distances(1:5), polyval(Coeffs, Distances(1:5)), 'r', 'LineWidth', 1.5);
    
%    plot(Distances(end), DistanceUnDirMeans{i}(:,end), 'ko', 'MarkerSize', 4, 'LineWidth', 0.25);
%    errorbar(Distances(end)-10, nanmean(DistanceUnDirMeans{i}(:,end)), nanstd(DistanceUnDirMeans{i}(:,end))/sqrt(length(find(~isnan(DistanceUnDirMeans{i}(:,end))))), 'ks-', 'LineWidth', 1.5, 'MarkerSize', 8);
    axis tight;
    Temp = axis;
    Temp = [-10 Distances(end-1)+5 1.15*Temp(3) 1.15*Temp(4)];
    axis(Temp);
    
    set(gca, 'XTick', [0 20 60 110 165 195], 'XTickLabel', {'0' '20' '60' '110' '165' 'UN'});
    clear MeanSEM;
    ylabel(YLabelStrings{i});
    set(gca, 'YTick', [-3:1:3]);
    
    plot(Temp(1:2), [0 0], 'k--');
    
    % Now to do the correlation and put the r and p values in the corner
    [CorrR, CorrP] = corr(Distances(reshape(repmat(1:1:5, length(BirdNames), 1), length(BirdNames)*5, 1))', reshape(FeatureDiffs{i}(:,1:5), length(BirdNames)*5, 1), 'rows', 'complete', 'type', 'Pearson');
    text(130, Temp(4)/1.08, {['r = ', num2str(CorrR)]; ['p = ', num2str(CorrP)]});
end

clear MeanSEM;
p(4,1).select();
hold on;
plot(Distances(1:5), FFCVDiffs(:,1:5)', 'ko-', 'LineWidth', 0.25, 'MarkerSize', 4);
for j = 1:5,
    if (~isempty(find(~isnan(FFCVDiffs(:,j)))))
        MeanSEM(j,:) = [nanmean(FFCVDiffs(:,j)) nanstd(FFCVDiffs(:,j))/sqrt(length(find(~isnan(FFCVDiffs(:,j)))))];
    else
        MeanSEM(j,:) = [NaN NaN];
    end
end
errorbar(Distances(1:5)+1, MeanSEM(:,1), MeanSEM(:,2), 'rs', 'LineWidth', 1.5, 'MarkerSize', 8);
% Fit a line to the means and plot the fit line
Coeffs = polyfit(DistancesMatrix(find(~isnan(FFCVDiffs(:,1:5)))), FFCVDiffs(find(~isnan(FFCVDiffs(:,1:5)))), 1);
plot(Distances(1:5), polyval(Coeffs, Distances(1:5)), 'r', 'LineWidth', 1.5);
    
% plot(Distances(end), AllBirdDistanceFFUnDirCVs(:,end), 'ko', 'MarkerSize', 4, 'LineWidth', 0.25);
% errorbar(Distances(end)-10, nanmean(AllBirdDistanceFFUnDirCVs(:,end)), nanstd(AllBirdDistanceFFUnDirCVs(:,end))/sqrt(length(find(~isnan(AllBirdDistanceFFUnDirCVs(:,end))))), 'ks-', 'LineWidth', 1.5, 'MarkerSize', 8);
axis tight;
Temp = axis;
Temp = [-10 Distances(end-1)+5 0 1.25*Temp(4)];
axis(Temp);
plot(Temp(1:2), [0 0], 'k--');
set(gca, 'XTick', [0 20 60 110 165 195], 'XTickLabel', {'0' '20' '60' '110' '165' 'UN'});
clear MeanSEM;
xlabel('Distance from female (cm)');
ylabel('CV of FF');

% Now to do the correlation and put the r and p values in the corner
[CorrR, CorrP] = corr(Distances(reshape(repmat(1:1:5, length(BirdNames), 1), length(BirdNames)*5, 1))', reshape(FFCVDiffs(:,1:5), length(BirdNames)*5, 1), 'rows', 'complete', 'type', 'Pearson');
text(130, Temp(4)/1.08, {['r = ', num2str(CorrR)]; ['p = ', num2str(CorrP)]});
 
p.fontsize = 12;
p.marginleft = 20;
p.margintop = 15;
set(gcf, 'ReSize', 'off');
set(gcf, 'PaperPositionMode', 'auto');
print(fullfile(FinalFigureDir, ['GroupDataZScores_OnlyDir_AbsDiffsFromUndir.', LabelString, '.eps']), '-depsc2', '-r600');
print(fullfile(FinalFigureDir, ['GroupDataZScores_OnlyDir_AbsDiffsFromUndir.', LabelString, '.png']), '-dpng', '-r600');    


disp('Done with stats and plots');
