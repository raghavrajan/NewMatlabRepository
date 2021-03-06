function [L0Dir_Undir_PValue] = Harini_AmplitudeZScores_Plots_OnlyDir(DirBout_Stats, DirUndirBout_Stats, UndirBout_Stats, StructString, MinTrialNo, Distances, LabelString, BirdNames)

FinalFigureDir = '/home/raghav/StudentRelated/Harini/Manuscript/results';

% Now to do the stats
% Now, we want to normalize all of the features across all the data for
% each bird. Then plot the normalized changes for directed songs for all
% features for a given bird in the same plot.
% Normalization is just a z-score across all data for that bird and that
% feature

Colours = 'rbkcmg';

MinTrialNo = 3;

for i = 1:length(DirBout_Stats)
    disp(['Bird #', num2str(i), ': ', BirdNames{i}]);
    % GroupData has 3 columns, the first one has the data corresponding to
    % the variable that that I'm looking at like TotalINNumber, etc. 2nd
    % column has the index of the distance at which the recordings were
    % done (0 for L0, 1 for L1 ... 6 for UN) and 3rd column has 1 for Dir,
    % 2 for Dir/Undir - DUN - and 3 for Undir
    if (~isfield(UndirBout_Stats{i}{end}, 'LogAmplitude_SyllLabel'))
        continue;
    end
    for SyllIndex = 1:length(UndirBout_Stats{i}{end}.LogAmplitude_SyllLabel),
        LogAmplitudeGroupData{i}{SyllIndex} = [];
        SyllDurationGroupData{i}{SyllIndex} = [];
        for j = 1:length(DirBout_Stats{i}),
            if (~isempty(DirBout_Stats{i}{j}))
                LogAmplitudeGroupData{i}{SyllIndex} = [LogAmplitudeGroupData{i}{SyllIndex}; [DirBout_Stats{i}{j}.LogAmplitude{SyllIndex}(:) ones(size(DirBout_Stats{i}{j}.LogAmplitude{SyllIndex}(:)))*j ones(size(DirBout_Stats{i}{j}.LogAmplitude{SyllIndex}(:)))*1]];
                SyllDurationGroupData{i}{SyllIndex} = [SyllDurationGroupData{i}{SyllIndex}; [DirBout_Stats{i}{j}.SyllDuration{SyllIndex}(:) ones(size(DirBout_Stats{i}{j}.SyllDuration{SyllIndex}(:)))*j ones(size(DirBout_Stats{i}{j}.SyllDuration{SyllIndex}(:)))*1]];
            end
        end
    
        for j = 1:length(DirUndirBout_Stats{i}),
            if (~isempty(DirUndirBout_Stats{i}{j}))
                LogAmplitudeGroupData{i}{SyllIndex} = [LogAmplitudeGroupData{i}{SyllIndex}; [DirUndirBout_Stats{i}{j}.LogAmplitude{SyllIndex}(:) ones(size(DirUndirBout_Stats{i}{j}.LogAmplitude{SyllIndex}(:)))*j ones(size(DirUndirBout_Stats{i}{j}.LogAmplitude{SyllIndex}(:)))*2]];
                SyllDurationGroupData{i}{SyllIndex} = [SyllDurationGroupData{i}{SyllIndex}; [DirUndirBout_Stats{i}{j}.SyllDuration{SyllIndex}(:) ones(size(DirUndirBout_Stats{i}{j}.SyllDuration{SyllIndex}(:)))*j ones(size(DirUndirBout_Stats{i}{j}.SyllDuration{SyllIndex}(:)))*1]];
            end
        end
        
        for j = 1:length(UndirBout_Stats{i}),
            if (~isempty(UndirBout_Stats{i}{j}))
                LogAmplitudeGroupData{i}{SyllIndex} = [LogAmplitudeGroupData{i}{SyllIndex}; [UndirBout_Stats{i}{j}.LogAmplitude{SyllIndex}(:) ones(size(UndirBout_Stats{i}{j}.LogAmplitude{SyllIndex}(:)))*j ones(size(UndirBout_Stats{i}{j}.LogAmplitude{SyllIndex}(:)))*3]];
                SyllDurationGroupData{i}{SyllIndex} = [SyllDurationGroupData{i}{SyllIndex}; [UndirBout_Stats{i}{j}.SyllDuration{SyllIndex}(:) ones(size(UndirBout_Stats{i}{j}.SyllDuration{SyllIndex}(:)))*j ones(size(UndirBout_Stats{i}{j}.SyllDuration{SyllIndex}(:)))*3]];
            end
        end

        % Remove NaN values
        NaNIndices = find(isnan(LogAmplitudeGroupData{i}{SyllIndex}(:,1)));
        if (~isempty(NaNIndices))
            LogAmplitudeGroupData{i}{SyllIndex}(NaNIndices,:) = [];
            SyllDurationGroupData{i}{SyllIndex}(NaNIndices,:) = [];
        end
        
        % Remove outliers
        % First find outliers based on duration and then based on amplitude
        LowerOutlierThreshold = prctile(SyllDurationGroupData{i}{SyllIndex}(:,1), 25) - 3*iqr(SyllDurationGroupData{i}{SyllIndex}(:,1));
        UpperOutlierThreshold = prctile(SyllDurationGroupData{i}{SyllIndex}(:,1), 75) + 3*iqr(SyllDurationGroupData{i}{SyllIndex}(:,1));
        OutlierIndices = find((SyllDurationGroupData{i}{SyllIndex}(:,1) < LowerOutlierThreshold) | (SyllDurationGroupData{i}{SyllIndex}(:,1) > UpperOutlierThreshold));
        disp(['Removed ', num2str(length(OutlierIndices)), ' for ', BirdNames{i}, ': Syll #', num2str(SyllIndex), ' based on syllable duration']);
        LogAmplitudeGroupData{i}{SyllIndex}(OutlierIndices,:) = [];
        
        LowerOutlierThreshold = prctile(LogAmplitudeGroupData{i}{SyllIndex}(:,1), 25) - 3*iqr(LogAmplitudeGroupData{i}{SyllIndex}(:,1));
        UpperOutlierThreshold = prctile(LogAmplitudeGroupData{i}{SyllIndex}(:,1), 75) + 3*iqr(LogAmplitudeGroupData{i}{SyllIndex}(:,1));
        OutlierIndices = find((LogAmplitudeGroupData{i}{SyllIndex}(:,1) < LowerOutlierThreshold) | (LogAmplitudeGroupData{i}{SyllIndex}(:,1) > UpperOutlierThreshold));
        disp(['Removed ', num2str(length(OutlierIndices)), ' for ', BirdNames{i}, ': Syll #', num2str(SyllIndex), ' based on syllable amplitude']);
        LogAmplitudeGroupData{i}{SyllIndex}(OutlierIndices,:) = [];
        
        LogAmplitudeGroupData{i}{SyllIndex}(:,1) = zscore(LogAmplitudeGroupData{i}{SyllIndex}(:,1));
    
    end
end

DistanceLogAmplitudeDirMeans = [];
DistanceLogAmplitudeDirUnDirMeans = [];
DistanceLogAmplitudeUnDirMeans = [];

Index = 0;

for i = 1:length(BirdNames),
    Index = i;
    for SyllIndex = 1:length(LogAmplitudeGroupData{i}),
        % First for directed songs
        for j = 1:6,
            DirSongIndices = find((LogAmplitudeGroupData{i}{SyllIndex}(:,3) == 1) & (LogAmplitudeGroupData{i}{SyllIndex}(:,2) == j));
            if (length(DirSongIndices) >= MinTrialNo)
                DistanceLogAmplitudeDirMeans{i}(SyllIndex, j) = mean(LogAmplitudeGroupData{i}{SyllIndex}(DirSongIndices,1));
                DistanceLogAmplitudeDirSEMs{i}(SyllIndex, j) = std(LogAmplitudeGroupData{i}{SyllIndex}(DirSongIndices,1))/sqrt(length(DirSongIndices));
            else
                DistanceLogAmplitudeDirMeans{i}(SyllIndex, j) = NaN;
                DistanceLogAmplitudeDirSEMs{i}(SyllIndex,j) = NaN;
            end

            DirUnDirSongIndices = find((LogAmplitudeGroupData{i}{SyllIndex}(:,3) == 2) & (LogAmplitudeGroupData{i}{SyllIndex}(:,2) == j));
            if (length(DirUnDirSongIndices) >= MinTrialNo)
                DistanceLogAmplitudeDirUnDirMeans{i}(SyllIndex, j) = mean(LogAmplitudeGroupData{i}{SyllIndex}(DirUnDirSongIndices,1));
                DistanceLogAmplitudeDirUnDirSEMs{i}(SyllIndex, j) = std(LogAmplitudeGroupData{i}{SyllIndex}(DirUnDirSongIndices,1))/sqrt(length(DirUnDirSongIndices));
            else
                DistanceLogAmplitudeDirUnDirMeans{i}(SyllIndex, j) = NaN;
                DistanceLogAmplitudeDirUnDirSEMs{i}(SyllIndex,j) = NaN;
            end
            
            UnDirSongIndices = find((LogAmplitudeGroupData{i}{SyllIndex}(:,3) == 3) & (LogAmplitudeGroupData{i}{SyllIndex}(:,2) == j));
            if (length(UnDirSongIndices) >= MinTrialNo)
                DistanceLogAmplitudeUnDirMeans{i}(SyllIndex, j) = mean(LogAmplitudeGroupData{i}{SyllIndex}(UnDirSongIndices,1));
                DistanceLogAmplitudeUnDirSEMs{i}(SyllIndex, j) = std(LogAmplitudeGroupData{i}{SyllIndex}(UnDirSongIndices,1))/sqrt(length(UnDirSongIndices));
            else
                DistanceLogAmplitudeUnDirMeans{i}(SyllIndex, j) = NaN;
                DistanceLogAmplitudeUnDirSEMs{i}(SyllIndex,j) = NaN;
            end
        end
    end
    if (~isempty(LogAmplitudeGroupData{i}))
        if (size(DistanceLogAmplitudeDirMeans{i}, 1) > 1)
            AllBirdDistanceLogAmplitudeDirMeans(i,:) = mean(DistanceLogAmplitudeDirMeans{i});
            AllBirdDistanceLogAmplitudeUnDirMeans(i,:) = mean(DistanceLogAmplitudeUnDirMeans{i});
        else
            AllBirdDistanceLogAmplitudeDirMeans(i,:) = DistanceLogAmplitudeDirMeans{i};
            AllBirdDistanceLogAmplitudeUnDirMeans(i,:) = DistanceLogAmplitudeUnDirMeans{i};
        end
    else
        AllBirdDistanceLogAmplitudeDirMeans(i,:) = ones(1,6)*NaN;
        AllBirdDistanceLogAmplitudeUnDirMeans(i,:) = ones(1,6)*NaN;
    end
end

YLabelStrings = {'# of INs (zscore)', '# of motifs/bout (zscore)', 'First motif duration (zscore)'};

DistancesMatrix = repmat(Distances(1:5), size(AllBirdDistanceLogAmplitudeDirMeans,1), 1);

figure;
set(gcf, 'Color', 'w');
set(gcf, 'Position', [680 162 440 250]);
p = panel();
p.pack(1,1);
p.de.margin = 15;

clear MeanSEM;
p(1,1).select();
hold on;
plot(Distances(1:5), AllBirdDistanceLogAmplitudeDirMeans(:,1:5)', 'ko-', 'LineWidth', 0.25, 'MarkerSize', 4);
for j = 1:5,
    if (~isempty(find(~isnan(AllBirdDistanceLogAmplitudeDirMeans(:,j)))))
        MeanSEM(j,:) = [nanmean(AllBirdDistanceLogAmplitudeDirMeans(:,j)) nanstd(AllBirdDistanceLogAmplitudeDirMeans(:,j))/sqrt(length(find(~isnan(AllBirdDistanceLogAmplitudeDirMeans(:,j)))))];
    else
        MeanSEM(j,:) = [NaN NaN];
    end
end
errorbar(Distances(1:5)+7, MeanSEM(:,1), MeanSEM(:,2), 'rs', 'LineWidth', 1.5, 'MarkerSize', 8);
% Fit a line to the means and plot the fit line
Coeffs = polyfit(DistancesMatrix(find(~isnan(AllBirdDistanceLogAmplitudeDirMeans(:,1:5)))), AllBirdDistanceLogAmplitudeDirMeans(find(~isnan(AllBirdDistanceLogAmplitudeDirMeans(:,1:5)))), 1);
plot(Distances(1:5), polyval(Coeffs, Distances(1:5)), 'r', 'LineWidth', 1.5);
TempCoeffs = robustfit(DistancesMatrix(find(~isnan(AllBirdDistanceLogAmplitudeDirMeans(:,1:5)))), AllBirdDistanceLogAmplitudeDirMeans(find(~isnan(AllBirdDistanceLogAmplitudeDirMeans(:,1:5)))));
plot(Distances(1:5), polyval(flipud(TempCoeffs), Distances(1:5)), 'b', 'LineWidth', 1.5);

plot(Distances(end), AllBirdDistanceLogAmplitudeUnDirMeans(:,end), 'ko', 'MarkerSize', 4, 'LineWidth', 0.25);
errorbar(Distances(end)-10, nanmean(AllBirdDistanceLogAmplitudeUnDirMeans(:,end)), nanstd(AllBirdDistanceLogAmplitudeUnDirMeans(:,end))/sqrt(length(find(~isnan(AllBirdDistanceLogAmplitudeUnDirMeans(:,end))))), 'ks-', 'LineWidth', 1.5, 'MarkerSize', 8);
axis tight;
Temp = axis;
Temp = [-10 Distances(end)+5 0 1.25*Temp(4)];
axis(Temp);
plot(Temp(1:2), ones(1,2)*nanmean(AllBirdDistanceLogAmplitudeUnDirMeans(:,end)), 'k--');
set(gca, 'XTick', [0 20 60 110 165 195], 'XTickLabel', {'0' '20' '60' '110' '165' 'UN'});
clear MeanSEM;
xlabel('Distance from female (cm)');
ylabel('Mean of LogAmplitude');

% Now to do the correlation and put the r and p values in the corner
[CorrR, CorrP] = corr(Distances(reshape(repmat(1:1:5, length(BirdNames), 1), length(BirdNames)*5, 1))', reshape(AllBirdDistanceLogAmplitudeDirMeans(:,1:5), length(BirdNames)*5, 1), 'rows', 'complete', 'type', 'Pearson');
text(130, Temp(4)/1.08, {['r = ', num2str(CorrR)]; ['p = ', num2str(CorrP)]});

p.fontsize = 12;
p.marginleft = 25;
p.margintop = 15;
set(gcf, 'ReSize', 'off');
set(gcf, 'PaperPositionMode', 'auto');
print(fullfile(FinalFigureDir, ['LogAmplitudeZScores_OnlyDir.', LabelString, '.eps']), '-depsc2', '-r600');
print(fullfile(FinalFigureDir, ['LogAmplitudeZScores_OnlyDir.', LabelString, '.png']), '-dpng', '-r600');    

% Now plot a version where I plot the absolute differences between the
% values and undirected songs in the absence of the female. This should
% help to show that it does become more like undirected song

LogAmplitudeMeanDiffs = ((AllBirdDistanceLogAmplitudeDirMeans - repmat(AllBirdDistanceLogAmplitudeUnDirMeans(:,end), 1, size(AllBirdDistanceLogAmplitudeDirMeans, 2)))./repmat(AllBirdDistanceLogAmplitudeUnDirMeans(:,end), 1, size(AllBirdDistanceLogAmplitudeDirMeans, 2))).^2;

figure;
set(gcf, 'Color', 'w');
set(gcf, 'Position', [680 162 440 250]);
p = panel();
p.pack(1,1);
p.de.margin = 15;

clear MeanSEM;
p(1,1).select();
hold on;
plot(Distances(1:5), LogAmplitudeMeanDiffs(:,1:5)', 'ko-', 'LineWidth', 0.25, 'MarkerSize', 4);
for j = 1:5,
    if (~isempty(find(~isnan(LogAmplitudeMeanDiffs(:,j)))))
        MeanSEM(j,:) = [nanmean(LogAmplitudeMeanDiffs(:,j)) nanstd(LogAmplitudeMeanDiffs(:,j))/sqrt(length(find(~isnan(LogAmplitudeMeanDiffs(:,j)))))];
    else
        MeanSEM(j,:) = [NaN NaN];
    end
end
errorbar(Distances(1:5)+1, MeanSEM(:,1), MeanSEM(:,2), 'rs', 'LineWidth', 1.5, 'MarkerSize', 8);
% Fit a line to the means and plot the fit line
Coeffs = polyfit(DistancesMatrix(find(~isnan(LogAmplitudeMeanDiffs(:,1:5)))), LogAmplitudeMeanDiffs(find(~isnan(LogAmplitudeMeanDiffs(:,1:5)))), 1);
plot(Distances(1:5), polyval(Coeffs, Distances(1:5)), 'r', 'LineWidth', 1.5);
TempCoeffs = robustfit(DistancesMatrix(find(~isnan(LogAmplitudeMeanDiffs(:,1:5)))), LogAmplitudeMeanDiffs(find(~isnan(LogAmplitudeMeanDiffs(:,1:5)))));
plot(Distances(1:5), polyval(flipud(TempCoeffs), Distances(1:5)), 'b', 'LineWidth', 1.5);

% plot(Distances(end), AllBirdDistanceLogAmplitudeUnDirMeans(:,end), 'ko', 'MarkerSize', 4, 'LineWidth', 0.25);
% errorbar(Distances(end)-10, nanmean(AllBirdDistanceLogAmplitudeUnDirMeans(:,end)), nanstd(AllBirdDistanceLogAmplitudeUnDirMeans(:,end))/sqrt(length(find(~isnan(AllBirdDistanceLogAmplitudeUnDirMeans(:,end))))), 'ks-', 'LineWidth', 1.5, 'MarkerSize', 8);
axis tight;
Temp = axis;
Temp = [-10 Distances(end-1)+5 0 1.25*Temp(4)];
axis(Temp);
plot(Temp(1:2), [0 0], 'k--');
set(gca, 'XTick', [0 20 60 110 165 195], 'XTickLabel', {'0' '20' '60' '110' '165' 'UN'});
clear MeanSEM;
xlabel('Distance from female (cm)');
ylabel('Mean of LogAmplitude');

% Now to do the correlation and put the r and p values in the corner
[CorrR, CorrP] = corr(Distances(reshape(repmat(1:1:5, length(BirdNames), 1), length(BirdNames)*5, 1))', reshape(LogAmplitudeMeanDiffs(:,1:5), length(BirdNames)*5, 1), 'rows', 'complete', 'type', 'Pearson');
text(130, Temp(4)/1.08, {['r = ', num2str(CorrR)]; ['p = ', num2str(CorrP)]});
 
p.fontsize = 12;
p.marginleft = 20;
p.margintop = 15;
set(gcf, 'ReSize', 'off');
set(gcf, 'PaperPositionMode', 'auto');
print(fullfile(FinalFigureDir, ['LogAmpitudeZScores_OnlyDir_AbsDiffsFromUndir.', LabelString, '.eps']), '-depsc2', '-r600');
print(fullfile(FinalFigureDir, ['LogAmpitudeZScores_OnlyDir_AbsDiffsFromUndir.', LabelString, '.png']), '-dpng', '-r600');    

disp('Done with stats and plots');
