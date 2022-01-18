function [L0Dir_Undir_PValue] = Harini_AllFeatureStats_Plots_LogAmpI(DirBout_Stats, DirUndirBout_Stats, UndirBout_Stats, StructString, MinTrialNo, Distances, LabelString, BirdNames)

FinalFigureDir = '/home/raghav/StudentRelated/Harini/Manuscript/results';

Colours = 'rbkcmg';

MinTrialNo = 10;

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
        for j = 1:length(DirBout_Stats{i}),
            if (~isempty(DirBout_Stats{i}{j}))
                LogAmplitudeGroupData{i}{SyllIndex} = [LogAmplitudeGroupData{i}{SyllIndex}; [DirBout_Stats{i}{j}.LogAmplitude{SyllIndex}(:) ones(size(DirBout_Stats{i}{j}.LogAmplitude{SyllIndex}(:)))*j ones(size(DirBout_Stats{i}{j}.LogAmplitude{SyllIndex}(:)))*1]];
            end
        end
    
        for j = 1:length(DirUndirBout_Stats{i}),
            if (~isempty(DirUndirBout_Stats{i}{j}))
                LogAmplitudeGroupData{i}{SyllIndex} = [LogAmplitudeGroupData{i}{SyllIndex}; [DirUndirBout_Stats{i}{j}.LogAmplitude{SyllIndex}(:) ones(size(DirUndirBout_Stats{i}{j}.LogAmplitude{SyllIndex}(:)))*j ones(size(DirUndirBout_Stats{i}{j}.LogAmplitude{SyllIndex}(:)))*2]];
            end
        end
        
        for j = 1:length(UndirBout_Stats{i}),
            if (~isempty(UndirBout_Stats{i}{j}))
                LogAmplitudeGroupData{i}{SyllIndex} = [LogAmplitudeGroupData{i}{SyllIndex}; [UndirBout_Stats{i}{j}.LogAmplitude{SyllIndex}(:) ones(size(UndirBout_Stats{i}{j}.LogAmplitude{SyllIndex}(:)))*j ones(size(UndirBout_Stats{i}{j}.LogAmplitude{SyllIndex}(:)))*3]];
            end
        end

        % Remove NaN values
        NaNIndices = find(isnan(LogAmplitudeGroupData{i}{SyllIndex}(:,1)));
        if (~isempty(NaNIndices))
            LogAmplitudeGroupData{i}{SyllIndex}(NaNIndices,:) = [];
        end
    
        LogAmplitudeGroupData{i}{SyllIndex}(:,1) = zscore(LogAmplitudeGroupData{i}{SyllIndex}(:,1));
    
    end
end

DistanceLogAmplitudeDirMeans = [];
DistanceLogAmplitudeDirUnDirMeans = [];
DistanceLogAmplitudeUnDirMeans = [];
DistanceLogAmplitudeDirCVs = [];
DistanceLogAmplitudeDirUnDirCVs = [];
DistanceLogAmplitudeUnDirCVs = [];

Index = 0;

AllBirdIndividualSyllableDirMeans = [];
AllBirdIndividualSyllableUnDirMeans = [];

for i = 1:length(BirdNames),
    Index = i;
    for SyllIndex = 1:length(LogAmplitudeGroupData{i}),
        % First for directed songs
        for j = 1:6,
            DirSongIndices = find((LogAmplitudeGroupData{i}{SyllIndex}(:,3) == 1) & (LogAmplitudeGroupData{i}{SyllIndex}(:,2) == j));
            if (length(DirSongIndices) >= MinTrialNo)
                DistanceLogAmplitudeDirMeans{i}(SyllIndex, j) = mean(LogAmplitudeGroupData{i}{SyllIndex}(DirSongIndices,1));
                DistanceLogAmplitudeDirCVs{i}(SyllIndex, j) = std(LogAmplitudeGroupData{i}{SyllIndex}(DirSongIndices,1))/mean(LogAmplitudeGroupData{i}{SyllIndex}(DirSongIndices,1));
                DistanceLogAmplitudeDirSEMs{i}(SyllIndex, j) = std(LogAmplitudeGroupData{i}{SyllIndex}(DirSongIndices,1))/sqrt(length(DirSongIndices));
            else
                DistanceLogAmplitudeDirMeans{i}(SyllIndex, j) = NaN;
                DistanceLogAmplitudeDirCVs{i}(SyllIndex, j) = NaN;
                DistanceLogAmplitudeDirSEMs{i}(SyllIndex,j) = NaN;
            end

            DirUnDirSongIndices = find((LogAmplitudeGroupData{i}{SyllIndex}(:,3) == 2) & (LogAmplitudeGroupData{i}{SyllIndex}(:,2) == j));
            if (length(DirUnDirSongIndices) >= MinTrialNo)
                DistanceLogAmplitudeDirUnDirMeans{i}(SyllIndex, j) = mean(LogAmplitudeGroupData{i}{SyllIndex}(DirUnDirSongIndices,1));
                DistanceLogAmplitudeDirUnDirCVs{i}(SyllIndex, j) = std(LogAmplitudeGroupData{i}{SyllIndex}(DirUnDirSongIndices,1))/mean(LogAmplitudeGroupData{i}{SyllIndex}(DirUnDirSongIndices,1));
                DistanceLogAmplitudeDirUnDirSEMs{i}(SyllIndex, j) = std(LogAmplitudeGroupData{i}{SyllIndex}(DirUnDirSongIndices,1))/sqrt(length(DirUnDirSongIndices));
            else
                DistanceLogAmplitudeDirUnDirMeans{i}(SyllIndex, j) = NaN;
                DistanceLogAmplitudeDirUnDirCVs{i}(SyllIndex, j) = NaN;
                DistanceLogAmplitudeDirUnDirSEMs{i}(SyllIndex,j) = NaN;
            end
            
            UnDirSongIndices = find((LogAmplitudeGroupData{i}{SyllIndex}(:,3) == 3) & (LogAmplitudeGroupData{i}{SyllIndex}(:,2) == j));
            if (length(UnDirSongIndices) >= MinTrialNo)
                DistanceLogAmplitudeUnDirMeans{i}(SyllIndex, j) = mean(LogAmplitudeGroupData{i}{SyllIndex}(UnDirSongIndices,1));
                DistanceLogAmplitudeUnDirCVs{i}(SyllIndex, j) = std(LogAmplitudeGroupData{i}{SyllIndex}(UnDirSongIndices,1))/mean(LogAmplitudeGroupData{i}{SyllIndex}(UnDirSongIndices,1));
                DistanceLogAmplitudeUnDirSEMs{i}(SyllIndex, j) = std(LogAmplitudeGroupData{i}{SyllIndex}(UnDirSongIndices,1))/sqrt(length(UnDirSongIndices));
            else
                DistanceLogAmplitudeUnDirMeans{i}(SyllIndex, j) = NaN;
                DistanceLogAmplitudeUnDirCVs{i}(SyllIndex, j) = NaN;
                DistanceLogAmplitudeUnDirSEMs{i}(SyllIndex,j) = NaN;
            end
        end
    end
    if (~isempty(LogAmplitudeGroupData{i}))
        AllBirdIndividualSyllableDirMeans = [AllBirdIndividualSyllableDirMeans; DistanceLogAmplitudeDirMeans{i}];
        AllBirdIndividualSyllableUnDirMeans = [AllBirdIndividualSyllableUnDirMeans; DistanceLogAmplitudeUnDirMeans{i}];
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

figure;
set(gcf, 'Color', 'w');
set(gcf, 'Position', [680 162 750 350]);
p = panel();
p.pack(1,2);
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
errorbar(Distances(1:5)+1, MeanSEM(:,1), MeanSEM(:,2), 'rs-', 'LineWidth', 1.5, 'MarkerSize', 8);
plot(Distances(end), AllBirdDistanceLogAmplitudeUnDirMeans(:,end), 'ko', 'MarkerSize', 4, 'LineWidth', 0.25);
errorbar(Distances(end)-10, nanmean(AllBirdDistanceLogAmplitudeUnDirMeans(:,end)), nanstd(AllBirdDistanceLogAmplitudeUnDirMeans(:,end))/sqrt(length(find(~isnan(AllBirdDistanceLogAmplitudeUnDirMeans(:,end))))), 'ks-', 'LineWidth', 1.5, 'MarkerSize', 8);
axis tight;
Temp1 = axis;

set(gca, 'XTick', [0 20 60 110 165 195], 'XTickLabel', {'0' '20' '60' '110' '165' 'UN'});
clear MeanSEM;
xlabel('Distance from female (cm)');
ylabel('Zscore of LogAmplitude');
% set(gca, 'YTick', [-3:1:3]);

p(1,2).select();
hold on;
plot(Distances(1:5), AllBirdDistanceLogAmplitudeUnDirMeans(:,1:5)', 'ko-', 'LineWidth', 0.25, 'MarkerSize', 4);
for j = 1:5,
    if (~isempty(find(~isnan(AllBirdDistanceLogAmplitudeUnDirMeans(:,j)))))
        MeanSEM(j,:) = [nanmean(AllBirdDistanceLogAmplitudeUnDirMeans(:,j)) nanstd(AllBirdDistanceLogAmplitudeUnDirMeans(:,j))/sqrt(length(find(~isnan(AllBirdDistanceLogAmplitudeUnDirMeans(:,j)))))];
    else
        MeanSEM(j,:) = [NaN NaN];
    end
end
errorbar(Distances(1:5)+1, MeanSEM(:,1), MeanSEM(:,2), 'bs-', 'LineWidth', 1.5, 'MarkerSize', 8);
plot(Distances(end), AllBirdDistanceLogAmplitudeUnDirMeans(:,end), 'ko', 'MarkerSize', 4, 'LineWidth', 0.25);
errorbar(Distances(end)-10, nanmean(AllBirdDistanceLogAmplitudeUnDirMeans(:,end)), nanstd(AllBirdDistanceLogAmplitudeUnDirMeans(:,end))/sqrt(length(find(~isnan(AllBirdDistanceLogAmplitudeUnDirMeans(:,end))))), 'ks-', 'LineWidth', 1.5, 'MarkerSize', 8);
axis tight;
Temp2 = axis;
Temp = [-10 Distances(end)+5 1.15*min(Temp1(3), Temp2(3)) 1.25*max(Temp1(4), Temp2(4))];
axis(Temp);
plot(Temp(1:2), [0 0], 'k--');
set(gca, 'XTick', [0 20 60 110 165 195], 'XTickLabel', {'0' '20' '60' '110' '165' 'UN'});
clear MeanSEM;
xlabel('Distance from female (cm)');
set(gca, 'YTick', [-3:1:3]);

p(1,1).select();
axis(Temp);
plot(Temp(1:2), [0 0], 'k--');

% Now to do the correlation and put the r and p values in the corner
[CorrR, CorrP] = corr(Distances(reshape(repmat(1:1:5, 7, 1), 7*5, 1))', reshape(AllBirdDistanceLogAmplitudeDirMeans(:,1:5), 7*5, 1), 'rows', 'complete', 'type', 'Pearson');
text(130, Temp(4)/1.08, {['r = ', num2str(CorrR)]; ['p = ', num2str(CorrP)]});

p(1,2).select();
% Now to do the correlation and put the r and p values in the corner
[CorrR, CorrP] = corr(Distances(reshape(repmat(1:1:5, 7, 1), 7*5, 1))', reshape(AllBirdDistanceLogAmplitudeUnDirMeans(:,1:5), 7*5, 1), 'rows', 'complete', 'type', 'Pearson');
text(130, Temp(4)/1.08, {['r = ', num2str(CorrR)]; ['p = ', num2str(CorrP)]});

p.fontsize = 12;
p.marginleft = 20;
p.margintop = 15;
% set(gcf, 'ReSize', 'off');
set(gcf, 'PaperPositionMode', 'auto');
print(fullfile(FinalFigureDir, ['Fig.6.AmplitudeGroupDataZScores.', LabelString, '.eps']), '-depsc2', '-r600');
print(fullfile(FinalFigureDir, ['Fig.6.AmplitudeGroupDataZScores.', LabelString, '.png']), '-dpng', '-r600');    


figure;
set(gcf, 'Color', 'w');
set(gcf, 'Position', [680 162 750 350]);
p = panel();
p.pack(1,2);
p.de.margin = 15;

% Flip data depending on whether dir is greater or undir is greater
for i = 1:size(AllBirdDistanceLogAmplitudeDirMeans),
    if (AllBirdDistanceLogAmplitudeDirMeans(i,1) < AllBirdDistanceLogAmplitudeUnDirMeans(i,end))
        AllBirdDistanceLogAmplitudeDirMeans(i,:) = -AllBirdDistanceLogAmplitudeDirMeans(i,:);
        AllBirdDistanceLogAmplitudeUnDirMeans(i,:) = -AllBirdDistanceLogAmplitudeUnDirMeans(i,:);
    end
end

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
errorbar(Distances(1:5)+1, MeanSEM(:,1), MeanSEM(:,2), 'rs-', 'LineWidth', 1.5, 'MarkerSize', 8);
plot(Distances(end), AllBirdDistanceLogAmplitudeUnDirMeans(:,end), 'ko', 'MarkerSize', 4, 'LineWidth', 0.25);
errorbar(Distances(end)-10, nanmean(AllBirdDistanceLogAmplitudeUnDirMeans(:,end)), nanstd(AllBirdDistanceLogAmplitudeUnDirMeans(:,end))/sqrt(length(find(~isnan(AllBirdDistanceLogAmplitudeUnDirMeans(:,end))))), 'ks-', 'LineWidth', 1.5, 'MarkerSize', 8);
axis tight;
Temp1 = axis;

set(gca, 'XTick', [0 20 60 110 165 195], 'XTickLabel', {'0' '20' '60' '110' '165' 'UN'});
clear MeanSEM;
xlabel('Distance from female (cm)');
ylabel('Zscore of LogAmplitude');
% set(gca, 'YTick', [-3:1:3]);

p(1,2).select();
hold on;
plot(Distances(1:5), AllBirdDistanceLogAmplitudeUnDirMeans(:,1:5)', 'ko-', 'LineWidth', 0.25, 'MarkerSize', 4);
for j = 1:5,
    if (~isempty(find(~isnan(AllBirdDistanceLogAmplitudeUnDirMeans(:,j)))))
        MeanSEM(j,:) = [nanmean(AllBirdDistanceLogAmplitudeUnDirMeans(:,j)) nanstd(AllBirdDistanceLogAmplitudeUnDirMeans(:,j))/sqrt(length(find(~isnan(AllBirdDistanceLogAmplitudeUnDirMeans(:,j)))))];
    else
        MeanSEM(j,:) = [NaN NaN];
    end
end
errorbar(Distances(1:5)+1, MeanSEM(:,1), MeanSEM(:,2), 'bs-', 'LineWidth', 1.5, 'MarkerSize', 8);
plot(Distances(end), AllBirdDistanceLogAmplitudeUnDirMeans(:,end), 'ko', 'MarkerSize', 4, 'LineWidth', 0.25);
errorbar(Distances(end)-10, nanmean(AllBirdDistanceLogAmplitudeUnDirMeans(:,end)), nanstd(AllBirdDistanceLogAmplitudeUnDirMeans(:,end))/sqrt(length(find(~isnan(AllBirdDistanceLogAmplitudeUnDirMeans(:,end))))), 'ks-', 'LineWidth', 1.5, 'MarkerSize', 8);
axis tight;
Temp2 = axis;
Temp = [-10 Distances(end)+5 1.15*min(Temp1(3), Temp2(3)) 1.25*max(Temp1(4), Temp2(4))];
axis(Temp);
plot(Temp(1:2), [0 0], 'k--');
set(gca, 'XTick', [0 20 60 110 165 195], 'XTickLabel', {'0' '20' '60' '110' '165' 'UN'});
clear MeanSEM;
xlabel('Distance from female (cm)');
set(gca, 'YTick', [-3:1:3]);

p(1,1).select();
axis(Temp);
plot(Temp(1:2), [0 0], 'k--');

% Now to do the correlation and put the r and p values in the corner
[CorrR, CorrP] = corr(Distances(reshape(repmat(1:1:5, 7, 1), 7*5, 1))', reshape(AllBirdDistanceLogAmplitudeDirMeans(:,1:5), 7*5, 1), 'rows', 'complete', 'type', 'Pearson');
text(130, Temp(4)/1.08, {['r = ', num2str(CorrR)]; ['p = ', num2str(CorrP)]});

p(1,2).select();
% Now to do the correlation and put the r and p values in the corner
[CorrR, CorrP] = corr(Distances(reshape(repmat(1:1:5, 7, 1), 7*5, 1))', reshape(AllBirdDistanceLogAmplitudeUnDirMeans(:,1:5), 7*5, 1), 'rows', 'complete', 'type', 'Pearson');
text(130, Temp(4)/1.08, {['r = ', num2str(CorrR)]; ['p = ', num2str(CorrP)]});

p.fontsize = 12;
p.marginleft = 20;
p.margintop = 15;
% set(gcf, 'ReSize', 'off');
set(gcf, 'PaperPositionMode', 'auto');
print(fullfile(FinalFigureDir, ['Fig.6.AmplitudeGroupDataZScores.Flipped.', LabelString, '.eps']), '-depsc2', '-r600');
print(fullfile(FinalFigureDir, ['Fig.6.AmplitudeGroupDataZScores.Flipped.', LabelString, '.png']), '-dpng', '-r600');    


figure;
set(gcf, 'Color', 'w');
set(gcf, 'Position', [680 162 750 350]);
p = panel();
p.pack(1,2);
p.de.margin = 15;

% Flip data depending on whether dir is greater or undir is greater
for i = 1:size(AllBirdIndividualSyllableDirMeans),
    if (AllBirdIndividualSyllableDirMeans(i,1) < AllBirdIndividualSyllableUnDirMeans(i,end))
        AllBirdIndividualSyllableDirMeans(i,:) = -AllBirdIndividualSyllableDirMeans(i,:);
        AllBirdIndividualSyllableUnDirMeans(i,:) = -AllBirdIndividualSyllableUnDirMeans(i,:);
    end
end

clear MeanSEM;
p(1,1).select();
hold on;
plot(Distances(1:5), AllBirdIndividualSyllableDirMeans(:,1:5)', 'ko-', 'LineWidth', 0.25, 'MarkerSize', 4);
for j = 1:5,
    if (~isempty(find(~isnan(AllBirdIndividualSyllableDirMeans(:,j)))))
        MeanSEM(j,:) = [nanmean(AllBirdIndividualSyllableDirMeans(:,j)) nanstd(AllBirdIndividualSyllableDirMeans(:,j))/sqrt(length(find(~isnan(AllBirdIndividualSyllableDirMeans(:,j)))))];
    else
        MeanSEM(j,:) = [NaN NaN];
    end
end
errorbar(Distances(1:5)+1, MeanSEM(:,1), MeanSEM(:,2), 'rs-', 'LineWidth', 1.5, 'MarkerSize', 8);
plot(Distances(end), AllBirdIndividualSyllableUnDirMeans(:,end), 'ko', 'MarkerSize', 4, 'LineWidth', 0.25);
errorbar(Distances(end)-10, nanmean(AllBirdIndividualSyllableUnDirMeans(:,end)), nanstd(AllBirdIndividualSyllableUnDirMeans(:,end))/sqrt(length(find(~isnan(AllBirdIndividualSyllableUnDirMeans(:,end))))), 'ks-', 'LineWidth', 1.5, 'MarkerSize', 8);
axis tight;
Temp1 = axis;

set(gca, 'XTick', [0 20 60 110 165 195], 'XTickLabel', {'0' '20' '60' '110' '165' 'UN'});
clear MeanSEM;
xlabel('Distance from female (cm)');
ylabel('Zscore of LogAmplitude');
% set(gca, 'YTick', [-3:1:3]);

p(1,2).select();
hold on;
plot(Distances(1:5), AllBirdIndividualSyllableUnDirMeans(:,1:5)', 'ko-', 'LineWidth', 0.25, 'MarkerSize', 4);
for j = 1:5,
    if (~isempty(find(~isnan(AllBirdIndividualSyllableUnDirMeans(:,j)))))
        MeanSEM(j,:) = [nanmean(AllBirdIndividualSyllableUnDirMeans(:,j)) nanstd(AllBirdIndividualSyllableUnDirMeans(:,j))/sqrt(length(find(~isnan(AllBirdIndividualSyllableUnDirMeans(:,j)))))];
    else
        MeanSEM(j,:) = [NaN NaN];
    end
end
errorbar(Distances(1:5)+1, MeanSEM(:,1), MeanSEM(:,2), 'bs-', 'LineWidth', 1.5, 'MarkerSize', 8);
plot(Distances(end), AllBirdIndividualSyllableUnDirMeans(:,end), 'ko', 'MarkerSize', 4, 'LineWidth', 0.25);
errorbar(Distances(end)-10, nanmean(AllBirdIndividualSyllableUnDirMeans(:,end)), nanstd(AllBirdIndividualSyllableUnDirMeans(:,end))/sqrt(length(find(~isnan(AllBirdIndividualSyllableUnDirMeans(:,end))))), 'ks-', 'LineWidth', 1.5, 'MarkerSize', 8);
axis tight;
Temp2 = axis;
Temp = [-10 Distances(end)+5 1.15*min(Temp1(3), Temp2(3)) 1.25*max(Temp1(4), Temp2(4))];
axis(Temp);
plot(Temp(1:2), [0 0], 'k--');
set(gca, 'XTick', [0 20 60 110 165 195], 'XTickLabel', {'0' '20' '60' '110' '165' 'UN'});
clear MeanSEM;
xlabel('Distance from female (cm)');
set(gca, 'YTick', [-3:1:3]);

p(1,1).select();
axis(Temp);
plot(Temp(1:2), [0 0], 'k--');

% Now to do the correlation and put the r and p values in the corner
[CorrR, CorrP] = corr(Distances(reshape(repmat(1:1:5, 16, 1), 16*5, 1))', reshape(AllBirdIndividualSyllableDirMeans(:,1:5), 16*5, 1), 'rows', 'complete', 'type', 'Pearson');
text(130, Temp(4)/1.08, {['r = ', num2str(CorrR)]; ['p = ', num2str(CorrP)]});

p(1,2).select();
% Now to do the correlation and put the r and p values in the corner
[CorrR, CorrP] = corr(Distances(reshape(repmat(1:1:5, 16, 1), 16*5, 1))', reshape(AllBirdIndividualSyllableUnDirMeans(:,1:5), 16*5, 1), 'rows', 'complete', 'type', 'Pearson');
text(130, Temp(4)/1.08, {['r = ', num2str(CorrR)]; ['p = ', num2str(CorrP)]});

p.fontsize = 12;
p.marginleft = 20;
p.margintop = 15;
% set(gcf, 'ReSize', 'off');
set(gcf, 'PaperPositionMode', 'auto');
print(fullfile(FinalFigureDir, ['Fig.6.AmplitudeIndividualSyllDataZScores.Flipped.', LabelString, '.eps']), '-depsc2', '-r600');
print(fullfile(FinalFigureDir, ['Fig.6.AmplitudeIndividualSyllDataZScores.Flipped.', LabelString, '.png']), '-dpng', '-r600');    


disp('Done with stats and plots');
