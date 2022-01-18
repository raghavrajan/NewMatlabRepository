function [L0Dir_Undir_PValue] = Harini_INNum_MotifNumDur_L0Dir_UnDir_Diffs(DirBout_Stats, DirUndirBout_Stats, UndirBout_Stats, Distances, BirdNames)

DistanceCode = {'L0' 'L1' 'L2' 'L3' 'L4' 'UN'};
FinalFigureDir = '/home/raghav/StudentRelated/Harini/Manuscript/results';

% We want to do the analysis per bird and for the group as a whole. The
% things that we will analyze are (1) # of INs at the beginning based on
% the 500ms criterion, (2) # of complete motifs in the bout, (3) all and
% first motif duration, (4), CV of FF and (5) amplitude of individual
% motif syllables only for HF mic birds

FeaturesToConsider = {'TotalINNumber_500ms' 'CompleteMotifNumber'};
MinTrialNo = 3;
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
    
    % Now find and remove any distance that doesn't have 3 bouts of 
    % directed song
    for j = 1:max(GroupData{i}(:,end-1)),
        DirIndicesForDistance = find((GroupData{i}(:,end-1) == j) & (GroupData{i}(:,end) == 1));
        if ((length(DirIndicesForDistance) > 0) & (length(DirIndicesForDistance) < MinTrialNo))
            GroupData{i}(DirIndicesForDistance,:) = [];
            disp(['Removed distance ', DistanceCode{j}, ' for ', BirdNames{i}, ' as there were only ', num2str(length(DirIndicesForDistance)), ' bouts']);
        end
    end
    
    L0DirIndices = find((GroupData{i}(:,end) == 1) & (GroupData{i}(:,end-1) == 1));
    UnDirIndices = find((GroupData{i}(:,end-1) == 6) & (GroupData{i}(:,end) == 3));
    
    for Feat = 1:length(FeaturesToConsider),
        IndividualL0Dir_Undir_PValues(i,Feat) = kruskalwallis([GroupData{i}(L0DirIndices, Feat); GroupData{i}(UnDirIndices, Feat)], [ones(length(L0DirIndices),1); ones(length(UnDirIndices), 1)*2], 'off');
        GroupMeans{Feat}(i,:) = [mean(GroupData{i}(L0DirIndices, Feat)) mean(GroupData{i}(UnDirIndices, Feat))];
    end
end

% Now to do the same for all motif duration

MinTrialNo = 7; % for all motifs and for FF and amplitude

for i = 1:length(DirBout_Stats)
    disp(['Bird #', num2str(i), ': ', BirdNames{i}]);
    AllMotifDur_GroupData{i} = [];
    % GroupData has 3 columns, the first one has the data corresponding to
    % the variable that that I'm looking at like TotalINNumber, etc. 2nd
    % column has the index of the distance at which the recordings were
    % done (0 for L0, 1 for L1 ... 6 for UN) and 3rd column has 1 for Dir,
    % 2 for Dir/Undir - DUN - and 3 for Undir
    
    for j = 1:length(DirBout_Stats{i}),
        if (~isempty(DirBout_Stats{i}{j}))
            AllMotifDur_GroupData{i} = [AllMotifDur_GroupData{i}; [DirBout_Stats{i}{j}.AllMotifDuration(:) ones(size(DirBout_Stats{i}{j}.AllMotifDuration(:)))*j ones(size(DirBout_Stats{i}{j}.AllMotifDuration(:)))*1]];
        end
    end
    
    for j = 1:length(DirUndirBout_Stats{i}),
        if (~isempty(DirUndirBout_Stats{i}{j}))
            AllMotifDur_GroupData{i} = [AllMotifDur_GroupData{i}; [DirUndirBout_Stats{i}{j}.AllMotifDuration(:) ones(size(DirUndirBout_Stats{i}{j}.AllMotifDuration(:)))*j ones(size(DirUndirBout_Stats{i}{j}.AllMotifDuration(:)))*2]];
        end
    end
    
    for j = 1:length(UndirBout_Stats{i}),
        if (~isempty(UndirBout_Stats{i}{j}))
            AllMotifDur_GroupData{i} = [AllMotifDur_GroupData{i}; [UndirBout_Stats{i}{j}.AllMotifDuration(:) ones(size(UndirBout_Stats{i}{j}.AllMotifDuration(:)))*j ones(size(UndirBout_Stats{i}{j}.AllMotifDuration(:)))*3]];
        end
    end

    % Remove NaN values
    NaNIndices = find(isnan(AllMotifDur_GroupData{i}(:,1)));
    if (~isempty(NaNIndices))
        AllMotifDur_GroupData{i}(NaNIndices,:) = [];
    end
    
    % Remove outliers based on the criterion that data is either 3*IQR 
    % above 75th percentile or below 25th percentile
    OutlierThreshold = [(prctile(AllMotifDur_GroupData{i}(:,1), 25) - 3*iqr(AllMotifDur_GroupData{i}(:,1))) (prctile(AllMotifDur_GroupData{i}(:,1), 75) + 3*iqr(AllMotifDur_GroupData{i}(:,1)))];
    OutlierIndices = find((AllMotifDur_GroupData{i}(:,1) < OutlierThreshold(1)) | (AllMotifDur_GroupData{i}(:,1) > OutlierThreshold(2)));
    if (~isempty(OutlierIndices))
        AllMotifDur_GroupData{i}(OutlierIndices, :) = [];
        disp(['Removed ', num2str(length(OutlierIndices)), ' from all motif durations:', BirdNames{i}]);
    end
    
    % Now find and remove any distance that doesn't have 3 bouts of 
    % directed song
    for j = 1:max(AllMotifDur_GroupData{i}(:,end-1)),
        DirIndicesForDistance = find((AllMotifDur_GroupData{i}(:,end-1) == j) & (AllMotifDur_GroupData{i}(:,end) == 1));
        if ((length(DirIndicesForDistance) > 0) & (length(DirIndicesForDistance) < MinTrialNo))
            AllMotifDur_GroupData{i}(DirIndicesForDistance,:) = [];
            disp(['Removed distance ', DistanceCode{j}, ' for ', BirdNames{i}, ' as there were only ', num2str(length(DirIndicesForDistance)), ' motifs - all motif durations']);
        end
    end
    
    L0DirIndices = find((AllMotifDur_GroupData{i}(:,end) == 1) & (AllMotifDur_GroupData{i}(:,end-1) == 1));
    UnDirIndices = find((AllMotifDur_GroupData{i}(:,end-1) == 6) & (AllMotifDur_GroupData{i}(:,end) == 3));
    
    [Hyp, PValue] = ttest2(AllMotifDur_GroupData{i}(L0DirIndices, 1), AllMotifDur_GroupData{i}(UnDirIndices, 1));
    IndividualL0Dir_Undir_PValues(i,3) = PValue;
    GroupMeans{3}(i,:) = [mean(AllMotifDur_GroupData{i}(L0DirIndices, 1)) mean(AllMotifDur_GroupData{i}(UnDirIndices, 1))];
end

YLabelString = {'Mean # of INs' 'Mean # of motifs/bout' 'Mean motif duration (msec)'};
figure;
p = panel();
p.pack('h', {1/3 1/3 1/3});

for i = 1:length(GroupMeans),
    p(i).select();
    hold on;
    for j = 1:size(GroupMeans{i},2),
        GroupMeanBar{i}(j) = bar(j, mean(GroupMeans{i}(:,j)));
        set(GroupMeanBar{i}(j), 'FaceColor', 'none');
        errorbar(j, mean(GroupMeans{i}(:,j)), std(GroupMeans{i}(:,j))/sqrt(length(GroupMeans{i}(:,j))), 'k.');
    end
    for j = 1:size(GroupMeans{i},1),
        if (IndividualL0Dir_Undir_PValues(j,i) < 0.05)
            MarkerFill = 'k';
        else
            MarkerFill = 'none';
        end
        plot([1.15 1.85], GroupMeans{i}(j,:), 'ko-', 'MarkerFaceColor', MarkerFill, 'MarkerSize', 6);
    end
    if (i < length(GroupMeans))
        axis([0.5 2.5 0 1.1*max(GroupMeans{i}(:))]);
    else
        axis([0.5 2.5 0.99*min(GroupMeans{i}(:)) 1.1*max(GroupMeans{i}(:))]);
    end
    GroupPValue(i) = signrank(GroupMeans{i}(:,1), GroupMeans{i}(:,2));
    sigstar({[1 2]}, GroupPValue(i));
    set(gca, 'XTick', [1 2], 'XTickLabel', {'L0 Songs' 'UNDIR'}, 'XTickLabelRotation', 45);
    ylabel(YLabelString{i});
    SigBirds = find(IndividualL0Dir_Undir_PValues(:,i) < 0.05);
    L0MoreThanUndir = intersect(SigBirds, find(GroupMeans{i}(:,1) > GroupMeans{i}(:,2)));
    disp([YLabelString{i}, ': # of individual birds with sig diffs = ', num2str(length(SigBirds)), ' / ', num2str(size(GroupMeans{i}, 1)), ' birds; # with L0 > Undir = ', num2str(length(L0MoreThanUndir)), '; group p-value = ', num2str(GroupPValue(i)), ' Wilcoxon sign-rank test']); 
end

p.fontsize = 12;
% p.fontname = 'Times';
% p.margintop = 10;
p.marginleft = 20;
% p.marginright = 20;
p.de.margin = 20;
p.marginbottom = 20;

set(gcf, 'Color', 'w');
set(gcf, 'Units', 'inches');
set(gcf, 'Position', [1 1 8.3 3.5]);
set(gcf, 'PaperPositionMode', 'auto');
print(fullfile(FinalFigureDir, ['L0Dir_Undir_Diffs_INNum_MotifNum_MotifDur.eps']), '-depsc2', '-r600');
print(fullfile(FinalFigureDir, ['L0Dir_Undir_Diffs_INNum_MotifNum_MotifDur.png']), '-dpng', '-r600');

disp('Done with stats and plots');
