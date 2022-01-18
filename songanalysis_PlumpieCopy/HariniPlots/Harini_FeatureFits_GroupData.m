function [L0Dir_Undir_PValue] = Harini_FeatureFits_GroupData(DirBout_Stats, DirUndirBout_Stats, UndirBout_Stats, StructString, MinTrialNo, Distances, LabelString, BirdNames)

FinalFigureDir = '/home/raghav/StudentRelated/Harini/Manuscript/results';
DistanceCode = {'L0' 'L1' 'L2' 'L3' 'L4' 'UN'};

% We want to do the analysis per bird and for the group as a whole. The
% things that we will analyze are (1) # of INs at the beginning based on
% the 500ms criterion, (2) # of complete motifs in the bout, (3) all and
% first motif duration, (4), CV of FF and (5) amplitude of individual
% motif syllables only for HF mic birds

FeaturesToConsider = {'TotalINNumber_500ms' 'CompleteMotifNumber' 'FirstMotifDuration'};
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
    
    % Remove outliers only for first motif duration based on the criterion
    % that data is either 3*IQR above 75th percentile or below 25th
    % percentile
    MotifDurationFeatures = find(cellfun(@length, strfind(FeaturesToConsider, 'MotifDuration')));
    for FeatToRemoveOutliers = MotifDurationFeatures(:)',
        OutlierThreshold = [(prctile(GroupData{i}(:,FeatToRemoveOutliers), 25) - 3*iqr(GroupData{i}(:,FeatToRemoveOutliers))) (prctile(GroupData{i}(:,FeatToRemoveOutliers), 75) + 3*iqr(GroupData{i}(:,FeatToRemoveOutliers)))];
        OutlierIndices = find((GroupData{i}(:,FeatToRemoveOutliers) < OutlierThreshold(1)) | (GroupData{i}(:,FeatToRemoveOutliers) > OutlierThreshold(2)));
        if (~isempty(OutlierIndices))
            GroupData{i}(OutlierIndices, FeatToRemoveOutliers) = NaN;
            disp(['Removed ', num2str(length(OutlierIndices)), ' from ', FeaturesToConsider{FeatToRemoveOutliers}, ':', BirdNames{i}]);
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
end

% Now get means and std for different distances
for i = 1:length(FeaturesToConsider),
    DirDistanceMeans{i} = ones(length(DirBout_Stats), length(DirBout_Stats{1}))*NaN;
    DirDistanceSEMs{i} = ones(length(DirBout_Stats), length(DirBout_Stats{1}))*NaN;
    DirDistanceSTDs{i} = ones(length(DirBout_Stats), length(DirBout_Stats{1}))*NaN;
    UndirMeans{i} = ones(length(DirBout_Stats), 1)*NaN;
    UndirSEMs{i} = ones(length(DirBout_Stats), 1)*NaN;
end

for i = 1:length(GroupData),
    for Feat = 1:length(FeaturesToConsider),
        for j = 1:length(DirBout_Stats{1}),
            DirIndices = find((GroupData{i}(:,end-1) == j) & (GroupData{i}(:,end) == 1));
            DirDistanceMeans{Feat}(i,j) = mean(GroupData{i}(DirIndices,Feat));
            DirDistanceSEMs{Feat}(i,j) = std(GroupData{i}(DirIndices,Feat));
            DirDistanceSTDs{Feat}(i,j) = std(GroupData{i}(DirIndices,Feat));
            if (j == 6)
                UndirIndices = find((GroupData{i}(:,end-1) == j) & (GroupData{i}(:,end) == 3));
                UndirMeans{Feat}(i,j) = mean(GroupData{i}(UndirIndices,Feat));
                UndirSEMs{Feat}(i,j) = std(GroupData{i}(UndirIndices,Feat));
            end
        end
    end
end

% Now do the group data plots, normalize with respect to directed song at
% L0.

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
end

% Now get means and std for different distances
AllMotifDur_DirDistanceMeans = ones(length(DirBout_Stats), length(DirBout_Stats{1}))*NaN;
AllMotifDur_DirDistanceSEMs = ones(length(DirBout_Stats), length(DirBout_Stats{1}))*NaN;
AllMotifDur_DirDistanceSTDs = ones(length(DirBout_Stats), length(DirBout_Stats{1}))*NaN;
AllMotifDur_UndirMeans = ones(length(DirBout_Stats), 1)*NaN;
AllMotifDur_UndirSEMs = ones(length(DirBout_Stats), 1)*NaN;

for i = 1:length(AllMotifDur_GroupData),
    for j = 1:length(DirBout_Stats{1}),
        DirIndices = find((AllMotifDur_GroupData{i}(:,end-1) == j) & (AllMotifDur_GroupData{i}(:,end) == 1));
        AllMotifDur_DirDistanceMeans(i,j) = mean(AllMotifDur_GroupData{i}(DirIndices,1));
        AllMotifDur_DirDistanceSEMs(i,j) = std(AllMotifDur_GroupData{i}(DirIndices,1));
        AllMotifDur_DirDistanceSTDs(i,j) = std(AllMotifDur_GroupData{i}(DirIndices,1));
        if (j == 6)
            UndirIndices = find((AllMotifDur_GroupData{i}(:,end-1) == j) & (AllMotifDur_GroupData{i}(:,end) == 3));
            AllMotifDur_UndirMeans(i,j) = mean(AllMotifDur_GroupData{i}(UndirIndices,1));
            AllMotifDur_UndirSEMs(i,j) = std(AllMotifDur_GroupData{i}(UndirIndices,1));
        end
    end
end

MinTrialNo = 7;
% Now to do the log amplitude calculated similar to Kao et al.
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
                SyllDurationGroupData{i}{SyllIndex} = [SyllDurationGroupData{i}{SyllIndex}; [DirUndirBout_Stats{i}{j}.SyllDuration{SyllIndex}(:) ones(size(DirUndirBout_Stats{i}{j}.SyllDuration{SyllIndex}(:)))*j ones(size(DirUndirBout_Stats{i}{j}.SyllDuration{SyllIndex}(:)))*2]];
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
        
        % Remove outliers based on the criterion that data is either 3*IQR 
        % above 75th percentile or below 25th percentile for syllable
        % duration and then for amplitude
        OutlierThreshold = [(prctile(SyllDurationGroupData{i}{SyllIndex}(:,1), 25) - 3*iqr(SyllDurationGroupData{i}{SyllIndex}(:,1))) (prctile(SyllDurationGroupData{i}{SyllIndex}(:,1), 75) + 3*iqr(SyllDurationGroupData{i}{SyllIndex}(:,1)))];
        OutlierIndices = find((SyllDurationGroupData{i}{SyllIndex}(:,1) < OutlierThreshold(1)) | (SyllDurationGroupData{i}{SyllIndex}(:,1) > OutlierThreshold(2)));
        if (~isempty(OutlierIndices))
            SyllDurationGroupData{i}{SyllIndex}(OutlierIndices, :) = [];
            LogAmplitudeGroupData{i}{SyllIndex}(OutlierIndices, :) = [];
            disp(['Removed ', num2str(length(OutlierIndices)), ' from syllable ', UndirBout_Stats{i}{end}.LogAmplitude_SyllLabel(SyllIndex), ' based on syll duration:', BirdNames{i}]);
        end
    
        % Now find and remove any distance that doesn't have 3 bouts of 
        % directed song
        for j = 1:max(LogAmplitudeGroupData{i}{SyllIndex}(:,end-1)),
            DirIndicesForDistance = find((LogAmplitudeGroupData{i}{SyllIndex}(:,end-1) == j) & (LogAmplitudeGroupData{i}{SyllIndex}(:,end) == 1));
            if ((length(DirIndicesForDistance) > 0) & (length(DirIndicesForDistance) < MinTrialNo))
                LogAmplitudeGroupData{i}{SyllIndex}(DirIndicesForDistance,:) = [];
                disp(['Removed distance ', DistanceCode{j}, ' for ', BirdNames{i}, ' as there were only ', num2str(length(DirIndicesForDistance)), ' syllables for syllable ', UndirBout_Stats{i}{end}.LogAmplitude_SyllLabel(SyllIndex), ' - log ampitude']);
            end
        end
    end
end

% Now get means and std for different distances
AllLogAmplitude_DirDistanceMeans = ones(length(DirBout_Stats), length(DirBout_Stats{1}))*NaN;
AllLogAmplitude_DirDistanceSEMs = ones(length(DirBout_Stats), length(DirBout_Stats{1}))*NaN;
AllLogAmplitude_DirDistanceSTDs = ones(length(DirBout_Stats), length(DirBout_Stats{1}))*NaN;
AllLogAmplitude_UndirMeans = ones(length(DirBout_Stats), 1)*NaN;
AllLogAmplitude_UndirSEMs = ones(length(DirBout_Stats), 1)*NaN;

for i = 1:length(LogAmplitudeGroupData),
    TempAmplitudeData = [];
    TempAmplitudeUnDirData = [];
    for j = 1:length(DirBout_Stats{1}),
        for k = 1:length(LogAmplitudeGroupData{i}),
            DirIndices = find((LogAmplitudeGroupData{i}{k}(:,end-1) == j) & (LogAmplitudeGroupData{i}{k}(:,end) == 1));
            TempAmplitudeData(k,j) = mean(LogAmplitudeGroupData{i}{k}(DirIndices,1));
        end
        if (j == 6)
            UnDirIndices = find((LogAmplitudeGroupData{i}{k}(:,end-1) == j) & (LogAmplitudeGroupData{i}{k}(:,end) == 3));
            TempAmplitudeUnDirData(k) = mean(LogAmplitudeGroupData{i}{k}(UnDirIndices,1));
        end
    end
    if (k > 1)
        AllLogAmplitude_DirDistanceMeans(i,:) = mean(TempAmplitudeData);
        AllLogAmplitude_UndirMeans(i,:) = mean(TempAmplitudeUnDirData);
    else
        AllLogAmplitude_DirDistanceMeans(i,:) = TempAmplitudeData;
        AllLogAmplitude_UndirMeans(i,:) = TempAmplitudeUnDirData;
    end
end

MinTrialNo = 10;
%Now to do the FF calculations
for i = 1:length(DirBout_Stats)
    disp(['Bird #', num2str(i), ': ', BirdNames{i}]);
%     GroupData has 3 columns, the first one has the data corresponding to
%     the variable that that I'm looking at like TotalINNumber, etc. 2nd
%     column has the index of the distance at which the recordings were
%     done (0 for L0, 1 for L1 ... 6 for UN) and 3rd column has 1 for Dir,
%     2 for Dir/Undir - DUN - and 3 for Undir
    if (~isfield(UndirBout_Stats{i}{end}, 'FF_SyllLabel'))
        continue;
    end
    clear SyllDurationGroupData;
    
    for SyllIndex = 1:length(UndirBout_Stats{i}{end}.FF_SyllLabel),
        FFGroupData{i}{SyllIndex} = [];
        SyllDurationGroupData{i}{SyllIndex} = [];
        for j = 1:length(DirBout_Stats{i}),
            if (~isempty(DirBout_Stats{i}{j}))
                FFGroupData{i}{SyllIndex} = [FFGroupData{i}{SyllIndex}; [DirBout_Stats{i}{j}.FF{SyllIndex}(:) ones(size(DirBout_Stats{i}{j}.FF{SyllIndex}(:)))*j ones(size(DirBout_Stats{i}{j}.FF{SyllIndex}(:)))*1]];
                SyllDurationGroupData{i}{SyllIndex} = [SyllDurationGroupData{i}{SyllIndex}; [DirBout_Stats{i}{j}.FF_SyllDuration{SyllIndex}(:) ones(size(DirBout_Stats{i}{j}.FF_SyllDuration{SyllIndex}(:)))*j ones(size(DirBout_Stats{i}{j}.FF_SyllDuration{SyllIndex}(:)))*1]];
            end
        end
    
        for j = 1:length(DirUndirBout_Stats{i}),
            if (~isempty(DirUndirBout_Stats{i}{j}))
                FFGroupData{i}{SyllIndex} = [FFGroupData{i}{SyllIndex}; [DirUndirBout_Stats{i}{j}.FF{SyllIndex}(:) ones(size(DirUndirBout_Stats{i}{j}.FF{SyllIndex}(:)))*j ones(size(DirUndirBout_Stats{i}{j}.FF{SyllIndex}(:)))*2]];
                SyllDurationGroupData{i}{SyllIndex} = [SyllDurationGroupData{i}{SyllIndex}; [DirUndirBout_Stats{i}{j}.FF_SyllDuration{SyllIndex}(:) ones(size(DirUndirBout_Stats{i}{j}.FF_SyllDuration{SyllIndex}(:)))*j ones(size(DirUndirBout_Stats{i}{j}.FF_SyllDuration{SyllIndex}(:)))*2]];
            end
        end
        
        for j = 1:length(UndirBout_Stats{i}),
            if (~isempty(UndirBout_Stats{i}{j}))
                FFGroupData{i}{SyllIndex} = [FFGroupData{i}{SyllIndex}; [UndirBout_Stats{i}{j}.FF{SyllIndex}(:) ones(size(UndirBout_Stats{i}{j}.FF{SyllIndex}(:)))*j ones(size(UndirBout_Stats{i}{j}.FF{SyllIndex}(:)))*3]];
                SyllDurationGroupData{i}{SyllIndex} = [SyllDurationGroupData{i}{SyllIndex}; [UndirBout_Stats{i}{j}.FF_SyllDuration{SyllIndex}(:) ones(size(UndirBout_Stats{i}{j}.FF_SyllDuration{SyllIndex}(:)))*j ones(size(UndirBout_Stats{i}{j}.FF_SyllDuration{SyllIndex}(:)))*3]];
            end
        end

        Remove NaN values
        NaNIndices = find(isnan(FFGroupData{i}{SyllIndex}(:,1)));
        if (~isempty(NaNIndices))
            FFGroupData{i}{SyllIndex}(NaNIndices,:) = [];
            SyllDurationGroupData{i}{SyllIndex}(NaNIndices,:) = [];
        end
        
        % Remove outliers based on the criterion that data is either 3*IQR 
        % above 75th percentile or below 25th percentile for syllable
        % duration and then for amplitude
        OutlierThreshold = [(prctile(SyllDurationGroupData{i}{SyllIndex}(:,1), 25) - 2*iqr(SyllDurationGroupData{i}{SyllIndex}(:,1))) (prctile(SyllDurationGroupData{i}{SyllIndex}(:,1), 75) + 2*iqr(SyllDurationGroupData{i}{SyllIndex}(:,1)))];
        OutlierIndices = find((SyllDurationGroupData{i}{SyllIndex}(:,1) < OutlierThreshold(1)) | (SyllDurationGroupData{i}{SyllIndex}(:,1) > OutlierThreshold(2)));
        if (~isempty(OutlierIndices))
            SyllDurationGroupData{i}{SyllIndex}(OutlierIndices, :) = [];
            FFGroupData{i}{SyllIndex}(OutlierIndices, :) = [];
            disp(['Removed ', num2str(length(OutlierIndices)), ' from syllable ', UndirBout_Stats{i}{end}.FF_SyllLabel(SyllIndex), ' based on syll duration:', BirdNames{i}]);
        end
    
        % Remove outliers based on the criterion that data is either 3*IQR 
        % above 75th percentile or below 25th percentile for syllable
        % duration and then for amplitude
        OutlierThreshold = [(prctile(FFGroupData{i}{SyllIndex}(:,1), 25) - 2*iqr(FFGroupData{i}{SyllIndex}(:,1))) (prctile(FFGroupData{i}{SyllIndex}(:,1), 75) + 2*iqr(FFGroupData{i}{SyllIndex}(:,1)))];
        OutlierIndices = find((FFGroupData{i}{SyllIndex}(:,1) < OutlierThreshold(1)) | (FFGroupData{i}{SyllIndex}(:,1) > OutlierThreshold(2)));
        if (~isempty(OutlierIndices))
            SyllDurationGroupData{i}{SyllIndex}(OutlierIndices, :) = [];
            FFGroupData{i}{SyllIndex}(OutlierIndices, :) = [];
            disp(['Removed ', num2str(length(OutlierIndices)), ' from syllable ', UndirBout_Stats{i}{end}.FF_SyllLabel(SyllIndex), ' based on syll FF:', BirdNames{i}]);
        end
        
%         Now find and remove any distance that doesn't have 3 bouts of 
%         directed song
        for j = 1:max(FFGroupData{i}{SyllIndex}(:,end-1)),
            DirIndicesForDistance = find((FFGroupData{i}{SyllIndex}(:,end-1) == j) & (FFGroupData{i}{SyllIndex}(:,end) == 1));
            if ((length(DirIndicesForDistance) > 0) & (length(DirIndicesForDistance) < MinTrialNo))
                FFGroupData{i}{SyllIndex}(DirIndicesForDistance,:) = [];
                disp(['Removed distance ', DistanceCode{j}, ' for ', BirdNames{i}, ' as there were only ', num2str(length(DirIndicesForDistance)), ' syllables for syllable ', UndirBout_Stats{i}{end}.FF_SyllLabel(SyllIndex), ' - log ampitude']);
            end
        end
    
        CV = @(x)std(x)/mean(x);
        
        FFCVGroupData{i}{SyllIndex} = [];
        for j = 1:max(FFGroupData{i}{SyllIndex}(:,end-1)),
            DirIndicesForDistance = find((FFGroupData{i}{SyllIndex}(:,end-1) == j) & (FFGroupData{i}{SyllIndex}(:,end) == 1));
            if (~isempty(DirIndicesForDistance))
                FFCVGroupData{i}{SyllIndex} = [FFCVGroupData{i}{SyllIndex}; [jackknife(CV, FFGroupData{i}{SyllIndex}(DirIndicesForDistance,1)) ones(length(DirIndicesForDistance),1)*j ones(length(DirIndicesForDistance),1)*1]];
            end
        end
        
        FFUndirCVGroupData{i}{SyllIndex} = [];
        UnDirIndicesForDistance = find((FFGroupData{i}{SyllIndex}(:,end-1) == 6) & (FFGroupData{i}{SyllIndex}(:,end) == 3));
        if (~isempty(UnDirIndicesForDistance))
            FFUndirCVGroupData{i}{SyllIndex} = [FFUndirCVGroupData{i}{SyllIndex}; [jackknife(CV, FFGroupData{i}{SyllIndex}(UnDirIndicesForDistance,1)) ones(length(UnDirIndicesForDistance),1)*6 ones(length(UnDirIndicesForDistance),1)*3]];
        end
    end
end

disp('Done with stats and plots');
