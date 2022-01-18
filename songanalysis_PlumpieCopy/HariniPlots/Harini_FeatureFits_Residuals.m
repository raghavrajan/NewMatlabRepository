function [L0Dir_Undir_PValue] = Harini_FeatureFits_Residuals(DirBout_Stats, DirUndirBout_Stats, UndirBout_Stats, StructString, MinTrialNo, Distances, LabelString, BirdNames)

FinalFigureDir = '/home/raghav/StudentRelated/Harini/Manuscript/results';
DistanceCode = {'L0' 'L1' 'L2' 'L3' 'L4' 'UN'};

% We want to do the analysis per bird and for the group as a whole. The
% things that we will analyze are (1) # of INs at the beginning based on
% the 500ms criterion, (2) # of complete motifs in the bout, (3) all and
% first motif duration, (4), CV of FF and (5) amplitude of individual
% motif syllables only for HF mic birds

% FeaturesToConsider = {'TotalINNumber_500ms' 'CompleteMotifNumber' 'FirstMotifDuration'};
% MinTrialNo = 3;
% for i = 1:length(DirBout_Stats)
%     disp(['Bird #', num2str(i), ': ', BirdNames{i}]);
%     GroupData{i} = [];
%     % GroupData has 3 columns, the first one has the data corresponding to
%     % the variable that that I'm looking at like TotalINNumber, etc. 2nd
%     % column has the index of the distance at which the recordings were
%     % done (0 for L0, 1 for L1 ... 6 for UN) and 3rd column has 1 for Dir,
%     % 2 for Dir/Undir - DUN - and 3 for Undir
%     
%     for j = 1:length(DirBout_Stats{i}),
%         if (~isempty(DirBout_Stats{i}{j}))
%             TempGroupData = [];
%             for k = 1:length(FeaturesToConsider),
%                 if (k == length(FeaturesToConsider))
%                     TempGroupData = [TempGroupData [eval(['DirBout_Stats{i}{j}.', FeaturesToConsider{k}, '(:)']) ones(size(eval(['DirBout_Stats{i}{j}.', FeaturesToConsider{k}, '(:)'])))*j ones(size(eval(['DirBout_Stats{i}{j}.', FeaturesToConsider{k}, '(:)'])))*1]];
%                 else
%                     TempGroupData = [TempGroupData eval(['DirBout_Stats{i}{j}.', FeaturesToConsider{k}, '(:)'])];
%                 end
%             end
%             GroupData{i} = [GroupData{i}; TempGroupData];
%         end
%     end
%     
%     for j = 1:length(DirUndirBout_Stats{i}),
%         if (~isempty(DirUndirBout_Stats{i}{j}))
%             TempGroupData = [];
%             for k = 1:length(FeaturesToConsider),
%                 if (k == length(FeaturesToConsider))
%                     TempGroupData = [TempGroupData [eval(['DirUndirBout_Stats{i}{j}.', FeaturesToConsider{k}, '(:)']) ones(size(eval(['DirUndirBout_Stats{i}{j}.', FeaturesToConsider{k}, '(:)'])))*j ones(size(eval(['DirUndirBout_Stats{i}{j}.', FeaturesToConsider{k}, '(:)'])))*2]];
%                 else
%                     TempGroupData = [TempGroupData eval(['DirUndirBout_Stats{i}{j}.', FeaturesToConsider{k}, '(:)'])];
%                 end
%             end
%             GroupData{i} = [GroupData{i}; TempGroupData];
%         end
%     end
%     
%     for j = 1:length(UndirBout_Stats{i}),
%         if (~isempty(UndirBout_Stats{i}{j}))
%             TempGroupData = [];
%             for k = 1:length(FeaturesToConsider),
%                 if (k == length(FeaturesToConsider))
%                     TempGroupData = [TempGroupData [eval(['UndirBout_Stats{i}{j}.', FeaturesToConsider{k}, '(:)']) ones(size(eval(['UndirBout_Stats{i}{j}.', FeaturesToConsider{k}, '(:)'])))*j ones(size(eval(['UndirBout_Stats{i}{j}.', FeaturesToConsider{k}, '(:)'])))*3]];
%                 else
%                     TempGroupData = [TempGroupData eval(['UndirBout_Stats{i}{j}.', FeaturesToConsider{k}, '(:)'])];
%                 end
%             end
%             GroupData{i} = [GroupData{i}; TempGroupData];
%         end
%     end
% 
%     % Remove NaN values
%     for k = 1:length(FeaturesToConsider),
%         NaNIndices = find(isnan(GroupData{i}(:,k)));
%         if (~isempty(NaNIndices))
%             GroupData{i}(NaNIndices,:) = [];
%         end
%     end
%     
%     % Remove outliers only for first motif duration based on the criterion
%     % that data is either 3*IQR above 75th percentile or below 25th
%     % percentile
%     MotifDurationFeatures = find(cellfun(@length, strfind(FeaturesToConsider, 'MotifDuration')));
%     for FeatToRemoveOutliers = MotifDurationFeatures(:)',
%         OutlierThreshold = [(prctile(GroupData{i}(:,FeatToRemoveOutliers), 25) - 3*iqr(GroupData{i}(:,FeatToRemoveOutliers))) (prctile(GroupData{i}(:,FeatToRemoveOutliers), 75) + 3*iqr(GroupData{i}(:,FeatToRemoveOutliers)))];
%         OutlierIndices = find((GroupData{i}(:,FeatToRemoveOutliers) < OutlierThreshold(1)) | (GroupData{i}(:,FeatToRemoveOutliers) > OutlierThreshold(2)));
%         if (~isempty(OutlierIndices))
%             GroupData{i}(OutlierIndices, FeatToRemoveOutliers) = NaN;
%             disp(['Removed ', num2str(length(OutlierIndices)), ' from ', FeaturesToConsider{FeatToRemoveOutliers}, ':', BirdNames{i}]);
%         end
%     end
%     
%     % Now find and remove any distance that doesn't have 3 bouts of 
%     % directed song
%     for j = 1:max(GroupData{i}(:,end-1)),
%         DirIndicesForDistance = find((GroupData{i}(:,end-1) == j) & (GroupData{i}(:,end) == 1));
%         if ((length(DirIndicesForDistance) > 0) & (length(DirIndicesForDistance) < MinTrialNo))
%             GroupData{i}(DirIndicesForDistance,:) = [];
%             disp(['Removed distance ', DistanceCode{j}, ' for ', BirdNames{i}, ' as there were only ', num2str(length(DirIndicesForDistance)), ' bouts']);
%         end
%     end
%     
%     UnNormalizedGroupData{i} = GroupData{i};
%     DirIndices = find(GroupData{i}(:,end) == 1);
%     % First plot the group data
%     figure;
%     for Feat = 1:length(FeaturesToConsider),
%         subplot(length(FeaturesToConsider),3,(Feat-1)*length(FeaturesToConsider) + 1);
%         plot(Distances(GroupData{i}(DirIndices, 4)) + (rand(1, length(DirIndices))*4 - 2), GroupData{i}(DirIndices, Feat), 'ko');
%         hold on;
%     
%         % Now plot the mean, SEM and the best fit line with robust regression
%         for Location = 1:max(GroupData{i}(DirIndices,4)),
%             LocIndices = intersect(DirIndices, find(GroupData{i}(:,4) == Location));
%             if (~isempty(LocIndices))
%                 if (length(LocIndices) == 1)
%                     plot(Distances(Location)+3, GroupData{i}(LocIndices, Feat), 'ks', 'MarkerSize', 6, 'MarkerFaceColor', 'k');
%                 else
%                     errorbar(Distances(Location)+3, nanmean(GroupData{i}(LocIndices, Feat)), nanstd(GroupData{i}(LocIndices, Feat))/sqrt(length(LocIndices)), 'ks', 'MarkerSize', 6, 'MarkerFaceColor', 'k');
%                 end
%             end
%         end
%         
%         % Now best fit line
%         RobustCoeffs = robustfit(Distances(GroupData{i}(DirIndices,4)), GroupData{i}(DirIndices, Feat));
%         plot(unique(Distances(GroupData{i}(DirIndices,4))), polyval(flipud(RobustCoeffs), unique(Distances(GroupData{i}(DirIndices,4)))), 'k', 'LineWidth', 1.5);
%         
%         [FixedX_BootStrapCoeffs, Random_BootStrapCoeffs, RandomProb, Rsq] = BootstrapTestRegressionSignificance(RobustCoeffs, Distances(GroupData{i}(DirIndices,4))', GroupData{i}(DirIndices, Feat), 'linear');
%         % [Rsq, F, Prob] = CalculateGoodnessofLinearFit(flipud(RobustCoeffs), Distances(GroupData{i}(DirIndices,4))', GroupData{i}(DirIndices, Feat));
%         axis tight;
%         Temp = axis;
%         Temp = [-5 Temp(2)+5 0.99*Temp(3) 1.2*Temp(4)];
%         axis(Temp);
%         text(Temp(1)+2, Temp(4)/1.15, {['Rsq = ', num2str(Rsq)]; ['p = ', num2str(RandomProb(end))]});
%         text(Temp(1)+2, Temp(4)/1.05, {['Slope = ', num2str(FixedX_BootStrapCoeffs(1,end)), ' : ', num2str(FixedX_BootStrapCoeffs(2,end)), ' : ', num2str(FixedX_BootStrapCoeffs(3,end))]});
%         ylabel(FeaturesToConsider{Feat})
%         if (Feat == length(FeaturesToConsider))
%             xlabel('Distance (cm)');
%         end
%         if (Feat == 1)
%             title([BirdNames{i}, 'Fit to feature value']);
%         end
%         
%         % Now calculate residuals
%         YPred = (Distances(GroupData{i}(DirIndices,4)) * RobustCoeffs(2)) + RobustCoeffs(1);
%         Error = GroupData{i}(DirIndices, Feat) - YPred(:);
%         
%         % Now plot residuals as a function of x-value with mean and
%         % errorbars
%         subplot(length(FeaturesToConsider),3,(Feat-1)*length(FeaturesToConsider) + 2);
%         hold on;
%         plot(Distances(GroupData{i}(DirIndices, 4)), Error, 'ko');
%         for Location = 1:max(GroupData{i}(DirIndices,4)),
%             LocIndices = intersect(DirIndices, find(GroupData{i}(:,4) == Location));
%             if (~isempty(LocIndices))
%                 if (length(LocIndices) == 1)
%                     plot(Distances(Location)+3, Error(LocIndices), 'ks', 'MarkerSize', 6, 'MarkerFaceColor', 'k');
%                 else
%                     errorbar(Distances(Location)+3, mean(Error(LocIndices)), std(Error(LocIndices))/sqrt(length(LocIndices)), 'ks', 'MarkerSize', 6, 'MarkerFaceColor', 'k');
%                 end
%             end
%         end
%         axis tight;
%         Temp = axis;
%         Temp = [-5 Temp(2)+5 1.02*Temp(3) 1.02*Temp(4)];
%         axis(Temp);
%         plot(Temp(1:2), zeros(1,2), 'k--');
%         
%         if (Feat == 1)
%             title('Distance vs. Residuals');
%         end
%         
%         if (Feat == length(FeaturesToConsider))
%             xlabel('Distance (cm)');
%         end
%         ylabel('Residual');
%         
%         subplot(length(FeaturesToConsider),3,(Feat-1)*length(FeaturesToConsider) + 3);
%         hold on;
%         Edges = linspace(min(Error)*0.95, max(Error)*1.05, 20);
%         plot(Edges, histc(Error, Edges), 'ko-');
%         if (Feat == 1)
%             title('Residual distribution');
%         end
%         ylabel('#');
%         if (Feat == length(FeaturesToConsider))
%             xlabel('Residual');
%         end
%     end
%     
%     % GroupData{i}(:,1:length(FeaturesToConsider)) = zscore(GroupData{i}(:,1:length(FeaturesToConsider)));
% end
% 
% % Now to do the same for all motif duration
% 
% MinTrialNo = 7; % for all motifs and for FF and amplitude
% 
% for i = 1:length(DirBout_Stats)
%     disp(['Bird #', num2str(i), ': ', BirdNames{i}]);
%     AllMotifDur_GroupData{i} = [];
%     % GroupData has 3 columns, the first one has the data corresponding to
%     % the variable that that I'm looking at like TotalINNumber, etc. 2nd
%     % column has the index of the distance at which the recordings were
%     % done (0 for L0, 1 for L1 ... 6 for UN) and 3rd column has 1 for Dir,
%     % 2 for Dir/Undir - DUN - and 3 for Undir
%     
%     for j = 1:length(DirBout_Stats{i}),
%         if (~isempty(DirBout_Stats{i}{j}))
%             AllMotifDur_GroupData{i} = [AllMotifDur_GroupData{i}; [DirBout_Stats{i}{j}.AllMotifDuration(:) ones(size(DirBout_Stats{i}{j}.AllMotifDuration(:)))*j ones(size(DirBout_Stats{i}{j}.AllMotifDuration(:)))*1]];
%         end
%     end
%     
%     for j = 1:length(DirUndirBout_Stats{i}),
%         if (~isempty(DirUndirBout_Stats{i}{j}))
%             AllMotifDur_GroupData{i} = [AllMotifDur_GroupData{i}; [DirUndirBout_Stats{i}{j}.AllMotifDuration(:) ones(size(DirUndirBout_Stats{i}{j}.AllMotifDuration(:)))*j ones(size(DirUndirBout_Stats{i}{j}.AllMotifDuration(:)))*2]];
%         end
%     end
%     
%     for j = 1:length(UndirBout_Stats{i}),
%         if (~isempty(UndirBout_Stats{i}{j}))
%             AllMotifDur_GroupData{i} = [AllMotifDur_GroupData{i}; [UndirBout_Stats{i}{j}.AllMotifDuration(:) ones(size(UndirBout_Stats{i}{j}.AllMotifDuration(:)))*j ones(size(UndirBout_Stats{i}{j}.AllMotifDuration(:)))*3]];
%         end
%     end
% 
%     % Remove NaN values
%     NaNIndices = find(isnan(AllMotifDur_GroupData{i}(:,1)));
%     if (~isempty(NaNIndices))
%         AllMotifDur_GroupData{i}(NaNIndices,:) = [];
%     end
%     
%     % Remove outliers based on the criterion that data is either 3*IQR 
%     % above 75th percentile or below 25th percentile
%     OutlierThreshold = [(prctile(AllMotifDur_GroupData{i}(:,1), 25) - 3*iqr(AllMotifDur_GroupData{i}(:,1))) (prctile(AllMotifDur_GroupData{i}(:,1), 75) + 3*iqr(AllMotifDur_GroupData{i}(:,1)))];
%     OutlierIndices = find((AllMotifDur_GroupData{i}(:,1) < OutlierThreshold(1)) | (AllMotifDur_GroupData{i}(:,1) > OutlierThreshold(2)));
%     if (~isempty(OutlierIndices))
%         AllMotifDur_GroupData{i}(OutlierIndices, :) = [];
%         disp(['Removed ', num2str(length(OutlierIndices)), ' from all motif durations:', BirdNames{i}]);
%     end
%     
%     % Now find and remove any distance that doesn't have 3 bouts of 
%     % directed song
%     for j = 1:max(AllMotifDur_GroupData{i}(:,end-1)),
%         DirIndicesForDistance = find((AllMotifDur_GroupData{i}(:,end-1) == j) & (AllMotifDur_GroupData{i}(:,end) == 1));
%         if ((length(DirIndicesForDistance) > 0) & (length(DirIndicesForDistance) < MinTrialNo))
%             AllMotifDur_GroupData{i}(DirIndicesForDistance,:) = [];
%             disp(['Removed distance ', DistanceCode{j}, ' for ', BirdNames{i}, ' as there were only ', num2str(length(DirIndicesForDistance)), ' motifs - all motif durations']);
%         end
%     end
%     
%     UnNormalizedGroupData{i} = AllMotifDur_GroupData{i};
%     DirIndices = find(AllMotifDur_GroupData{i}(:,end) == 1);
%     % First plot the group data
%     figure;
%     subplot(1,3,1);
%     plot(Distances(AllMotifDur_GroupData{i}(DirIndices, end-1)) + (rand(1, length(DirIndices))*4 - 2), AllMotifDur_GroupData{i}(DirIndices, 1), 'ko');
%     hold on;
% 
%     % Now plot the mean, SEM and the best fit line with robust regression
%     for Location = 1:max(AllMotifDur_GroupData{i}(DirIndices,end-1)),
%         LocIndices = intersect(DirIndices, find(AllMotifDur_GroupData{i}(:,end-1) == Location));
%         if (~isempty(LocIndices))
%             if (length(LocIndices) == 1)
%                 plot(Distances(Location)+3, AllMotifDur_GroupData{i}(LocIndices, 1), 'ks', 'MarkerSize', 6, 'MarkerFaceColor', 'k');
%             else
%                 errorbar(Distances(Location)+3, nanmean(AllMotifDur_GroupData{i}(LocIndices, 1)), nanstd(AllMotifDur_GroupData{i}(LocIndices, 1))/sqrt(length(LocIndices)), 'ks', 'MarkerSize', 6, 'MarkerFaceColor', 'k');
%             end
%         end
%     end
% 
%     % Now best fit line
%     RobustCoeffs = robustfit(Distances(AllMotifDur_GroupData{i}(DirIndices,end-1)), AllMotifDur_GroupData{i}(DirIndices, 1));
%     plot(unique(Distances(AllMotifDur_GroupData{i}(DirIndices,end-1))), polyval(flipud(RobustCoeffs), unique(Distances(AllMotifDur_GroupData{i}(DirIndices,end-1)))), 'k', 'LineWidth', 1.5);
% 
%     [FixedX_BootStrapCoeffs, Random_BootStrapCoeffs, RandomProb, Rsq] = BootstrapTestRegressionSignificance(RobustCoeffs, Distances(AllMotifDur_GroupData{i}(DirIndices,end-1))', AllMotifDur_GroupData{i}(DirIndices, 1), 'linear');
%     axis tight;
%     Temp = axis;
%     Temp = [-5 Temp(2)+5 0.99*Temp(3) 1.2*Temp(4)];
%     axis(Temp);
%     text(Temp(1)+2, Temp(4)/1.15, {['Rsq = ', num2str(Rsq)]; ['p = ', num2str(RandomProb(end))]});
%     text(Temp(1)+2, Temp(4)/1.05, {['Slope = ', num2str(FixedX_BootStrapCoeffs(1,end)), ' : ', num2str(FixedX_BootStrapCoeffs(2,end)), ' : ', num2str(FixedX_BootStrapCoeffs(3,end))]});
%     ylabel('All motif durations')
%     xlabel('Distance (cm)');
%     title([BirdNames{i}, 'Fit to feature value']);
% 
%     % Now calculate residuals
%     YPred = (Distances(AllMotifDur_GroupData{i}(DirIndices,end-1)) * RobustCoeffs(2)) + RobustCoeffs(1);
%     Error = AllMotifDur_GroupData{i}(DirIndices, 1) - YPred(:);
% 
%     % Now plot residuals as a function of x-value with mean and
%     % errorbars
%     subplot(1,3,2);
%     hold on;
%     plot(Distances(AllMotifDur_GroupData{i}(DirIndices, end-1)), Error, 'ko');
%     for Location = 1:max(AllMotifDur_GroupData{i}(DirIndices,end-1)),
%         LocIndices = intersect(DirIndices, find(AllMotifDur_GroupData{i}(:,end-1) == Location));
%         if (~isempty(LocIndices))
%             if (length(LocIndices) == 1)
%                 plot(Distances(Location)+3, Error(LocIndices), 'ks', 'MarkerSize', 6, 'MarkerFaceColor', 'k');
%             else
%                 errorbar(Distances(Location)+3, mean(Error(LocIndices)), std(Error(LocIndices))/sqrt(length(LocIndices)), 'ks', 'MarkerSize', 6, 'MarkerFaceColor', 'k');
%             end
%         end
%     end
%     axis tight;
%     Temp = axis;
%     Temp = [-5 Temp(2)+5 1.02*Temp(3) 1.02*Temp(4)];
%     axis(Temp);
%     plot(Temp(1:2), zeros(1,2), 'k--');
%     title('Distance vs. Residuals');
%     xlabel('Distance (cm)');
%     ylabel('Residual');
% 
%     subplot(1,3,3);
%     hold on;
%     Edges = linspace(min(Error)*0.95, max(Error)*1.05, 20);
%     plot(Edges, histc(Error, Edges), 'ko-');
%     title('Residual distribution');
%     ylabel('#');
%     xlabel('Residual');
% end
% 
% MinTrialNo = 7;
% % Now to do the log amplitude calculated similar to Kao et al.
% for i = 1:length(DirBout_Stats)
%     disp(['Bird #', num2str(i), ': ', BirdNames{i}]);
%     % GroupData has 3 columns, the first one has the data corresponding to
%     % the variable that that I'm looking at like TotalINNumber, etc. 2nd
%     % column has the index of the distance at which the recordings were
%     % done (0 for L0, 1 for L1 ... 6 for UN) and 3rd column has 1 for Dir,
%     % 2 for Dir/Undir - DUN - and 3 for Undir
%     if (~isfield(UndirBout_Stats{i}{end}, 'LogAmplitude_SyllLabel'))
%         continue;
%     end
%     for SyllIndex = 1:length(UndirBout_Stats{i}{end}.LogAmplitude_SyllLabel),
%         LogAmplitudeGroupData{i}{SyllIndex} = [];
%         SyllDurationGroupData{i}{SyllIndex} = [];
%         for j = 1:length(DirBout_Stats{i}),
%             if (~isempty(DirBout_Stats{i}{j}))
%                 LogAmplitudeGroupData{i}{SyllIndex} = [LogAmplitudeGroupData{i}{SyllIndex}; [DirBout_Stats{i}{j}.LogAmplitude{SyllIndex}(:) ones(size(DirBout_Stats{i}{j}.LogAmplitude{SyllIndex}(:)))*j ones(size(DirBout_Stats{i}{j}.LogAmplitude{SyllIndex}(:)))*1]];
%                 SyllDurationGroupData{i}{SyllIndex} = [SyllDurationGroupData{i}{SyllIndex}; [DirBout_Stats{i}{j}.SyllDuration{SyllIndex}(:) ones(size(DirBout_Stats{i}{j}.SyllDuration{SyllIndex}(:)))*j ones(size(DirBout_Stats{i}{j}.SyllDuration{SyllIndex}(:)))*1]];
%             end
%         end
%     
%         for j = 1:length(DirUndirBout_Stats{i}),
%             if (~isempty(DirUndirBout_Stats{i}{j}))
%                 LogAmplitudeGroupData{i}{SyllIndex} = [LogAmplitudeGroupData{i}{SyllIndex}; [DirUndirBout_Stats{i}{j}.LogAmplitude{SyllIndex}(:) ones(size(DirUndirBout_Stats{i}{j}.LogAmplitude{SyllIndex}(:)))*j ones(size(DirUndirBout_Stats{i}{j}.LogAmplitude{SyllIndex}(:)))*2]];
%                 SyllDurationGroupData{i}{SyllIndex} = [SyllDurationGroupData{i}{SyllIndex}; [DirUndirBout_Stats{i}{j}.SyllDuration{SyllIndex}(:) ones(size(DirUndirBout_Stats{i}{j}.SyllDuration{SyllIndex}(:)))*j ones(size(DirUndirBout_Stats{i}{j}.SyllDuration{SyllIndex}(:)))*2]];
%             end
%         end
%         
%         for j = 1:length(UndirBout_Stats{i}),
%             if (~isempty(UndirBout_Stats{i}{j}))
%                 LogAmplitudeGroupData{i}{SyllIndex} = [LogAmplitudeGroupData{i}{SyllIndex}; [UndirBout_Stats{i}{j}.LogAmplitude{SyllIndex}(:) ones(size(UndirBout_Stats{i}{j}.LogAmplitude{SyllIndex}(:)))*j ones(size(UndirBout_Stats{i}{j}.LogAmplitude{SyllIndex}(:)))*3]];
%                 SyllDurationGroupData{i}{SyllIndex} = [SyllDurationGroupData{i}{SyllIndex}; [UndirBout_Stats{i}{j}.SyllDuration{SyllIndex}(:) ones(size(UndirBout_Stats{i}{j}.SyllDuration{SyllIndex}(:)))*j ones(size(UndirBout_Stats{i}{j}.SyllDuration{SyllIndex}(:)))*3]];
%             end
%         end
% 
%         % Remove NaN values
%         NaNIndices = find(isnan(LogAmplitudeGroupData{i}{SyllIndex}(:,1)));
%         if (~isempty(NaNIndices))
%             LogAmplitudeGroupData{i}{SyllIndex}(NaNIndices,:) = [];
%             SyllDurationGroupData{i}{SyllIndex}(NaNIndices,:) = [];
%         end
%         
%         % Remove outliers based on the criterion that data is either 3*IQR 
%         % above 75th percentile or below 25th percentile for syllable
%         % duration and then for amplitude
%         OutlierThreshold = [(prctile(SyllDurationGroupData{i}{SyllIndex}(:,1), 25) - 3*iqr(SyllDurationGroupData{i}{SyllIndex}(:,1))) (prctile(SyllDurationGroupData{i}{SyllIndex}(:,1), 75) + 3*iqr(SyllDurationGroupData{i}{SyllIndex}(:,1)))];
%         OutlierIndices = find((SyllDurationGroupData{i}{SyllIndex}(:,1) < OutlierThreshold(1)) | (SyllDurationGroupData{i}{SyllIndex}(:,1) > OutlierThreshold(2)));
%         if (~isempty(OutlierIndices))
%             SyllDurationGroupData{i}{SyllIndex}(OutlierIndices, :) = [];
%             LogAmplitudeGroupData{i}{SyllIndex}(OutlierIndices, :) = [];
%             disp(['Removed ', num2str(length(OutlierIndices)), ' from syllable ', UndirBout_Stats{i}{end}.LogAmplitude_SyllLabel(SyllIndex), ' based on syll duration:', BirdNames{i}]);
%         end
%     
%         % Now find and remove any distance that doesn't have 3 bouts of 
%         % directed song
%         for j = 1:max(LogAmplitudeGroupData{i}{SyllIndex}(:,end-1)),
%             DirIndicesForDistance = find((LogAmplitudeGroupData{i}{SyllIndex}(:,end-1) == j) & (LogAmplitudeGroupData{i}{SyllIndex}(:,end) == 1));
%             if ((length(DirIndicesForDistance) > 0) & (length(DirIndicesForDistance) < MinTrialNo))
%                 LogAmplitudeGroupData{i}{SyllIndex}(DirIndicesForDistance,:) = [];
%                 disp(['Removed distance ', DistanceCode{j}, ' for ', BirdNames{i}, ' as there were only ', num2str(length(DirIndicesForDistance)), ' syllables for syllable ', UndirBout_Stats{i}{end}.LogAmplitude_SyllLabel(SyllIndex), ' - log ampitude']);
%             end
%         end
%     
%         DirIndices = find(LogAmplitudeGroupData{i}{SyllIndex}(:,end) == 1);
%         % First plot the group data
%         figure;
%         subplot(1,3,1);
%         plot(Distances(LogAmplitudeGroupData{i}{SyllIndex}(DirIndices, end-1)) + (rand(1, length(DirIndices))*4 - 2), LogAmplitudeGroupData{i}{SyllIndex}(DirIndices, 1), 'ko');
%         hold on;
% 
%         % Now plot the mean, SEM and the best fit line with robust regression
%         for Location = 1:max(LogAmplitudeGroupData{i}{SyllIndex}(DirIndices,end-1)),
%             LocIndices = intersect(DirIndices, find(LogAmplitudeGroupData{i}{SyllIndex}(:,end-1) == Location));
%             if (~isempty(LocIndices))
%                 if (length(LocIndices) == 1)
%                     plot(Distances(Location)+3, LogAmplitudeGroupData{i}{SyllIndex}(LocIndices, 1), 'ks', 'MarkerSize', 6, 'MarkerFaceColor', 'k');
%                 else
%                     errorbar(Distances(Location)+3, nanmean(LogAmplitudeGroupData{i}{SyllIndex}(LocIndices, 1)), nanstd(LogAmplitudeGroupData{i}{SyllIndex}(LocIndices, 1))/sqrt(length(LocIndices)), 'ks', 'MarkerSize', 6, 'MarkerFaceColor', 'k');
%                 end
%             end
%         end
% 
%         % Now best fit line
%         RobustCoeffs = robustfit(Distances(LogAmplitudeGroupData{i}{SyllIndex}(DirIndices,end-1)), LogAmplitudeGroupData{i}{SyllIndex}(DirIndices, 1));
%         plot(unique(Distances(LogAmplitudeGroupData{i}{SyllIndex}(DirIndices,end-1))), polyval(flipud(RobustCoeffs), unique(Distances(LogAmplitudeGroupData{i}{SyllIndex}(DirIndices,end-1)))), 'k', 'LineWidth', 1.5);
% 
%         [FixedX_BootStrapCoeffs, Random_BootStrapCoeffs, RandomProb, Rsq] = BootstrapTestRegressionSignificance(RobustCoeffs, Distances(LogAmplitudeGroupData{i}{SyllIndex}(DirIndices,end-1))', LogAmplitudeGroupData{i}{SyllIndex}(DirIndices, 1), 'linear');
%         axis tight;
%         Temp = axis;
%         Temp = [-5 Temp(2)+5 0.99*Temp(3) 1.2*Temp(4)];
%         axis(Temp);
%         text(Temp(1)+2, Temp(4)/1.15, {['Rsq = ', num2str(Rsq)]; ['p = ', num2str(RandomProb(end))]});
%         text(Temp(1)+2, Temp(4)/1.05, {['Slope = ', num2str(FixedX_BootStrapCoeffs(1,end)), ' : ', num2str(FixedX_BootStrapCoeffs(2,end)), ' : ', num2str(FixedX_BootStrapCoeffs(3,end))]});
%         ylabel('Syllable ampitude')
%         xlabel('Distance (cm)');
%         title([BirdNames{i}, 'Fit to feature value for syllable ', UndirBout_Stats{i}{end}.LogAmplitude_SyllLabel(SyllIndex)]);
% 
%         % Now calculate residuals
%         YPred = (Distances(LogAmplitudeGroupData{i}{SyllIndex}(DirIndices,end-1)) * RobustCoeffs(2)) + RobustCoeffs(1);
%         Error = LogAmplitudeGroupData{i}{SyllIndex}(DirIndices, 1) - YPred(:);
% 
%         % Now plot residuals as a function of x-value with mean and
%         % errorbars
%         subplot(1,3,2);
%         hold on;
%         plot(Distances(LogAmplitudeGroupData{i}{SyllIndex}(DirIndices, end-1)), Error, 'ko');
%         for Location = 1:max(LogAmplitudeGroupData{i}{SyllIndex}(DirIndices,end-1)),
%             LocIndices = intersect(DirIndices, find(LogAmplitudeGroupData{i}{SyllIndex}(:,end-1) == Location));
%             if (~isempty(LocIndices))
%                 if (length(LocIndices) == 1)
%                     plot(Distances(Location)+3, Error(LocIndices), 'ks', 'MarkerSize', 6, 'MarkerFaceColor', 'k');
%                 else
%                     errorbar(Distances(Location)+3, mean(Error(LocIndices)), std(Error(LocIndices))/sqrt(length(LocIndices)), 'ks', 'MarkerSize', 6, 'MarkerFaceColor', 'k');
%                 end
%             end
%         end
%         axis tight;
%         Temp = axis;
%         Temp = [-5 Temp(2)+5 1.02*Temp(3) 1.02*Temp(4)];
%         axis(Temp);
%         plot(Temp(1:2), zeros(1,2), 'k--');
%         title('Distance vs. Residuals');
%         xlabel('Distance (cm)');
%         ylabel('Residual');
% 
%         subplot(1,3,3);
%         hold on;
%         Edges = linspace(min(Error)*0.95, max(Error)*1.05, 20);
%         plot(Edges, histc(Error, Edges), 'ko-');
%         title('Residual distribution');
%         ylabel('#');
%         xlabel('Residual');
% 
%     end
% end

MinTrialNo = 10;
% Now to do the FF calculations
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

        % Remove NaN values
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
        
        % Now find and remove any distance that doesn't have 3 bouts of 
        % directed song
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
        
        DirIndices = find(FFCVGroupData{i}{SyllIndex}(:,end) == 1);
        % First plot the group data
        figure;
        subplot(1,3,1);
%        plot(Distances(FFCVGroupData{i}{SyllIndex}(DirIndices, end-1)) + (rand(1, length(DirIndices))*4 - 2), FFCVGroupData{i}{SyllIndex}(DirIndices, 1), 'ko');
        hold on;

        % Now plot the mean, SEM and the best fit line with robust regression
%         for Location = 1:max(FFCVGroupData{i}{SyllIndex}(DirIndices,end-1)),
%             LocIndices = int  ersect(DirIndices, find(FFCVGroupData{i}{SyllIndex}(:,end-1) == Location));
%             if (~isempty(LocIndices))
%                 if (length(LocIndices) == 1)
%                     plot(Distances(Location)+3, FFCVGroupData{i}{SyllIndex}(LocIndices, 1), 'ks', 'MarkerSize', 6, 'MarkerFaceColor', 'k');
%                 else
%                     errorbar(Distances(Location)+3, nanmean(FFCVGroupData{i}{SyllIndex}(LocIndices, 1)), nanstd(FFCVGroupData{i}{SyllIndex}(LocIndices, 1))/sqrt(length(LocIndices)), 'ks', 'MarkerSize', 6, 'MarkerFaceColor', 'k');
%                 end
%             end
%         end


        DirIndices = find(FFGroupData{i}{SyllIndex}(:,end) == 1);
        for Location = 1:max(FFGroupData{i}{SyllIndex}(DirIndices,end-1)),
            LocIndices = intersect(DirIndices, find(FFGroupData{i}{SyllIndex}(:,end-1) == Location));
            if (~isempty(LocIndices))
                plot(Distances(Location), nanstd(FFGroupData{i}{SyllIndex}(LocIndices, 1))/nanmean(FFGroupData{i}{SyllIndex}(LocIndices, 1)), 'ks', 'MarkerSize', 6, 'MarkerFaceColor', 'k');
            end
        end
        UndirIndices = find((FFGroupData{i}{SyllIndex}(:,end) == 3) & ((FFGroupData{i}{SyllIndex}(:,end-1) == 6)));
        plot(Distances(6), nanstd(FFGroupData{i}{SyllIndex}(UndirIndices, 1))/nanmean(FFGroupData{i}{SyllIndex}(UndirIndices, 1)), 'ks', 'MarkerSize', 6, 'MarkerFaceColor', 'k');
%         % Now best fit line
%         RobustCoeffs = robustfit(Distances(FFCVGroupData{i}{SyllIndex}(DirIndices,end-1)), FFCVGroupData{i}{SyllIndex}(DirIndices, 1));
%         plot(unique(Distances(FFCVGroupData{i}{SyllIndex}(DirIndices,end-1))), polyval(flipud(RobustCoeffs), unique(Distances(FFCVGroupData{i}{SyllIndex}(DirIndices,end-1)))), 'k', 'LineWidth', 1.5);

%         [FixedX_BootStrapCoeffs, Random_BootStrapCoeffs, RandomProb, Rsq] = BootstrapTestRegressionSignificance(RobustCoeffs, Distances(FFCVGroupData{i}{SyllIndex}(DirIndices,end-1))', FFCVGroupData{i}{SyllIndex}(DirIndices, 1), 'linear');
        axis tight;
        Temp = axis;
        Temp = [-5 Temp(2)+5 0.99*Temp(3) 1.2*Temp(4)];
        axis(Temp);
%         text(Temp(1)+2, Temp(4)/1.15, {['Rsq = ', num2str(Rsq)]; ['p = ', num2str(RandomProb(end))]});
%         text(Temp(1)+2, Temp(4)/1.05, {['Slope = ', num2str(FixedX_BootStrapCoeffs(1,end)), ' : ', num2str(FixedX_BootStrapCoeffs(2,end)), ' : ', num2str(FixedX_BootStrapCoeffs(3,end))]});
        ylabel('Syllable FF CV')
        xlabel('Distance (cm)');
        title([BirdNames{i}, 'Fit to feature value for syllable ', UndirBout_Stats{i}{end}.FF_SyllLabel(SyllIndex)]);

%         % Now calculate residuals
%         YPred = (Distances(FFCVGroupData{i}{SyllIndex}(DirIndices,end-1)) * RobustCoeffs(2)) + RobustCoeffs(1);
%         Error = FFCVGroupData{i}{SyllIndex}(DirIndices, 1) - YPred(:);
% 
%         % Now plot residuals as a function of x-value with mean and
%         % errorbars
%         subplot(1,3,2);
%         hold on;
%         plot(Distances(FFCVGroupData{i}{SyllIndex}(DirIndices, end-1)), Error, 'ko');
%         for Location = 1:max(FFCVGroupData{i}{SyllIndex}(DirIndices,end-1)),
%             LocIndices = intersect(DirIndices, find(FFCVGroupData{i}{SyllIndex}(:,end-1) == Location));
%             if (~isempty(LocIndices))
%                 if (length(LocIndices) == 1)
%                     plot(Distances(Location)+3, Error(LocIndices), 'ks', 'MarkerSize', 6, 'MarkerFaceColor', 'k');
%                 else
%                     errorbar(Distances(Location)+3, mean(Error(LocIndices)), std(Error(LocIndices))/sqrt(length(LocIndices)), 'ks', 'MarkerSize', 6, 'MarkerFaceColor', 'k');
%                 end
%             end
%         end
%         axis tight;
%         Temp = axis;
%         Temp = [-5 Temp(2)+5 1.02*Temp(3) 1.02*Temp(4)];
%         axis(Temp);
%         plot(Temp(1:2), zeros(1,2), 'k--');
%         title('Distance vs. Residuals');
%         xlabel('Distance (cm)');
%         ylabel('Residual');
% 
%         subplot(1,3,3);
%         hold on;
%         Edges = linspace(min(Error)*0.95, max(Error)*1.05, 20);
%         plot(Edges, histc(Error, Edges), 'ko-');
%         title('Residual distribution');
%         ylabel('#');
%         xlabel('Residual');

    end
end

disp('Done with stats and plots');
