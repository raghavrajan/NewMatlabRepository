function [L0Dir_Undir_PValue] = Harini_FeatureFits_INNum_MotifNumDur_Time(DirBout_Stats, DirUndirBout_Stats, UndirBout_Stats, Distances, BirdNames)

DistanceCode = {'L0' 'L1' 'L2' 'L3' 'L4' 'UN'};
FinalFigureDir = '/home/raghav/StudentRelated/Harini/Manuscript/results';
YLabelString = {'# of INs' '# of motifs/bout' 'Motif duration (msec)'};

% We want to do the analysis per bird and for the group as a whole. The
% things that we will analyze are (1) # of INs at the beginning based on
% the 500ms criterion, (2) # of complete motifs in the bout, (3) all and
% first motif duration, (4), CV of FF and (5) amplitude of individual
% motif syllables only for HF mic birds

for i = 1:length(DirBout_Stats),
    figure(i);
    FigPanel(i).p = panel();
    FigPanel(i).p.pack('h', {1/3 1/3 1/3});
end

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
    
    L0DirIndices = find((GroupData{i}(:,end-1) == 1) & (GroupData{i}(:,end) == 1));
    DirIndices = find(GroupData{i}(:,end) == 1);
    UnDirIndices = find((GroupData{i}(:,end-1) == 6) & (GroupData{i}(:,end) == 3));
    
    for Feat = 1:length(FeaturesToConsider),
        FigPanel(i).p(Feat).select();
        hold on;
    
        % Now plot the mean, SEM and the best fit line with robust regression
        for Location = 1:max(GroupData{i}(DirIndices,end-1)),
            LocIndices = intersect(DirIndices, find(GroupData{i}(:,end-1) == Location));
            if (~isempty(LocIndices))
                errorbar(Distances(Location), mean(GroupData{i}(LocIndices, Feat)), std(GroupData{i}(LocIndices, Feat))/sqrt(length(LocIndices)), 'ks', 'MarkerSize', 6, 'MarkerFaceColor', 'k');
            end
        end
        errorbar(Distances(Location)+30, mean(GroupData{i}(UnDirIndices, Feat)), std(GroupData{i}(UnDirIndices, Feat))/sqrt(length(UnDirIndices)), 'ks', 'MarkerSize', 6, 'MarkerFaceColor', 'k');
        
        % Now best fit line
        RobustCoeffs = robustfit(Distances(GroupData{i}(DirIndices,end-1)), GroupData{i}(DirIndices, Feat));
        [FixedX_BootStrapCoeffs, Random_BootStrapCoeffs, RandomProb, Rsq(Feat,i)] = BootstrapTestRegressionSignificance(RobustCoeffs, Distances(GroupData{i}(DirIndices,end-1))', GroupData{i}(DirIndices, Feat), 'linear');
        
        FitCoeffs{Feat}(i,:) = RobustCoeffs;
        BootStrapCI{Feat}{i} = FixedX_BootStrapCoeffs;
    
        if (RobustCoeffs(2) < 0)
            if ((FixedX_BootStrapCoeffs(1,2) < 0) && (FixedX_BootStrapCoeffs(3,2) < 0))
                PlotLineColor = 'b';
                SlopeSignificance(Feat,i) = -1;
            else
                PlotLineColor = 'k';
                SlopeSignificance(Feat,i) = 0;
            end
        else
            if (RobustCoeffs(2) > 0)
                if ((FixedX_BootStrapCoeffs(1,2) > 0) && (FixedX_BootStrapCoeffs(3,2) > 0))
                    PlotLineColor = 'r';
                    SlopeSignificance(Feat,i) = 1;
                else
                    PlotLineColor = 'k';
                    SlopeSignificance(Feat,i) = 0;
                end
            else
                PlotLineColor = 'k';
                SlopeSignificance(Feat,i) = 0;
            end
        end
        plot(unique(Distances(GroupData{i}(DirIndices,end-1))), polyval(flipud(RobustCoeffs), unique(Distances(GroupData{i}(DirIndices,end-1)))), PlotLineColor, 'LineWidth', 1.5);
        patch([unique(Distances(GroupData{i}(DirIndices,end-1))) fliplr(unique(Distances(GroupData{i}(DirIndices,end-1))))], [polyval(fliplr(FixedX_BootStrapCoeffs(1,:)), unique(Distances(GroupData{i}(DirIndices,end-1)))) fliplr(polyval(fliplr(FixedX_BootStrapCoeffs(3,:)), unique(Distances(GroupData{i}(DirIndices,end-1)))))], PlotLineColor, 'EdgeColor', 'none', 'FaceAlpha', 0.2);
        
        
        % [Rsq, F, Prob] = CalculateGoodnessofLinearFit(flipud(RobustCoeffs), Distances(GroupData{i}(DirIndices,4))', GroupData{i}(DirIndices, Feat));
        axis tight;
        Temp = axis;
        Temp = [-5 Temp(2)+5 0.99*Temp(3) 1.2*Temp(4)];
        axis(Temp);
        plot([-5 Temp(2)+5], ones(1,2)*RobustCoeffs(1), 'k--');
        text(Temp(2)-23, Temp(3), '//', 'FontSize', 16)
        %text(Temp(1)+2, Temp(4)/1.15, {['Rsq = ', num2str(Rsq)]; ['p = ', num2str(RandomProb(end))]});
        %text(Temp(1)+2, Temp(4)/1.05, {['Slope = ', num2str(FixedX_BootStrapCoeffs(1,end)), ' : ', num2str(FixedX_BootStrapCoeffs(2,end)), ' : ', num2str(FixedX_BootStrapCoeffs(3,end))]});
        ylabel(YLabelString{Feat})
        if (Feat == length(FeaturesToConsider))
            xlabel('Distance (cm)');
        end
        if (Feat == 2)
            title(BirdNames{i});
        end
        L0Dir_UnDir_Means{Feat}(i,:) = [mean(GroupData{i}(L0DirIndices, Feat)) mean(GroupData{i}(UnDirIndices, Feat))];
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
    
    L0DirIndices = find((AllMotifDur_GroupData{i}(:,end-1) == 1) & (AllMotifDur_GroupData{i}(:,end) == 1));
    DirIndices = find(AllMotifDur_GroupData{i}(:,end) == 1);
    UnDirIndices = find((AllMotifDur_GroupData{i}(:,end-1) == 6) & (AllMotifDur_GroupData{i}(:,end) == 3));
    % First plot the group data
    FigPanel(i).p(length(FeaturesToConsider)+1).select();
    hold on;

    % Now plot the mean, SEM and the best fit line with robust regression
    for Location = 1:max(AllMotifDur_GroupData{i}(DirIndices,end-1)),
        LocIndices = intersect(DirIndices, find(AllMotifDur_GroupData{i}(:,end-1) == Location));
        if (~isempty(LocIndices))
            errorbar(Distances(Location), mean(AllMotifDur_GroupData{i}(LocIndices, 1)), std(AllMotifDur_GroupData{i}(LocIndices, 1))/sqrt(length(LocIndices)), 'ks', 'MarkerSize', 6, 'MarkerFaceColor', 'k');
        end
    end

    errorbar(Distances(Location)+30, mean(AllMotifDur_GroupData{i}(UnDirIndices, 1)), std(AllMotifDur_GroupData{i}(UnDirIndices, 1))/sqrt(length(UnDirIndices)), 'ks', 'MarkerSize', 6, 'MarkerFaceColor', 'k');
    % Now best fit line
    RobustCoeffs = robustfit(Distances(AllMotifDur_GroupData{i}(DirIndices,end-1)), AllMotifDur_GroupData{i}(DirIndices, 1));
    [FixedX_BootStrapCoeffs, Random_BootStrapCoeffs, RandomProb, Rsq(3,i)] = BootstrapTestRegressionSignificance(RobustCoeffs, Distances(AllMotifDur_GroupData{i}(DirIndices,end-1))', AllMotifDur_GroupData{i}(DirIndices, 1), 'linear');
    
    FitCoeffs{3}(i,:) = RobustCoeffs;
    BootStrapCI{3}{i} = FixedX_BootStrapCoeffs;
    
    if (RobustCoeffs(2) < 0)
       if ((FixedX_BootStrapCoeffs(1,2) < 0) && (FixedX_BootStrapCoeffs(3,2) < 0))
            PlotLineColor = 'b';
            SlopeSignificance(3,i) = -1; % for decrease
        else
            PlotLineColor = 'k';
            SlopeSignificance(3,i) = 0;
        end
    else
        if (RobustCoeffs(2) > 0)
            if ((FixedX_BootStrapCoeffs(1,2) > 0) && (FixedX_BootStrapCoeffs(3,2) > 0))
                PlotLineColor = 'r';
                SlopeSignificance(3,i) = 1; % for increase
            else
                PlotLineColor = 'k';
                SlopeSignificance(3,i) = 0;
            end
        else
            PlotLineColor = 'k';
            SlopeSignificance(3,i) = 0;
        end
    end
    
    plot(unique(Distances(AllMotifDur_GroupData{i}(DirIndices,end-1))), polyval(flipud(RobustCoeffs), unique(Distances(AllMotifDur_GroupData{i}(DirIndices,end-1)))), PlotLineColor, 'LineWidth', 1.5);
    patch([unique(Distances(AllMotifDur_GroupData{i}(DirIndices,end-1))) fliplr(unique(Distances(AllMotifDur_GroupData{i}(DirIndices,end-1))))], [polyval(fliplr(FixedX_BootStrapCoeffs(1,:)), unique(Distances(AllMotifDur_GroupData{i}(DirIndices,end-1)))) fliplr(polyval(fliplr(FixedX_BootStrapCoeffs(3,:)), unique(Distances(AllMotifDur_GroupData{i}(DirIndices,end-1)))))], PlotLineColor, 'EdgeColor', 'none', 'FaceAlpha', 0.2);
        
    axis tight;
    Temp = axis;
    Temp = [-5 Temp(2)+5 0.99*Temp(3) 1.01*Temp(4)];
    axis(Temp);
    plot([-5 Temp(2)+5], ones(1,2)*RobustCoeffs(1), 'k--');
    text(Temp(2)-23, Temp(3), '//', 'FontSize', 16)
    %text(Temp(1)+2, Temp(4)/1.15, {['Rsq = ', num2str(Rsq)]; ['p = ', num2str(RandomProb(end))]});
    %text(Temp(1)+2, Temp(4)/1.05, {['Slope = ', num2str(FixedX_BootStrapCoeffs(1,end)), ' : ', num2str(FixedX_BootStrapCoeffs(2,end)), ' : ', num2str(FixedX_BootStrapCoeffs(3,end))]});
    ylabel(YLabelString{3});
    if (Feat == length(FeaturesToConsider))
        xlabel('Distance (cm)');
    end
    L0Dir_UnDir_Means{3}(i,:) = [mean(AllMotifDur_GroupData{i}(L0DirIndices, 1)) mean(AllMotifDur_GroupData{i}(UnDirIndices, 1))];
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

for i = 1:length(DirBout_Stats),
    FigPanel(i).p.fontsize = 12;
    FigPanel(i).p.marginleft = 20;
    FigPanel(i).p.de.margin = 20;

    figure(i);
    set(gcf, 'Color', 'w');
    set(gcf, 'Units', 'inches');
    set(gcf, 'Position', [1 1 8.3 3.5]);
    set(gcf, 'PaperPositionMode', 'auto');
    print(fullfile(FinalFigureDir, [BirdNames{i}, 'FeatureFits_INNum_MotifNum_MotifDur_Time.eps']), '-depsc2', '-r600');
    print(fullfile(FinalFigureDir, [BirdNames{i}, 'FeatureFits_INNum_MotifNum_MotifDur_Time.png']), '-dpng', '-r600');
end

close all;
BirdsToPlot = {'brwn49org68' 'org05blue12'};

for BirdIndex = 1:length(BirdsToPlot),
    i = strmatch(BirdsToPlot{BirdIndex}, BirdNames, 'exact');
    BirdsToPlotIndices(BirdIndex) = i;
end

figure;
p = panel();
p.pack('h', {1/3 1/3 1/3});

Symbols = 'sod^v><+hx*p';
% Now to plot the group data
for i = 1:length(FeaturesToConsider)+1,
    if (i > length(FeaturesToConsider))
        DirectedDistanceMeans = AllMotifDur_DirDistanceMeans(:,1:5)./repmat(AllMotifDur_DirDistanceMeans(:,1), 1, 5);
        UndirectedDistanceMeans = AllMotifDur_UndirMeans(:,end)./(AllMotifDur_DirDistanceMeans(:,1));
    else
        DirectedDistanceMeans = DirDistanceMeans{i}(:,1:5)./repmat(DirDistanceMeans{i}(:,1), 1, 5);
        UndirectedDistanceMeans = UndirMeans{i}(:,end)./(DirDistanceMeans{i}(:,1));
    end
    
    p(i).select();
    hold on;
    for j = 1:size(DirectedDistanceMeans,1),
        if (~isempty(find(BirdsToPlotIndices == j)))
            plot(Distances(1:5), DirectedDistanceMeans(j,:), ['k-', Symbols(j)], 'Color', [0.25 0.25 0.25]);
            plot(Distances(5) + 25, UndirectedDistanceMeans(j), ['k-', Symbols(j)], 'Color', [0.25 0.25 0.25]);
        else
            plot(Distances(1:5), DirectedDistanceMeans(j,:), ['k-', Symbols(j)], 'Color', [0.75 0.75 0.75]);
            plot(Distances(5) + 25, UndirectedDistanceMeans(j), ['k-', Symbols(j)], 'Color', [0.75 0.75 0.75]);
        end
    end
    
    for j = 1:size(DirectedDistanceMeans,2),
        NonNaNValues = find(~isnan(DirectedDistanceMeans(:,j)));
        errorbar(Distances(j), mean(DirectedDistanceMeans(NonNaNValues,j)), std(DirectedDistanceMeans(NonNaNValues,j))/sqrt(length(NonNaNValues)), 'ks-', 'MarkerSize', 6, 'MarkerFaceColor', 'k');
    end
    
    errorbar(Distances(5)+25, mean(UndirectedDistanceMeans), std(UndirectedDistanceMeans)/sqrt(length(UndirectedDistanceMeans)), 'ks-', 'MarkerSize', 6, 'MarkerFaceColor', 'k');
    [GroupDataRobustCoeffs, RandomPosPValue, SamePosPValue] = MissingValuesPermutationTest(repmat(Distances(1:5), size(DirectedDistanceMeans,1),1), DirectedDistanceMeans);
    if (SamePosPValue(2) < 0.05)
        if (GroupDataRobustCoeffs(2) < 0)
            plot(Distances(1:5), polyval(flipud(GroupDataRobustCoeffs(:)), Distances(1:5)), 'b', 'LineWidth', 2);
        else
            plot(Distances(1:5), polyval(flipud(GroupDataRobustCoeffs(:)), Distances(1:5)), 'r', 'LineWidth', 2);
        end
    else
        plot(Distances(1:5), polyval(flipud(GroupDataRobustCoeffs(:)), Distances(1:5)), 'k', 'LineWidth', 2);
    end
    switch i
        case 1
            ylabel('Normalized # of INs');
        case 2
            ylabel('Normalized # of motifs/bout');
        case 3
            ylabel('Normalized motif duration');
    end
    
    if (i == 2)
        xlabel('Distance from female (cm)');
    end
    DistanceXLabelString = [];
    DistanceXValues = [];
    for Dist = 1:length(Distances)-1,
        DistanceXValues(Dist) = Distances(Dist);
        DistanceXLabelString{Dist} = num2str(Distances(Dist));
    end
    DistanceXValues(end+1) = DistanceXValues(end) + 25;
    DistanceXLabelString{end+1} = 'UN';
    set(gca, 'XTick', DistanceXValues, 'XTickLabel', DistanceXLabelString, 'XTickLabelRotation', 45);
    axis tight;
    Temp = axis;
    Temp = [-5 Temp(2)+5 0.98*Temp(3) 1.02*Temp(4)];
    axis(Temp)
    plot(Temp(1:2), ones(1,2), 'k--');
    text(110, Temp(4), ['p = ', num2str(round(SamePosPValue(2)*100)/100)]);
    text(170, Temp(3), '//', 'FontSize', 12)
end
   
p.fontsize = 12;
p.de.margin = 20;
p.marginleft = 25;
p.margintop = 10;
p.marginbottom = 30;
p.marginright = 10;

set(gcf, 'Color', 'w');
set(gcf, 'Units', 'inches');
set(gcf, 'Position', [10 1.4 8.3 3.1]);
set(gcf, 'PaperPositionMode', 'auto');
print(fullfile(FinalFigureDir, ['AllBirds_Summary_INNum_MotifNum_MotifDur_Time.eps']), '-depsc2', '-r600');
print(fullfile(FinalFigureDir, ['AllBirds_Summary_INNum_MotifNum_MotifDur_Time.png']), '-dpng', '-r600');

close all;
figure;
p = panel();
p.pack('h', {1/3 1/3 1/3});

% Now to plot the slopes contingent on whether mean during L0 Dir >
% or < Undir
for Feat = 1:length(L0Dir_UnDir_Means),
    p(Feat).select();
    hold on;
    L0MoreThanUndirIndices = find(L0Dir_UnDir_Means{Feat}(:,1) > L0Dir_UnDir_Means{Feat}(:,2));
    SlopeSigBar(Feat,1) = bar(1, mean(FitCoeffs{Feat}(L0MoreThanUndirIndices, 2)));
    set(SlopeSigBar(Feat,1), 'FaceColor', 'none');
    errorbar(1, mean(FitCoeffs{Feat}(L0MoreThanUndirIndices, 2)), std(FitCoeffs{Feat}(L0MoreThanUndirIndices, 2))/sqrt(length(L0MoreThanUndirIndices)), 'k', 'MarkerFaceColor', 'k');
    for i = L0MoreThanUndirIndices(:)',
        switch SlopeSignificance(Feat,i)
            case -1
                MarkerColorSymbol = 'b';
                MarkerSymbolFaceColor = 'b';
                
            case 0
                MarkerColorSymbol = 'k';
                MarkerSymbolFaceColor = 'none';
                
            case 1
                MarkerColorSymbol = 'r';
                MarkerSymbolFaceColor = 'r';
        end
        plot(1.15, FitCoeffs{Feat}(i,2), [MarkerColorSymbol, 'o'], 'MarkerSize', 4, 'MarkerFaceColor', MarkerSymbolFaceColor);
    end
    L0LessThanUndirIndices = find(L0Dir_UnDir_Means{Feat}(:,1) < L0Dir_UnDir_Means{Feat}(:,2));
    SlopeSigBar(Feat,1) = bar(2, mean(FitCoeffs{Feat}(L0LessThanUndirIndices, 2)));
    set(SlopeSigBar(Feat,1), 'FaceColor', 'none');
    errorbar(2, mean(FitCoeffs{Feat}(L0LessThanUndirIndices, 2)), std(FitCoeffs{Feat}(L0LessThanUndirIndices, 2))/sqrt(length(L0LessThanUndirIndices)), 'k', 'MarkerFaceColor', 'k');
    for i = L0LessThanUndirIndices(:)',
        switch SlopeSignificance(Feat,i)
            case -1
                MarkerColorSymbol = 'b';
                MarkerSymbolFaceColor = 'b';
                
            case 0
                MarkerColorSymbol = 'k';
                MarkerSymbolFaceColor = 'none';
                
            case 1
                MarkerColorSymbol = 'r';
                MarkerSymbolFaceColor = 'r';
        end
        plot(1.85, FitCoeffs{Feat}(i,2), [MarkerColorSymbol, 'o'], 'MarkerSize', 4, 'MarkerFaceColor', MarkerSymbolFaceColor);
    end
    axis tight;
    Temp = axis;
    Temp = [0.5 2.5 1.2*Temp(3) 1.2*Temp(4)];
    axis(Temp);
    % plot(Temp(1:2), zeros(2,1), 'k--');
    set(gca, 'XTick', [1 2], 'XTickLabel', {'L0 DIR > UNDIR' 'L0 DIR < UNDIR'}, 'XTickLabelRotation', 45);
    if (Feat == 1)
        ylabel('Slope');
    end
    switch (Feat)
        case 1
            title('# of INs');
        case 2
            title('# of motifs/bout');
        case 3
            title('Motif duration');
    end
end

p.fontsize = 12;
p.de.margin = 20;
p.marginleft = 25;
p.margintop = 10;
p.marginbottom = 30;
p.marginright = 10;

set(gcf, 'Color', 'w');
set(gcf, 'Units', 'inches');
set(gcf, 'Position', [10 1.4 8.3 3.1]);
set(gcf, 'PaperPositionMode', 'auto');
print(fullfile(FinalFigureDir, ['AllBirds_Slope_Summary_INNum_MotifNum_MotifDur_Time.eps']), '-depsc2', '-r600');
print(fullfile(FinalFigureDir, ['AllBirds_Slope_Summary_INNum_MotifNum_MotifDur_Time.png']), '-dpng', '-r600');
disp('Done with stats and plots');
