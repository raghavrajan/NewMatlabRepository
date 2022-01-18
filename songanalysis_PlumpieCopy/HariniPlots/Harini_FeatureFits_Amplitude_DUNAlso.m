function [L0Dir_Undir_PValue] = Harini_FeatureFits_Amplitude(DirBout_Stats, DirUndirBout_Stats, UndirBout_Stats, Distances, BirdNames)

FinalFigureDir = '/home/raghav/StudentRelated/Harini/Manuscript/results';
DistanceCode = {'L0' 'L1' 'L2' 'L3' 'L4' 'UN'};

% We want to do the analysis per bird and for the group as a whole. The
% things that we will analyze are (1) # of INs at the beginning based on
% the 500ms criterion, (2) # of complete motifs in the bout, (3) all and
% first motif duration, (4), CV of FF and (5) amplitude of individual
% motif syllables only for HF mic birds

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
    figure(i);
    FigPanel(i).p = panel();
    FigPanel(i).p.pack(1, length(UndirBout_Stats{i}{end}.LogAmplitude_SyllLabel));
    
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
            disp(['Removed ', num2str(length(OutlierIndices)), '(', num2str(100*length(OutlierIndices)/length(SyllDurationGroupData{i}{SyllIndex})), '%) from syllable ', UndirBout_Stats{i}{end}.LogAmplitude_SyllLabel(SyllIndex), ' based on syll duration:', BirdNames{i}]);
            SyllDurationGroupData{i}{SyllIndex}(OutlierIndices, :) = [];
            LogAmplitudeGroupData{i}{SyllIndex}(OutlierIndices, :) = [];
        end
    
        % Remove outliers based on the criterion that data is either 3*IQR 
        % above 75th percentile or below 25th percentile for syllable
        % duration and then for amplitude
        OutlierThreshold = [(prctile(LogAmplitudeGroupData{i}{SyllIndex}(:,1), 25) - 3*iqr(LogAmplitudeGroupData{i}{SyllIndex}(:,1))) (prctile(LogAmplitudeGroupData{i}{SyllIndex}(:,1), 75) + 3*iqr(LogAmplitudeGroupData{i}{SyllIndex}(:,1)))];
        OutlierIndices = find((LogAmplitudeGroupData{i}{SyllIndex}(:,1) < OutlierThreshold(1)) | (LogAmplitudeGroupData{i}{SyllIndex}(:,1) > OutlierThreshold(2)));
        if (~isempty(OutlierIndices))
            disp(['Removed ', num2str(length(OutlierIndices)), '(', num2str(100*length(OutlierIndices)/length(SyllDurationGroupData{i}{SyllIndex})), '%) from syllable ', UndirBout_Stats{i}{end}.LogAmplitude_SyllLabel(SyllIndex), ' based on syll amplitude:', BirdNames{i}]);
            SyllDurationGroupData{i}{SyllIndex}(OutlierIndices, :) = [];
            LogAmplitudeGroupData{i}{SyllIndex}(OutlierIndices, :) = [];
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
    
        DirIndices = find(LogAmplitudeGroupData{i}{SyllIndex}(:,end) == 1);
        UnDirIndices = find((LogAmplitudeGroupData{i}{SyllIndex}(:,end-1) == 6) & (LogAmplitudeGroupData{i}{SyllIndex}(:,end) == 3));
        % First plot the group data
        FigPanel(i).p(1,SyllIndex).select();
        hold on;

        % Now plot the mean, SEM and the best fit line with robust regression
        for Location = 1:max(LogAmplitudeGroupData{i}{SyllIndex}(DirIndices,end-1)),
            LocIndices = intersect(DirIndices, find(LogAmplitudeGroupData{i}{SyllIndex}(:,end-1) == Location));
            if (~isempty(LocIndices))
                errorbar(Distances(Location), mean(LogAmplitudeGroupData{i}{SyllIndex}(LocIndices, 1)), std(LogAmplitudeGroupData{i}{SyllIndex}(LocIndices, 1))/sqrt(length(LocIndices)), 'ks', 'MarkerSize', 6, 'MarkerFaceColor', 'k');
            end
        end
        errorbar(Distances(Location)+30, mean(LogAmplitudeGroupData{i}{SyllIndex}(UnDirIndices, 1)), std(LogAmplitudeGroupData{i}{SyllIndex}(UnDirIndices, 1))/sqrt(length(UnDirIndices)), 'ks', 'MarkerSize', 6, 'MarkerFaceColor', 'k');
        
        % Now best fit line
        RobustCoeffs = robustfit(Distances(LogAmplitudeGroupData{i}{SyllIndex}(DirIndices,end-1)), LogAmplitudeGroupData{i}{SyllIndex}(DirIndices, 1));
        FitCoeffs{i}{SyllIndex} = RobustCoeffs;
        
        [FixedX_BootStrapCoeffs, Random_BootStrapCoeffs, RandomProb, Rsq] = BootstrapTestRegressionSignificance(RobustCoeffs, Distances(LogAmplitudeGroupData{i}{SyllIndex}(DirIndices,end-1))', LogAmplitudeGroupData{i}{SyllIndex}(DirIndices, 1), 'linear');
        BootStrapCI{i}{SyllIndex} = FixedX_BootStrapCoeffs;
        
        if (RobustCoeffs(2) < 0)
            if ((FixedX_BootStrapCoeffs(1,2) < 0) && (FixedX_BootStrapCoeffs(3,2) < 0))
                PlotLineColor = 'b';
                SlopeSignificance{i}(SyllIndex) = -1;
            else
                PlotLineColor = 'k';
                SlopeSignificance{i}(SyllIndex) = 0;
            end
        else
            if (RobustCoeffs(2) > 0)
                if ((FixedX_BootStrapCoeffs(1,2) > 0) && (FixedX_BootStrapCoeffs(3,2) > 0))
                    PlotLineColor = 'r';
                    SlopeSignificance{i}(SyllIndex) = 1;
                else
                    PlotLineColor = 'k';
                    SlopeSignificance{i}(SyllIndex) = 0;
                end
            else
                PlotLineColor = 'k';
                SlopeSignificance{i}(SyllIndex) = 0;
            end
        end
        plot(unique(Distances(LogAmplitudeGroupData{i}{SyllIndex}(DirIndices,end-1))), polyval(flipud(RobustCoeffs), unique(Distances(LogAmplitudeGroupData{i}{SyllIndex}(DirIndices,end-1)))), PlotLineColor, 'LineWidth', 1.5);
        patch([unique(Distances(LogAmplitudeGroupData{i}{SyllIndex}(DirIndices,end-1))) fliplr(unique(Distances(LogAmplitudeGroupData{i}{SyllIndex}(DirIndices,end-1))))], [polyval(fliplr(FixedX_BootStrapCoeffs(1,:)), unique(Distances(LogAmplitudeGroupData{i}{SyllIndex}(DirIndices,end-1)))) fliplr(polyval(fliplr(FixedX_BootStrapCoeffs(3,:)), unique(Distances(LogAmplitudeGroupData{i}{SyllIndex}(DirIndices,end-1)))))], PlotLineColor, 'EdgeColor', 'none', 'FaceAlpha', 0.2);
        
        axis tight;
        Temp = axis;
        Temp = [-5 Temp(2)+5 0.99*Temp(3) 1.01*Temp(4)];
        text(Temp(2)-23, Temp(3), '//', 'FontSize', 16);
        plot(Temp(1:2), ones(1,2)*RobustCoeffs(1), 'k--');
        axis(Temp);
        ylabel('Syllable ampitude (au)')
        xlabel('Distance (cm)');

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
    print(fullfile(FinalFigureDir, [BirdNames{i}, 'FeatureFits_Amplitude.eps']), '-depsc2', '-r600');
    print(fullfile(FinalFigureDir, [BirdNames{i}, 'FeatureFits_Amplitude.png']), '-dpng', '-r600');
end

BirdDistanceDirMeans = ones(length(LogAmplitudeGroupData), length(Distances))*NaN;
BirdUnDirMeans = ones(length(LogAmplitudeGroupData), 1)*NaN;

SyllableDistanceDirMeans = ones(sum(cellfun(@length, LogAmplitudeGroupData)), length(Distances))*NaN;
SyllableUnDirMeans = ones(sum(cellfun(@length, LogAmplitudeGroupData)), 1)*NaN;
SyllableBirdIndices = ones(sum(cellfun(@length, LogAmplitudeGroupData)), 1)*NaN;

UnNormalizedSyllableDistanceDirMeans = ones(sum(cellfun(@length, LogAmplitudeGroupData)), length(Distances))*NaN;
UnNormalizedSyllableUnDirMeans = ones(sum(cellfun(@length, LogAmplitudeGroupData)), 1)*NaN;

DirL0_UnDir_Stats = ones(sum(cellfun(@length, LogAmplitudeGroupData)), 1)*NaN;

SyllIndex = 0;
for i = 1:length(LogAmplitudeGroupData),
    TempDistanceDirMeans = ones(length(LogAmplitudeGroupData{i}), length(Distances))*NaN;
    TempDistanceUnDirMeans = ones(length(LogAmplitudeGroupData{i}),1)*NaN;
    TempDirL0_UnDir_Stats = ones(length(LogAmplitudeGroupData{i}),1)*NaN;
    for j = 1:length(LogAmplitudeGroupData{i}),
        DirIndices = find(LogAmplitudeGroupData{i}{j}(:,end) == 1);
        DirL0Indices = find((LogAmplitudeGroupData{i}{j}(:,end) == 1) & (LogAmplitudeGroupData{i}{j}(:,end-1) == 1));
        for k = 1:max(LogAmplitudeGroupData{i}{j}(DirIndices,end-1)),
            DirLocIndices = find((LogAmplitudeGroupData{i}{j}(:,end) == 1) & (LogAmplitudeGroupData{i}{j}(:,end-1) == k));
            TempDistanceDirMeans(j,k) = mean(LogAmplitudeGroupData{i}{j}(DirLocIndices,1));
        end
        UnDirIndices = find((LogAmplitudeGroupData{i}{j}(:,end) == 3) & (LogAmplitudeGroupData{i}{j}(:,end-1) == 6));
        [TempHypothesis, TempDirL0_UnDir_Stats(j)] = ttest2(LogAmplitudeGroupData{i}{j}(DirL0Indices,1), LogAmplitudeGroupData{i}{j}(UnDirIndices));
        TempDistanceUnDirMeans(j) = mean(LogAmplitudeGroupData{i}{j}(UnDirIndices,1));
    end
    DirL0_UnDir_Stats(SyllIndex+1:SyllIndex+length(TempDirL0_UnDir_Stats)) = TempDirL0_UnDir_Stats;
    
    UnNormalizedSyllableDistanceDirMeans(SyllIndex+1:SyllIndex+size(TempDistanceDirMeans,1), :) = TempDistanceDirMeans;
    UnNormalizedSyllableUnDirMeans(SyllIndex+1:SyllIndex+size(TempDistanceDirMeans,1), 1) = TempDistanceUnDirMeans;
    
    % Normalize each syllable mean to it's directed L0 value
    TempDistanceUnDirMeans = TempDistanceUnDirMeans./TempDistanceDirMeans(:,1);
    TempDistanceDirMeans = TempDistanceDirMeans./(repmat(TempDistanceDirMeans(:,1), 1, length(Distances)));
    
    SyllableDistanceDirMeans(SyllIndex+1:SyllIndex+size(TempDistanceDirMeans,1), :) = TempDistanceDirMeans;
    SyllableUnDirMeans(SyllIndex+1:SyllIndex+size(TempDistanceDirMeans,1), 1) = TempDistanceUnDirMeans;
    SyllableBirdIndices(SyllIndex+1:SyllIndex+size(TempDistanceDirMeans,1), 1) = i;
    SyllIndex = SyllIndex + size(TempDistanceDirMeans, 1);
    
    if (length(LogAmplitudeGroupData{i}) > 1)
        BirdDistanceDirMeans(i,:) = mean(TempDistanceDirMeans);
    else
        BirdDistanceDirMeans(i,:) = TempDistanceDirMeans;
    end
    BirdUnDirMeans(i) = mean(TempDistanceUnDirMeans);
end

% Plot with data averaged within each bird across syllables
BirdsToPlot = {'brwn89red54' 'brwn99red65'};
BirdTitles = {'Bird #12' 'Bird #11'};

for BirdIndex = 1:length(BirdsToPlot),
    i = strmatch(BirdsToPlot{BirdIndex}, BirdNames, 'exact');
    BirdsToPlotIndices(BirdIndex) = i;
end

close all;
figure;
p = panel();
p.pack('h', {1/3 1/3 1/3});

Symbols = 'sod^v><+hx*p';
p(1).select();
hold on;
UniqueBirdIndices = unique(SyllableBirdIndices);
for i = UniqueBirdIndices(:)',
    Syllables = find(SyllableBirdIndices == i)
    for j = Syllables(:)',
        if (DirL0_UnDir_Stats(j) < 0.05)
            if (UnNormalizedSyllableDistanceDirMeans(j,1) > UnNormalizedSyllableUnDirMeans(j))
                MarkerSymbolFaceColor = 'b';
                SymbolColor = 'b';
            else
                MarkerSymbolFaceColor = 'r';
                SymbolColor = 'r';
            end
        else
            MarkerSymbolFaceColor = 'none';
            SymbolColor = 'k';
        end
        plot(i, UnNormalizedSyllableDistanceDirMeans(j,1)/UnNormalizedSyllableUnDirMeans(j), [SymbolColor, Symbols(i)], 'MarkerFaceColor', MarkerSymbolFaceColor);
    end
end
axis tight;
Temp = axis;
Temp = [0.75 max(UniqueBirdIndices)+0.25 0.96*Temp(3) 1.02*Temp(4)];
axis(Temp);
plot(Temp(1:2), [1 1], 'k--');
ylabel('L0 songs / Undirected');
BirdFinalIndices = [10 11 12 13 5];

for i = 1:length(UniqueBirdIndices),
    BirdIndexXLabelString{i} = ['Bird #', num2str(BirdFinalIndices(i))];
end
set(gca, 'XTick', 1:1:max(UniqueBirdIndices), 'XTickLabel', BirdIndexXLabelString, 'XTickLabelRotation', 45);

% Now to plot the slopes contingent on whether mean during L0 Dir >
% or < Undir
p(2).select();
hold on;

L0MoreThanUndirIndices = intersect(find(~isnan(BirdUnDirMeans)), find(BirdUnDirMeans < 1));
L0MoreSlopes = [];
for i = L0MoreThanUndirIndices(:)',
    for j = 1:length(FitCoeffs{i}),
        L0MoreSlopes(end+1) = FitCoeffs{i}{j}(2);
        switch (SlopeSignificance{i}(j))
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
        plot(1.15, FitCoeffs{i}{j}(2), [MarkerColorSymbol, Symbols(i)], 'MarkerSize', 4, 'MarkerFaceColor', MarkerSymbolFaceColor);    
    end
end

SlopeSigBar(1) = bar(1, mean(L0MoreSlopes));
set(SlopeSigBar(1), 'FaceColor', 'none');
errorbar(1, mean(L0MoreSlopes), std(L0MoreSlopes)/sqrt(length(L0MoreSlopes)), 'k', 'MarkerFaceColor', 'k');

L0LessThanUndirIndices = intersect(find(~isnan(BirdUnDirMeans)), find(BirdUnDirMeans > 1));
L0LessSlopes = [];
for i = L0LessThanUndirIndices(:)',
    for j = 1:length(FitCoeffs{i}),
        L0LessSlopes(end+1) = FitCoeffs{i}{j}(2);
        switch (SlopeSignificance{i}(j))
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
        plot(1.85, FitCoeffs{i}{j}(2), [MarkerColorSymbol, Symbols(i)], 'MarkerSize', 4, 'MarkerFaceColor', MarkerSymbolFaceColor);    
    end
end

SlopeSigBar(2) = bar(2, mean(L0LessSlopes));
set(SlopeSigBar(2), 'FaceColor', 'none');
errorbar(2, mean(L0LessSlopes), std(L0LessSlopes)/sqrt(length(L0LessSlopes)), 'k', 'MarkerFaceColor', 'k');

axis tight;
Temp = axis;
Temp = [0.5 2.5 1.2*Temp(3) 1.2*Temp(4)];
axis(Temp);
% plot(Temp(1:2), zeros(2,1), 'k--');
set(gca, 'XTick', [1 2], 'XTickLabel', {'L0 Songs > UNDIR' 'L0 Songs < UNDIR'}, 'XTickLabelRotation', 45);
ylabel('Slope');

% Now plot group data for all birds relative to dir L0

p(3).select();
hold on;
DirectedDistanceMeans = BirdDistanceDirMeans(:,1:5);
UndirectedDistanceMeans = BirdUnDirMeans;
% Now to invert the birds for which undir is greater than L0 Dir
InvertBirds = find(UndirectedDistanceMeans > 1);
for i = InvertBirds(:)',
    DirectedDistanceMeans(i,:) = 1./DirectedDistanceMeans(i,:);
    UndirectedDistanceMeans(i) = 1/UndirectedDistanceMeans(i);
end

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

errorbar(Distances(5)+25, nanmean(UndirectedDistanceMeans), nanstd(UndirectedDistanceMeans)/sqrt(length(find(~isnan(UndirectedDistanceMeans)))), 'ks-', 'MarkerSize', 6, 'MarkerFaceColor', 'k');
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
ylabel({'Normalized mean'; 'syllable amplitude'})
xlabel('Distance from female (cm)');

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
Temp = [-5 Temp(2)+5 0.96*Temp(3) 1.02*Temp(4)];
axis(Temp)
plot(Temp(1:2), ones(1,2), 'k--');
text(125, Temp(4), ['p = ', num2str(SamePosPValue(2))]);
text(170, Temp(3), '//', 'FontSize', 12)

p.fontsize = 12;
p.de.margin = 20;
p.marginleft = 25;
p.margintop = 10;
p.marginbottom = 35;
p.marginright = 10;

set(gcf, 'Color', 'w');
set(gcf, 'Units', 'inches');
set(gcf, 'Position', [10 1.4 8.3 3.1]);
set(gcf, 'PaperPositionMode', 'auto');
print(fullfile(FinalFigureDir, ['DUNAlso_HFMicBirds_Example_Summary_Amplitude.eps']), '-depsc2', '-r600');
print(fullfile(FinalFigureDir, ['DUNAlso_HFMicBirds_Example_Summary_Amplitude.png']), '-dpng', '-r600');


% Now plot group data for all syllables relative to dir L0
close all;
figure;
p = panel();
p.pack({1});
p(1).select();
hold on;
DirectedDistanceMeans = SyllableDistanceDirMeans(:,1:5);
UndirectedDistanceMeans = SyllableUnDirMeans;

% Now to invert the birds for which undir is greater than L0 Dir
InvertBirds = find(UndirectedDistanceMeans > 1);
for i = InvertBirds(:)',
    DirectedDistanceMeans(i,:) = 1./DirectedDistanceMeans(i,:);
    UndirectedDistanceMeans(i) = 1/UndirectedDistanceMeans(i);
end

for j = 1:size(DirectedDistanceMeans,1),
    if (~isempty(find(BirdsToPlotIndices == SyllableBirdIndices(j))))
        plot(Distances(1:5), DirectedDistanceMeans(j,:), ['k-', Symbols(SyllableBirdIndices(j))], 'Color', [0.25 0.25 0.25]);
        plot(Distances(5) + 25, UndirectedDistanceMeans(j), ['k-', Symbols(SyllableBirdIndices(j))], 'Color', [0.25 0.25 0.25]);
    else
        plot(Distances(1:5), DirectedDistanceMeans(j,:), ['k-', Symbols(SyllableBirdIndices(j))], 'Color', [0.75 0.75 0.75]);
        plot(Distances(5) + 25, UndirectedDistanceMeans(j), ['k-', Symbols(SyllableBirdIndices(j))], 'Color', [0.75 0.75 0.75]);
    end
end

for j = 1:size(DirectedDistanceMeans,2),
    NonNaNValues = find(~isnan(DirectedDistanceMeans(:,j)));
    errorbar(Distances(j), mean(DirectedDistanceMeans(NonNaNValues,j)), std(DirectedDistanceMeans(NonNaNValues,j))/sqrt(length(NonNaNValues)), 'ks-', 'MarkerSize', 6, 'MarkerFaceColor', 'k');
end

errorbar(Distances(5)+25, nanmean(UndirectedDistanceMeans), nanstd(UndirectedDistanceMeans)/sqrt(length(find(~isnan(UndirectedDistanceMeans)))), 'ks-', 'MarkerSize', 6, 'MarkerFaceColor', 'k');
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
ylabel({'Normalized mean'; 'syllable amplitude'})
xlabel('Distance from female (cm)');

DistanceXLabelString = [];
DistanceXValues = [];
for Dist = 1:length(Distances)-1,
    DistanceXValues(Dist) = Distances(Dist);
    DistanceXLabelString{Dist} = num2str(Distances(Dist));
end
DistanceXValues(end+1) = DistanceXValues(end) + 25;
DistanceXLabelString{end+1} = 'UN';
set(gca, 'XTick', DistanceXValues, 'XTickLabel', DistanceXLabelString);
axis tight;
Temp = axis;
Temp = [-5 Temp(2)+5 0.98*Temp(3) 1.02*Temp(4)];
axis(Temp)
plot(Temp(1:2), ones(1,2), 'k--');
text(125, Temp(4), ['p = ', num2str(SamePosPValue(2))]);
text(170, Temp(3), '//', 'FontSize', 12)

p.fontsize = 12;
p.de.margin = 20;
p.marginleft = 25;
p.margintop = 10;
p.marginbottom = 30;
p.marginright = 10;

set(gcf, 'Color', 'w');
set(gcf, 'Units', 'inches');
set(gcf, 'Position', [10 1.4 4 3]);
set(gcf, 'PaperPositionMode', 'auto');
print(fullfile(FinalFigureDir, ['DUNAlso_HFMicBirds_SyllableWise_Summary_Amplitude.eps']), '-depsc2', '-r600');
print(fullfile(FinalFigureDir, ['DUNAlso_HFMicBirds_SyllableWise_Summary_Amplitude.png']), '-dpng', '-r600');

disp('Done with stats and plots');
