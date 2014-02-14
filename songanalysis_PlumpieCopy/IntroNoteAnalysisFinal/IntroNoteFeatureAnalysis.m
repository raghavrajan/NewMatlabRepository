function [INFAnalysisResults] = IntroNoteFeatureAnalysis(IntroNoteResults, BoutType, PlotFeatCols, varargin)

% Using a system for AllINFeatLabels that keeps track of whether the intro
% note was the first, last or a middle intro note. The way I do this is by
% having 3 boolean flags for first, middle and last intro note. For
% instance if there was only one intro note, it would have the flags 1 0 1
% to indicate that is the first, it is also the last.

if (nargin == 3)
    MinNumber = 10;
else
    MinNumber = varargin{1};
end

FeatureCols = [1 2 3 4];

AllINFeats = [];
AllINFeatLabels = [];
AllINFeatNoofINs = [];
AllINTimeToMotifCommonEnd = [];
AllINTimeToMotifCommonStart = [];
for i = 1:length(IntroNoteResults.NoofINs),
    TempFeats = IntroNoteResults.BoutDetails(i).Feats(IntroNoteResults.INs{i},FeatureCols);
    if (IntroNoteResults.NoofINs(i) > 0)
        AllINTimeToMotifCommonEnd = [AllINTimeToMotifCommonEnd; -(IntroNoteResults.BoutDetails(i).onsets(IntroNoteResults.MotifStartIndex(i)) - (IntroNoteResults.BoutDetails(i).onsets(IntroNoteResults.INs{i})))];
        AllINTimeToMotifCommonStart = [AllINTimeToMotifCommonStart; -(IntroNoteResults.BoutDetails(i).onsets(IntroNoteResults.INs{i}(1)) - (IntroNoteResults.BoutDetails(i).onsets(IntroNoteResults.INs{i})))];
    end
        
    AllINFeats = [AllINFeats; TempFeats];
    AllINFeatNoofINs = [AllINFeatNoofINs; ones(size(TempFeats,1),1)*size(TempFeats,1)];
    if (~isempty(TempFeats))
        if (size(TempFeats,1) == 1)
            AllINFeatLabels = [AllINFeatLabels; [[1 0 1] [-size(TempFeats,1):1:-1]'] [1:1:size(TempFeats,1)]' ones(size(TempFeats,1),1)*-1000]; 
        else
            if (size(TempFeats,1) == 2)
                AllINFeatLabels = [AllINFeatLabels; [[1 0 0; 0 0 1] [-size(TempFeats,1):1:-1]'] [1:1:size(TempFeats,1)]' ones(size(TempFeats,1),1)*-1000];
            else
                AllINFeatLabels = [AllINFeatLabels; [[1 0 0; repmat([0 1 0], (size(TempFeats,1)-2), 1); [0 0 1]] [-size(TempFeats,1):1:-1]'] [1:1:size(TempFeats,1)]' ones(size(TempFeats,1),1)*-1000];
            end
        end
    end
end

BoutIndices = find(IntroNoteResults.WithinBoutNoofINs(:,1) > 0);
for i = 1:length(BoutIndices),
    AllINTimeToMotifCommonEnd = [AllINTimeToMotifCommonEnd; -(IntroNoteResults.BoutDetails(IntroNoteResults.WithinBoutINBoutIndices(BoutIndices(i))).onsets(IntroNoteResults.WithinBoutNoofINs(BoutIndices(i),4)) - IntroNoteResults.BoutDetails(IntroNoteResults.WithinBoutINBoutIndices(BoutIndices(i))).onsets(IntroNoteResults.WithinBoutINs{BoutIndices(i)}))];
    AllINTimeToMotifCommonStart = [AllINTimeToMotifCommonStart; -(IntroNoteResults.BoutDetails(IntroNoteResults.WithinBoutINBoutIndices(BoutIndices(i))).onsets(IntroNoteResults.WithinBoutINs{BoutIndices(i)}(1)) - IntroNoteResults.BoutDetails(IntroNoteResults.WithinBoutINBoutIndices(BoutIndices(i))).onsets(IntroNoteResults.WithinBoutINs{BoutIndices(i)}))];
    
    TempFeats = IntroNoteResults.BoutDetails(IntroNoteResults.WithinBoutINBoutIndices(BoutIndices(i))).Feats(IntroNoteResults.WithinBoutINs{BoutIndices(i)},FeatureCols);
    AllINFeats = [AllINFeats; TempFeats];
    AllINFeatNoofINs = [AllINFeatNoofINs; ones(size(TempFeats,1),1)*size(TempFeats,1)];
    if (~isempty(TempFeats))
        if (size(TempFeats,1) == 1)
            AllINFeatLabels = [AllINFeatLabels; [[1 0 1] [-size(TempFeats,1):1:-1]'] [1:1:size(TempFeats,1)]' ones(size(TempFeats,1),1)*IntroNoteResults.WithinBoutNoofINs(BoutIndices(i),5)]; 
        else
            if (size(TempFeats,1) == 2)
                AllINFeatLabels = [AllINFeatLabels; [[1 0 0; 0 0 1] [-size(TempFeats,1):1:-1]'] [1:1:size(TempFeats,1)]' ones(size(TempFeats,1),1)*IntroNoteResults.WithinBoutNoofINs(BoutIndices(i),5)];
            else
                AllINFeatLabels = [AllINFeatLabels; [[1 0 0; repmat([0 1 0], (size(TempFeats,1)-2), 1); [0 0 1]] [-size(TempFeats,1):1:-1]'] [1:1:size(TempFeats,1)]' ones(size(TempFeats,1),1)*IntroNoteResults.WithinBoutNoofINs(BoutIndices(i),5)];
            end
        end
    end
end

MeanAllINFeats = mean(AllINFeats);
MADAllINFeats = std(AllINFeats);

%NormAllINFeats = (AllINFeats - repmat(MeanAllINFeats, size(AllINFeats, 1), 1))./repmat(MADAllINFeats, size(AllINFeats, 1), 1);
%NormAllINFeats = AllINFeats./repmat(MADAllINFeats, size(AllINFeats, 1), 1);
NormAllINFeats = AllINFeats;

MaxINs = max([max(IntroNoteResults.NoofINs) max(IntroNoteResults.WithinBoutNoofINs(:,1))]);

for i = 1:MaxINs,
    INFeats{i} = [];
end

% switch (BoutType)
%     case 'Beginning'
%         for i = 1:MaxINs,
%             INFeats{i} = [];
%             for j = 1:length(IntroNoteResults.NoofINs),
%                 if (IntroNoteResults.NoofINs(j) == i)
%                     INFeats{i} = [INFeats{i}; [IntroNoteResults.BoutDetails(j).Feats([IntroNoteResults.INs{j}],FeatureCols) (i:-1:1)']];
%                 end
%             end
%         end
% 
%     case 'Within'
%         BoutIndices = find(IntroNoteResults.WithinBoutNoofINs(:,1) > 0);
%         for i = 1:MaxINs,
%             if (i > length(INFeats))
%                 INFeats{i} = [];
%             end
%             for j = 1:length(BoutIndices),
%                 if (IntroNoteResults.WithinBoutNoofINs(BoutIndices(j),1) == i)
%                     INFeats{i} = [INFeats{i}; [IntroNoteResults.BoutDetails(IntroNoteResults.WithinBoutINBoutIndices(BoutIndices(j))).Feats(IntroNoteResults.WithinBoutINs{BoutIndices(j)},FeatureCols) (i:-1:1)']];
%                 end
%             end
%         end
%         
%     case 'All'
%                 for i = 1:MaxINs,
%             INFeats{i} = [];
%             for j = 1:length(IntroNoteResults.NoofINs),
%                 if (IntroNoteResults.NoofINs(j) == i)
%                     INFeats{i} = [INFeats{i}; [IntroNoteResults.BoutDetails(j).Feats([IntroNoteResults.INs{j}],FeatureCols) (i:-1:1)']];
%                 end
%             end
%         end
%         BoutIndices = find(IntroNoteResults.WithinBoutNoofINs(:,1) > 0);
%         for i = 1:MaxINs,
%             if (i > length(INFeats))
%                 INFeats{i} = [];
%             end
%             for j = 1:length(BoutIndices),
%                 if (IntroNoteResults.WithinBoutNoofINs(BoutIndices(j),1) == i)
%                     INFeats{i} = [INFeats{i}; [IntroNoteResults.BoutDetails(IntroNoteResults.WithinBoutINBoutIndices(BoutIndices(j))).Feats(IntroNoteResults.WithinBoutINs{BoutIndices(j)},FeatureCols) (i:-1:1)']];
%                 end
%             end
%         end
% end
% 
% for i = 1:MaxINs,
%     if (~isempty(INFeats{i}))
%         NormINFeats{i} = (INFeats{i}(:,1:end-1) - repmat(MeanAllINFeats, size(INFeats{i},1), 1))./repmat(MADAllINFeats, size(INFeats{i},1),1);
%         % NormINFeats{i} = INFeats{i}(:,1:end-1)./repmat(MADAllINFeats, size(INFeats{i},1),1);
%         % NormINFeats{i} = INFeats{i}(:,1:end-1);
%     end
% end

for i = 1:MaxINs,
    NumTrials(i) = length(find(AllINFeatNoofINs == i))/i;
end
[MaxTrials, MaxTrialIndex] = max(NumTrials(2:end));
MaxTrialIndex = MaxTrialIndex + 1;

ValidTrials = [];
NumValidTrials = [];
ValidIntroNoteNos = [];
SeqIndex = 1;
for i = 1:MaxINs,
    Indices = find(AllINFeatNoofINs == i);
    if (length(Indices) >= MinNumber*i)
        ValidTrials = [ValidTrials; Indices];
        NumValidTrials = [NumValidTrials; length(ValidTrials)];
        ValidIntroNoteNos = [ValidIntroNoteNos; i];
    end
end

INFAnalysisResults.ValidTrials = ValidTrials;
INFAnalysisResults.NumValidTrials = NumValidTrials;
INFAnalysisResults.ValidIntroNoteNos = ValidIntroNoteNos;

NumReps = 100;
% [INFAnalysisResults.RandomOrderSequences] = CalculateRandomOrderingAcousticFeatures(AllINFeatLabels, AllINFeats, AllINFeatNoofINs, ValidIntroNoteNos, NumReps);

INFAnalysisResults.CommonLast.SimilarPositionDistances = [];
INFAnalysisResults.CommonLast.OtherPositionDistances = [];

INNo = 1;
for i = ValidIntroNoteNos',
    Indices1 = find(AllINFeatNoofINs == i);
    for j = min(AllINFeatLabels(Indices1,4)):1:max(AllINFeatLabels(Indices1,4)),
        TempSamePosDist = [];
        TempOtherPosDist = [];
        INGroup1 = find((AllINFeatLabels(:,4) == j) & (AllINFeatNoofINs == i));
        for k = ValidIntroNoteNos',
            Indices2 = find(AllINFeatNoofINs == k);
            for m = min(AllINFeatLabels(Indices2,4)):1:max(AllINFeatLabels(Indices2,4)),
                INGroup2 = find((AllINFeatLabels(:,4) == m) & (AllINFeatNoofINs == k));
                if ((j == m) && (i == k))
                    continue;
                else
                    Distance = pdist2(mean(AllINFeats(INGroup1,:)), mean(AllINFeats(INGroup2,:)), 'mahalanobis', cov(AllINFeats(INGroup2,:)));
                    if (j == m)
                        TempSamePosDist = [TempSamePosDist; Distance];
                        INFAnalysisResults.CommonLast.SimilarPositionDistances = [INFAnalysisResults.CommonLast.SimilarPositionDistances; Distance];
                    else
                        TempOtherPosDist = [TempOtherPosDist; Distance];
                        INFAnalysisResults.CommonLast.OtherPositionDistances = [INFAnalysisResults.CommonLast.OtherPositionDistances; Distance];
                    end
                end
            end
        end
        if (~isempty(TempSamePosDist) & ~isempty(TempOtherPosDist))
            INFAnalysisResults.CommonLastIndividualDistances(INNo, :) = [mean(TempSamePosDist) mean(TempOtherPosDist) i j];
            INNo = INNo + 1;
        end
    end
end

INFAnalysisResults.CommonFirst.SimilarPositionDistances = [];
INFAnalysisResults.CommonFirst.OtherPositionDistances = [];

INNo = 1;
for i = ValidIntroNoteNos',
    Indices1 = find(AllINFeatNoofINs == i);
    for j = min(AllINFeatLabels(Indices1,5)):1:max(AllINFeatLabels(Indices1,5)),
        TempSamePosDist = [];
        TempOtherPosDist = [];
        INGroup1 = find((AllINFeatLabels(:,5) == j) & (AllINFeatNoofINs == i));
        for k = ValidIntroNoteNos',
            Indices2 = find(AllINFeatNoofINs == k);
            for m = min(AllINFeatLabels(Indices2,5)):1:max(AllINFeatLabels(Indices2,5)),
                INGroup2 = find((AllINFeatLabels(:,5) == m) & (AllINFeatNoofINs == k));
                if ((j == m) && (i == k))
                    continue;
                else
                    Distance = pdist2(mean(AllINFeats(INGroup1,:)), mean(AllINFeats(INGroup2,:)), 'mahalanobis', cov(AllINFeats(INGroup2,:)));
                    if (j == m)
                        TempSamePosDist = [TempSamePosDist; Distance];
                        INFAnalysisResults.CommonFirst.SimilarPositionDistances = [INFAnalysisResults.CommonFirst.SimilarPositionDistances; Distance];
                    else
                        TempOtherPosDist = [TempOtherPosDist; Distance];
                        INFAnalysisResults.CommonFirst.OtherPositionDistances = [INFAnalysisResults.CommonFirst.OtherPositionDistances; Distance];
                    end
                end
            end
        end
        if (~isempty(TempSamePosDist) & ~isempty(TempOtherPosDist))
            INFAnalysisResults.CommonFirstIndividualDistances(INNo, :) = [mean(TempSamePosDist) mean(TempOtherPosDist) i j];
            INNo = INNo + 1;
        end
    end
end


for i = min(AllINFeatLabels(:,4)):1:max(AllINFeatLabels(:,4)),
    Indices = find(AllINFeatLabels(ValidTrials,4) == i);
    Indices = ValidTrials(Indices);
    CommonLastPositionIndices{max(AllINFeatLabels(:,4))-i+1} = Indices;
    MeanINValues(max(AllINFeatLabels(:,4))-i+1,:) = mean(NormAllINFeats(Indices,:));
    CommonLastMeanINValues(max(AllINFeatLabels(:,4))-i+1,:) = mean(NormAllINFeats(Indices,:));
end

for i = min(AllINFeatLabels(:,5)):1:max(AllINFeatLabels(:,5)),
    Indices = find(AllINFeatLabels(ValidTrials,5) == i);
    Indices = ValidTrials(Indices);
    CommonFirstPositionIndices{i} = Indices;
    CommonFirstMeanINValues(i,:) = mean(NormAllINFeats(Indices,:));
end

Colours = ['brkgmcy'];
Symbols = ['sdvp'];
FinalAxis = 1.02*[min(NormAllINFeats(:,PlotFeatCols(1))) max(NormAllINFeats(:,PlotFeatCols(1))) min(NormAllINFeats(:,PlotFeatCols(2))) max(NormAllINFeats(:,PlotFeatCols(2)))];

INFAnalysisResults.PosBasedMeanDistanceToLast = [];

SeqIndex = 1;
for i = 1:MaxINs,
    RemovalIndices = [];
    Indices = find(AllINFeatNoofINs == i);
%     for j = 1:i:length(Indices),
%         if (~isempty(find(abs(NormAllINFeats(Indices(j:j+i-1),:)) > 40000)))
%             RemovalIndices = [RemovalIndices; [j:j+i-1]'];
%         end
%     end
    if (~isempty(RemovalIndices))
        Indices(RemovalIndices) = [];
    end
    
    if (length(Indices) >= MinNumber*i)
        figure;
        set(gcf, 'Color', 'w');
        set(gcf, 'Position', [108 317 240 300]);
        hold on;
        for j = 1:i,
            plot(NormAllINFeats(Indices(j:i:end),PlotFeatCols(1)), NormAllINFeats(Indices(j:i:end),PlotFeatCols(2)), [Colours(mod(i-j, length(Colours))+1), '+'], 'MarkerSize', 4, 'MarkerFaceColor', [Colours(mod(i-j,length(Colours))+1)]);
            hold on;
            PlotConfidenceEllipse(NormAllINFeats(Indices(j:i:end),PlotFeatCols(1:2)), Colours(mod(i-j,length(Colours))+1), 1);
        end
        
        for j = 1:i,
            %MeanFeatValues = mean(NormAllINFeats(Indices(j:i:end),:));
            %plot(MeanFeatValues(PlotFeatCols(1)), MeanFeatValues(PlotFeatCols(2)), ['w', '+'], 'MarkerSize', 10, 'LineWidth', 3);
            %plot(MeanFeatValues(PlotFeatCols(1)), MeanFeatValues(PlotFeatCols(2)), [Colours(mod(i-j, length(Colours))+1), '+'], 'MarkerSize', 8, 'LineWidth', 2);
            
            %plot(MeanINValues(i-j+1,PlotFeatCols(1)), MeanINValues(i-j+1,PlotFeatCols(2)), ['w', 'x'], 'MarkerSize', 6, 'LineWidth', 3);
            %plot(MeanINValues(i-j+1,PlotFeatCols(1)), MeanINValues(i-j+1,PlotFeatCols(2)), [Colours(mod(i-j, length(Colours))+1), 'x'], 'MarkerSize', 4, 'LineWidth', 2);
        end
%         plot(MeanLastINValues(PlotFeatCols(1)), MeanLastINValues(PlotFeatCols(2)), 'wx', 'MarkerSize', 12, 'LineWidth', 3);
%         plot(MeanLastINValues(PlotFeatCols(1)), MeanLastINValues(PlotFeatCols(2)), 'gx', 'MarkerSize', 10, 'LineWidth', 2);
        SeqIndex = SeqIndex + 1;
        axis(FinalAxis);
        set(gca, 'FontSize', 10, 'FontName', 'Arial');
        set(gca, 'Box', 'off');
    end
end

set(gcf, 'Color', 'w');

% Now plot considering a common first IN reference scheme

Colours = ['krbgmcy'];
for i = 1:MaxINs,
    Indices = find(AllINFeatNoofINs == i);

    if (length(Indices) >= MinNumber*i)
        figure;
        set(gcf, 'Color', 'w');
        set(gcf, 'Position', [108 317 240 300]);
        hold on;
        for j = 1:i,
            plot(NormAllINFeats(Indices(j:i:end),PlotFeatCols(1)), NormAllINFeats(Indices(j:i:end),PlotFeatCols(2)), [Colours(mod(j-1, length(Colours)) +1), '.'], 'MarkerSize', 4, 'MarkerFaceColor', [Colours(mod(j-1,length(Colours)) + 1)]);
            hold on;
            PlotConfidenceEllipse(NormAllINFeats(Indices(j:i:end),PlotFeatCols(1:2)), Colours(mod(j-1,length(Colours))+1), 1);
        end
        
        for j = 1:i,
            MeanFeatValues = mean(NormAllINFeats(Indices(j:i:end),:));
            %plot(MeanFeatValues(PlotFeatCols(1)), MeanFeatValues(PlotFeatCols(2)), ['w', '+'], 'MarkerSize', 10, 'LineWidth', 3);
            %plot(MeanFeatValues(PlotFeatCols(1)), MeanFeatValues(PlotFeatCols(2)), [Colours(mod(j, length(Colours))), '+'], 'MarkerSize', 8, 'LineWidth', 2);
            
            %plot(CommonFirstMeanINValues(j,PlotFeatCols(1)), CommonFirstMeanINValues(j,PlotFeatCols(2)), ['w', 'x'], 'MarkerSize', 6, 'LineWidth', 3);
            %plot(CommonFirstMeanINValues(j,PlotFeatCols(1)), CommonFirstMeanINValues(j,PlotFeatCols(2)), [Colours(mod(j, length(Colours))), 'x'], 'MarkerSize', 4, 'LineWidth', 2);
        end
%         plot(MeanLastINValues(PlotFeatCols(1)), MeanLastINValues(PlotFeatCols(2)), 'wx', 'MarkerSize', 12, 'LineWidth', 3);
%         plot(MeanLastINValues(PlotFeatCols(1)), MeanLastINValues(PlotFeatCols(2)), 'gx', 'MarkerSize', 10, 'LineWidth', 2);
        SeqIndex = SeqIndex + 1;
        axis(FinalAxis);
        set(gca, 'FontSize', 10, 'FontName', 'Arial');
        set(gca, 'Box', 'off');
    end
end

set(gcf, 'Color', 'w');


INFAnalysisResults.AllINFeats = AllINFeats;
INFAnalysisResults.NormAllINFeats = NormAllINFeats;
INFAnalysisResults.AllINFeatLabels = AllINFeatLabels;
INFAnalysisResults.AllINFeatNoofINs = AllINFeatNoofINs;

% INFAnalysisResults.INFeats = INFeats;
% INFAnalysisResults.NormINFeats = NormINFeats;

INFAnalysisResults.MeanAllINFeats = MeanAllINFeats;
INFAnalysisResults.MADAllINFeats = MADAllINFeats;

INFAnalysisResults.PosBasedDistanceToLast = [];
INFAnalysisResults.PosBasedVar = [];
LastINIndices = find(AllINFeatLabels(ValidTrials,4) == -1);
LastINIndices = ValidTrials(LastINIndices);
INFAnalysisResults.LastINIndices = LastINIndices;

for i = min(AllINFeatLabels(ValidTrials,4)):1:max(AllINFeatLabels(ValidTrials,4));
    Indices = find(AllINFeatLabels(ValidTrials,4) == i);
    Indices = ValidTrials(Indices);
    if (length(Indices) >= MinNumber)
        INFAnalysisResults.PosBasedDistanceToLast = [INFAnalysisResults.PosBasedDistanceToLast; [i mean(pdist2(NormAllINFeats(Indices,:), mean(NormAllINFeats(LastINIndices,:)), 'mahalanobis', cov(NormAllINFeats(LastINIndices,:))))]];
        INFAnalysisResults.PosBasedVar = [INFAnalysisResults.PosBasedVar; [i det(cov(NormAllINFeats(Indices,:)))]];
    end
end

FirstINIndices = find(AllINFeatLabels(ValidTrials,1) == 1);
FirstINIndices = ValidTrials(FirstINIndices);
INFAnalysisResults.FirstINIndices = FirstINIndices;

INFAnalysisResults.DistanceToLast = pdist2(NormAllINFeats, mean(NormAllINFeats(LastINIndices,:)), 'mahalanobis', cov(NormAllINFeats(LastINIndices,:)));
INFAnalysisResults.FirstINsDistanceToLast = [AllINFeatNoofINs(FirstINIndices) pdist2(NormAllINFeats(FirstINIndices,:), mean(NormAllINFeats(LastINIndices,:)), 'mahalanobis', cov(NormAllINFeats(LastINIndices,:)))];
disp('Finished feature analysis');
