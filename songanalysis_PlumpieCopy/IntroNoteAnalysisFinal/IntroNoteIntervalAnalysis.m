function [INIAnalysisResults] = IntroNoteIntervalAnalysis(IntroNoteResults, varargin)

% Using a system for AllINFeatLabels that keeps track of whether the intro
% note was the first, last or a middle intro note. The way I do this is by
% having 3 boolean flags for first, middle and last intro note. For
% instance if there was only one intro note, it would have the flags 1 0 1
% to indicate that is the first, it is also the last.

if (nargin == 1)
    MinNumber = 10;
else
    MinNumber = varargin{1};
end

AllINFeats = [];
AllINFeatRatios = [];
AllINFeatLabels = [];
AllINFeatNoofINs = [];
for i = 1:length(IntroNoteResults.NoofINs),
    if (IntroNoteResults.NoofINs(i) == 0)
        continue;
    end
    TempFeats = IntroNoteResults.BoutDetails(i).onsets(IntroNoteResults.INs{i} + 1) - IntroNoteResults.BoutDetails(i).offsets(IntroNoteResults.INs{i}); 
    AllINFeats = [AllINFeats; TempFeats];
    AllINFeatRatios = [AllINFeatRatios; TempFeats(2:end)./TempFeats(1:end-1)];
    AllINFeatNoofINs = [AllINFeatNoofINs; ones(size(TempFeats,1),1)*(size(TempFeats,1))];
    if (size(TempFeats,1) == 1)
        AllINFeatLabels = [AllINFeatLabels; [[0 0 0] [-(size(TempFeats,1)-1):1:0]' [1:1:size(TempFeats)]']]; 
    else
        if (size(TempFeats,1) == 2)
            AllINFeatLabels = [AllINFeatLabels; [[1 0 1; 0 0 0] [-(size(TempFeats,1)-1):1:0]' [1:1:size(TempFeats)]']]; 
        else
            if (size(TempFeats,1) == 3)
                AllINFeatLabels = [AllINFeatLabels; [[1 0 0; 0 0 1; 0 0 0] [-(size(TempFeats,1)-1):1:0]' [1:1:size(TempFeats)]']]; 
            else
                AllINFeatLabels = [AllINFeatLabels; [[1 0 0; repmat([0 1 0], (size(TempFeats,1)-3), 1); 0 0 1; 0 0 0] [-(size(TempFeats,1)-1):1:0]' [1:1:size(TempFeats)]']];
            end
        end
    end
end

BoutIndices = find(IntroNoteResults.WithinBoutNoofINs(:,1) > 0);
for i = 1:length(BoutIndices),
    TempFeats = IntroNoteResults.BoutDetails(IntroNoteResults.WithinBoutINBoutIndices(BoutIndices(i))).onsets(IntroNoteResults.WithinBoutINs{BoutIndices(i)}+1) - IntroNoteResults.BoutDetails(IntroNoteResults.WithinBoutINBoutIndices(BoutIndices(i))).offsets(IntroNoteResults.WithinBoutINs{BoutIndices(i)});
    AllINFeats = [AllINFeats; TempFeats];
    AllINFeatRatios = [AllINFeatRatios; TempFeats(2:end)./TempFeats(1:end-1)];
    AllINFeatNoofINs = [AllINFeatNoofINs; ones(size(TempFeats,1),1)*(size(TempFeats,1))];
    if (size(TempFeats,1) == 1)
        AllINFeatLabels = [AllINFeatLabels; [[0 0 0] [-(size(TempFeats,1)-1):1:0]' [1:1:size(TempFeats)]']]; 
    else
        if (size(TempFeats,1) == 2)
            AllINFeatLabels = [AllINFeatLabels; [[1 0 1; 0 0 0] [-(size(TempFeats,1)-1):1:0]' [1:1:size(TempFeats)]']]; 
        else
            if (size(TempFeats,1) == 3)
                AllINFeatLabels = [AllINFeatLabels; [[1 0 0; 0 0 1; 0 0 0] [-(size(TempFeats,1)-1):1:0]' [1:1:size(TempFeats)]']]; 
            else
                AllINFeatLabels = [AllINFeatLabels; [[1 0 0; repmat([0 1 0], (size(TempFeats,1)-3), 1); 0 0 1; 0 0 0] [-(size(TempFeats,1)-1):1:0]' [1:1:size(TempFeats)]']];
            end
        end
    end
end

AllINFeats = AllINFeats * 1000;
% MeanAllINFeats = mean(AllINFeats);
% MADAllINFeats = std(AllINFeats);
% 
% NormAllINFeats = (AllINFeats - repmat(MeanAllINFeats, size(AllINFeats, 1), 1))./repmat(MADAllINFeats, size(AllINFeats, 1), 1);

MaxINs = max([max(IntroNoteResults.NoofINs) max(IntroNoteResults.WithinBoutNoofINs(:,1))]);

% for i = 1:MaxINs,
%     INInts{i} = [];
%     for j = 1:length(IntroNoteResults.NoofINs),
%         if (IntroNoteResults.NoofINs(j) == i)
%             INInts{i} = [INInts{i}; [(IntroNoteResults.BoutDetails(j).onsets([IntroNoteResults.INs{j}+1]) - IntroNoteResults.BoutDetails(j).offsets([IntroNoteResults.INs{j}]))]'];
%         end
%     end
% end
% 
% BoutIndices = find(IntroNoteResults.WithinBoutNoofINs(:,1) > 0);
% for i = 1:MaxINs,
%     if (i > length(INInts))
%         INInts{i} = [];
%     end
%     for j = 1:length(BoutIndices),
%         if (IntroNoteResults.WithinBoutNoofINs(BoutIndices(j),1) == i)
%             INInts{i} = [INInts{i}; [(IntroNoteResults.BoutDetails(IntroNoteResults.WithinBoutINBoutIndices(BoutIndices(j))).onsets(IntroNoteResults.WithinBoutINs{BoutIndices(j)}+1) - IntroNoteResults.BoutDetails(IntroNoteResults.WithinBoutINBoutIndices(BoutIndices(j))).offsets(IntroNoteResults.WithinBoutINs{BoutIndices(j)}))]'];
%         end
%     end
% end

ValidTrials = [];
NumValidTrials = [];
ValidIntroNoteNos = [];

for i = 1:MaxINs,
    Indices = find(AllINFeatNoofINs == i);
    if (length(Indices) >= MinNumber*i)
        ValidTrials = [ValidTrials; Indices];
        NumValidTrials = [NumValidTrials; length(Indices)];
        ValidIntroNoteNos = [ValidIntroNoteNos; i];
    end
end

for i = min(AllINFeatLabels(ValidTrials,4)):1:max(AllINFeatLabels(ValidTrials,4)),
    Indices = find(AllINFeatLabels(ValidTrials,4) == i);
    Indices = ValidTrials(Indices);
    INIAnalysisResults.CommonLastMedians(max(AllINFeatLabels(ValidTrials,4))-i+1,:) = [i median(AllINFeats(Indices))];
    INIAnalysisResults.CommonLastIQRs(max(AllINFeatLabels(ValidTrials,4))-i+1,:) = [i iqr(AllINFeats(Indices))];
    INIAnalysisResults.CommonLastCVs(max(AllINFeatLabels(ValidTrials,4))-i+1,:) = [i std(AllINFeats(Indices))/mean(AllINFeats(Indices))];
end

INIAnalysisResults.CommonFirstMedians = [];
INIAnalysisResults.CommonFirstIQRs = [];

for i = min(AllINFeatLabels(ValidTrials,5)):1:max(AllINFeatLabels(ValidTrials,5)),
    Indices = find(AllINFeatLabels(ValidTrials,5) == i);
    Indices = ValidTrials(Indices);

    INIAnalysisResults.CommonFirstMedians(i,:) = [i median(AllINFeats(Indices))];
    INIAnalysisResults.CommonFirstIQRs(i,:) = [i iqr(AllINFeats(Indices))];
    INIAnalysisResults.CommonFirstCVs(i,:) = [i std(AllINFeats(Indices))/mean(AllINFeats(Indices))];
end

% Figure for intervals for each sequence length based on position relative
% to common last IN
figure;
set(gcf, 'Position', [108 317 240 300]);
hold on;
Colours = ['brcgmky'];
Symbols = ['osdvp'];

AllInts = ones(MaxINs, MaxINs);
SeqIndex = 1;
for i = 1:MaxINs,
    Indices = find(AllINFeatNoofINs == i);
    TempINFeats = [];
        
    if (length(Indices) >= MinNumber*i)
        %figure;
        %set(gcf, 'Position', [108 317 240 300]);
        %set(gcf, 'Color', 'w');
        for j = 1:i,
            TempINFeats = [TempINFeats AllINFeats(Indices(j:i:end))];
            if (MaxINs > 5)
                Offset = -1/4.5 + (i-1)*1/9;
            else
                Offset = -1/3.5 + (i-1)*1/7;
            end
            %errorbar((j-i)+(Offset), mean(AllINFeats(Indices(j:i:end))), std(AllINFeats(Indices(j:i:end))), [Colours(mod(i-1,length(Colours))+1), Symbols(mod(floor(i/length(Colours)), length(Symbols)) + 1), '-'], 'MarkerSize', 5); 
        end
        % errorbar([-(i-1):1:0]+Offset, mean(TempINFeats), std(TempINFeats), [Colours(mod(i-1, length(Colours)) + 1), 'o-'], 'MarkerSize', 5);
        SortedTempINFeats = sort(TempINFeats, 1);
        if (MaxINs > 5)
            errorbar([-(i-1):1:0]*2+Offset, median(TempINFeats), median(TempINFeats)-SortedTempINFeats(round(0.25*size(TempINFeats,1)),:), SortedTempINFeats(round(0.75*size(TempINFeats,1)),:)-median(TempINFeats), [Colours(mod(i-1, length(Colours)) + 1), 'o-'], 'MarkerSize', 5);
        else
            errorbar([-(i-1):1:0]+Offset, median(TempINFeats), median(TempINFeats)-SortedTempINFeats(round(0.25*size(TempINFeats,1)),:), SortedTempINFeats(round(0.75*size(TempINFeats,1)),:)-median(TempINFeats), [Colours(mod(i-1, length(Colours)) + 1), 'o-'], 'MarkerSize', 5);
        end
        INIAnalysisResults.IndividualSeqMedians{SeqIndex} = [[-(i-1):1:0]' median(TempINFeats)' flipud(INIAnalysisResults.CommonLastMedians(1:i,2)) INIAnalysisResults.CommonFirstMedians(1:i,2) [1:1:size(TempINFeats,2)]'];
        SeqIndex = SeqIndex + 1;
        AllInts(end-i+1:end,i) = mean(TempINFeats);
        AllInts(1:end-i,i) = NaN;
        %errorbar(-(i-1):1:0, mean(TempINFeats), std(TempINFeats), 'ks-', 'MarkerSize', 6, 'MarkerFaceColor', 'k');
        set(gca, 'FontSize', 10, 'FontName', 'Arial');
        set(gca, 'Box', 'off');
        set(gca, 'XTickLabel', []);
    end
end
%plot(INIAnalysisResults.CommonLastMedians(:,1), INIAnalysisResults.CommonLastMedians(:,2), 'ks', 'MarkerSize', 6);
set(gca, 'Box', 'off');
set(gcf, 'Color', 'w');

% Figure for intervals for each sequence length based on position relative
% to common first IN
figure;
set(gcf, 'Position', [108 317 240 300]);
hold on;
Colours = ['brcgmky'];
Symbols = ['osdvp'];

AllInts = ones(MaxINs, MaxINs);
for i = 1:MaxINs,
    Indices = find(AllINFeatNoofINs == i);
    TempINFeats = [];
        
    if (length(Indices) >= MinNumber*i)
        %figure;
        %set(gcf, 'Position', [108 317 240 300]);
        %set(gcf, 'Color', 'w');
        for j = 1:i,
            TempINFeats = [TempINFeats AllINFeats(Indices(j:i:end))];
            Offset = +1/3.5 - (i-1)*1/7;
            %errorbar((j-i)+(Offset), mean(AllINFeats(Indices(j:i:end))), std(AllINFeats(Indices(j:i:end))), [Colours(mod(i-1,length(Colours))+1), Symbols(mod(floor(i/length(Colours)), length(Symbols)) + 1), '-'], 'MarkerSize', 5); 
        end
        % errorbar([-(i-1):1:0]+Offset, mean(TempINFeats), std(TempINFeats), [Colours(mod(i-1, length(Colours)) + 1), 'o-'], 'MarkerSize', 5);
        SortedTempINFeats = sort(TempINFeats, 1);
        errorbar([1:1:size(TempINFeats,2)]+Offset, median(TempINFeats), median(TempINFeats)-SortedTempINFeats(round(0.25*size(TempINFeats,1)),:), SortedTempINFeats(round(0.75*size(TempINFeats,1)),:)-median(TempINFeats), [Colours(mod(i-1, length(Colours)) + 1), 'o-'], 'MarkerSize', 5);
        AllInts(end-i+1:end,i) = mean(TempINFeats);
        AllInts(1:end-i,i) = NaN;
        %errorbar(-(i-1):1:0, mean(TempINFeats), std(TempINFeats), 'ks-', 'MarkerSize', 6, 'MarkerFaceColor', 'k');
        set(gca, 'FontSize', 10, 'FontName', 'Arial');
        set(gca, 'Box', 'off');
        set(gca, 'XTickLabel', []);
    end
end
%plot(INIAnalysisResults.CommonFirstMedians(:,1), INIAnalysisResults.CommonFirstMedians(:,2), 'ks', 'MarkerSize', 6);

NumReps = 100;
INIAnalysisResults.RandomOrderSequences = CalculateRandomOrdering(AllINFeatLabels, AllINFeats, AllINFeatNoofINs, ValidIntroNoteNos, NumReps);

SeqIndex = 1;
INIAnalysisResults.CommonLast.SimilarPositionDistances = [];
INIAnalysisResults.CommonLast.OtherPositionDistances = [];
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
                    Distance = abs(median(AllINFeats(INGroup1)) - median(AllINFeats(INGroup2)));
                    if (j == m)
                        TempSamePosDist = [TempSamePosDist; Distance];
                        INIAnalysisResults.CommonLast.SimilarPositionDistances = [INIAnalysisResults.CommonLast.SimilarPositionDistances; Distance];
                    else
                        TempOtherPosDist = [TempOtherPosDist; Distance];
                        INIAnalysisResults.CommonLast.OtherPositionDistances = [INIAnalysisResults.CommonLast.OtherPositionDistances; Distance];
                    end
                end
            end
        end
        if (~isempty(TempSamePosDist) & ~isempty(TempOtherPosDist))
            INIAnalysisResults.CommonLastIndividualDistances(SeqIndex, :) = [mean(TempSamePosDist) mean(TempOtherPosDist) i j];
            SeqIndex = SeqIndex + 1;
        end
    end
end
   
SeqIndex = 1;
INIAnalysisResults.CommonFirst.SimilarPositionDistances = [];
INIAnalysisResults.CommonFirst.OtherPositionDistances = [];
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
                    Distance = abs(median(AllINFeats(INGroup1)) - median(AllINFeats(INGroup2)));
                    if (j == m)
                        TempSamePosDist = [TempSamePosDist; Distance];
                        INIAnalysisResults.CommonFirst.SimilarPositionDistances = [INIAnalysisResults.CommonFirst.SimilarPositionDistances; Distance];
                    else
                        TempOtherPosDist = [TempOtherPosDist; Distance];
                        INIAnalysisResults.CommonFirst.OtherPositionDistances = [INIAnalysisResults.CommonFirst.OtherPositionDistances; Distance];
                    end
                end
            end
        end
        
        if (~isempty(TempSamePosDist) && ~isempty(TempOtherPosDist))
            if ((length(TempSamePosDist) > 1) && (length(TempOtherPosDist) > 1))
                INIAnalysisResults.CommonFirstIndividualDistances(SeqIndex, :) = [mean(TempSamePosDist) mean(TempOtherPosDist) std(TempSamePosDist) std(TempOtherPosDist) i j];
            else
                if (length(TempSamePosDist) > 1)
                    INIAnalysisResults.CommonFirstIndividualDistances(SeqIndex, :) = [mean(TempSamePosDist) mean(TempOtherPosDist) std(TempSamePosDist) NaN i j];
                else
                    if (length(TempOtherPosDist) > 1)
                        INIAnalysisResults.CommonFirstIndividualDistances(SeqIndex, :) = [mean(TempSamePosDist) mean(TempOtherPosDist) NaN std(TempOtherPosDist) i j];
                    else
                        INIAnalysisResults.CommonFirstIndividualDistances(SeqIndex, :) = [mean(TempSamePosDist) mean(TempOtherPosDist) NaN NaN i j];
                    end
                end
            end
            SeqIndex = SeqIndex + 1;
        end
    end
end

% figure;
% set(gcf, 'Position', [108 109 240 125]);
% hold on;
Colours = ['brkgmcy'];
Symbols = ['sdvp'];
INIAnalysisResults.MeanInts = [];
INIAnalysisResults.STDInts = [];
INIAnalysisResults.IntCV = [];
for i = max(AllINFeatLabels(:,4)):-1:min(AllINFeatLabels(:,4)),
    Indices = find(AllINFeatLabels(:,4) == i);
    if (length(Indices) >= MinNumber)
        % errorbar(i, mean(AllINFeats(Indices)), std(AllINFeats(Indices)), [Colours(mod(abs(i),length(Colours))+1), Symbols(mod(floor(abs(i)/length(Colours)), length(Symbols)) + 1), '-'], 'MarkerSize', 6, 'MarkerFaceColor', Colours(mod(abs(i),length(Colours)) + 1), 'LineWidth', 1); 
        INIAnalysisResults.MeanInts = [INIAnalysisResults.MeanInts; [i mean(AllINFeats(Indices))]];
        INIAnalysisResults.STDInts = [INIAnalysisResults.STDInts; [i std(AllINFeats(Indices))]];
        INIAnalysisResults.IntCV = [INIAnalysisResults.IntCV; [i std(AllINFeats(Indices))./mean(AllINFeats(Indices))]];
    end
    %errorbar(INIAnalysisResults.MeanInts(:,1), INIAnalysisResults.MeanInts(:,2), INIAnalysisResults.STDInts(:,2), 'ks-', 'MarkerSize', 3); 
end

set(gca, 'Box', 'off');
set(gcf, 'Color', 'w');

figure;
set(gcf, 'Position', [108 317 240 300]);
hold on;
Colours = ['brcgmky'];
Symbols = ['osdvp'];

INIAnalysisResults.AllSeqIQRs = [];
SeqIndex = 1;
for i = 1:MaxINs,
    Indices = find(AllINFeatNoofINs == i);
    TempINFeats = [];
    
    if (length(Indices) >= MinNumber*i)
        %figure;
        %set(gcf, 'Position', [108 317 240 300]);
        %set(gcf, 'Color', 'w');
        for j = 1:i,
            TempINFeats = [TempINFeats AllINFeats(Indices(j:i:end))];
            Offset = -1/3.5 + (i-1)*1/7;
            %plot((j-i)+(Offset), std(AllINFeats(Indices(j:i:end)))./mean(AllINFeats(Indices(j:i:end))), [Colours(mod(i-1,length(Colours))+1), Symbols(mod(floor(i/length(Colours)), length(Symbols)) + 1), '-'], 'MarkerSize', 5); 
            %plot((j-i)+(Offset), iqr(AllINFeats(Indices(j:i:end))), [Colours(mod(i-1,length(Colours))+1), Symbols(mod(floor(i/length(Colours)), length(Symbols)) + 1), '-'], 'MarkerSize', 5); 
        end
        plot([-(i-1):1:0]+Offset, iqr(TempINFeats), [Colours(mod(i-1, length(Colours)) + 1), 'o-'], 'MarkerSize', 5);
        INIAnalysisResults.AllSeqIQRs = [INIAnalysisResults.AllSeqIQRs; [[-(i-1):1:0]' iqr(TempINFeats)']];
        INIAnalysisResults.IndividualSeqIQRs{SeqIndex} = [[-(i-1):1:0]' iqr(TempINFeats)' flipud(INIAnalysisResults.CommonLastIQRs(1:i,2)) INIAnalysisResults.CommonFirstIQRs(1:i,2) [1:1:size(TempINFeats,2)]'];
        SeqIndex = SeqIndex + 1;
        set(gca, 'FontSize', 10, 'FontName', 'Arial');
        set(gca, 'Box', 'off');
        set(gca, 'XTickLabel', []);
    end
end

% figure;
% set(gcf, 'Position', [108 109 240 125]);
% hold on;
Colours = ['brkgmcy'];
Symbols = ['sdvp'];
INIAnalysisResults.MeanInts = [];
INIAnalysisResults.MedianInts = [];
INIAnalysisResults.STDInts = [];
INIAnalysisResults.IntCV = [];
INIAnalysisResults.IntIQR = [];
INIAnalysisResults.IntIQRNormalised = [];
for i = max(AllINFeatLabels(:,4)):-1:min(AllINFeatLabels(:,4)),
    Indices = find(AllINFeatLabels(:,4) == i);
    if (length(Indices) >= MinNumber)
        % errorbar(i, mean(AllINFeats(Indices)), std(AllINFeats(Indices)), [Colours(mod(abs(i),length(Colours))+1), Symbols(mod(floor(abs(i)/length(Colours)), length(Symbols)) + 1), '-'], 'MarkerSize', 6, 'MarkerFaceColor', Colours(mod(abs(i),length(Colours)) + 1), 'LineWidth', 1); 
        INIAnalysisResults.MeanInts = [INIAnalysisResults.MeanInts; [i mean(AllINFeats(Indices))]];
        INIAnalysisResults.MedianInts = [INIAnalysisResults.MedianInts; [i median(AllINFeats(Indices))]];
        INIAnalysisResults.STDInts = [INIAnalysisResults.STDInts; [i std(AllINFeats(Indices))]];
        INIAnalysisResults.IntCV = [INIAnalysisResults.IntCV; [i std(AllINFeats(Indices))./mean(AllINFeats(Indices))]];
        INIAnalysisResults.IntIQR = [INIAnalysisResults.IntIQR; [i iqr(AllINFeats(Indices))]];
        TempINSortedFeats = sort(AllINFeats(Indices));
        INIAnalysisResults.IntIQRNormalised = [INIAnalysisResults.IntIQRNormalised; [i (TempINSortedFeats(round(0.75*length(TempINSortedFeats))) - TempINSortedFeats(round(0.25*length(TempINSortedFeats))))/(TempINSortedFeats(round(0.75*length(TempINSortedFeats))) + TempINSortedFeats(round(0.25*length(TempINSortedFeats))))]];
    end
    %plot(INIAnalysisResults.MedianInts(:,1), INIAnalysisResults.IntIQR(:,2), 'ks-', 'MarkerSize', 3); 
end

set(gca, 'Position', [0.13 0.24 0.7750 0.7]);
set(gca, 'Box', 'off');
set(gca, 'FontSize', 10, 'FontName', 'Arial');
set(gcf, 'Color', 'w');

INIAnalysisResults.AllINFeats = AllINFeats;
INIAnalysisResults.AllINFeatLabels = AllINFeatLabels;
INIAnalysisResults.AllINFeatNoofINs = AllINFeatNoofINs;
INIAnalysisResults.AllINFeatRatios = AllINFeatRatios;
INIAnalysisResults.ValidTrials = ValidTrials;


disp('Finished feature analysis');
