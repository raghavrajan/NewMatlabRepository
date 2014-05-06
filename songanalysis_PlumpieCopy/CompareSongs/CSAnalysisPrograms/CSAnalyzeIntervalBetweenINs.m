function [] = CSAnalyzeIntervalBetweenINs(CSData)

% Algorithm
% For each day and for each bout, first get the first motif syllable. The
% syllables before this first motif syllable can be counted as Total
% Syllables before first motif syllable. In addition, the number of INs in
% these syllables can be counted for #of INs.

MotifSyllArray = cellstr(char(ones(length(CSData.AllLabels), 1)*double(CSData.MotifSyllLabels)));
INArray = cellstr(char(ones(length(CSData.AllLabels), 1)*double(CSData.INLabels)));
MotifInitiationSyllArray = cellstr(char(ones(length(CSData.AllLabels), 1)*double(CSData.MotifInitiationSyllLabels)));
    
% using cellfun so that i iterate over each element of the cell array.
% To use cellfun, all of the other inputs also have to be in the form
% of cell arrays of the same length - so the previous three lines
% convert file type, data dir and output dir - common parameters for
% all of the files into cell arrays
    
[INs, Motifs, Bouts] = cellfun(@CSIdentifyINs, CSData.AllLabels', MotifSyllArray, INArray, 'UniformOutput', 0);

for i = 1:length(INs),
    MinINPos(i) = min(sum(INs{i}.PosFromLast));
end

MinINPos = min(MinINPos);

figure;
hold on;
for i = 1:CSData.NoofDays,
    for j = MinINPos:-1,
        PosFromLast = sum(INs{i}.PosFromLast);
        Indices = find(PosFromLast == j);
        Intervals = CSData.AllOnsets{i}(INs{i}.Indices(Indices) + 1) - CSData.AllOffsets{i}(INs{i}.Indices(Indices));
        MedianInterval((j + abs(MinINPos) + 1), i) = median(Intervals);
        IQRInterval((j + abs(MinINPos) + 1), i) = iqr(Intervals);
        Percentile_25((j + abs(MinINPos) + 1), i) = prctile(Intervals, 25);
        Percentile_75((j + abs(MinINPos) + 1), i) = prctile(Intervals, 75);
        NumIntervals{j + abs(MinINPos) + 1,i} = ['(', num2str(length(Indices)), ')'];
    end
end
BarPlotHandle = bar(MedianInterval);
hold on;
colormap('gray');
for j = 1:size(MedianInterval, 2),
    XVal = get(get(BarPlotHandle(j), 'children'), 'xData');
    errorbar(mean(XVal,1), MedianInterval(:,j), MedianInterval(:,j) - Percentile_25(:,j), Percentile_75(:,j) - MedianInterval(:,j), 'k.', 'MarkerSize', 2);
    text(XVal(1,:), Percentile_75(:,j)*1.05, NumIntervals(:,j), 'FontSize', 8);
end
axis tight;
Temp = axis;
axis([(Temp(1) - 0.1) (Temp(2) + 0.1) 0 1.1*Temp(4)]);
set(gca, 'XTick', 1:1:abs(MinINPos), 'XTickLabel', (MinINPos:1:-1), 'FontSize', 16);
ylabel('Median interval (msec)', 'FontSize', 16);
title('Intervals between successive INs as a function of position of IN relative to common last IN', 'FontSize', 16);
set(gcf, 'Position', [102 262 850 425]);

disp('Finished analyzing intervals between INs');