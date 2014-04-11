function [] = CSAnalyseINFrequency(CSData)

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

figure;

for i = 1:length(INs),
    MeanNumINs(1,i) = mean(INs{i}.NumINs(Motifs{i}.BoutBeginningMotifs));
    STDNumINs(1,i) = std(INs{i}.NumINs(Motifs{i}.BoutBeginningMotifs));
    
    MeanNumINs(2,i) = mean(INs{i}.NumINs(Motifs{i}.WithinBoutMotifs));
    STDNumINs(2,i) = std(INs{i}.NumINs(Motifs{i}.WithinBoutMotifs));
end
BarPlotHandle = bar(MeanNumINs);
hold on;
colormap('gray');
for j = 1:size(MeanNumINs, 2),
    XVal = get(get(BarPlotHandle(j), 'children'), 'xData');
    errorbar(mean(XVal,1), MeanNumINs(:,j), STDNumINs(:,j), 'k.', 'MarkerSize', 2);
end
axis tight;
Temp = axis;
axis([(Temp(1) - 0.1) (Temp(2) + 0.1) 0 1.1*Temp(4)]);
set(gca, 'XTick', [1 2], 'XTickLabel', [{'Bout Beginning'}; {'Within Bouts'}], 'FontSize', 16);
title('Number of INs', 'FontSize', 16);
ylabel('Number of INs', 'FontSize', 16);

disp('Finished plotting IN frequencies');
