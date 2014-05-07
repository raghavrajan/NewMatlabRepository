function [] = CSPlotINDistribution(CSData)

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

BoutBeginningFig = figure;
WithinBoutFig = figure;

for i = 1:length(INs),
    MaxINs(i) = max(INs{i}.NumINs(Motifs{i}.BoutBeginningMotifs));
end

MaxINs = max(MaxINs) + 1;
XVals = 0:1:MaxINs;

Colors = 'rgbcmk';
Symbols = 'o+>';

for i = 1:length(INs),
    BoutBeginningINDist(i,:) = histc(INs{i}.NumINs(Motifs{i}.BoutBeginningMotifs), XVals);
    BoutBeginningINDist(i,:) = BoutBeginningINDist(i,:)/sum(BoutBeginningINDist(i,:));
    figure(BoutBeginningFig);
    hold on;
    plot(XVals, BoutBeginningINDist(i,:)*100, [Colors(mod(i-1, length(Colors)) + 1), Symbols(ceil(i/length(Colors))), '-'], 'LineWidth', 2);
    
    WithinBoutINDist(i,:) = histc(INs{i}.NumINs(Motifs{i}.WithinBoutMotifs), XVals);
    WithinBoutINDist(i,:) = WithinBoutINDist(i,:)/sum(WithinBoutINDist(i,:));
    figure(WithinBoutFig);
    hold on;
    plot(XVals, WithinBoutINDist(i,:)*100, [Colors(mod(i-1, length(Colors)) + 1), Symbols(ceil(i/length(Colors))), '-'], 'LineWidth', 2);
end

figure(BoutBeginningFig);
axis([XVals(1) XVals(end) 0 max(BoutBeginningINDist(:))*105]);
title('Bout beginnning', 'FontSize', 16);
xlabel('Number of INs', 'FontSize', 16);
ylabel('% of trials', 'FontSize', 16);
set(gca, 'FontSize', 16);

figure(WithinBoutFig);
axis([XVals(1) XVals(end) 0 max(WithinBoutINDist(:))*105]);
title('Within bouts', 'FontSize', 16);
xlabel('Number of INs', 'FontSize', 16);
ylabel('% of trials', 'FontSize', 16);
set(gca, 'FontSize', 16);

disp('Finished plotting IN frequencies');
