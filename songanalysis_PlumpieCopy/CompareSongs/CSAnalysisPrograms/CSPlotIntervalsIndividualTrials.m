function [] = CSPlotIntervalsIndividualTrials(CSData)

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
for i = 1:CSData.NoofDays,
    subplot(CSData.NoofDays, 1, i);
    hold on;
    for j = 1:length(INs{i}.Starts),
        Indices = INs{i}.Starts(j):INs{i}.Ends(j);
        Intervals = CSData.AllOnsets{i}(Indices + 1) - CSData.AllOffsets{i}(Indices);
        plot(-(length(Indices)):1:-1, Intervals, 'ko-', 'MarkerSize', 2);
    end
end

for i = 1:CSData.NoofDays,
    subplot(CSData.NoofDays, 1, i);
    set(gca, 'YScale', 'log');
    axis tight;
    Temp = axis;
    Temp(1) = Temp(1) - 0.5;
    Temp(2) = Temp(2) + 0.5;
    Temp(3) = 10;
    Temp(4) = Temp(4)*1.02
end


figure;
hold on;
for i = 1:CSData.NoofDays,
    hold on;
    for j = 1:length(INs{i}.Starts),
        Indices = INs{i}.Starts(j):INs{i}.Ends(j);
        Intervals = CSData.AllOnsets{i}(Indices + 1) - CSData.AllOffsets{i}(Indices);
        plot(-(length(Indices)):1:-1, Intervals, [Colors(i), 'o-'], 'MarkerSize', 2);
    end
end

for i = 1:CSData.NoofDays,
    set(gca, 'YScale', 'log');
    axis tight;
    Temp = axis;
    Temp(1) = Temp(1) - 0.5;
    Temp(2) = Temp(2) + 0.5;
    Temp(3) = 10;
    Temp(4) = Temp(4)*1.02
end


disp('Finished analyzing intervals between INs');