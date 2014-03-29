function [] = CSAnalyseINFrequency(CSData)

% Algorithm
% For each day and for each bout, first get the first motif syllable. The
% syllables before this first motif syllable can be counted as Total
% Syllables before first motif syllable. In addition, the number of INs in
% these syllables can be counted for #of INs.

for i = 1:CSData.NoofDays,
    MotifSyllArray = cellstr(char(ones(length(CSData.AllLabels{i}), 1)*double(CSData.MotifSyllLabels)));
    INArray = cellstr(char(ones(length(CSData.AllLabels{i}), 1)*double(CSData.INLabels)));
    MotifInitiationSyllArray = cellstr(char(ones(length(CSData.AllLabels{i}), 1)*double(CSData.MotifInitiationSyllLabels)));
    
    % using cellfun so that i iterate over each element of the cell array.
    % To use cellfun, all of the other inputs also have to be in the form
    % of cell arrays of the same length - so the previous three lines
    % convert file type, data dir and output dir - common parameters for
    % all of the files into cell arrays
    
    INs{i} = cellfun(@CSIdentifyINs, CSData.AllLabels{i}', MotifSyllArray, INArray, MotifInitiationSyllArray, 'UniformOutput', 0);
end

disp('Finished plotting IN frequencies');
