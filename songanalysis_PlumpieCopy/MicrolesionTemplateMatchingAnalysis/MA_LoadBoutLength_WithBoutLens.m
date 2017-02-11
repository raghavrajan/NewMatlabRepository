function [Results] = MA_LoadBoutLength_WithBoutLens(SongFile, OutputDir, Label, SongFileNo)

if (~exist(fullfile(OutputDir, [SongFile, '.', Label, '.TempMatch.mat']), 'file'))
    Results = [];
    return;
end

Temp = load(fullfile(OutputDir, [SongFile, '.', Label, '.TempMatch.mat']));

Results = [];
for i = 1:length(Temp.Bout),
    Results = [Results; [i Temp.Bout{i}.BoutOnset Temp.Bout{i}.BoutOffset (Temp.Bout{i}.BoutOffset - Temp.Bout{i}.BoutOnset) SongFileNo]];
end       
