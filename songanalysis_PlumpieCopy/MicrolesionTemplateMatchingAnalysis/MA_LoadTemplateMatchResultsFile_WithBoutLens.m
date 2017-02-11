function [Results] = MA_LoadTemplateMatchResultsFile_WithBoutLens(SongFile, OutputDir, Threshold, Label, SongFileNo)

SmoothingKernelLen = 3;
SmoothingKernel = ones(SmoothingKernelLen,1)/SmoothingKernelLen;

if (~exist(fullfile(OutputDir, [SongFile, '.', Label, '.TempMatch.mat']), 'file'))
    Results = [];
    return;
end

Temp = load(fullfile(OutputDir, [SongFile, '.', Label, '.TempMatch.mat']));

Results = [];
for i = 1:length(Temp.Bout),
    [Pks, Locs] = findpeaks(conv(Temp.Bout{i}.MaxBoutSeqMatch, SmoothingKernel, 'same'), 'MINPEAKHEIGHT', Threshold);
    Pks = Pks(:);
    Locs = Locs(:);
    
    Results = [Results; [Pks (Temp.Bout{i}.T(Locs)' + Temp.Bout{i}.BoutOnset/1000) ones(size(Pks))*SongFileNo]];
end       
