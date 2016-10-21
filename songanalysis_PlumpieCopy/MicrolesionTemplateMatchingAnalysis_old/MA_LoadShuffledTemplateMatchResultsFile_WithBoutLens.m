function [Results] = MA_LoadShuffledTemplateMatchResultsFile_WithBoutLens(SongFile, OutputDir, Label)

SmoothingKernelLen = 3;
SmoothingKernel = ones(SmoothingKernelLen,1)/SmoothingKernelLen;

Files = dir(fullfile(OutputDir, [SongFile, '.', Label, '.*.TempMatch.mat']));

if isempty(Files)
    Results = [];
    return;
end

Results = [];
for i = 1:length(Files),
    Temp = load(fullfile(OutputDir, Files(i).name));
    Results = [Results; max(conv(Temp.Bout.MaxBoutSeqMatch, SmoothingKernel, 'same'))];
end       
