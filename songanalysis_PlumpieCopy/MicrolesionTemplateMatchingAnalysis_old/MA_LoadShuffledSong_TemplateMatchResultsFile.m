function [Results] = MA_LoadShuffledSong_TemplateMatchResultsFile(SongFile, OutputDir, Label)

cd(OutputDir);

SmoothingKernelLen = 3;
SmoothingKernel = ones(SmoothingKernelLen,1)/SmoothingKernelLen;

Files = dir([OutputDir, SongFile, '.', Label, '.*.TempMatch.mat']);
for i = 1:length(Files),
    Temp = load(Files(i).name);
    Results{i} = conv(Temp.Bout.MaxBoutSeqMatch, SmoothingKernel, 'same');
    Results{i} = Results{i}(:);
end
