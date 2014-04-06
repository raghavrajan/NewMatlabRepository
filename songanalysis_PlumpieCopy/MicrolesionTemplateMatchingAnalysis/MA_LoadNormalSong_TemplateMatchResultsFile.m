function [Results] = MA_LoadNormalSong_TemplateMatchResultsFile(SongFile, OutputDir, Label)

cd(OutputDir);

SmoothingKernelLen = 3;
SmoothingKernel = ones(SmoothingKernelLen,1)/SmoothingKernelLen;

Files = dir([OutputDir, SongFile, '.', Label, '.TempMatch.mat']);
Temp = load(Files(1).name);
Results = max(conv(Temp.Bout.MaxBoutSeqMatch, SmoothingKernel, 'same'));
