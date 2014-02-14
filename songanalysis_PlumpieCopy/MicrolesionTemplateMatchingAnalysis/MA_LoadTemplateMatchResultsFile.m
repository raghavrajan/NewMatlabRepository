function [Results] = MA_LoadTemplateMatchResultsFile(SongFile, OutputDir, Threshold, Label, SongFileNo)

SmoothingKernelLen = 3;
SmoothingKernel = ones(SmoothingKernelLen,1)/SmoothingKernelLen;

Temp = load([OutputDir, SongFile, '.', Label, '.TempMatch.mat']);

[Pks, Locs] = findpeaks(conv(Temp.Bout.MaxBoutSeqMatch, SmoothingKernel, 'same'), 'MINPEAKHEIGHT', Threshold);
Pks = Pks(:);
Locs = Locs(:);
        
Results = [Pks Temp.Bout.T(Locs)' ones(size(Pks))*SongFileNo];
        
