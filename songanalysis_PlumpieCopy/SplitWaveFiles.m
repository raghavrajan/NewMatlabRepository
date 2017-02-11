function [] = SplitWaveFiles(DirName, FileIdentifierString, MaxFileDur)

cd(DirName);
mkdir('SplitWavFiles');
Files = dir(['*', FileIdentifierString, '*']);

for i = 1:length(Files),
    [Data, Fs] = wavread(Files(i).name);
    FileLen = length(Data)/Fs;
    NumToSplit = ceil(FileLen/MaxFileDur);
    fprintf(['Splitting file ', Files(i).name, ':']);
    cd('SplitWavFiles');
    for j = 1:NumToSplit,
        fprintf('>');
        StartIndex = ((j-1)*MaxFileDur*Fs) + 1;
        EndIndex = j*MaxFileDur*Fs;
        if (EndIndex > length(Data))
            EndIndex = length(Data);
        end
        wavwrite(Data(StartIndex:EndIndex), Fs, [Files(i).name, '.', num2str(j), '.wav']);
    end
    fprintf('\n');
    cd(DirName);
end
        