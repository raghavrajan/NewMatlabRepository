function [] = BatchSongSeqComparison(DirectoryName, FileType, Motif, OutputFileName)

cd(DirectoryName);

if (ispc)
    if (DirectoryName(end) ~= '\')
        DirectoryName(end+1) = '\';
    end
else
    if (DirectoryName(end) ~= '/')
        DirectoryName(end+1) = '/';
    end
end

x = 1:1:size(Motif,2);

LogFile = fopen([OutputFileName, '.log'], 'w');
ResultsFile = fopen([OutputFileName, '.txt'], 'w');

WarpIndex = 1;
for i = [0 5 10 15],
    xx = linspace(x(1), x(end), (size(Motif,2) * (1 + i/100)));
    clear WMotif;
    for j = 1:size(Motif, 1);
        WMotif(j,:) = spline(x, Motif(j,:), xx);
    end
    DSeqMatchValues{WarpIndex} = CalculateMotifSeqMatch(DirectoryName, 'DirSongFileList.txt', FileType, WMotif, 'off', LogFile, ResultsFile);
    USeqMatchValues{WarpIndex} = CalculateMotifSeqMatch(DirectoryName, 'UnDirSongFileList.txt', FileType, WMotif, 'off', LogFile, ResultsFile);
    WarpIndex = WarpIndex + 1;
end

save(OutputFileName, 'DSeqMatchValues', 'USeqMatchValues');

cd ../../12-10-2007/RawDataFiles
DirectoryName = pwd;

if (ispc)
    if (DirectoryName(end) ~= '\')
        DirectoryName(end+1) = '\';
    end
else
    if (DirectoryName(end) ~= '/')
        DirectoryName(end+1) = '/';
    end
end

x = 1:1:size(Motif,2);

WarpIndex = 1;
for i = [0 5 10 15],
    xx = linspace(x(1), x(end), (size(Motif,2) * (1 + i/100)));
    clear WMotif;
    for j = 1:size(Motif, 1);
        WMotif(j,:) = spline(x, Motif(j,:), xx);
    end
    DSeqMatchValues{WarpIndex} = CalculateMotifSeqMatch(DirectoryName, 'DirSongFileList.txt', FileType, WMotif, 'off');
    USeqMatchValues{WarpIndex} = CalculateMotifSeqMatch(DirectoryName, 'UnDirSongFileList.txt', FileType, WMotif, 'off');
    WarpIndex = WarpIndex + 1;
end

save(OutputFileName, 'DSeqMatchValues', 'USeqMatchValues');