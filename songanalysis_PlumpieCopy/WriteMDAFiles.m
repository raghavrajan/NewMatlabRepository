function [] = WriteMDAFiles(DataDir, FileList, FileType, ChanNos)

Fid = fopen(FileList, 'r');
Files = textscan(Fid, '%s', 'DeLimiter', '\n');
Files = Files{1};
fclose(Fid);

for i = 1:length(Files),
    fprintf('%i >> ', i);
    [RawData, Fs] = GetData(DataDir, Files{i}, FileType, 1);
    FileLen(i) = length(RawData)/Fs;
end

FilesRemaining = Files;
FilesRemainingLen = FileLen;

WritingFileIndex = 0;
Index = 0;
while (~isempty(FilesRemaining))
    CumulativeFileLen = cumsum(FilesRemainingLen);
    FileIndex = find(CumulativeFileLen < (800), 1, 'last');
    Data = [];
    for i = 1:FileIndex,
        WritingFileIndex = WritingFileIndex + 1;
        fprintf('%i >> ', WritingFileIndex);
        TempData = [];
        for j = ChanNos(:)',
            [RawData, Fs] = GetData(DataDir, FilesRemaining{i}, FileType, j);
            TempData = [TempData; RawData(:)'];
        end
        Data = [Data TempData];
    end
    Index = Index + 1;
    Fid = fopen([FileList, '.', num2str(Index), '_txt'], 'w');
    for i = 1:FileIndex,
        fprintf(Fid, '%s\n', FilesRemaining{i});
    end
    fclose(Fid);
    FilesRemaining(1:FileIndex) = [];
    FilesRemainingLen(1:FileIndex) = [];
    writemda32(Data, [FileList, '.', num2str(Index), '_txt.mda']);
end
