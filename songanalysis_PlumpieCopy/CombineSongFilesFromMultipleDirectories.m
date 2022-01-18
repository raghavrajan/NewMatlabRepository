function [] = CombineSongFilesFromMultipleDirectories(DirList, FileLists, OutputFileName, varargin)

if (nargin > 3)
    CSVFile = varargin{1};
    Fid = fopen(CSVFile, 'r');
    Text = textscan(Fid, '%s', 'DeLimiter', '\n');
    Text = Text{1};
    fclose(Fid);
    for i = 1:length(Text),
        TempText = textscan(Text{i}, '%s', 'DeLimiter', ',');
        DirList{i} = TempText{1}{1};
        FileLists{i} = TempText{1}{2};
    end
end

PresentDir = pwd;

[OutputFileDir, OutputFileNameOnly, OutputFileExt] = fileparts(OutputFileName);
if (isempty(OutputFileDir))
    OutputFileName = fullfile(PresentDir, OutputFileName);
end

OutputFid = fopen(OutputFileName, 'w');

for i = 1:length(DirList),
    Fid = fopen(fullfile(DirList{i}, FileLists{i}), 'r');
    Files = textscan(Fid, '%s', 'DeLimiter', '\n');
    Files = Files{1};
    for j = 1:length(Files),
        fprintf(OutputFid, '%s\n', fullfile(DirList{i}, Files{j}));
    end
    fclose(Fid);
end
fclose(OutputFid);
    