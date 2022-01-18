function [] = WriteTimesWithSyllLabels(DataDir, FileList, FileType, FileTimesFile, Mode)

% This is a script to write out all the syllable labels and the time when
% the syllable occured.
% Inputs:
%   DataDir - directory where raw data is
%   FileList - list of song files to be considered
%   FileType - type of file
%   FileTimesFile - a file with the times of each file
%   Mode - either absolute or relative - relative to the first syllable
%   time


%% First read in all the file names and get the times of the files
Fid = fopen(fullfile(DataDir, FileList), 'r');
SongFiles = textscan(Fid, '%s', 'DeLimiter', '\n');
SongFiles = SongFiles{1};
fclose(Fid);

Fid = fopen(fullfile(DataDir, FileTimesFile), 'r');
FileTimesStrings = textscan(Fid, '%s', 'DeLimiter', '\n');
FileTimesStrings = FileTimesStrings{1};
fclose(Fid);

for i = 1:length(SongFiles),
    FileTimeStringIndex = find(contains(FileTimesStrings, SongFiles{i}));
    if (~isempty(FileTimeStringIndex))
        TempIndex = strfind(FileTimesStrings{FileTimeStringIndex}, 'Duration');
        FileTimes(i) = 1000 * (str2double(FileTimesStrings{FileTimeStringIndex}(TempIndex-15:TempIndex-14))*3600 + str2double(FileTimesStrings{FileTimeStringIndex}(TempIndex-12:TempIndex-11))*60 + str2double(FileTimesStrings{FileTimeStringIndex}(TempIndex-9:TempIndex-8)));
    end
end

if ~isempty(find(strcmp('relative', Mode)))
    FileTimes = FileTimes - FileTimes(1);
end

%% Now read in the notes files and write out a csv file with separator ;
% The file should have the syllable label and the time of occurence
Fid = fopen(fullfile(DataDir, [FileList, '.Labels.Times.csv']), 'w');
% Write header with ; as separator
fprintf(Fid, '#;Syll Label;Syll Onset Time (ms); Syll Offset Time (ms)\n');
SyllIndex = 0;
for i = 1:length(SongFiles),
    Temp = load(fullfile(DataDir, 'ASSLNoteFiles', [SongFiles{i}, '.not.mat']));
    for j = 1:length(Temp.onsets),
        SyllIndex = SyllIndex + 1;
        fprintf(Fid, '%i;%c;%8.3f;%8.3f\n', SyllIndex, Temp.labels(j), FileTimes(i) + (Temp.onsets(j)), FileTimes(i) + (Temp.offsets(j)));
    end
end
disp('Finished');