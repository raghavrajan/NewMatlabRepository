function [] = WriteMaxMatchFileList(BirdName, Date, MaxBouts, BoutStatsFile, Condition)

PresentDir = pwd;

clear MaxMatches SongFileNames;
Temp = textread(['/home/raghav/BoutStatistics/', BoutStatsFile], '%s', 'delimiter', '\n');
MaxMatchLines = find(cellfun(@length, strfind(Temp, 'Maximum')));
FileNames = find(cellfun(@length, strfind(Temp, 'Songfile name')));
MatchValIndex = 1;
for i = 1:length(MaxMatchLines),
    [s, s1, s2, s3, Max1, s4, s5, Max2] = strread(Temp{MaxMatchLines(i)}, '%s %s %s %s %f %s %s %f');
    MaxMatches(MatchValIndex) = Max1;
    SongFileIndex = find(FileNames < MaxMatchLines(i), 1, 'last');
    Strings = strread(Temp{FileNames(SongFileIndex)}, '%s');
    SongFileNames{MatchValIndex} = Strings{end};
    MatchValIndex = MatchValIndex + 1;
    clear Strings Max1;
end

[SortedVals, SortedIndices] = sort(MaxMatches);
Fid = fopen(['/net/doupe3/raghav/HVC_Microlesions/', BirdName, '/', BirdName, '_', Date, '_', Condition, '_maxMatchFiles.txt'], 'w');
if (length(SortedIndices) < MaxBouts)
    for i = 1:length(SortedIndices),
        fprintf(Fid, '%s\n', SongFileNames{SortedIndices(end - i + 1)});
    end
else
    for i = 1:MaxBouts,
        fprintf(Fid, '%s\n', SongFileNames{SortedIndices(end - (i - 1))});
    end
end
fclose(Fid);

cd(PresentDir);