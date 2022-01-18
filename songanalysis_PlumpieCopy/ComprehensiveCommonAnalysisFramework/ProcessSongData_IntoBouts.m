function [BirdParameters, Flag] = ProcessSongData_IntoBouts(CSVTextFile, InterboutInterval, SavedDataDir)

%% First step parse the CSV text file and load up filenames and note data
% First get details from the CSV text file
disp(['Getting header data from CSV file ', CSVTextFile, ' ...']);
[HeaderLine, BirdDetails] = LSINA_GetDetailsFromCSVFile(CSVTextFile);

% Now parse all the lines into the appropriate variables based on the
% header line
disp('Getting data from CSV file ...');
[BirdParameters] = LSINA_ParseHeaderBirdData(HeaderLine, BirdDetails);

% Now put in the input inter-bout interval
for i = 1:length(BirdParameters),
    BirdParameters(i).Interboutinterval = InterboutInterval;
end

% Now for each of the birds, load up all the filenames
disp('Loading up filenames ...');
for i = 1:length(BirdParameters),
    FileNames_FileName = fullfile(SavedDataDir, [BirdParameters(i).BirdName, '.', BirdParameters(i).UndirSongFileList, '.', num2str(BirdParameters(i).Interboutinterval), '.FileNames.mat']);
    if (exist(FileNames_FileName, 'file'))
        load(FileNames_FileName, 'SongFileNames');
        BirdParameters(i).SongFileNames = SongFileNames;
    else
        [BirdParameters(i).SongFileNames] = LSINA_GetDataFileNames(BirdParameters(i));
        SongFileNames = BirdParameters(i).SongFileNames;
        save(FileNames_FileName, 'SongFileNames');
    end
end

% Now check if notes files have unequal numbers of labels, onsets and
% offsets
disp('Checking note data ...');
Flag = zeros(length(BirdParameters), 1);
for i = 1:length(BirdParameters),
    if (exist(BirdParameters(i).DataDirectory, 'dir'))
        if (isfield(BirdParameters(i), 'UndirSongFileList'))
            Flag(i) = CheckLengthsOnsetsOffsetsLabels(BirdParameters(i).DataDirectory, BirdParameters(i).UndirSongFileList);
        else
            Flag(i) = CheckLengthsOnsetsOffsetsLabels(BirdParameters(i).DataDirectory, BirdParameters(i).SongFileList);
        end
    end
end

% If there are files with unequal numbers, then don't run the script - just
% return
if (sum(Flag) > 0)
    return;
else
    disp('All notes files are ok');
end


% Now load up the note files and the length of each file
disp('Loading up note data ...');
for i = 1:length(BirdParameters),
    Notes_FileName = fullfile(SavedDataDir, [BirdParameters(i).BirdName, '.', BirdParameters(i).UndirSongFileList, '.', num2str(BirdParameters(i).Interboutinterval), '.Notes.mat']);
    if (exist(Notes_FileName, 'file'))
        load(Notes_FileName, 'NoteInfo', 'FileLen');
        BirdParameters(i).NoteInfo = NoteInfo;
        BirdParameters(i).FileLen = FileLen;
    else
        [BirdParameters(i).NoteInfo, BirdParameters(i).FileLen] = LSINA_LoadNoteFileInfo(BirdParameters(i));
        NoteInfo = BirdParameters(i).NoteInfo;
        FileLen = BirdParameters(i).FileLen;
        save(Notes_FileName, 'NoteInfo', 'FileLen');
    end
end

%% Now divided up data into bouts using the inter-bout interval specified
disp('Identifying bouts ...');
for i = 1:length(BirdParameters),
    fprintf('%d >> ', i);
    BoutInfo_FileName = fullfile(SavedDataDir, [BirdParameters(i).BirdName, '.', BirdParameters(i).UndirSongFileList, '.', num2str(BirdParameters(i).Interboutinterval), '.BoutInfo.mat']);
    if (exist(BoutInfo_FileName, 'file'))
        load(BoutInfo_FileName, 'Bouts', 'BoutDetails');
        BirdParameters(i).Bouts = Bouts;
        BirdParameters(i).BoutDetails = BoutDetails;
    else
        [BirdParameters(i).Bouts, Gaps, BirdParameters(i).BoutDetails] = LSINA_GetBoutInfo(BirdParameters(i));
        Bouts = BirdParameters(i).Bouts;
        BoutDetails = BirdParameters(i).BoutDetails;
        save(BoutInfo_FileName, 'Bouts', 'BoutDetails');
    end
    fprintf('%d >> ', i);
end
fprintf('\n');

% Now make a BoutDetails field for the structure that has the labels,
% onsets and offsets for bouts, depending on whether it is continuous or
% not

% Now to write the bout information into a text file - I want to write bout
% #, filename, bout onset and offset in sec to a text file
disp('Writing bout info to text file ...');
for i = 1:length(BirdParameters),
    fprintf('%d >> ', i);
    Fid = fopen(fullfile(SavedDataDir, [BirdParameters(i).BirdName, '.', BirdParameters(i).UndirSongFileList, '.', num2str(BirdParameters(i).Interboutinterval), '.BoutInfo.txt']), 'w');
    fprintf(Fid, 'Bout #\tFileName\tBout onset (sec)\tBout offset (sec)\n');
    for j = 1:size(BirdParameters(i).Bouts, 1),
        if ((BirdParameters(i).Bouts(j,7) == 1) && (BirdParameters(i).Bouts(j,8) > 0) && (BirdParameters(i).Bouts(j,9) > 1))
            fprintf(Fid, '%d\t%s\t%3.5f\t%3.5f\n', j, BirdParameters(i).SongFileNames{BirdParameters(i).Bouts(j,3)}, BirdParameters(i).Bouts(j,5)/1000, BirdParameters(i).Bouts(j,6)/1000);
        end
    end
    fclose(Fid);
end
fprintf('\n');

%% Now for all data, get the time of each bout
for i = 1:length(BirdParameters),
    % For this, I need to first get the times for each of the files
    for j = 1:length(BirdParameters(i).SongFileNames),
        switch BirdParameters(i).FileType
            case 'wav'
                ExtensionIndex = strfind(BirdParameters(i).SongFileNames{j}, '.wav');
                
            case 'okrank'
                ExtensionIndex = length(BirdParameters(i).SongFileNames{j}) + 1;
                
            case 'obs'
                ExtensionIndex = strfind(BirdParameters(i).SongFileNames{j}, '.cbin');
        end
        FileTimeString = BirdParameters(i).SongFileNames{j}((ExtensionIndex - 6):(ExtensionIndex - 1)); % string from filename in the format hhmmss (h for hour, m for minutes, s for seconds)
        BirdParameters(i).FileTime(j) = (str2double(FileTimeString(1:2))) + (str2double(FileTimeString(3:4))/60) + (str2double(FileTimeString(5:6))/3600); % in hours\
    end
    BirdParameters(i).FileTime = BirdParameters(i).FileTime(:);
    BirdParameters(i).BoutOnsetTimes_WithinDay = BirdParameters(i).FileTime(BirdParameters(i).Bouts(:,3)) + (BirdParameters(i).Bouts(:,5)/(1000 * 3600)); % in hours
end

% Now for all data, I have to get the recording day index

for i = 1:length(BirdParameters),
    BirdNames{i} = BirdParameters(i).BirdName;
    BirdParameters(i).RecordingDayIndex = NaN;
    BirdParameters(i).SessionNumWithinDay = NaN;
    BirdParameters(i).OverallSessionNum = NaN;
end
UniqueBirds = unique(BirdNames);

for i = 1:length(UniqueBirds),
    Matches = find(strcmp(UniqueBirds{i}, BirdNames));
    for j = 1:length(Matches),
        Date{i}{j} = BirdParameters(Matches(j)).DataLabel;
        DateNumber{i}(j) = datenum(Date{i}{j}, 'ddmmyy');
        SessionOnsetTime{i}(j) = BirdParameters(Matches(j)).FileTime(1);
    end
    for j = 1:length(Matches),
        BirdParameters(Matches(j)).RecordingDayIndex = DateNumber{i}(j) - min(DateNumber{i}) + 1;
    end
    UniqueRecordingDayNumbers = unique(DateNumber{i});
    for j = 1:length(UniqueRecordingDayNumbers),
        DayMatches = find(DateNumber{i} == UniqueRecordingDayNumbers(j));
        DaySessionOnsetTimes = SessionOnsetTime{i}(DayMatches);
        [SortedVals, SortedIndices] = sort(DaySessionOnsetTimes);
        WithinDaySessionIndices{i}(DayMatches(SortedIndices)) = 1:1:length(SortedIndices);
    end
    for j = 1:length(Matches),
        BirdParameters(Matches(j)).SessionNumWithinDay = WithinDaySessionIndices{i}(j);
    end
    [SortedMatches, SortedMatchIndices] = sortrows([DateNumber{i}(:) SessionOnsetTime{i}(:)]);
    for j = 1:length(SortedMatchIndices),
        BirdParameters(Matches(SortedMatchIndices(j))).OverallSessionNum = j;
    end
end


disp('Finished Analysis');