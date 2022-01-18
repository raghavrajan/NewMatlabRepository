function [] = Harini_MonotonePlaybackAnalysis(CSVTextFile)

BackgroundNoiseLevel = 46.6; % in dB

[BirdParameters, Flag] = ProcessSongData(CSVTextFile, 2000);

% Correct song amplitude by subtracting the background noise level
% mentioned abvoe
for i = 1:length(BirdParameters),
    AmplitudeColIndex = strmatch(BirdParameters(i).SAPFeat_FieldNames, 'LogAmplitude', 'exact');
    BirdParameters(i).SAPFeatsMatrix(:,AmplitudeColIndex) = BirdParameters(i).SAPFeatsMatrix(:,AmplitudeColIndex) + 70 - BackgroundNoiseLevel;
end

if (sum(Flag) > 0)
    return;
end

Colours = 'rgbcmk';
Symbols = '+o<sd';

% === Some things that are common across all data =========================
% Get unique BirdNames to get an idea of how many birds there are
for i = 1:length(BirdParameters),
    BirdNames{i} = BirdParameters(i).BirdName;
end
CurrentIndex = 1;
while (CurrentIndex < length(BirdNames)),
    ToRemove = [];
    for j = (CurrentIndex+1):length(BirdNames),
        if (strcmp(BirdNames{CurrentIndex}, BirdNames{j}))
            ToRemove(end+1) = j;
        end
    end
    BirdNames(ToRemove) = [];
    CurrentIndex = CurrentIndex + 1;
end

% For each unique bird, get the indices of that bird's data from the main
% structure BirdParameters

for i = 1:length(BirdNames),
    BirdNameIndices{i} = [];
    for j = 1:length(BirdParameters),
        if (strcmp(BirdParameters(j).BirdName, BirdNames{i}))
            BirdNameIndices{i}(end+1) = j;
        end
    end
end

% =============== File times and syllable times across entire day =========
% Now to get the times of all syllables
% First we need to get file times and this will automatically give syllable
% times
for i = 1:length(BirdParameters),
    for j = 1:length(BirdParameters(i).SongFileNames),
        FileTimeString = BirdParameters(i).SongFileNames{j}(end-9:end-4);
        BirdParameters(i).FileTime(j) = str2double(FileTimeString(1:2)) + str2double(FileTimeString(3:4))/60 + str2double(FileTimeString(5:6))/3600;
    end
    BirdParameters(i).FileTime = BirdParameters(i).FileTime(:);
    % To get time of each syllable, I will take the onset time of the file
    % and add sylllable onset time to it - file onset time is in hours and
    % syllable onset time is ms - have to convert properly
    BirdParameters(i).SyllableTimes = BirdParameters(i).FileTime(BirdParameters(i).SyllableData(:,2)) + BirdParameters(i).SyllableData(:,4)/(1000*3600);
    BirdParameters(i).SyllableTimes = BirdParameters(i).SyllableTimes(:);
    % The day is divided into 2 hour chunks - 6-8; 8-10; 10-12; 12-14;
    % 14-16; 16-18; 18-20;
    % And I have to put the syllable times into one of these categories
    % So will floor all syllable time values and then subtract the starting
    % time which is 6 and divide by the duration of each chunk - in this
    % case 2. Will floor this and should get an integer which corresponds
    % to the category in which the time falls
    BirdParameters(i).SyllableTimeCategories = floor((floor(BirdParameters(i).SyllableTimes) - 6)/2);
    BirdParameters(i).SyllableTimeCategories = BirdParameters(i).SyllableTimeCategories(:);
end

% ========== Dates of recording ===========================================

% First get the dates from data label - Data label has ddmmyy format for
% the date as a string
for i = 1:length(BirdParameters),
    BirdParameters(i).Day = str2double(BirdParameters(i).DataLabel(1:2));
    BirdParameters(i).Month = str2double(BirdParameters(i).DataLabel(3:4));
    BirdParameters(i).Year = str2double(BirdParameters(i).DataLabel(5:6));
end

% Now to sort the data for each bird separately, first based on year, then
% month, then day
for i = 1:length(BirdNames),
    Dates = [[BirdParameters(BirdNameIndices{i}).Year]' [BirdParameters(BirdNameIndices{i}).Month]' [BirdParameters(BirdNameIndices{i}).Day]'];
    [SortedVals, SortedIndices] = sortrows(Dates);
    SortedBirdNameIndices{i} = BirdNameIndices{i}(SortedIndices);
    BirdNameIndices{i} = SortedBirdNameIndices{i};
end

% I will also make a new structure called IndividualBirds that will now
% have each individual birds data. Within this, there will be a field
% called SortedBirdParameters that will now contain all the relevant data
% for that bird
for i = 1:length(BirdNames),
    IndividualBirds(i).SortedBirdParameters = BirdParameters(SortedBirdNameIndices{i});
end

% Now to go through the sorted list and get individual days and number the
% days in a column from 1 onwards for each bird

for i = 1:length(IndividualBirds),
    Index = 1;
    DayIndex = 1;
    TempBirdParameters = IndividualBirds(i).SortedBirdParameters;
    while (Index <= length(TempBirdParameters))
        MatchingDays = find(([TempBirdParameters.Day] == TempBirdParameters(Index).Day) & ([TempBirdParameters.Month] == TempBirdParameters(Index).Month) & ([TempBirdParameters.Year] == TempBirdParameters(Index).Year));
        for j = MatchingDays(:)',
            IndividualBirds(i).SortedBirdParameters(j).RecordingDayIndex = DayIndex;
        end
        IndividualBirds(i).RecordingDays{DayIndex} = TempBirdParameters(Index).DataLabel;
        IndividualBirds(i).RecordingDates(DayIndex,:) = [TempBirdParameters(Index).Day TempBirdParameters(Index).Month TempBirdParameters(Index).Year];
        DayIndex = DayIndex + 1;
        Index = MatchingDays(end) + 1;
    end
end

% ======= One last thing - put all the syllables together in one large
% vector with all the conditions too, separately for each bird ============

for i = 1:length(IndividualBirds),
    IndividualBirds(i).AllSyllableData = [];
    IndividualBirds(i).AllSyllableFeatValues = [];
    IndividualBirds(i).AllRecordingDayIndices = [];
    IndividualBirds(i).AllSyllableTimeCategories = [];
    IndividualBirds(i).SongFileNames = [];
    IndividualBirds(i).FileTime = [];
    IndividualBirds(i).Bouts = [];
    IndividualBirds(i).BoutTimesRelativeToStart = [];
    TotalNumSongFiles = 0;
    TotalNumSyllables = 0;
    
    for j = 1:length(IndividualBirds(i).SortedBirdParameters),
        TempBouts = IndividualBirds(i).SortedBirdParameters(j).Bouts;
        % The above matrix has information about bouts. The first and
        % second column have the syllable #s for the onset and offset of a
        % bout. This has to incremented by the total number of syllables
        % each time I go to a new filelist. The syllables #s are for within
        % a bout. Therefore, to increment this properly, I have to first
        % make this a cumulative # of syllables across files
        % The 3rd and 4th column have the onset and offset file #s for the
        % bouts and this has to be incremented by the total number of song
        % files in the filelist
        
        % But before I do this will get times for all bout relative to
        % female introduction time
        % Now also add the bout time relative to female introduction time
        TempBoutTimesRelativeToStart = (TempBouts(:,5)/(1000*3600)) + IndividualBirds(i).SortedBirdParameters(j).FileTime(IndividualBirds(i).SortedBirdParameters(j).Bouts(:,3));
        % The first part (bout onset time) in the right-side of above
        % expression is in ms and has to be converted to hours which is the
        % unit for the file time and female introduction times
        
        FemaleIntroductionTime = str2double(IndividualBirds(i).SortedBirdParameters(j).Femaleintroductiontime(1:2)) + str2double(IndividualBirds(i).SortedBirdParameters(j).Femaleintroductiontime(3:4))/60 + str2double(IndividualBirds(i).SortedBirdParameters(j).Femaleintroductiontime(5:6))/3600; 
        TempBoutTimesRelativeToStart = TempBoutTimesRelativeToStart(:) - FemaleIntroductionTime;
        
        IndividualBirds(i).BoutTimesRelativeToStart = [IndividualBirds(i).BoutTimesRelativeToStart(:); TempBoutTimesRelativeToStart];
        
        for k = 1:length(IndividualBirds(i).SortedBirdParameters(j).SongFileNames),
            IndividualBirds(i).SongFileNames{end+1} = fullfile(IndividualBirds(i).SortedBirdParameters(j).DataDirectory, IndividualBirds(i).SortedBirdParameters(j).SongFileNames{k});
            OtherFileBouts = find(TempBouts(:,3) > k);
            NumSylls = length(IndividualBirds(i).SortedBirdParameters(j).NoteInfo{k}.labels);
            TempBouts(OtherFileBouts,1:2) = TempBouts(OtherFileBouts,1:2) + NumSylls;
        end
        
        % Add the file times for each file into this now
        IndividualBirds(i).FileTime = [IndividualBirds(i).FileTime(:); IndividualBirds(i).SortedBirdParameters(j).FileTime(:)];
        
        TempSyllableData = IndividualBirds(i).SortedBirdParameters(j).SyllableData;
        % The second and third column in the above matrix has onset and
        % offset file #s for each syllable. Since I'm putting together data
        % across many different filelists, I'm going to keep a counter that
        % keeps track of total # of files and then add it to the second and
        % third column of this matrix
        
        TempSyllableData(:,2:3) = TempSyllableData(:,2:3) + TotalNumSongFiles;
        IndividualBirds(i).AllSyllableData = [IndividualBirds(i).AllSyllableData; TempSyllableData];

        TempBouts(:,1:2) = TempBouts(:,1:2) + TotalNumSyllables;
        TempBouts(:,3:4) = TempBouts(:,3:4) + TotalNumSongFiles;
        
        IndividualBirds(i).Bouts = [IndividualBirds(i).Bouts; TempBouts];
        
        % Now to increase the value of the counter of total number of song
        % files and syllables
        TotalNumSongFiles = TotalNumSongFiles + length(IndividualBirds(i).SortedBirdParameters(j).SongFileNames);
        TotalNumSyllables = TotalNumSyllables + size(TempSyllableData,1);
        
        IndividualBirds(i).AllSyllableFeatValues = [IndividualBirds(i).AllSyllableFeatValues; IndividualBirds(i).SortedBirdParameters(j).SAPFeatsMatrix];
        IndividualBirds(i).AllRecordingDayIndices = [IndividualBirds(i).AllRecordingDayIndices; ones(size(IndividualBirds(i).SortedBirdParameters(j).SyllableData,1), 1)*IndividualBirds(i).SortedBirdParameters(j).RecordingDayIndex];
        IndividualBirds(i).AllSyllableTimeCategories = [IndividualBirds(i).AllSyllableTimeCategories; IndividualBirds(i).SortedBirdParameters(j).SyllableTimeCategories];
    end
end

AmplitudeColIndex = strmatch(IndividualBirds(1).SortedBirdParameters(1).SAPFeat_FieldNames, 'LogAmplitude', 'exact');
AmplitudeColIndex = 2;
for i = 1:length(IndividualBirds),
    figure;
    hold on;
    plot(IndividualBirds(i).AllRecordingDayIndices, IndividualBirds(i).AllSyllableFeatValues(:,AmplitudeColIndex), 'ko', 'MarkerSize', 4);
    for j = 1:length(IndividualBirds(i).RecordingDays),
        DayIndices = find((IndividualBirds(i).AllRecordingDayIndices == j) & (~isnan(IndividualBirds(i).AllSyllableFeatValues(:, AmplitudeColIndex))));
        errorbar(j+0.1,nanmean(IndividualBirds(i).AllSyllableFeatValues(DayIndices, AmplitudeColIndex)), nanstd(IndividualBirds(i).AllSyllableFeatValues(DayIndices, AmplitudeColIndex))/sqrt(length(DayIndices)), 'ks-', 'MarkerSize', 8, 'LineWidth', 1.5);
    end
    set(gca, 'XTick', 1:1:length(IndividualBirds(i).RecordingDays), 'XTickLabel', IndividualBirds(i).RecordingDays, 'XTickLabelRotation', 45);
    title([BirdNames{i}, ': Monotone playback']);
end
disp('Finished plotting');