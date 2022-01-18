function [] = HariniPlots_ScriptFile(CSVTextFile, PlotOption, BirdOption, varargin)

% ================================= Help ==================================
% 'Plot Option Strings:');
%     disp('1. DayWiseFeaturePlot - plots the feature for each session in chronological order of days');
%     disp('2. ConditionWiseFeaturePlot - plots the feature for each session in order of conditions');
%     disp('3. TimeOfDayPlot - plots the feature for all sessions based on time of day');
%     disp('4. ChangeWithinSameDay - plots the change within sessions that occured on the same day');
%     disp('5. AverageAcrossConditions - plots the average feature values for each condition averaged over all sessions of a particular condition - plots for individual birds and the group');
%     disp('6. VariabilityForDifferentGroupings - plots the variability after grouping the data into different categories');
%     disp('7. BoutLength - plots the bout length for different sessions in chronological order');
%     disp('8. NumMotifs - plots the number of motifs for different sessions in chronological order');
%     disp('9. BoutStatistics - plots the number of INs, number of motifs, bout length, motif duration, etc. ');
%     disp('10. ProportionDirectedSongsAtEachLocation - plots the proportion of directed songs at each location');
%     disp('11. PlotDirSongTimes - plots the times of directed songs at each location relative to female introduction time');
%     disp('12. PlotINNumber_DirOnlySongs - plots the mean number of INs at each location for only directed songs');
%     disp('13. PlotMotifNumber_DirOnlySongs - plots the mean number of motifs at each location for only directed songs');
%     disp('14. PlotRawWaveforms_LogAmplitudes - plots the raw waveform and
%     log amplitude');
%
%     disp('Bird Option Strings:');
%     disp('All - for all birds');
%     disp('Otherwise, specific bird name for plots for only that bird');
% end
  
BackgroundNoiseLevel = 46.6; % in dB

if (strcmp(PlotOption, 'help'))
    disp('Usage: ');
    eval(['dbtype HariniPlots_ScriptFile 1']);
    disp('Plot Option Strings:');
    disp('1. DayWisePlot - plots the feature for each session in chronological order of days');
    disp('2. ConditionWisePlot - plots the feature for each session in order of conditions');
    disp('3. TimeOfDayPlot - plots the feature for all sessions based on time of day');
    disp('4. ChangeWithinSameDay - plots the change within sessions that occured on the same day');
    disp('5. AverageAcrossConditions - plots the average feature values for each condition averaged over all sessions of a particular condition - plots for individual birds and the group');
    disp('6. VariabilityForDifferentGroupings - plots the variability after grouping the data into different categories');
    disp('7. BoutLength - plots the bout length for different sessions in chronological order');
    disp('8. NumMotifs - plots the number of motifs for different sessions in chronological order');
    disp('9. BoutStatistics - plots the number of INs, number of motifs, bout length, motif duration, etc. ');
    disp('10. ProportionDirectedSongsAtEachLocation - plots the proportion of directed songs at each location');
    disp('11. PlotDirSongTimes - plots the times of directed songs at each location relative to female introduction time');
    disp('12. PlotINNumber_DirOnlySongs - plots the mean number of INs at each location for only directed songs');
    disp('13. PlotMotifNumber_DirOnlySongs - plots the mean number of motifs at each location for only directed songs');
    disp('14. PlotRawWaveforms_LogAmplitudes - plots the raw waveform and log amplitude');
    disp('Bird Option Strings:');
    disp('All - for all birds');
    disp('Otherwise, specific bird name for plots for only that bird');
    
    disp('Varargin - add video scoring csv file if present');
    return;
end

[BirdParameters, Flag] = ProcessSongData(CSVTextFile, 1000);

% One more pre-processing step. I currently have a list of bouts and I also
% have a separate list of syllables that is compiled across all files. I
% have to have a set of indices that allow me to know which bout each
% syllable in the list of syllables belongs to
for i = 1:length(BirdParameters),
    BirdParameters(i).SyllableListBoutNum = zeros(size(BirdParameters(i).SyllableData,1), 1);
    TotalSyllNo = 0;
    for j = 1:length(BirdParameters(i).SongFileNames),
        BoutIndices = find(BirdParameters(i).Bouts(:,3) == j);
        for k = 1:length(BoutIndices),
            BirdParameters(i).SyllableListBoutNum((TotalSyllNo + BirdParameters(i).Bouts(BoutIndices(k), 1)):(TotalSyllNo + BirdParameters(i).Bouts(BoutIndices(k), 2))) = BoutIndices(k);
        end
        TotalSyllNo = TotalSyllNo + length(BirdParameters(i).NoteInfo{j}.labels);
    end
end

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

if (nargin > 3)
    VideoScoringCSVTextFile = varargin{1};
    % First get details from the CSV text file
    disp('Getting header data from CSV file ...');  
    [HeaderLine, BirdDetails] = LSINA_GetDetailsFromCSVFile(VideoScoringCSVTextFile);

    % Now parse all the lines into the appropriate variables based on the
    % header line
    disp('Getting data from CSV file ...');
    [VideoScoringData] = LSINA_ParseHeaderBirdData(HeaderLine, BirdDetails);
else
    VideoScoringData = [];
end

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

% ============= Conditions ================================================
% Now to get the different conditions and convert this to an array of 1s
% and 0s for each bird
for i = 1:length(BirdNames),
    TempBirdParameters = IndividualBirds(i).SortedBirdParameters;
    Conditions = [];
    for j = 1:length(TempBirdParameters),
        Conditions{j} = TempBirdParameters(j).Condition;
    end
    
    ToRemove = [];
    CurrentIndex = 1;
    while (CurrentIndex < length(Conditions)),
        ToRemove = [];
        for j = (CurrentIndex+1):length(Conditions),
            if (strcmp(Conditions{CurrentIndex}, Conditions{j}))
                ToRemove(end+1) = j;
            end
        end
        Conditions(ToRemove) = [];
        CurrentIndex = CurrentIndex + 1;
    end
    
    IndividualBirds(i).Conditions = Conditions;
    IndividualBirds(i).Conditions = mat2cell(sortrows(reshape(cell2mat(IndividualBirds(i).Conditions), length(IndividualBirds(i).Conditions{1}), length(IndividualBirds(i).Conditions))'), ones(length(IndividualBirds(i).Conditions), 1), length(IndividualBirds(i).Conditions{1}));
    
    for j = 1:length(IndividualBirds(i).SortedBirdParameters),
        for k = 1:length(IndividualBirds(i).Conditions),
            if (strcmp(IndividualBirds(i).SortedBirdParameters(j).Condition, IndividualBirds(i).Conditions{k}))
                switch (IndividualBirds(i).Conditions{k})
                    case 'L0'
                        IndividualBirds(i).SortedBirdParameters(j).ConditionIndices = 1;
                        
                    case 'L1'
                        IndividualBirds(i).SortedBirdParameters(j).ConditionIndices = 2;
                        
                    case 'L2'
                        IndividualBirds(i).SortedBirdParameters(j).ConditionIndices = 3;
                        
                    case 'L3'
                        IndividualBirds(i).SortedBirdParameters(j).ConditionIndices = 4;    
                        
                    case 'L4'
                        IndividualBirds(i).SortedBirdParameters(j).ConditionIndices = 5;    
                        
                    case 'UN'
                        IndividualBirds(i).SortedBirdParameters(j).ConditionIndices = 6;
                        
                    otherwise
                        IndividualBirds(i).SortedBirdParameters(j).ConditionIndices = 1;
                end
                break;
            end
        end
    end
end

% ============= Microphone ================================================
% Now to get the different microphone conditions and convert this to an array of 1s
% and 0s for each bird

for i = 1:length(BirdNames),
    TempBirdParameters = IndividualBirds(i).SortedBirdParameters;
    Microphone = [];
    for j = 1:length(TempBirdParameters),
        Microphone{j} = TempBirdParameters(j).Microphone;
    end
    
    ToRemove = [];
    CurrentIndex = 1;
    while (CurrentIndex < length(Microphone)),
        ToRemove = [];
        for j = (CurrentIndex+1):length(Microphone),
            if (strcmp(Microphone{CurrentIndex}, Microphone{j}))
                ToRemove(end+1) = j;
            end
        end
        Microphone(ToRemove) = [];
        CurrentIndex = CurrentIndex + 1;
    end
    
    IndividualBirds(i).Microphone = Microphone;
    
    for j = 1:length(IndividualBirds(i).SortedBirdParameters),
        for k = 1:length(IndividualBirds(i).Microphone),
            if (strcmp(IndividualBirds(i).SortedBirdParameters(j).Microphone, IndividualBirds(i).Microphone{k}))
                IndividualBirds(i).SortedBirdParameters(j).MicrophoneIndices = k;
                break;
            end
        end
    end
end

% ======= One last thing - put all the syllables together in one large
% vector with all the conditions too, separately for each bird ============

for i = 1:length(IndividualBirds),
    IndividualBirds(i).AllFFTLogAmplitudes = [];
    IndividualBirds(i).AllSyllableLogAmpStatus = [];
    IndividualBirds(i).AllSyllableLogAmpBaselineAmpValue = [];
    IndividualBirds(i).AllSyllableData = [];
    IndividualBirds(i).AllSyllableFeatValues = [];
    IndividualBirds(i).AllSyllableLogAmplitudeKao = [];
    IndividualBirds(i).AllSyllableLogAmplitudeRMS = [];
    IndividualBirds(i).AllSyllableFFAutocorrKao = [];
    IndividualBirds(i).AllConditionIndices = [];
    IndividualBirds(i).AllMicrophoneIndices = [];
    IndividualBirds(i).AllRecordingDayIndices = [];
    IndividualBirds(i).AllSyllableTimeCategories = [];
    IndividualBirds(i).SongFileNames = [];
    IndividualBirds(i).FileTime = [];
    IndividualBirds(i).Bouts = [];
    IndividualBirds(i).BoutTimesRelativeToStart = [];
    IndividualBirds(i).AllSyllableTemplateMatchValues = [];
    IndividualBirds(i).AllSyllableListBoutNum = [];
    TotalNumSongFiles = 0;
    TotalNumSyllables = 0;
    TotalBoutNum = 0;
    
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
        
        for k = 1:length(IndividualBirds(i).SortedBirdParameters(j).FFTLogAmplitudes),
            IndividualBirds(i).AllFFTLogAmplitudes = [IndividualBirds(i).AllFFTLogAmplitudes; IndividualBirds(i).SortedBirdParameters(j).FFTLogAmplitudes{k}(:)];
            IndividualBirds(i).AllSyllableLogAmpStatus = [IndividualBirds(i).AllSyllableLogAmpStatus; IndividualBirds(i).SortedBirdParameters(j).SyllableStatus{k}(:)];
            IndividualBirds(i).AllSyllableLogAmpBaselineAmpValue = [IndividualBirds(i).AllSyllableLogAmpBaselineAmpValue; IndividualBirds(i).SortedBirdParameters(j).BaselineAmpValue{k}(:)];
        end
%        IndividualBirds(i).AllSyllableTemplateMatchValues = [IndividualBirds(i).AllSyllableTemplateMatchValues; cell2mat(IndividualBirds(i).SortedBirdParameters(j).TemplateMatchValues)'];
        IndividualBirds(i).AllSyllableFeatValues = [IndividualBirds(i).AllSyllableFeatValues; IndividualBirds(i).SortedBirdParameters(j).SAPFeatsMatrix];
        IndividualBirds(i).AllSyllableLogAmplitudeKao = [IndividualBirds(i).AllSyllableLogAmplitudeKao; IndividualBirds(i).SortedBirdParameters(j).AmplitudeKao(:)];
        IndividualBirds(i).AllSyllableLogAmplitudeRMS = [IndividualBirds(i).AllSyllableLogAmplitudeRMS; IndividualBirds(i).SortedBirdParameters(j).AmplitudeRMS(:)];
        IndividualBirds(i).AllSyllableFFAutocorrKao = [IndividualBirds(i).AllSyllableFFAutocorrKao; IndividualBirds(i).SortedBirdParameters(j).FFAutocorrKao(:)];
        IndividualBirds(i).AllConditionIndices = [IndividualBirds(i).AllConditionIndices; ones(size(IndividualBirds(i).SortedBirdParameters(j).SyllableData,1), 1)*IndividualBirds(i).SortedBirdParameters(j).ConditionIndices];
        IndividualBirds(i).AllMicrophoneIndices = [IndividualBirds(i).AllMicrophoneIndices; ones(size(IndividualBirds(i).SortedBirdParameters(j).SyllableData,1), 1)*IndividualBirds(i).SortedBirdParameters(j).MicrophoneIndices];
        IndividualBirds(i).AllRecordingDayIndices = [IndividualBirds(i).AllRecordingDayIndices; ones(size(IndividualBirds(i).SortedBirdParameters(j).SyllableData,1), 1)*IndividualBirds(i).SortedBirdParameters(j).RecordingDayIndex];
        IndividualBirds(i).AllSyllableTimeCategories = [IndividualBirds(i).AllSyllableTimeCategories; IndividualBirds(i).SortedBirdParameters(j).SyllableTimeCategories];
        IndividualBirds(i).AllSyllableListBoutNum = [IndividualBirds(i).AllSyllableListBoutNum; IndividualBirds(i).SortedBirdParameters(j).SyllableListBoutNum(:)+TotalBoutNum];
        TotalBoutNum = TotalBoutNum + size(IndividualBirds(i).SortedBirdParameters(j).Bouts,1);
    end
end

% ============= Bout statistics ===========================================
% In this section, I will calculate bout statistics for each bout
% Things that will be calculated
% 1. Bout length
% 2. # of motifs
% 3. Whether the bout has only INs (1) or INs + calls (0) at the beginning
% 4. First motif duration
% 5. All motif durations
% 6. Mean motif duration
% 7. # of INs
% 8. Total # of syllables (INs + calls)
% 9. Average FF within a bout
% 10. FF of first occurence of the syllable
% 11. All FFs within a bout
% 12. Time in the bout when the first motif syllable appears
% 13. CV of FF within a bout
% 14. Average FM within a bout
% 15. Average Mean Freq within a bout
% All of these will be calculated only for bouts that have motif syllables
% in them and have 2s silence before and after the bout

disp('Now calculating bout statistics ...');

% First to list out all the things that we are calculating
for i = 1:length(IndividualBirds),
    IndividualBirds(i).BoutStatisticsColumnNames{1} = 'BoutIndex'; % index of the bout
    IndividualBirds(i).BoutStatisticsColumnNames{2} = 'BoutLength'; % bout length
    IndividualBirds(i).BoutStatisticsColumnNames{3} = 'FirstMotifSyllTime'; % Time of occurence of first motif syllable within the bout
    IndividualBirds(i).BoutStatisticsColumnNames{4} = 'AllINs'; % Num INs before onset of first motif - all INs irrespective of interval between them and irrespective of whether they're the last set before motif onset
    IndividualBirds(i).BoutStatisticsColumnNames{5} = 'NumINs'; % Last set of consecutive INs with intervals < 500ms 
    IndividualBirds(i).BoutStatisticsColumnNames{6} = 'AllSylls'; % Total number of syllables before first motif onset irrespective of interval between them
    IndividualBirds(i).BoutStatisticsColumnNames{7} = 'NumSylls'; % Last set of syllables before motif onset with intervals < 500ms
    IndividualBirds(i).BoutStatisticsColumnNames{8} = 'FirstMotifSyllTime_GapsLessThan500ms'; % Considers only the last set of consecutive syllables with < 500ms gaps and calculates first motif syllable onset time from the onset of the first of these consecutive syllables
    IndividualBirds(i).BoutStatisticsColumnNames{9} = 'OnlyINBoutOrNot'; % whether a bout has only INs (1) before motif onset or has INs and calls/other syllables (0)
    IndividualBirds(i).BoutStatisticsColumnNames{10} = 'NumMotifs'; % number of motifs in the bout
    IndividualBirds(i).BoutStatisticsColumnNames{11} = 'FirstMotifDur'; % first motif duration calculated as time between onset of first syllable and offset of last syllable
    IndividualBirds(i).BoutStatisticsColumnNames{12} = 'AverageMotifDur'; % average motif duration calculated as above
    IndividualBirds(i).BoutStatisticsColumnNames{13} = 'FirstMotifDurOnset'; % first motif duration calculated as time between onest of first syllable and onset of last syllable
    IndividualBirds(i).BoutStatisticsColumnNames{14} = 'AverageMotifDurOnset'; % average motif duration calculated as time between onset of first syllable and onset of last syllable
    IndividualBirds(i).BoutStatisticsColumnNames{15} = 'AverageSyllMeanFreq'; % average of mean frequency for all sylls - calculated as average of average of mean freq for a given syll type - for all syll types
    IndividualBirds(i).BoutStatisticsColumnNames{16} = 'AverageSyllMeanFM'; % average of mean FM for all sylls - calculated as average of average of mean FM for a given syll type - for all syll types
    IndividualBirds(i).BoutStatisticsColumnNames{17} = 'AverageSyllMeanLogAmplitude'; % average of mean log amplitude for all sylls - calculated as average of average of mean log amplitude for a given syll type - for all syll types
    IndividualBirds(i).BoutStatisticsColumnNames{18} = 'AverageSyllFF'; % average of mean FF for all sylls - calculated as average of average of mean FF for a given syll type - for all syll types
    IndividualBirds(i).BoutStatisticsColumnNames{19} = 'CVSyllFF'; % average of CV of FF for all sylls - calculated as average of CV of FF for a given syll type - for all syll types
    IndividualBirds(i).BoutStatisticsColumnNames{20} = 'Microphone'; % Microphone with which this was recorded
    IndividualBirds(i).BoutStatisticsColumnNames{21} = 'RecordingDay'; % Recording Day index for the corresponding bout
    IndividualBirds(i).BoutStatisticsColumnNames{22} = 'DateLabel'; % Date label for the bout
    IndividualBirds(i).BoutStatisticsColumnNames{23} = 'Condition'; % Condition in which this was recorded
end

for i = 1:length(IndividualBirds),
    disp(['Bird #', num2str(i), ' ...']);
    ValidSongBouts = find((IndividualBirds(i).Bouts(:,7) == 1) & (IndividualBirds(i).Bouts(:,8) > 0) & (IndividualBirds(i).Bouts(:,9) > 1));
    Index = 1;
    INLabels = IndividualBirds(i).SortedBirdParameters(1).INLabels;
    MotifLabels = IndividualBirds(i).SortedBirdParameters(1).MotifLabels;
    CommonMotifs = IndividualBirds(i).SortedBirdParameters(1).CommonMotifs;
    IndividualBirds(i).UniqueSyllLabels = unique(char(IndividualBirds(i).AllSyllableData(:,1)));
    
    % Get all the syllables for which we have measured FF
    BirdParameters_FieldNames = fieldnames(IndividualBirds(i).SortedBirdParameters(1));
    FF_Fields = find(cellfun(@length, strfind(BirdParameters_FieldNames, 'FFSyllable')));
    
    FF_Index = 1;
    for j = 1:length(FF_Fields),
        FFParameters_String = eval(['IndividualBirds(', num2str(i), ').SortedBirdParameters(', num2str(1), ').', BirdParameters_FieldNames{FF_Fields(j)}]);
        ColonIndex = find(FFParameters_String == ':');
        if (~isempty(ColonIndex))
            IndividualBirds(i).FF_SyllLabels(FF_Index) = FFParameters_String(ColonIndex - 1);
            FF_Index = FF_Index + 1;
        end
    end
    
    for j = ValidSongBouts(:)',
        BoutLabels = char(IndividualBirds(i).AllSyllableData(IndividualBirds(i).Bouts(j,1):IndividualBirds(i).Bouts(j,2),1));
        BoutOnsets = IndividualBirds(i).AllSyllableData(IndividualBirds(i).Bouts(j,1):IndividualBirds(i).Bouts(j,2),4);
        BoutOffsets = IndividualBirds(i).AllSyllableData(IndividualBirds(i).Bouts(j,1):IndividualBirds(i).Bouts(j,2),5);
        BoutSAPFeats = IndividualBirds(i).AllSyllableFeatValues(IndividualBirds(i).Bouts(j,1):IndividualBirds(i).Bouts(j,2),:);
        BoutCondition = IndividualBirds(i).AllConditionIndices(IndividualBirds(i).Bouts(j,1):IndividualBirds(i).Bouts(j,2),:);
        BoutMicrophone = IndividualBirds(i).AllMicrophoneIndices(IndividualBirds(i).Bouts(j,1):IndividualBirds(i).Bouts(j,2),:);
        BoutRecordingDay = IndividualBirds(i).AllRecordingDayIndices(IndividualBirds(i).Bouts(j,1):IndividualBirds(i).Bouts(j,2),:);
        TempBoutDateLabel = cell2mat((IndividualBirds(i).RecordingDays(BoutRecordingDay))');
        BoutDateLabel =[];
        for h = 1:size(TempBoutDateLabel,1),
            ConcatBoutDateLabel = strcat(TempBoutDateLabel(h,1:end));
            BoutDateLabel = [BoutDateLabel; str2double(ConcatBoutDateLabel)];
            clear ConcatBoutDateLabel;
        end
        TempBoutDateLabel =[];
        % Now, first we need to check if this bout has only INs at the
        % beginning. If it does not, then skip this bout.
        
        % Now, first check if there is atleast one full motif in the bout,
        % if not, then skip the bout
        if (isempty(strfind(BoutLabels(:)', CommonMotifs{1})))
            disp(['Skipped bout without 1 motif: ', IndividualBirds(i).SongFileNames{IndividualBirds(i).Bouts(j,3)}]);
            continue;
        end
        
        % Now, also skip any bout that starts with a motif from the first
        % syllable as this is likely to be a case of previous syllables not
        % being labelled.
        TempMotifs = strfind(BoutLabels(:)', CommonMotifs{1});
        if (TempMotifs(1) == 1)
            disp(['Skipped bout starts with the motif: ', IndividualBirds(i).SongFileNames{IndividualBirds(i).Bouts(j,3)}]);
            continue;
        end
        
        % Add another two variable to connect whether these bouts are
        % directed or not and whether female responded or not - this will
        % come the video scoring data if present
        if (~isempty(VideoScoringData))
            [Dir, SongFileName, Ext] = fileparts(IndividualBirds(i).SongFileNames{IndividualBirds(i).Bouts(j,3)});
            SongFileName = [SongFileName, Ext];
            MatchFlag = 0;
            for SongFileIndex = 1:length(VideoScoringData),
                if (strcmp(VideoScoringData(SongFileIndex).FileName, SongFileName))
                    if (abs((VideoScoringData(SongFileIndex).BoutOnset * 1000) - BoutOnsets(1)) <= 15)
                        % Basically above two conditions checked that the
                        % filename is the same and the absolute difference
                        % in bout onset time for the bout and the data from
                        % the video scoring file do not differ by more than
                        % 15 ms, then it is the correct bout
                        MatchFlag = 1;
                        IndividualBirds(i).BoutCategorisation{Index} = VideoScoringData(SongFileIndex).DirecteDUNdirected;
                        IndividualBirds(i).BoutFemaleResponse{Index} = VideoScoringData(SongFileIndex).NRorR;
                        break;
                    end
                end
            end
            if (MatchFlag == 0)
                IndividualBirds(i).BoutCategorisation{Index} = 'NA';
                IndividualBirds(i).BoutFemaleResponse{Index} = 'NA';
            end
        else
            IndividualBirds(i).BoutCategorisation{Index} = 'NA';
            IndividualBirds(i).BoutFemaleResponse{Index} = 'NA';
        end
        
        % Now I also have to add the time of the bout relative to 0 which
        % is the time at which female was introduced
        % For this, it will be (BoutOnsetTime + FileTime -
        % Femaleintroductiontime)
        IndividualBirds(i).BoutStatisticsTimeRelativeToStart(Index) = IndividualBirds(i).BoutTimesRelativeToStart(j);
        
        ColumnIndex = find(strcmp('BoutIndex', IndividualBirds(i).BoutStatisticsColumnNames));
        IndividualBirds(i).BoutStatistics(Index, ColumnIndex) = j; 
        
        ColumnIndex = find(strcmp('BoutLength', IndividualBirds(i).BoutStatisticsColumnNames));
        IndividualBirds(i).BoutStatistics(Index, ColumnIndex) = BoutOffsets(end) - BoutOnsets(1); % bout length
 
        % This is checking for first motif syllable in the bout
%         for k = 1:length(BoutLabels),
%             if (~isempty(find(MotifLabels == BoutLabels(k))))
%                 FirstMotifSyllIndex = k;
%                 break;
%             end
%         end

        % This is checking for the first motif in the bout based on common
        % motifs
        TempMotifs = strfind(BoutLabels(:)', CommonMotifs{1});
        FirstMotifSyllIndex = TempMotifs(1);
        
        ColumnIndex = find(strcmp('FirstMotifSyllTime', IndividualBirds(i).BoutStatisticsColumnNames));
        if (k > length(BoutLabels))
            IndividualBirds(i).BoutStatistics(Index,ColumnIndex) = NaN; 
        else
            IndividualBirds(i).BoutStatistics(Index,ColumnIndex) = BoutOnsets(FirstMotifSyllIndex) - BoutOnsets(1); 
        end
        
        % Now to check # of introductory notes
        % Need to get only the last set of consecutive INs that don't have
        % more than 500ms between them as the number of INs - this is how
        % other papers count INs
        % Will also keep track of total # of INs
        
        NumINs = 0;
        INs = [];
        for k = 1:(FirstMotifSyllIndex-1),
            if (~isempty(find(INLabels == BoutLabels(k))))
                INs(end+1) = k;
                NumINs = NumINs + 1;
            end
        end
        
        ColumnIndex = find(strcmp('AllINs', IndividualBirds(i).BoutStatisticsColumnNames));
        IndividualBirds(i).BoutStatistics(Index,ColumnIndex) = NumINs; 
        
        
        % First to ensure that we only take the last set of consecutive INs
        Gaps = diff(INs);
        LongGaps = find(Gaps > 1);
        if (~isempty(LongGaps))
            INs = INs(LongGaps(end)+1:end);
        end
        
        % Now to check that the gap between these INs is not more than
        % 500ms
        Gaps = BoutOnsets(INs(2:end)) - BoutOffsets(INs(1:end-1));
        LongGaps = find(Gaps > 500);
        if (~isempty(LongGaps))
            INs = INs(LongGaps(end)+1:end);
        end
        NumINs = length(INs);
        
        ColumnIndex = find(strcmp('NumINs', IndividualBirds(i).BoutStatisticsColumnNames));
        IndividualBirds(i).BoutStatistics(Index,ColumnIndex) = NumINs; 
        
        ColumnIndex = find(strcmp('AllSylls', IndividualBirds(i).BoutStatisticsColumnNames));
        IndividualBirds(i).BoutStatistics(Index,ColumnIndex) = FirstMotifSyllIndex - 1; 
        
        % Now, for number of syllables, we're going to consider only the
        % last set of syllables that have < 500ms between them.
        
        % Now to check that the gap between these syllables is not more
        % than 500ms
        
        Gaps = BoutOnsets(2:(FirstMotifSyllIndex - 1)) - BoutOffsets(1:(FirstMotifSyllIndex - 2));
        LongGaps = find(Gaps > 500);
        if (~isempty(LongGaps))
            Sylls = (LongGaps(end) + 1):1:FirstMotifSyllIndex-1;
        else
            Sylls = 1:1:FirstMotifSyllIndex - 1;
        end
        NumSylls = length(Sylls);
        
        ColumnIndex = find(strcmp('NumSylls', IndividualBirds(i).BoutStatisticsColumnNames));
        IndividualBirds(i).BoutStatistics(Index,ColumnIndex) = NumSylls; 
        
        % Now to consider the first motif syll time based on the last set
        % of consecutive syllables that have gaps < 500ms. Will consider
        % the motif syllable onset time from the start of this sequence of
        % consecutive syllables
        
        ColumnIndex = find(strcmp('FirstMotifSyllTime_GapsLessThan500ms', IndividualBirds(i).BoutStatisticsColumnNames));
        IndividualBirds(i).BoutStatistics(Index,ColumnIndex) = BoutOnsets(FirstMotifSyllIndex) - BoutOnsets(Sylls(1));
        
        AllSyllIndex = find(strcmp('AllSylls', IndividualBirds(i).BoutStatisticsColumnNames));
        AllINIndex = find(strcmp('AllINs', IndividualBirds(i).BoutStatisticsColumnNames));
        
        ColumnIndex = find(strcmp('OnlyINBoutOrNot', IndividualBirds(i).BoutStatisticsColumnNames));
        if (IndividualBirds(i).BoutStatistics(Index,AllINIndex) == IndividualBirds(i).BoutStatistics(Index,AllSyllIndex))
            IndividualBirds(i).BoutStatistics(Index, ColumnIndex) = 1; % means the bout has only INs, since # of INs and # of syllables is equal
        else
            IndividualBirds(i).BoutStatistics(Index, ColumnIndex) = 0; % means the bout has both INs and calls (or other syllables) since # of INs not equal to # of syllables
        end
        
        
        % Assuming each bird has only one common motif
        ColumnIndex = find(strcmp('NumMotifs', IndividualBirds(i).BoutStatisticsColumnNames));
        IndividualBirds(i).BoutStatistics(Index,ColumnIndex) = length(strfind(BoutLabels(:)', CommonMotifs{1})); % # of motifs in the mout
                
        BoutMotifs = strfind(BoutLabels(:)', CommonMotifs{1});
        IndividualBirds(i).BoutMotifDurs{Index} = BoutOffsets(BoutMotifs - 1 + length(CommonMotifs{1})) - BoutOnsets(BoutMotifs); % measured as difference between offset of last syllable and onset of first syllable
        IndividualBirds(i).BoutMotifDursOnsets{Index} = BoutOnsets(BoutMotifs - 1 + length(CommonMotifs{1})) - BoutOnsets(BoutMotifs); % measured as difference between onset of last syllable and onset of first syllable
        
        if (isempty(BoutMotifs))
            ColumnIndex = find(strcmp('FirstMotifDur', IndividualBirds(i).BoutStatisticsColumnNames));
            IndividualBirds(i).BoutStatistics(Index,ColumnIndex) = NaN; % first motif duration
            
            ColumnIndex = find(strcmp('AverageMotifDur', IndividualBirds(i).BoutStatisticsColumnNames));
            IndividualBirds(i).BoutStatistics(Index,ColumnIndex) = NaN; % average motif duration
            
            ColumnIndex = find(strcmp('FirstMotifDurOnset', IndividualBirds(i).BoutStatisticsColumnNames));
            IndividualBirds(i).BoutStatistics(Index,ColumnIndex) = NaN; % first motif duration onsets
            
            ColumnIndex = find(strcmp('AverageMotifDurOnset', IndividualBirds(i).BoutStatisticsColumnNames));
            IndividualBirds(i).BoutStatistics(Index,ColumnIndex) = NaN; % average motif duration onsets
        else
            ColumnIndex = find(strcmp('FirstMotifDur', IndividualBirds(i).BoutStatisticsColumnNames));
            IndividualBirds(i).BoutStatistics(Index,ColumnIndex) = IndividualBirds(i).BoutMotifDurs{Index}(1); % first motif duration
            
            ColumnIndex = find(strcmp('AverageMotifDur', IndividualBirds(i).BoutStatisticsColumnNames));
            IndividualBirds(i).BoutStatistics(Index,ColumnIndex) = mean(IndividualBirds(i).BoutMotifDurs{Index}); % average motif duration
            
            ColumnIndex = find(strcmp('FirstMotifDurOnset', IndividualBirds(i).BoutStatisticsColumnNames));
            IndividualBirds(i).BoutStatistics(Index,ColumnIndex) = IndividualBirds(i).BoutMotifDursOnsets{Index}(1); % first motif duration with onsets
            
            ColumnIndex = find(strcmp('AverageMotifDurOnset', IndividualBirds(i).BoutStatisticsColumnNames));
            IndividualBirds(i).BoutStatistics(Index,ColumnIndex) = mean(IndividualBirds(i).BoutMotifDursOnsets{Index}); % average motif duration
        end
        
        ColumnIndex = find(strcmp('Condition', IndividualBirds(i).BoutStatisticsColumnNames));
        IndividualBirds(i).BoutStatistics(Index,ColumnIndex) = mean(BoutCondition); 
        
        ColumnIndex = find(strcmp('Microphone', IndividualBirds(i).BoutStatisticsColumnNames));
        IndividualBirds(i).BoutStatistics(Index,ColumnIndex) = mean(BoutMicrophone);
        
        ColumnIndex = find(strcmp('RecordingDay', IndividualBirds(i).BoutStatisticsColumnNames));
        IndividualBirds(i).BoutStatistics(Index,ColumnIndex) = mean(BoutRecordingDay);
        
        ColumnIndex = find(strcmp('DateLabel', IndividualBirds(i).BoutStatisticsColumnNames));
        IndividualBirds(i).BoutStatistics(Index,ColumnIndex) = mean(BoutDateLabel);
        clear BoutDateLabel;
        
        % Now to put all the FF and amplitudes of individual syllables
        % together
        % First FF - should do this for all the syllables that have FF boundaries set in the original csv file
        for k = 1:length(IndividualBirds(i).FF_SyllLabels),
            BoutSylls = find(BoutLabels(:) == IndividualBirds(i).FF_SyllLabels(k));
            if (~isempty(BoutSylls))
                ColumnIndex = find(strcmp(IndividualBirds(1).SortedBirdParameters(1).SAPFeat_FieldNames, 'MeanFundamentalFrequency'));
                IndividualBirds(i).FFSyll(k).MeanFF{Index} = BoutSAPFeats(BoutSylls, ColumnIndex)'; 
                
                ColumnIndex = find(strcmp(IndividualBirds(1).SortedBirdParameters(1).SAPFeat_FieldNames, 'MedianFundamentalFrequency'));
                IndividualBirds(i).FFSyll(k).MedianFF{Index} = BoutSAPFeats(BoutSylls, ColumnIndex)'; 
            else
                IndividualBirds(i).FFSyll(k).MeanFF{Index} = []; 
                
                IndividualBirds(i).FFSyll(k).MedianFF{Index} = []; 
            end
        end
        
        % Now for amplitude - should do this for all unique syllables
        for k = 1:length(IndividualBirds(i).UniqueSyllLabels),
            BoutSylls = find(BoutLabels(:) == IndividualBirds(i).UniqueSyllLabels(k));
            if (~isempty(BoutSylls))
                ColumnIndex = find(strcmp(IndividualBirds(1).SortedBirdParameters(1).SAPFeat_FieldNames, 'LogAmplitude'));
                IndividualBirds(i).AllSylls(k).LogAmplitude{Index} = BoutSAPFeats(BoutSylls, ColumnIndex)'; 
                
                ColumnIndex = find(strcmp(IndividualBirds(1).SortedBirdParameters(1).SAPFeat_FieldNames, 'FrequencyModulation'));
                IndividualBirds(i).AllSylls(k).FrequencyModulation{Index} = BoutSAPFeats(BoutSylls, ColumnIndex)'; 
                
                ColumnIndex = find(strcmp(IndividualBirds(1).SortedBirdParameters(1).SAPFeat_FieldNames, 'MeanFrequency'));
                IndividualBirds(i).AllSylls(k).MeanFrequency{Index} = BoutSAPFeats(BoutSylls, ColumnIndex)'; 
                IndividualBirds(i).AllSylls(k).SyllMicrophone{Index} = BoutMicrophone(BoutSylls);
            else
                IndividualBirds(i).AllSylls(k).LogAmplitude{Index} = []; 
                IndividualBirds(i).AllSylls(k).MeanFrequency{Index} = []; 
                IndividualBirds(i).AllSylls(k).FrequencyModulation{Index} = []; 
                IndividualBirds(i).AllSylls(k).SyllMicrophone{Index} = [];
            end
        end
        
        
        Index = Index + 1;
    end
    % Now I need to remove outliers for log amplitude, motif duration and
    % FF
    % After removing outliers for motif duration, I have to redo the motif
    % duration values
    % Outliers are removed based on the condition that these are > 3*IQR
    % above the 75th percentile of all values across conditions or < 3*IQR
    % below the 25th percentile of all values across conditions
    % First for FF
    disp('Removing outliers');
    
%     for j = 1:length(IndividualBirds(i).FFSyll),
%         % First for mean FF
%         SyllFFValues = cell2mat(IndividualBirds(i).FFSyll(j).MeanFF);
%         UpperOutlierThreshold = prctile(SyllFFValues, 75) + (3 * iqr(SyllFFValues));
%         LowerOutlierThreshold = prctile(SyllFFValues, 25) - (3 * iqr(SyllFFValues));
%         OutlierNum = 0;
%         for k = 1:length(IndividualBirds(i).FFSyll(j).MeanFF),
%             Indices = find((IndividualBirds(i).FFSyll(j).MeanFF{k} > UpperOutlierThreshold) | (IndividualBirds(i).FFSyll(j).MeanFF{k} < LowerOutlierThreshold));
%             if (~isempty(Indices))
%                 OutlierNum = OutlierNum + length(Indices);
%                 IndividualBirds(i).FFSyll(j).MeanFF{k}(Indices) = [];
%             end
%         end
%         
%         disp([BirdNames{i}, ': Removed ', num2str(OutlierNum), ' mean FF outliers out of a total of ', num2str(length(SyllFFValues)), ' syllable ', IndividualBirds(i).FF_SyllLabels(j), ' (', num2str(OutlierNum * 100 / length(SyllFFValues)), '%)']);
%         
%         % Next for median FF
%         SyllFFValues = cell2mat(IndividualBirds(i).FFSyll(j).MedianFF);
%         UpperOutlierThreshold = prctile(SyllFFValues, 75) + (3 * iqr(SyllFFValues));
%         LowerOutlierThreshold = prctile(SyllFFValues, 25) - (3 * iqr(SyllFFValues));
%         OutlierNum = 0;
%         for k = 1:length(IndividualBirds(i).FFSyll(j).MedianFF),
%             Indices = find((IndividualBirds(i).FFSyll(j).MedianFF{k} > UpperOutlierThreshold) | (IndividualBirds(i).FFSyll(j).MedianFF{k} < LowerOutlierThreshold));
%             if (~isempty(Indices))
%                 OutlierNum = OutlierNum + length(Indices);
%                 IndividualBirds(i).FFSyll(j).MedianFF{k}(Indices) = [];
%             end
%         end
%         disp([BirdNames{i}, ': Removed ', num2str(OutlierNum), ' median FF outliers out of a total of ', num2str(length(SyllFFValues)), ' syllable ', IndividualBirds(i).FF_SyllLabels(j), ' (', num2str(OutlierNum * 100 / length(SyllFFValues)), '%)']);
%     end
%     
%     % Next do the same for log amplitude values
%     for j = 1:length(IndividualBirds(i).AllSylls),
%         % First for mean FF
%         SyllFFValues = cell2mat(IndividualBirds(i).AllSylls(j).LogAmplitude);
%         UpperOutlierThreshold = prctile(SyllFFValues, 75) + (3 * iqr(SyllFFValues));
%         LowerOutlierThreshold = prctile(SyllFFValues, 25) - (3 * iqr(SyllFFValues));
%         OutlierNum = 0;
%         for k = 1:length(IndividualBirds(i).AllSylls(j).LogAmplitude),
%             Indices = find((IndividualBirds(i).AllSylls(j).LogAmplitude{k} > UpperOutlierThreshold) | (IndividualBirds(i).AllSylls(j).LogAmplitude{k} < LowerOutlierThreshold));
%             if (~isempty(Indices))
%                 OutlierNum = OutlierNum + length(Indices);
%                 IndividualBirds(i).AllSylls(j).LogAmplitude{k}(Indices) = [];
%             end
%         end
%         disp([BirdNames{i}, ': Removed ', num2str(OutlierNum), ' mean log ampitude outliers out of a total of ', num2str(length(SyllFFValues)), ' syllable ', IndividualBirds(i).UniqueSyllLabels(j), ' (', num2str(OutlierNum * 100 / length(SyllFFValues)), '%)']);
%     end
%     
%     % Now do the same for motif duration calculated in two ways
%     % First for motif durs calculated as offset of last syllable - onset of
%     % first syllable
%     MotifDurValues = cell2mat(IndividualBirds(i).BoutMotifDurs');
%     UpperOutlierThreshold = prctile(SyllFFValues, 75) + (3 * iqr(MotifDurValues));
%     LowerOutlierThreshold = prctile(SyllFFValues, 25) - (3 * iqr(MotifDurValues));
%     OutlierNum = 0;
%     for k = 1:length(IndividualBirds(i).BoutMotifDurs{k}),
%         Indices = find((IndividualBirds(i).BoutMotifDurs{k} > UpperOutlierThreshold) | (IndividualBirds(i).BoutMotifDurs{k} < LowerOutlierThreshold));
%         if (~isempty(Indices))
%             OutlierNum = OutlierNum + length(Indices);
%             IndividualBirds(i).BoutMotifDurs{k}(Indices) = [];
%         end
%     end
%     disp([BirdNames{i}, ': Removed ', num2str(OutlierNum), ' motif duration outliers out of a total of ', num2str(length(MotifDurValues)), ' (', num2str(OutlierNum * 100 / length(MotifDurValues)), '%)']);
% 
%     % Next for motif durations calculated as onset of last syllable - onset
%     % of first syllable
%     MotifDurValues = cell2mat(IndividualBirds(i).BoutMotifDursOnsets');
%     UpperOutlierThreshold = prctile(SyllFFValues, 75) + (3 * iqr(MotifDurValues));
%     LowerOutlierThreshold = prctile(SyllFFValues, 25) - (3 * iqr(MotifDurValues));
%     OutlierNum = 0;
%     for k = 1:length(IndividualBirds(i).BoutMotifDursOnsets{k}),
%         Indices = find((IndividualBirds(i).BoutMotifDursOnsets{k} > UpperOutlierThreshold) | (IndividualBirds(i).BoutMotifDursOnsets{k} < LowerOutlierThreshold));
%         if (~isempty(Indices))
%             OutlierNum = OutlierNum + length(Indices);
%             IndividualBirds(i).BoutMotifDursOnsets{k}(Indices) = [];
%         end
%     end
%     disp([BirdNames{i}, ': Removed ', num2str(OutlierNum), ' motif duration outliers based on onsets out of a total of ', num2str(length(MotifDurValues)), ' (', num2str(OutlierNum * 100 / length(MotifDurValues)), '%)']);
    
    % Now redo all the first motif duration values
    for Index = 1:length(IndividualBirds(i).BoutMotifDurs),
        if (isempty(IndividualBirds(i).BoutMotifDurs{Index}))
            ColumnIndex = find(strcmp('FirstMotifDur', IndividualBirds(i).BoutStatisticsColumnNames));
            IndividualBirds(i).BoutStatistics(Index,ColumnIndex) = NaN; % first motif duration
            
            ColumnIndex = find(strcmp('AverageMotifDur', IndividualBirds(i).BoutStatisticsColumnNames));
            IndividualBirds(i).BoutStatistics(Index,ColumnIndex) = NaN; % average motif duration
        else
            ColumnIndex = find(strcmp('FirstMotifDur', IndividualBirds(i).BoutStatisticsColumnNames));
            IndividualBirds(i).BoutStatistics(Index,ColumnIndex) = IndividualBirds(i).BoutMotifDurs{Index}(1); % first motif duration
            
            ColumnIndex = find(strcmp('AverageMotifDur', IndividualBirds(i).BoutStatisticsColumnNames));
            IndividualBirds(i).BoutStatistics(Index,ColumnIndex) = mean(IndividualBirds(i).BoutMotifDurs{Index}); % average motif duration
        end
        if (isempty(IndividualBirds(i).BoutMotifDursOnsets{Index}))
            ColumnIndex = find(strcmp('FirstMotifDurOnset', IndividualBirds(i).BoutStatisticsColumnNames));
            IndividualBirds(i).BoutStatistics(Index,ColumnIndex) = NaN; % first motif duration onsets
            
            ColumnIndex = find(strcmp('AverageMotifDurOnset', IndividualBirds(i).BoutStatisticsColumnNames));
            IndividualBirds(i).BoutStatistics(Index,ColumnIndex) = NaN; % average motif duration onsets
        else
            ColumnIndex = find(strcmp('FirstMotifDurOnset', IndividualBirds(i).BoutStatisticsColumnNames));
            IndividualBirds(i).BoutStatistics(Index,ColumnIndex) = IndividualBirds(i).BoutMotifDursOnsets{Index}(1); % first motif duration with onsets
            
            ColumnIndex = find(strcmp('AverageMotifDurOnset', IndividualBirds(i).BoutStatisticsColumnNames));
            IndividualBirds(i).BoutStatistics(Index,ColumnIndex) = mean(IndividualBirds(i).BoutMotifDursOnsets{Index}); % average motif duration
        end
    end
    % Now add the mean FF and CV FF for each bout only for syllables that
    % are both motif syllables and we have measured FF
    SyllsToUse = intersect(IndividualBirds(i).SortedBirdParameters(1).MotifLabels, IndividualBirds(i).FF_SyllLabels);
    for Index = 1:length(IndividualBirds(i).BoutMotifDurs),
        MeanFFVal = 0;
        CVFFVal = 0;
        MeanFreqVal = 0;
        MeanLogAmplitudeVal = 0;
        MeanFMVal = 0;
        
        BoutFlag = 1; % decide whether to use the bout or not for Mean FF and CV FF
        % decision is based on whether all of the above syllables are
        % present or not
        for k = 1:length(SyllsToUse),
            FFSyllIndex = find(IndividualBirds(i).FF_SyllLabels == SyllsToUse(k));
            if (isempty(IndividualBirds(i).FFSyll(FFSyllIndex).MeanFF{Index}))
                BoutFlag = 0;
                break;
            end
        end
        if (BoutFlag == 1)
            for k = 1:length(SyllsToUse),
                FFSyllIndex = find(IndividualBirds(i).FF_SyllLabels == SyllsToUse(k));
                MeanFFVal = MeanFFVal + mean(IndividualBirds(i).FFSyll(FFSyllIndex).MeanFF{Index});
                CVFFVal = CVFFVal + std(IndividualBirds(i).FFSyll(FFSyllIndex).MeanFF{Index})/mean(IndividualBirds(i).FFSyll(FFSyllIndex).MeanFF{Index});
                
                AllSyllIndex = find(IndividualBirds(i).UniqueSyllLabels == SyllsToUse(k));
                MeanFreqVal = MeanFreqVal + mean(IndividualBirds(i).AllSylls(AllSyllIndex).MeanFrequency{Index});
                MeanFMVal = MeanFMVal + mean(IndividualBirds(i).AllSylls(AllSyllIndex).FrequencyModulation{Index});
                MeanLogAmplitudeVal = MeanLogAmplitudeVal + mean(IndividualBirds(i).AllSylls(AllSyllIndex).LogAmplitude{Index});
            end
            ColumnIndex = find(strcmp('AverageSyllFF', IndividualBirds(i).BoutStatisticsColumnNames));
            IndividualBirds(i).BoutStatistics(Index,ColumnIndex) = MeanFFVal/length(SyllsToUse);
            ColumnIndex = find(strcmp('CVSyllFF', IndividualBirds(i).BoutStatisticsColumnNames));
            IndividualBirds(i).BoutStatistics(Index,ColumnIndex) = CVFFVal/length(SyllsToUse);
            
            ColumnIndex = find(strcmp('AverageSyllMeanFreq', IndividualBirds(i).BoutStatisticsColumnNames));
            IndividualBirds(i).BoutStatistics(Index,ColumnIndex) = MeanFreqVal/length(SyllsToUse);
            ColumnIndex = find(strcmp('AverageSyllMeanFM', IndividualBirds(i).BoutStatisticsColumnNames));
            IndividualBirds(i).BoutStatistics(Index,ColumnIndex) = MeanFMVal/length(SyllsToUse);
            ColumnIndex = find(strcmp('AverageSyllMeanLogAmplitude', IndividualBirds(i).BoutStatisticsColumnNames));
            IndividualBirds(i).BoutStatistics(Index,ColumnIndex) = MeanLogAmplitudeVal/length(SyllsToUse);
        else
            ColumnIndex = find(strcmp('AverageSyllFF', IndividualBirds(i).BoutStatisticsColumnNames));
            IndividualBirds(i).BoutStatistics(Index,ColumnIndex) = NaN;
            ColumnIndex = find(strcmp('CVSyllFF', IndividualBirds(i).BoutStatisticsColumnNames));
            IndividualBirds(i).BoutStatistics(Index,ColumnIndex) = NaN;
            
            ColumnIndex = find(strcmp('AverageSyllMeanFreq', IndividualBirds(i).BoutStatisticsColumnNames));
            IndividualBirds(i).BoutStatistics(Index,ColumnIndex) = NaN;
            ColumnIndex = find(strcmp('AverageSyllMeanFM', IndividualBirds(i).BoutStatisticsColumnNames));
            IndividualBirds(i).BoutStatistics(Index,ColumnIndex) = NaN;
            ColumnIndex = find(strcmp('AverageSyllMeanLogAmplitude', IndividualBirds(i).BoutStatisticsColumnNames));
            IndividualBirds(i).BoutStatistics(Index,ColumnIndex) = NaN;
        end 
    end  
            
end

% Now to add the actual distances for each of the conditions

for i = 1:length(IndividualBirds),
    for j = 1:length(IndividualBirds(i).Conditions),
        switch IndividualBirds(i).Conditions{j}
            case 'L0'
                IndividualBirds(i).ConditionDistances(j) = 0;
            
            case 'L1'
                IndividualBirds(i).ConditionDistances(j) = 20;
                
            case 'L2'
                IndividualBirds(i).ConditionDistances(j) = 60;
                
            case 'L3'
                IndividualBirds(i).ConditionDistances(j) = 110;
                
            case 'L4'
                IndividualBirds(i).ConditionDistances(j) = 165;
                
            case 'UN'
                IndividualBirds(i).ConditionDistances(j) = 200;
        end
    end
end
                
% ============ Now to make plots ==========================================
FeatureToPlot_ColNo = 1;
PlotType = 'mean';
DifferentMicrophones = [{'BP'} {'MM'} {'HF'}];

% Made one common way of plotting the results of any analysis. Basically,
% there is a script file called Harini_<<PlotOption>>.m and that is called
% to do the appropriate analysis. The first few lines are to check if the
% file exists - if not, then program returns with a message that such a
% program does not exist. 
% All scripts also take the same inputs to make things uniform

if (~isempty(which(['Harini_', PlotOption, '.m'])))
    eval(['Harini_', PlotOption, '(IndividualBirds, BirdNames, BirdOption, FeatureToPlot_ColNo, PlotType, DifferentMicrophones);']);
else
    disp(['Script Harini_', PlotOption, '.m does not exist']);
    return;
end



if (strcmp(PlotOption, 'VariabilityForDifferentGroupings'))
    % ======= Now to plot for all birds together ==============================
    % Plots to make
    % 1) Variability across different groups
    % a) Variability within a session
    % b) Variability within a day
    % c) Variability across consecutive days
    % d) Variability within a distance condition
    % e) Variability within a time category

    % Separately for the two types of microphones

    VariabilityTypes = [{'Within session'} {'Within day'} {'Across consecutive days'} {'Within condition'} {'Within a time category'}];
    DifferentMicrophones = [{'BP'} {'MM'} {'HF'}];
    
    for MicrophoneType = 1:length(DifferentMicrophones),
        clear WithinSessionVariability WithinDayVariability AcrossDayVariability WithinConditionVariability WithinTimeCategoryVariability;
        for i = 1:length(BirdNames),
%             if ~((strcmp(BirdOption, BirdNames{i})) || strcmp(BirdOption, 'All'))
%                 continue;
%             end
            MicrophoneIndex = find(cellfun(@length, strfind(IndividualBirds(i).Microphone, DifferentMicrophones{MicrophoneType})));
            if (isempty(MicrophoneIndex))
                continue;
            end
            % First within session variability
            for j = 1:length(IndividualBirds(i).SortedBirdParameters(1).MotifLabels),    
                WithinSessionVariability{i}{j} = [];
                for RecordingDayIndex = 1:length(IndividualBirds(i).RecordingDays),
                    for ConditionIndex = 1:length(IndividualBirds(i).Conditions),
                        MotifLabel = IndividualBirds(i).SortedBirdParameters(1).MotifLabels(j);
                        % Now find all indices that match microphone, condition,
                        % recording day and motif label criterion and not NaN
                        ValidIndices = find((char(IndividualBirds(i).AllSyllableData(:,1)) == MotifLabel) & (IndividualBirds(i).AllMicrophoneIndices == MicrophoneIndex) & (IndividualBirds(i).AllConditionIndices == ConditionIndex) & (IndividualBirds(i).AllRecordingDayIndices == RecordingDayIndex) & (~isnan(IndividualBirds(i).AllSyllableFeatValues(:,FeatureToPlot_ColNo))));
                        if (~isempty(ValidIndices))
                            WithinSessionVariability{i}{j}(end+1) = std(IndividualBirds(i).AllSyllableFeatValues(ValidIndices,FeatureToPlot_ColNo));
                        end
                    end
                end
            end

            % Next within a recording day variability
            for j = 1:length(IndividualBirds(i).SortedBirdParameters(1).MotifLabels),
                WithinDayVariability{i}{j} = [];
                for RecordingDayIndex = 1:length(IndividualBirds(i).RecordingDays),
                    MotifLabel = IndividualBirds(i).SortedBirdParameters(1).MotifLabels(j);
                    % Now find all indices that match microphone,
                    % recording day and motif label criterion and not NaN
                    ValidIndices = find((char(IndividualBirds(i).AllSyllableData(:,1)) == MotifLabel) & (IndividualBirds(i).AllMicrophoneIndices == MicrophoneIndex) & (IndividualBirds(i).AllRecordingDayIndices == RecordingDayIndex) & (~isnan(IndividualBirds(i).AllSyllableFeatValues(:,FeatureToPlot_ColNo))));
                    if (~isempty(ValidIndices))
                        WithinDayVariability{i}{j}(end+1) = std(IndividualBirds(i).AllSyllableFeatValues(ValidIndices,FeatureToPlot_ColNo));
                    end
                end
            end

            % Next varability across continuous recording days
            TempDays = datetime(fliplr(IndividualBirds(i).RecordingDates));
            for j = 1:length(IndividualBirds(i).SortedBirdParameters(1).MotifLabels),
                AcrossDayVariability{i}{j} = [];
                for RecordingDayIndex = 1:length(IndividualBirds(i).RecordingDays)-1,
                    if (diff(datenum(TempDays(RecordingDayIndex:RecordingDayIndex+1))) == 1)
                        MotifLabel = IndividualBirds(i).SortedBirdParameters(1).MotifLabels(j);
                        % Now find all indices that match microphone,
                        % recording day and the next day and motif label criterion and not NaN
                        ValidIndices = find((char(IndividualBirds(i).AllSyllableData(:,1)) == MotifLabel) & (IndividualBirds(i).AllMicrophoneIndices == MicrophoneIndex) & ((IndividualBirds(i).AllRecordingDayIndices == RecordingDayIndex) | (IndividualBirds(i).AllRecordingDayIndices == (RecordingDayIndex+1))) & (~isnan(IndividualBirds(i).AllSyllableFeatValues(:,FeatureToPlot_ColNo))));
                        if (~isempty(ValidIndices))
                            AcrossDayVariability{i}{j}(end+1) = std(IndividualBirds(i).AllSyllableFeatValues(ValidIndices,FeatureToPlot_ColNo));
                        end
                    end
                end
            end

            % Within a distance condition variability
            for j = 1:length(IndividualBirds(i).SortedBirdParameters(1).MotifLabels),
                WithinConditionVariability{i}{j} = [];
                for ConditionIndex = 1:length(IndividualBirds(i).Conditions),
                    MotifLabel = IndividualBirds(i).SortedBirdParameters(1).MotifLabels(j);
                    % Now find all indices that match microphone, condition,
                    % and motif label criterion and not NaN
                    ValidIndices = find((char(IndividualBirds(i).AllSyllableData(:,1)) == MotifLabel) & (IndividualBirds(i).AllMicrophoneIndices == MicrophoneIndex) & (IndividualBirds(i).AllConditionIndices == ConditionIndex) & (~isnan(IndividualBirds(i).AllSyllableFeatValues(:,FeatureToPlot_ColNo))));
                    if (~isempty(ValidIndices))
                        WithinConditionVariability{i}{j}(end+1) = std(IndividualBirds(i).AllSyllableFeatValues(ValidIndices,FeatureToPlot_ColNo));
                    end
                end
            end

            % Within a time category variability
            TimeCategories = unique(IndividualBirds(i).AllSyllableTimeCategories);
            for j = 1:length(IndividualBirds(i).SortedBirdParameters(1).MotifLabels),
                WithinTimeCategoryVariability{i}{j} = [];
                for TimeCategoryIndex = 1:length(TimeCategories),
                    MotifLabel = IndividualBirds(i).SortedBirdParameters(1).MotifLabels(j);
                    % Now find all indices that match microphone, time category
                    % and motif label criterion and not NaN
                    ValidIndices = find((char(IndividualBirds(i).AllSyllableData(:,1)) == MotifLabel) & (IndividualBirds(i).AllMicrophoneIndices == MicrophoneIndex) & (IndividualBirds(i).AllSyllableTimeCategories == TimeCategoryIndex) & (~isnan(IndividualBirds(i).AllSyllableFeatValues(:,FeatureToPlot_ColNo))));
                    if (~isempty(ValidIndices))
                        WithinTimeCategoryVariability{i}{j}(end+1) = std(IndividualBirds(i).AllSyllableFeatValues(ValidIndices,FeatureToPlot_ColNo));
                    end
                end
            end
        end
    
        % Now plot all of this as bar graphs
        % To pool all the data into one big matrix that has individual bird
        % syllables as the rows and the different types of variabilities as
        % columns. Two additional columns to specifiy the bird from which the data
        % came and the syllable # for that bird
        AllVariability = [];
        for i = 1:length(WithinSessionVariability),
            if (~isempty(WithinSessionVariability{i}))
                AllVariability = [AllVariability; [cellfun(@mean, WithinSessionVariability{i})' cellfun(@mean, WithinDayVariability{i})' cellfun(@mean, AcrossDayVariability{i})' cellfun(@mean, WithinConditionVariability{i})' cellfun(@mean, WithinTimeCategoryVariability{i})' (1:1:length(WithinSessionVariability{i}))' ones(length(WithinSessionVariability{i}),1)*i]]; 
            end
        end

        figure;
        hold on;
        Bars = bar(mean(AllVariability(:,1:5)));
        errorbar(mean(AllVariability(:,1:5)), std(AllVariability(:,1:5))/sqrt(size(AllVariability,1)), 'ks');
        set(Bars, 'FaceColor', 'none', 'EdgeColor', 'k');
        for i = 1:size(AllVariability,1),
            plot(AllVariability(i,1:5), [Colours(mod(AllVariability(i,end), length(Colours)) + 1), 'o-']);
        end
        set(gca, 'XTick', 1:1:length(VariabilityTypes), 'XTickLabel', VariabilityTypes, 'XTickLabelRotation', 45);
        ylabel(['Standard deviation of feature value (n=', num2str(size(AllVariability, 1)), ' syllables)']);
        title(['Variability in feature values based on different groupings - ', DifferentMicrophones(MicrophoneType)]);
    end    
end

% ========= Now to do stats on this =======================================

% % Trying n-way anova with 
% % Now to plot for each bird and each motif syllable within that bird
% for i = 1:length(BirdNames),
%     for j = 1:length(IndividualBirds(i).SortedBirdParameters(1).MotifLabels),
%         for k = 1:length(IndividualBirds(i).Microphone),
%             ValidIndices = find((IndividualBirds(i).AllMicrophoneIndices == k) & (char(IndividualBirds(i).AllSyllableData(:,1)) == IndividualBirds(i).SortedBirdParameters(1).MotifLabels(j)) & (~isnan(IndividualBirds(i).AllSyllableFeatValues(:,FeatureToPlot_ColNo))));
%             [p{i}{j}{k}, tabl{i}{j}{k}, stats{i}{j}{k}] = anovan(IndividualBirds(i).AllSyllableFeatValues(ValidIndices,FeatureToPlot_ColNo), {IndividualBirds(i).AllConditionIndices(ValidIndices(:)) IndividualBirds(i).AllRecordingDayIndices(ValidIndices(:)) IndividualBirds(i).AllSyllableTimeCategories(ValidIndices(:))}, 'display', 'off', 'varnames', {'Distance', 'Recording Day', 'Recording Time'});
%         end
%     end
% end

% close all;

% =========== Now to plot bout statistics =================================

if (strcmp(PlotOption, 'BoutStatistics'))
    % First for each bird separately, plot bout length, # of motifs, # of
    % INs, # of syllables, first motif duration, average motif duration,
    % first motif syllable time
    for i = 1:length(IndividualBirds),
        figure;
        for ConditionIndex = 1:max(IndividualBirds(i).BoutStatistics(:,end)),
            % THis is to use only bouts with INs and no calls
            % Indices = find((IndividualBirds(i).BoutStatistics(:,end) == ConditionIndex) & (IndividualBirds(i).BoutStatistics(:,7) > 0) & (IndividualBirds(i).BoutStatistics(:,6) > 0));
            
            % THis is to use all bouts irrespective of whether they have
            % calls, INs or both
            Indices = find((IndividualBirds(i).BoutStatistics(:,end) == ConditionIndex) & (IndividualBirds(i).BoutStatistics(:,7) > 0));

            % If there are no bouts with the criterion specified above,
            % then use all bout
            if (isempty(Indices))
                disp(['Empty indices: ', num2str(i)]);
                Indices = find((IndividualBirds(i).BoutStatistics(:,end) == ConditionIndex) & (IndividualBirds(i).BoutStatistics(:,7) > 0));
            end
            
            SubplotIndex = 0;
            for j = 1:length(IndividualBirds(i).BoutStatisticsColumnNames),
                switch (IndividualBirds(i).BoutStatisticsColumnNames{j})
                    case {'BoutIndex', 'OnlyINBoutOrNot', 'Condition'}
                        continue;
                        
                    otherwise
                        SubplotIndex = SubplotIndex + 1;
                        subplot(3,3,SubplotIndex); % Bout length
                        hold on;
                        errorbar(ConditionIndex-0.1, mean(IndividualBirds(i).BoutStatistics(Indices,j)), std(IndividualBirds(i).BoutStatistics(Indices,j))/sqrt(length(Indices)), 'bs');
%                        errorbar(ConditionIndex+0.1, median(IndividualBirds(i).BoutStatistics(Indices,j)), std(IndividualBirds(i).BoutStatistics(Indices,j))/sqrt(length(Indices)), 'rs');
%                         BootstrappedMeans = bootstrp(10000, @mean, IndividualBirds(i).BoutStatistics(Indices,j));
%                         plot(ConditionIndex, mean(IndividualBirds(i).BoutStatistics(Indices,j)), 'ks');
%                         plot(ones(1,2)*ConditionIndex, [prctile(BootstrappedMeans, 2.5) prctile(BootstrappedMeans, 97.5)], 'k');

                        % Now to plot out the directed bouts separately and
                        % then the female responding bouts separately and
                        % the directed bouts with female responding and not
                        % responding - 3 other categories
                        % directed, female responding
                        DirectedFemaleRespondingBoutIndices = Indices(find((strcmp(IndividualBirds(i).BoutFemaleResponse(Indices), 'R')) & (strcmp(IndividualBirds(i).BoutCategorisation(Indices), 'D'))));
                        if (~isempty(DirectedFemaleRespondingBoutIndices))
                            errorbar(ConditionIndex, mean(IndividualBirds(i).BoutStatistics(DirectedFemaleRespondingBoutIndices,j)), std(IndividualBirds(i).BoutStatistics(DirectedFemaleRespondingBoutIndices,j))/sqrt(length(DirectedFemaleRespondingBoutIndices)), 'rs');
                        end
                        
                        % directed, all but female responding
                        DirectedFemaleNotRespondingBoutIndices = Indices(find((~strcmp(IndividualBirds(i).BoutFemaleResponse(Indices), 'R')) & (strcmp(IndividualBirds(i).BoutCategorisation(Indices), 'D'))));
                        if (~isempty(DirectedFemaleNotRespondingBoutIndices))
                            errorbar(ConditionIndex+0.1, mean(IndividualBirds(i).BoutStatistics(DirectedFemaleNotRespondingBoutIndices,j)), std(IndividualBirds(i).BoutStatistics(DirectedFemaleNotRespondingBoutIndices,j))/sqrt(length(DirectedFemaleNotRespondingBoutIndices)), 'gs');
                        end
                        
                        % undirected, female responding
                        UnDirectedFemaleRespondingBoutIndices = Indices(find((strcmp(IndividualBirds(i).BoutFemaleResponse(Indices), 'R')) & (strcmp(IndividualBirds(i).BoutCategorisation(Indices), 'UN'))));
                        if (~isempty(UnDirectedFemaleRespondingBoutIndices))
                            errorbar(ConditionIndex+0.2, mean(IndividualBirds(i).BoutStatistics(UnDirectedFemaleRespondingBoutIndices,j)), std(IndividualBirds(i).BoutStatistics(UnDirectedFemaleRespondingBoutIndices,j))/sqrt(length(UnDirectedFemaleRespondingBoutIndices)), 'cs');
                        end
                        
                        % undirected, all but female responding
                        UnDirectedFemaleNotRespondingBoutIndices = Indices(find((~strcmp(IndividualBirds(i).BoutFemaleResponse(Indices), 'R')) & (strcmp(IndividualBirds(i).BoutCategorisation(Indices), 'UN'))));
                        if (~isempty(UnDirectedFemaleNotRespondingBoutIndices))
                            errorbar(ConditionIndex+0.3, mean(IndividualBirds(i).BoutStatistics(UnDirectedFemaleNotRespondingBoutIndices,j)), std(IndividualBirds(i).BoutStatistics(UnDirectedFemaleNotRespondingBoutIndices,j))/sqrt(length(UnDirectedFemaleNotRespondingBoutIndices)), 'ks');
                        end

                        
                        if (ConditionIndex == max(IndividualBirds(i).BoutStatistics(:,end)))
                            BootstrappedMeans = bootstrp(10000, @mean, IndividualBirds(i).BoutStatistics(Indices,j)); % Bootstrapped means for undirected and then will plot the 95% confidence intervals for undirected song
                            plot([0 5], ones(1,2)*prctile(BootstrappedMeans, 97.5), 'b--');
                            plot([0 5], ones(1,2)*prctile(BootstrappedMeans, 2.5), 'b--');
                            UndirCIs{j}(i,:) = [prctile(BootstrappedMeans, 2.5) prctile(BootstrappedMeans, 97.5)];
                            
                            BootstrappedMedians = bootstrp(10000, @median, IndividualBirds(i).BoutStatistics(Indices,j)); % Bootstrapped medians for undirected and then will plot the 95% confidence intervals for undirected song
%                             plot([0 5], ones(1,2)*prctile(BootstrappedMedians, 97.5), 'r--');
%                             plot([0 5], ones(1,2)*prctile(BootstrappedMedians, 2.5), 'r--');
                            UndirMedianCIs{j}(i,:) = [prctile(BootstrappedMedians, 2.5) prctile(BootstrappedMedians, 97.5)];
                        end
            
                        MeanFeatureValue{j}(i, ConditionIndex) = mean(IndividualBirds(i).BoutStatistics(Indices,j));
                        MedianFeatureValue{j}(i, ConditionIndex) = median(IndividualBirds(i).BoutStatistics(Indices,j));
                        
                        ylabel(IndividualBirds(i).BoutStatisticsColumnNames{j});
                        set(gca, 'XTick', 1:1:max(IndividualBirds(i).BoutStatistics(:,end)), 'XTickLabel', IndividualBirds(i).Conditions);
                end
            end
        end
    end
    
    % Now to work out significances
    for j = 1:length(IndividualBirds(1).BoutStatisticsColumnNames),
        figure;
        SubplotIndex = 0;
        for i = 1:length(IndividualBirds),
            switch (IndividualBirds(i).BoutStatisticsColumnNames{j})        
                case {'BoutIndex', 'OnlyINBoutOrNot', 'Condition'}
                    continue;
                otherwise
                    % First use kruskal-wallis to check if there is a
                    % difference
                    [p, tabl, stats] = kruskalwallis(IndividualBirds(i).BoutStatistics(:,j), IndividualBirds(i).BoutStatistics(:,end), 'off');
                    SignificanceMatrix{i}{j} = zeros(length(IndividualBirds(i).Conditions));
                    if (p < 0.05)
                        Comparisons = multcompare(stats, 'display', 'off');
                        SigIndices = find(Comparisons(:,end) < 0.05);
                        for k = 1:length(SigIndices),
                            SignificanceMatrix{i}{j}(Comparisons(SigIndices(k),1), Comparisons(SigIndices(k),2)) = 1;
                        end
                    end
                    SubplotIndex = SubplotIndex + 1;
                    subplot(3,3,SubplotIndex);
                    imagesc(SignificanceMatrix{i}{j});
                    set(gca, 'XTick', 1:1:max(IndividualBirds(i).BoutStatistics(:,end)), 'XTickLabel', IndividualBirds(i).Conditions);
                    set(gca, 'YTick', 1:1:max(IndividualBirds(i).BoutStatistics(:,end)), 'YTickLabel', IndividualBirds(i).Conditions);
                    title(IndividualBirds(i).BoutStatisticsColumnNames{j});
                    colormap('gray');
                    caxis([0 1]);
            end
        end
    end
    
    % Plot group summary plots for each of the features
       
    for j = 1:length(IndividualBirds(1).BoutStatisticsColumnNames),
        switch (IndividualBirds(i).BoutStatisticsColumnNames{j})        
            case {'BoutIndex', 'OnlyINBoutOrNot', 'Condition'}
                continue;
              
            otherwise
                figure;
                hold on;
                for i = 1:size(MeanFeatureValue{j},2),
                    Bar{j}(i) = bar(i, mean(MeanFeatureValue{j}(:,i)));
                    set(Bar{j}(i), 'FaceColor', 'none', 'EdgeColor', 'k', 'LineWidth', 2);
                    errorbar(i, mean(MeanFeatureValue{j}(:,i)), std(MeanFeatureValue{j}(:,i))/sqrt(size(MeanFeatureValue{j},1)), 'ks');
                end
                plot(MeanFeatureValue{j}', 'LineWidth', 0.5);
                ylabel(IndividualBirds(i).BoutStatisticsColumnNames{j});
                set(gca, 'XTick', 1:1:max(IndividualBirds(i).BoutStatistics(:,end)), 'XTickLabel', IndividualBirds(i).Conditions);
                title(IndividualBirds(i).BoutStatisticsColumnNames{j});
        end
    end
    
    
    % Plot group summary plots for each of the features relative to
    % undirected
       
    for j = 1:length(IndividualBirds(1).BoutStatisticsColumnNames),
        switch (IndividualBirds(i).BoutStatisticsColumnNames{j})        
            case {'BoutIndex', 'OnlyINBoutOrNot', 'Condition'}
                continue;
              
            otherwise
                TempMeanFeatureValue{j} = MeanFeatureValue{j};
                TempMeanFeatureValue{j} = TempMeanFeatureValue{j}./(repmat(TempMeanFeatureValue{j}(:,end), 1, size(TempMeanFeatureValue{j},2)));
                
                figure;
                hold on;
                for i = 1:size(TempMeanFeatureValue{j},2)-1,
                    errorbar(i, mean(TempMeanFeatureValue{j}(:,i)), std(TempMeanFeatureValue{j}(:,i))/sqrt(size(TempMeanFeatureValue{j},1)), 'ks');
                end
                plot(TempMeanFeatureValue{j}(:,1:end-1)', 'LineWidth', 0.5);
                ylabel(IndividualBirds(i).BoutStatisticsColumnNames{j});
                set(gca, 'XTick', 1:1:max(IndividualBirds(i).BoutStatistics(:,end))-1, 'XTickLabel', IndividualBirds(i).Conditions(1:end-1));
                title('Normalized to undirected');
        end
    end
    
    
    % Now to see which ones are significantly different from undirected
    % bootstrapped confidence intervals and plot this as a 1 (significant)
    % and 0 (not significant) plot
    for j = 1:length(IndividualBirds(1).BoutStatisticsColumnNames),
        switch (IndividualBirds(i).BoutStatisticsColumnNames{j})        
            case {'BoutIndex', 'OnlyINBoutOrNot', 'Condition'}
                continue;
                
            otherwise
                for i = 1:length(IndividualBirds),
                    % Basically this matrix has a -1 if the value is less
                    % than undir, 0 if it is not different from undir and
                    % +1 if it is greater than undir
                    GroupSignificanceMatrix{j}(i,:) = ((MeanFeatureValue{j}(i,:) < UndirCIs{j}(i,1)) * -1) + (MeanFeatureValue{j}(i,:) > UndirCIs{j}(i,2));
                end
        end
        figure;
        imagesc(GroupSignificanceMatrix{j});
        title(IndividualBirds(i).BoutStatisticsColumnNames{j});
        colormap('gray');
        colorbar;
    end
end


disp('Finished plotting');