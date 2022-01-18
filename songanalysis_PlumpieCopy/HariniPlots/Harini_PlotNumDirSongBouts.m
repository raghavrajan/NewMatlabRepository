function [] = Harini_PlotNumDirSongBouts(CSVTextFile, VideoScoringCSVTextFile, InterBoutInterval)

ContinuousFileTime = 30; % in seconds

[BirdParameters, Flag] = ProcessSongData(CSVTextFile, InterBoutInterval);

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

if (sum(Flag) > 0)
    return;
end

Colours = 'rgbcmk';
Symbols = '+o<sd';

% First get details from the CSV text file
disp('Getting header data from CSV file ...');  
[HeaderLine, BirdDetails] = LSINA_GetDetailsFromCSVFile(VideoScoringCSVTextFile);

% Now parse all the lines into the appropriate variables based on the
% header line
disp('Getting data from CSV file ...');
[VideoScoringData] = LSINA_ParseHeaderBirdData(HeaderLine, BirdDetails);

% Get all video scoring filenames into a cell array
for i = 1:length(VideoScoringData),
    VideoScoredBoutFileNames{i} = VideoScoringData(i).FileName;
end

% Now I need to get the first valid song bout and calculate the latency to
% the first song

% I also need the total number of song bouts. In this case, there are song
% bouts that sometimes spread over two files. In this case I need to check
% the time of the file and see whether two files are continuous or not.

% For this, I need to first get the times for each of the files
for i = 1:length(BirdParameters),
    for j = 1:length(BirdParameters(i).SongFileNames),
        ExtensionIndex = strfind(BirdParameters(i).SongFileNames{j}, '.wav');
        FileTimeString = BirdParameters(i).SongFileNames{j}((ExtensionIndex - 6):(ExtensionIndex - 1)); % string from filename in the format hhmmss (h for hour, m for minutes, s for seconds)
        BirdParameters(i).FileTime(j) = (str2double(FileTimeString(1:2))) + (str2double(FileTimeString(3:4))/60) + (str2double(FileTimeString(5:6))/3600); % in hours
    end
end


for i = 1:length(BirdParameters),
    FemaleIntroductionTime = (str2double(BirdParameters(i).Femaleintroductiontime(1:2))) + (str2double(BirdParameters(i).Femaleintroductiontime(3:4))/60) + (str2double(BirdParameters(i).Femaleintroductiontime(5:6))/3600); % in hours
    FemaleExitTime = (str2double(BirdParameters(i).Femaleexittime(1:2))) + (str2double(BirdParameters(i).Femaleexittime(3:4))/60) + (str2double(BirdParameters(i).Femaleexittime(5:6))/3600); % in hours
    SongBouts = find(BirdParameters(i).Bouts(:,7) == 1);
    if (~isempty(SongBouts))
        FirstSongBoutTime = BirdParameters(i).FileTime(BirdParameters(i).Bouts(SongBouts(1),3)) + (BirdParameters(i).Bouts(SongBouts(1),5)/(3600*1000)); % in hours
        for j = BirdParameters(i).Bouts(SongBouts(1),1):BirdParameters(i).Bouts(SongBouts(1),2),
            if (~isempty(find(BirdParameters(i).MotifLabels == char(BirdParameters(i).SyllableData(j,1)))))
                FirstMotifTime = BirdParameters(i).FileTime(BirdParameters(i).Bouts(SongBouts(1),3)) + (BirdParameters(i).SyllableData(j,4)/(3600*1000)); % in hours
                break;
            else
                FirstMotifTime = NaN;
            end
        end
        LatencyToFirstMotif(i) = (FirstMotifTime - FemaleIntroductionTime) * 60; % in minutes
        LatencyToFirstSong(i) = (FirstSongBoutTime - FemaleIntroductionTime) * 60; % in minutes
    else
        LatencyToFirstSong(i) = (FemaleExitTime - FemaleIntroductionTime) * 60; % in minutes
        LatencyToFirstMotif(i) = (FemaleExitTime - FemaleIntroductionTime) * 60; % in minutes
    end
    if (LatencyToFirstSong(i) < 0)
        disp([BirdParameters(i).SongFileNames{1}, ': female intro time = ', BirdParameters(i).Femaleintroductiontime, '; latency = ', num2str(LatencyToFirstSong(i)), ' mins or ', num2str(LatencyToFirstSong(i)*60), ' s']);
    end
end

LatencyToFirstSong(find((LatencyToFirstSong < 0) & (LatencyToFirstSong > -0.5))) = 0;
LatencyToFirstMotif(find((LatencyToFirstMotif < 0) & (LatencyToFirstMotif > -0.5))) = 0;

for i = 1:length(BirdParameters),
    BirdParameters(i).LatencyToFirstSong = LatencyToFirstSong(i);
    BirdParameters(i).LatencyToFirstMotif = LatencyToFirstMotif(i);
end

% For each bird, I also need a cumulative list of syllables for all bouts
for i = 1:length(BirdParameters),
    TotalSyllNum = 0;
    for j = 1:length(BirdParameters(i).SongFileNames),
        Matches = find(BirdParameters(i).Bouts(:,3) == j);
        if (~isempty(Matches))
            BirdParameters(i).BoutCumulativeSyllNum(Matches,:) = BirdParameters(i).Bouts(Matches,1:2) + TotalSyllNum;
            TotalSyllNum = BirdParameters(i).BoutCumulativeSyllNum(end,end);
        end
    end
end

% Get unique BirdNames to get an idea of how many birds there are
for i = 1:length(BirdParameters),
    BirdNames{i} = BirdParameters(i).BirdName;
    MicrophoneTypes{i} = BirdParameters(i).Microphone;
end

UniqueBirds = unique(BirdNames);

% Now to get the number of song bouts in total for each bird for each day
for i = 1:length(UniqueBirds),
    Matches = strmatch(UniqueBirds{i}, BirdNames);
    RecordingDates = [];
    for j = 1:length(Matches),
        Date{i}{j} = BirdParameters(Matches(j)).DataLabel;
        DateNumber{i}(j) = datenum(Date{i}{j}, 'ddmmyy');
        Condition{i}{j} = BirdParameters(Matches(j)).Condition;
    end
    [SortedVals, SortedIndices] = sort(DateNumber{i});
    Date{i} = Date{i}(SortedIndices);
    DateNumber{i} = DateNumber{i}(SortedIndices);
    Condition{i} = Condition{i}(SortedIndices);
    MatchIndices{i} = Matches(SortedIndices);
    
    for j = SortedIndices(:)',
        SkipNextBout = 0;
        ValidSongBouts = find(BirdParameters(Matches(j)).Bouts(:,7) == 1);
        BirdParameters(Matches(j)).BoutPreTime = [];
        BirdParameters(Matches(j)).BoutPostTime = [];
        BirdParameters(Matches(j)).BoutLength = [];
        BirdParameters(Matches(j)).NumMotifs = [];
        BirdParameters(Matches(j)).NumINs = [];
        BirdParameters(Matches(j)).MotifStartTime = [];
        BirdParameters(Matches(j)).FirstMotifDur = [];
        BirdParameters(Matches(j)).AllMotifDurs = [];
        BirdParameters(Matches(j)).BoutVideoScoringTag = [];
        BirdParameters(Matches(j)).BoutFemaleResponseTag = [];
        
        BirdParameters(Matches(j)).BoutContinuousWithNext = [];
        
        
        NumSongBouts{i}(j) = 0;
        BirdParameters(Matches(j)).ValidSongBouts = [];
        if (isempty(ValidSongBouts))
            continue;
        end
        
        if (length(ValidSongBouts) == 1)
            BirdParameters(Matches(j)).ValidSongBouts = ValidSongBouts;
            
            SongFileIndex = BirdParameters(Matches(j)).Bouts(ValidSongBouts,3);
            BoutOnsetTime = BirdParameters(Matches(j)).Bouts(ValidSongBouts,5)/1000; % in seconds
            BoutOffsetTime = BirdParameters(Matches(j)).Bouts(ValidSongBouts,6)/1000; % in seconds
            
            VideoScoreBoutIndex = strmatch(BirdParameters(Matches(j)).SongFileNames{SongFileIndex}, VideoScoredBoutFileNames);
            if (~isempty(VideoScoreBoutIndex))
                ChosenBout = NaN;
                for VideoBouts = VideoScoreBoutIndex(:)',
                    if (abs(VideoScoringData(VideoBouts).BoutOnset - BoutOnsetTime) < 0.1)
                        ChosenBout = VideoBouts;
                        break;
                    end
                end
                if (~isnan(ChosenBout))
                    BirdParameters(Matches(j)).BoutVideoScoringTag{end+1} = VideoScoringData(VideoBouts).DirecteDUNdirected;
                    BirdParameters(Matches(j)).BoutFemaleResponseTag{end+1} = VideoScoringData(VideoBouts).NRorR;
                else
                    BirdParameters(Matches(j)).BoutVideoScoringTag{end+1} = 'NA';
                    BirdParameters(Matches(j)).BoutFemaleResponseTag{end+1} = 'NA';
                end
            else
                BirdParameters(Matches(j)).BoutVideoScoringTag{end+1} = 'NA';
                BirdParameters(Matches(j)).BoutFemaleResponseTag{end+1} = 'NA';
            end

            BirdParameters(Matches(j)).BoutPreTime(end+1) = BoutOnsetTime;
            BirdParameters(Matches(j)).BoutPostTime(end+1) = ContinuousFileTime - BoutOffsetTime;
            
            BirdParameters(Matches(j)).BoutLength(end+1) = BoutOffsetTime - BoutOnsetTime;
            
            BoutLabels = char(BirdParameters(Matches(j)).SyllableData(BirdParameters(Matches(j)).BoutCumulativeSyllNum(ValidSongBouts,1):BirdParameters(Matches(j)).BoutCumulativeSyllNum(ValidSongBouts,2), 1));
            BoutOnsets = BirdParameters(Matches(j)).SyllableData(BirdParameters(Matches(j)).BoutCumulativeSyllNum(ValidSongBouts,1):BirdParameters(Matches(j)).BoutCumulativeSyllNum(ValidSongBouts,2), 4)/1000; % in seconds
            BoutOffsets = BirdParameters(Matches(j)).SyllableData(BirdParameters(Matches(j)).BoutCumulativeSyllNum(ValidSongBouts,1):BirdParameters(Matches(j)).BoutCumulativeSyllNum(ValidSongBouts,2), 5)/1000; % in seconds

            BoutMotifs = strfind(BoutLabels(:)', BirdParameters(Matches(j)).CommonMotifs{1});
            BirdParameters(Matches(j)).NumMotifs(end+1) = length(BoutMotifs);

            FirstMotifSyll = NaN;
            for BoutSylls = 1:length(BoutLabels),
                if (~isempty(find(BirdParameters(Matches(j)).MotifLabels == BoutLabels(BoutSylls))))
                    FirstMotifSyll = BoutSylls;
                    break;
                end
            end
            if (isempty(find(isnan(FirstMotifSyll))))
                BirdParameters(Matches(j)).MotifStartTime(end+1) = BoutOnsets(FirstMotifSyll) - BoutOnsets(1);
            else
                BirdParameters(Matches(j)).MotifStartTime(end+1) = NaN;
            end

            % Now to look in the syllables before the first motif syll and
            % see how many are INs and which ones are the last INs with >
            % 500ms between them.

            if (isempty(find(isnan(FirstMotifSyll))))
                INs = [];
                for INLabel = 1:FirstMotifSyll-1,
                    if (~isempty(find(BirdParameters(Matches(j)).INLabels == BoutLabels(INLabel))))
                        INs(end+1) = INLabel;
                    end
                end
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
                BirdParameters(Matches(j)).NumINs(end+1) = length(INs);
            else
                BirdParameters(Matches(j)).NumINs(end+1) = NaN;
            end



            if (~isempty(BoutMotifs))
                BirdParameters(Matches(j)).FirstMotifDur(end+1) = 1000*(BoutOffsets(BoutMotifs(1) + length(BirdParameters(Matches(j)).CommonMotifs{1}) - 1) - BoutOnsets(BoutMotifs(1)));
                BirdParameters(Matches(j)).AllMotifDurs{end+1} = 1000*(BoutOffsets(BoutMotifs + length(BirdParameters(Matches(j)).CommonMotifs{1}) - 1) - BoutOnsets(BoutMotifs));
            else
                BirdParameters(Matches(j)).FirstMotifDur(end+1) = NaN;
                BirdParameters(Matches(j)).AllMotifDurs{end+1} = NaN;
            end
            
            
            BirdParameters(Matches(j)).BoutContinuousWithNext(end+1) = 0;
            NumSongBouts{i}(j) = 1;
        else
            for k = ValidSongBouts(:)',
                if (SkipNextBout == 1)
                    SkipNextBout = 0;
                    continue;
                end
                
                SongFileIndex = BirdParameters(Matches(j)).Bouts(k,3);
                BoutOnsetTime = BirdParameters(Matches(j)).Bouts(k,5)/1000; % in seconds
                BoutOffsetTime = BirdParameters(Matches(j)).Bouts(k,6)/1000; % in seconds
                
                VideoScoreBoutIndex = strmatch(BirdParameters(Matches(j)).SongFileNames{SongFileIndex}, VideoScoredBoutFileNames);
                if (~isempty(VideoScoreBoutIndex))
                    ChosenBout = NaN;
                    for VideoBouts = VideoScoreBoutIndex(:)',
                        if (abs(VideoScoringData(VideoBouts).BoutOnset - BoutOnsetTime) < 0.1)
                            ChosenBout = VideoBouts;
                            break;
                        end
                    end
                    if (~isnan(ChosenBout))
                        BirdParameters(Matches(j)).BoutVideoScoringTag{end+1} = VideoScoringData(VideoBouts).DirecteDUNdirected;
                        BirdParameters(Matches(j)).BoutFemaleResponseTag{end+1} = VideoScoringData(VideoBouts).NRorR;
                    else
                        BirdParameters(Matches(j)).BoutVideoScoringTag{end+1} = 'NA';
                        BirdParameters(Matches(j)).BoutFemaleResponseTag{end+1} = 'NA';
                    end
                else
                    BirdParameters(Matches(j)).BoutVideoScoringTag{end+1} = 'NA';
                    BirdParameters(Matches(j)).BoutFemaleResponseTag{end+1} = 'NA';
                end
            
                BirdParameters(Matches(j)).BoutPreTime(end+1) = BoutOnsetTime;
                BirdParameters(Matches(j)).BoutPostTime(end+1) = ContinuousFileTime - BoutOffsetTime;
                
                BirdParameters(Matches(j)).BoutLength(end+1) = BoutOffsetTime - BoutOnsetTime;

                BoutLabels = char(BirdParameters(Matches(j)).SyllableData(BirdParameters(Matches(j)).BoutCumulativeSyllNum(k,1):BirdParameters(Matches(j)).BoutCumulativeSyllNum(k,2), 1));
                BoutOnsets = BirdParameters(Matches(j)).SyllableData(BirdParameters(Matches(j)).BoutCumulativeSyllNum(k,1):BirdParameters(Matches(j)).BoutCumulativeSyllNum(k,2), 4)/1000; % in seconds
                BoutOffsets = BirdParameters(Matches(j)).SyllableData(BirdParameters(Matches(j)).BoutCumulativeSyllNum(k,1):BirdParameters(Matches(j)).BoutCumulativeSyllNum(k,2), 5)/1000; % in seconds

                BoutMotifs = strfind(BoutLabels(:)', BirdParameters(Matches(j)).CommonMotifs{1});
                BirdParameters(Matches(j)).NumMotifs(end+1) = length(BoutMotifs);

                FirstMotifSyll = NaN;
                for BoutSylls = 1:length(BoutLabels),
                    if (~isempty(find(BirdParameters(Matches(j)).MotifLabels == BoutLabels(BoutSylls))))
                        FirstMotifSyll = BoutSylls;
                        break;
                    end
                end
                if (isempty(find(isnan(FirstMotifSyll))))
                    BirdParameters(Matches(j)).MotifStartTime(end+1) = BoutOnsets(FirstMotifSyll) - BoutOnsets(1);
                else
                    BirdParameters(Matches(j)).MotifStartTime(end+1) = NaN;
                end

                % Now to look in the syllables before the first motif syll and
                % see how many are INs and which ones are the last INs with >
                % 500ms between them.

                if (isempty(find(isnan(FirstMotifSyll))))
                    INs = [];
                    for INLabel = 1:FirstMotifSyll-1,
                        if (~isempty(find(BirdParameters(Matches(j)).INLabels == BoutLabels(INLabel))))
                            INs(end+1) = INLabel;
                        end
                    end
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
                    BirdParameters(Matches(j)).NumINs(end+1) = length(INs);
                else
                    BirdParameters(Matches(j)).NumINs(end+1) = NaN;
                end



                if (~isempty(BoutMotifs))
                    BirdParameters(Matches(j)).FirstMotifDur(end+1) = 1000*(BoutOffsets(BoutMotifs(1) + length(BirdParameters(Matches(j)).CommonMotifs{1}) - 1) - BoutOnsets(BoutMotifs(1)));
                    BirdParameters(Matches(j)).AllMotifDurs{end+1} = 1000*(BoutOffsets(BoutMotifs + length(BirdParameters(Matches(j)).CommonMotifs{1}) - 1) - BoutOnsets(BoutMotifs));
                else
                    BirdParameters(Matches(j)).FirstMotifDur(end+1) = NaN;
                    BirdParameters(Matches(j)).AllMotifDurs{end+1} = NaN;
                end
                
                if ((BoutOnsetTime >= InterBoutInterval/1000) && ((ContinuousFileTime - BoutOffsetTime) >= InterBoutInterval/1000))
                    NumSongBouts{i}(j) = NumSongBouts{i}(j) + 1;
                    BirdParameters(Matches(j)).ValidSongBouts(end+1) = k;
                    BirdParameters(Matches(j)).BoutContinuousWithNext(end+1) = 0;
                else
                    if ((ContinuousFileTime - BoutOffsetTime) >= InterBoutInterval/1000)
                        NumSongBouts{i}(j) = NumSongBouts{i}(j) + 1;
                        BirdParameters(Matches(j)).ValidSongBouts(end+1) = k;
                        BirdParameters(Matches(j)).BoutContinuousWithNext(end+1) = 0;
                    else
                        % First check whether there is a next file. If not,
                        % just include the bout and continue
                        if (k == ValidSongBouts(end))
                            NumSongBouts{i}(j) = NumSongBouts{i}(j) + 1;
                            BirdParameters(Matches(j)).ValidSongBouts(end+1) = k;
                            BirdParameters(Matches(j)).BoutContinuousWithNext(end+1) = 0;
                            continue;
                        end
                        
                        % Now to check whether next file in the
                        % chronological order is the next file based on
                        % time too. If no, then count as a bout.
                        % Now to check if next valid song bout is in
                        % the next file and is close to the start
                        CurrentBout = find(ValidSongBouts == k);
                        NextBout = ValidSongBouts(CurrentBout + 1);
                        if ((BirdParameters(Matches(j)).Bouts(NextBout,3) - BirdParameters(Matches(j)).Bouts(k,3)) == 1)
                            % Now if the time before the start of the
                            % next bout and the end of the current bout
                            % is less than 2s, then merge the two bouts
                            CurrentBoutOffsetTime = BirdParameters(Matches(j)).FileTime(SongFileIndex) + BoutOffsetTime/3600;
                            NextBoutOnsetTime = BirdParameters(Matches(j)).FileTime(BirdParameters(Matches(j)).Bouts(NextBout,3)) + BirdParameters(Matches(j)).Bouts(NextBout,5)/(1000*3600);
                            if ((NextBoutOnsetTime - CurrentBoutOffsetTime)*3600 < InterBoutInterval/1000)
                                NumSongBouts{i}(j) = NumSongBouts{i}(j) + 1;
                                BirdParameters(Matches(j)).ValidSongBouts(end+1) = k;
                                BirdParameters(Matches(j)).BoutContinuousWithNext(end+1) = 1;
                                SkipNextBout = 1;
                            else
                                NumSongBouts{i}(j) = NumSongBouts{i}(j) + 1;
                                BirdParameters(Matches(j)).ValidSongBouts(end+1) = k;
                                BirdParameters(Matches(j)).BoutContinuousWithNext(end+1) = 0;
                            end
                        else
                            NumSongBouts{i}(j) = NumSongBouts{i}(j) + 1;
                            BirdParameters(Matches(j)).ValidSongBouts(end+1) = k;
                            BirdParameters(Matches(j)).BoutContinuousWithNext(end+1) = 0;
                        end
                    end
                end
            end
        end
    end
end

% Now to plot data for each individual bird
for i = 1:length(UniqueBirds),
    BirdStats(i) = Harini_PlotIndividualBirdData(UniqueBirds{i}, BirdParameters(MatchIndices{i}), Date{i}, DateNumber{i}, Condition{i});
end

close all;

% Now to get number of INs, motifs, etc. for all bout types separately.

for i = 1:length(BirdStats),
    [FinalBirdStats(i)] = Harini_GetIndividualBirdNumINs_Motifs_All(BirdStats(i), InterBoutInterval, i);
end

% Now put in code here to make all data visualization plots
Harini_VisualizeData(FinalBirdStats);


% Now put data into a table and then display it.
DataSummaryTable = table;
DataSummaryTable.BirdNames = UniqueBirds';
DataSummaryTable.NumRecordingDays = [BirdStats.NumRecordingDays]';
DataSummaryTable.TotalNumDaysFromFirstToLast = [BirdStats.TotalNumDaysFromFirstToLast]';
DataSummaryTable.TotalNumSessions = [BirdStats.TotalNumSessions]';
DataSummaryTable.TotalNumSongBouts = [BirdStats.TotalNumSongBouts]';
DataSummaryTable.PercentContinuousSongBouts = [BirdStats.PercentContinousSongBouts]';
for i = 1:length(BirdStats),
    DataSummaryTable.PercentVideoScoredBouts(i) = 100 * (sum(BirdStats(i).NumBoutTypesPerCondition(:)))/(sum(BirdStats(i).NumAllSongBouts) - sum(BirdStats(i).NumContinuousSongBouts));
end

for i = 1:length(UniqueBirds),
    TempString = [];
    for j = 1:length(BirdStats(i).UniqueConditions),
        TempString = [TempString BirdStats(i).UniqueConditions{j}(:)'];
        if (j ~= length(BirdStats(i).UniqueConditions))
            TempString = [TempString ', '];
        end
    end
    DataSummaryTable.UniqueConditions{i} = TempString;
end
writetable(DataSummaryTable, 'TempDataSummaryTable.txt', 'DeLimiter', ';');

% Now to find out how many dir song bouts there are for each condition
% across all birds

ConditionSpecificBoutTables = table;
AllConditionTypes = {'L0' 'L1' 'L2' 'L3' 'L4' 'UN'};
BoutTypes = {'D' 'DUN' 'UN'};
for j = 1:length(AllConditionTypes),
    eval(['ConditionSpecificBoutTables.', AllConditionTypes{j}, ' = table;']);
end
for k = 1:length(BoutTypes),
    BoutTypesNumBouts{k} = [];
end
for i = 1:length(BirdStats),
    for j = 1:length(AllConditionTypes),
        Indices = strmatch(AllConditionTypes{j}, BirdStats(i).Condition, 'exact');
        for k = 1:length(BoutTypes),
            if (isempty(Indices))
                BoutTypesNumBouts{k}(i,j) = NaN;
                eval(['ConditionSpecificBoutTables.', AllConditionTypes{j}, '.', BoutTypes{k}, '(', num2str(i), ') = NaN;']);
            else
                BoutTypesNumBouts{k}(i,j) = 100 * sum(BirdStats(i).NumBoutTypes(Indices, k))/sum(sum(BirdStats(i).NumBoutTypes(Indices, :)));
                eval(['ConditionSpecificBoutTables.', AllConditionTypes{j}, '.', BoutTypes{k}, '(', num2str(i), ') = ', num2str(sum(BirdStats(i).NumBoutTypes(Indices, k))), ';']);
            end
        end
    end
end

% Now get data for latency to first song
for i = 1:length(AllConditionTypes),
    for j = 1:length(BirdStats),
        Indices = strmatch(AllConditionTypes{i}, BirdStats(j).UniqueConditions);
        if (~isempty(Indices))
            AllBirdsMedianLatencyToFirstSong(j,i) = median(BirdStats(j).LatencyToFirstSongPerCondition{Indices});
            AllBirdsIQRLatencyToFirstSong(j,i) = iqr(BirdStats(j).LatencyToFirstSongPerCondition{Indices});
        else
            AllBirdsMedianLatencyToFirstSong(j,i) = NaN;
            AllBirdsIQRLatencyToFirstSong(j,i) = NaN;
        end
    end
end

% Now get data for # INs, # motifs, first motif dur - only for 
for i = 1:length(AllConditionTypes),
    for j = 1:length(FinalBirdStats),
        Indices = strmatch(AllConditionTypes{i}, FinalBirdStats(j).UniqueConditions);
        if (~isempty(Indices))
            NumDSongBouts(j,i) = length(FinalBirdStats(j).ConditionSpecific_DSong_INs{Indices});
            if (length(FinalBirdStats(j).ConditionSpecific_DSong_INs{Indices}) >= 5)
                AllBirdsMeanNumINs(j,i) = mean(FinalBirdStats(j).ConditionSpecific_DSong_INs{Indices});
                AllBirdsMeanNumMotifs(j,i) = mean(FinalBirdStats(j).ConditionSpecific_DSong_NumMotifs{Indices});
                if (length(~isnan(FinalBirdStats(j).ConditionSpecific_DSong_FirstMotifDur{Indices})) >= 5)
                    AllBirdsMeanFirstMotifDur(j,i) = nanmean(FinalBirdStats(j).ConditionSpecific_DSong_FirstMotifDur{Indices});
                else
                    AllBirdsMeanFirstMotifDur(j,i) = NaN;
                end
            else
                AllBirdsMeanNumINs(j,i) = NaN;
                AllBirdsMeanNumMotifs(j,i) = NaN;
                AllBirdsMeanFirstMotifDur(j,i) = NaN;
            end
        else
            NumDSongBouts(j,i) = 0;
            AllBirdsMeanNumINs(j,i) = NaN;
            AllBirdsMeanNumMotifs(j,i) = NaN;
            AllBirdsMeanFirstMotifDur(j,i) = NaN;
        end
    end
end

FinalFigureDir = '/home/raghav/StudentRelated/Harini/Manuscript/JEBManuscript/Figures';

Distances = [0 20 60 110 165 185];
DistancesMatrix = repmat(Distances(:)', size(AllBirdsMedianLatencyToFirstSong,1),1);

% Temporary variables
TempAllBirdsMeanNumINs = AllBirdsMeanNumINs;
TempAllBirdsMeanNumMotifs = AllBirdsMeanNumMotifs;
TempAllBirdsMeanFirstMotifDur = AllBirdsMeanFirstMotifDur;

% % Normalize to undirected by subtracting IN number
% AllBirdsMeanNumINs = AllBirdsMeanNumINs - repmat(AllBirdsMeanNumINs(:,end), 1, size(AllBirdsMeanNumINs,2));
% AllBirdsMeanNumMotifs = AllBirdsMeanNumMotifs - repmat(AllBirdsMeanNumMotifs(:,end), 1, size(AllBirdsMeanNumMotifs,2));
% AllBirdsMeanFirstMotifDur = AllBirdsMeanFirstMotifDur - repmat(AllBirdsMeanFirstMotifDur(:,end), 1, size(AllBirdsMeanFirstMotifDur,2));

% Normalize to undirected by making it relative to undirected
AllBirdsMeanNumINs = AllBirdsMeanNumINs./repmat(AllBirdsMeanNumINs(:,end), 1, size(AllBirdsMeanNumINs,2));
AllBirdsMeanNumMotifs = AllBirdsMeanNumMotifs./repmat(AllBirdsMeanNumMotifs(:,end), 1, size(AllBirdsMeanNumMotifs,2));
AllBirdsMeanFirstMotifDur = AllBirdsMeanFirstMotifDur./repmat(AllBirdsMeanFirstMotifDur(:,end), 1, size(AllBirdsMeanFirstMotifDur,2));


for i = 1:length(Distances),
    if (i == length(Distances))
        DistanceString{i} = 'UN';
    else
        DistanceString{i} = num2str(Distances(i));
    end
end


% Now plot the mean IN number, mean motif # and mean motif duration plots
figure;
p = panel();
set(gcf, 'Color', 'w');
set(gcf, 'Position', [627 641 600 350]);
p.pack({1});
p(1).select();
hold on;
for i = 1:size(AllBirdsMeanNumINs, 1),
    NotNaNValues = find(~isnan(AllBirdsMeanNumINs(i,1:5)));
    plot(Distances(NotNaNValues), AllBirdsMeanNumINs(i, NotNaNValues), 'ko-', 'Color', [0.75 0.75 0.75]);
end
for i = 1:6,
    errorbar(Distances(i), nanmean(AllBirdsMeanNumINs(:,i)), nanstd(AllBirdsMeanNumINs(:,i))/sqrt(length(find(~isnan(AllBirdsMeanNumINs(:,i))))), 'ks', 'MarkerSize', 8, 'MarkerFaceColor', 'k');
end
p.fontsize = 14;
p.fontname = 'Arial';
set(gca, 'XTick', Distances, 'XTickLabel', DistanceString);
ylabel('Mean number of INs relative to UN');
axis tight;
Temp = axis;
Temp = [(Distances(1) - 5) (Distances(end)+5) 1.2*Temp(3) 1.2*Temp(4)];
axis(Temp);
plot(Temp(1:2), zeros(1,2), 'k--');
text(-33, 1.02*Temp(4), 'A', 'FontSize', 12, 'FontName', 'Arial');
text(172, Temp(3), '//', 'FontSize', 12, 'FontName', 'Arial');
[r, prob] = corr(DistancesMatrix(1:size(DistancesMatrix,1)*5)', AllBirdsMeanNumINs(1:size(DistancesMatrix,1)*5)', 'rows', 'complete');
text(60, 0.9*Temp(4), {['r = ', num2str(r)]; ['p = ', num2str(prob)]}, 'FontSize', 14, 'FontName', 'Arial');
xlabel('Distance from the female');
p.margintop = 10;
title('Normalized # of INs');
OutputDir = '/home/raghav/StudentRelated/Harini/Manuscript/JEBManuscript/Figures/';
set(gcf, 'PaperPositionMode', 'auto');
print(fullfile(OutputDir, ['NumberofINs.png']), '-dpng');

% Now plot the mean IN number, mean motif # and mean motif duration plots
figure;
p = panel();
set(gcf, 'Color', 'w');
set(gcf, 'Position', [627 641 600 500]);
p.pack({1});
p(1).select();
hold on;
for i = 1:size(AllBirdsMeanNumMotifs, 1),
    NotNaNValues = find(~isnan(AllBirdsMeanNumMotifs(i,1:5)));
    plot(Distances(NotNaNValues), AllBirdsMeanNumMotifs(i, NotNaNValues), 'ko-', 'Color', [0.75 0.75 0.75]);
end
for i = 1:6,
    errorbar(Distances(i), nanmean(AllBirdsMeanNumMotifs(:,i)), nanstd(AllBirdsMeanNumMotifs(:,i))/sqrt(length(find(~isnan(AllBirdsMeanNumMotifs(:,i))))), 'ks', 'MarkerSize', 8, 'MarkerFaceColor', 'k');
end
p.fontsize = 14;
p.fontname = 'Arial';
set(gca, 'XTick', Distances, 'XTickLabel', DistanceString);
ylabel('Mean number of motifs/bout relative to UN');
axis tight;
Temp = axis;
Temp = [(Distances(1) - 5) (Distances(end)+5) 1.2*Temp(3) 1.2*Temp(4)];
axis(Temp);
plot(Temp(1:2), zeros(1,2), 'k--');
text(-33, 1.02*Temp(4), 'A', 'FontSize', 12, 'FontName', 'Arial');
text(172, Temp(3), '//', 'FontSize', 12, 'FontName', 'Arial');
[r, prob] = corr(DistancesMatrix(1:size(DistancesMatrix,1)*5)', AllBirdsMeanNumMotifs(1:size(DistancesMatrix,1)*5)', 'rows', 'complete');
text(60, 0.9*Temp(4), {['r = ', num2str(r)]; ['p = ', num2str(prob)]}, 'FontSize', 14, 'FontName', 'Arial');
xlabel('Distance from the female');
p.margintop = 10;
title('Normalized # of motifs/bout');
OutputDir = '/home/raghav/StudentRelated/Harini/Manuscript/JEBManuscript/Figures/';
set(gcf, 'PaperPositionMode', 'auto');
print(fullfile(OutputDir, ['Numberofmotifs.png']), '-dpng');

% Now plot the mean IN number, mean motif # and mean motif duration plots
figure;
p = panel();
set(gcf, 'Color', 'w');
set(gcf, 'Position', [627 641 600 500]);
p.pack({1});
p(1).select();
hold on;
for i = 1:size(AllBirdsMeanFirstMotifDur, 1),
    NotNaNValues = find(~isnan(AllBirdsMeanFirstMotifDur(i,1:5)));
    plot(Distances(NotNaNValues), AllBirdsMeanFirstMotifDur(i, NotNaNValues), 'ko-', 'Color', [0.75 0.75 0.75]);
end
for i = 1:6,
    errorbar(Distances(i), nanmean(AllBirdsMeanFirstMotifDur(:,i)), nanstd(AllBirdsMeanFirstMotifDur(:,i))/sqrt(length(find(~isnan(AllBirdsMeanFirstMotifDur(:,i))))), 'ks', 'MarkerSize', 8, 'MarkerFaceColor', 'k');
end
p.fontsize = 14;
p.marginleft = 20;
p.fontname = 'Arial';
set(gca, 'XTick', Distances, 'XTickLabel', DistanceString);
ylabel('First motif dur (msec) relative to UN');
axis tight;
Temp = axis;
Temp = [(Distances(1) - 5) (Distances(end)+5) 1.2*Temp(3) 1.2*Temp(4)];
axis(Temp);
plot(Temp(1:2), zeros(1,2), 'k--');
text(-33, 1.02*Temp(4), 'A', 'FontSize', 12, 'FontName', 'Arial');
text(172, Temp(3), '//', 'FontSize', 12, 'FontName', 'Arial');
[r, prob] = corr(DistancesMatrix(1:size(DistancesMatrix,1)*5)', AllBirdsMeanFirstMotifDur(1:size(DistancesMatrix,1)*5)', 'rows', 'complete');
text(60, 0.9*Temp(4), {['r = ', num2str(r)]; ['p = ', num2str(prob)]}, 'FontSize', 14, 'FontName', 'Arial');
xlabel('Distance from the female');
p.margintop = 10;
title('Normalized first motif duration (msec)');
OutputDir = '/home/raghav/StudentRelated/Harini/Manuscript/JEBManuscript/Figures/';
set(gcf, 'PaperPositionMode', 'auto');
print(fullfile(OutputDir, ['FirstMotifDur.png']), '-dpng');

% Convert latency and iqr of latency to first songs to secs instead of mins
AllBirdsMedianLatencyToFirstSong = AllBirdsMedianLatencyToFirstSong * 60;
AllBirdsIQRLatencyToFirstSong = AllBirdsIQRLatencyToFirstSong * 60;

figure;
set(gcf, 'Color', 'w');
set(gcf, 'Units', 'centimeters');
set(gcf, 'Position', [1 16 9 21]);
p = panel();
p.pack({1/3 1/3 1/3});
% First plot median latency to first songs
p(1).select();
hold on;
for i = 1:size(AllBirdsMedianLatencyToFirstSong, 1),
    NotNaNValues = find(~isnan(AllBirdsMedianLatencyToFirstSong(i,:)));
    plot(Distances(NotNaNValues), AllBirdsMedianLatencyToFirstSong(i, NotNaNValues), 'ko-', 'Color', [0.75 0.75 0.75]);
end
for i = 1:6,
    errorbar(Distances(i), nanmean(AllBirdsMedianLatencyToFirstSong(:,i)), nanstd(AllBirdsMedianLatencyToFirstSong(:,i))/sqrt(length(find(~isnan(AllBirdsMedianLatencyToFirstSong(:,i))))), 'ks', 'MarkerSize', 8, 'MarkerFaceColor', 'k');
end
p.fontsize = 8;
p.fontname = 'Arial';
set(gca, 'XTick', Distances, 'XTickLabel', DistanceString);
ylabel('Median latency to first song (sec)');
set(gca, 'YScale', 'log');
axis tight;
Temp = axis;
Temp = [(Distances(1) - 5) (Distances(end)+5) 0.9*Temp(3) 1.2*Temp(4)];
axis(Temp);
text(-33, 1.02*Temp(4), 'A', 'FontSize', 12, 'FontName', 'Arial');
text(172, Temp(3), '//', 'FontSize', 12, 'FontName', 'Arial');
[r, prob] = corr(DistancesMatrix(1:size(DistancesMatrix,1)*5)', AllBirdsMedianLatencyToFirstSong(1:size(DistancesMatrix,1)*5)', 'rows', 'complete');
text(60, 10*Temp(3), {['r = ', num2str(r)]; ['p = ', num2str(prob)]}, 'FontSize', 8, 'FontName', 'Arial');

% Next plot iqr to first songs
p(2).select();
hold on;
for i = 1:size(AllBirdsIQRLatencyToFirstSong, 1),
    NotNaNValues = find(~isnan(AllBirdsIQRLatencyToFirstSong(i,:)));
    plot(Distances(NotNaNValues), AllBirdsIQRLatencyToFirstSong(i, NotNaNValues), 'ko-', 'Color', [0.75 0.75 0.75]);
end
for i = 1:6,
    errorbar(Distances(i), nanmean(AllBirdsIQRLatencyToFirstSong(:,i)), nanstd(AllBirdsIQRLatencyToFirstSong(:,i))/sqrt(length(find(~isnan(AllBirdsIQRLatencyToFirstSong(:,i))))), 'ks', 'MarkerSize', 8, 'MarkerFaceColor', 'k');
end
p.fontsize = 8;
p.fontname = 'Arial';
set(gca, 'XTick', Distances, 'XTickLabel', DistanceString);
ylabel('IQR of latency to first song (sec)');
set(gca, 'YScale', 'log');
axis tight;
Temp = axis;
Temp = [(Distances(1) - 5) (Distances(end)+5) 0.9*Temp(3) 1.2*Temp(4)];
axis(Temp);
text(-33, 1.02*Temp(4), 'B', 'FontSize', 12, 'FontName', 'Arial');
text(172, Temp(3), '//', 'FontSize', 12, 'FontName', 'Arial');
[r, prob] = corr(DistancesMatrix(1:size(DistancesMatrix,1)*5)', AllBirdsIQRLatencyToFirstSong(1:size(DistancesMatrix,1)*5)', 'rows', 'complete');
text(60, Temp(3)*10, {['r = ', num2str(r)]; ['p = ', num2str(prob)]}, 'FontSize', 8, 'FontName', 'Arial');

% Last plot number of directed songs
p(3).select();
hold on;
for i = 1:size(BoutTypesNumBouts{1}, 1),
    NotNaNValues = find(~isnan(BoutTypesNumBouts{1}(i,:)));
    plot(Distances(NotNaNValues), BoutTypesNumBouts{1}(i, NotNaNValues), 'ko-', 'Color', [0.75 0.75 0.75]);
end
for i = 1:6,
    errorbar(Distances(i), nanmean(BoutTypesNumBouts{1}(:,i)), nanstd(BoutTypesNumBouts{1}(:,i))/sqrt(length(find(~isnan(BoutTypesNumBouts{1}(:,i))))), 'ks', 'MarkerSize', 8, 'MarkerFaceColor', 'k');
end
p.fontsize = 8;
p.fontname = 'Arial';
xlabel('Distance from female cage (cm)');
set(gca, 'XTick', Distances, 'XTickLabel', DistanceString);
ylabel('Proportion of courtship songs (%)');
axis([(Distances(1) - 5) (Distances(end)+5) 0 100]);
text(-33, 103, 'C', 'FontSize', 12, 'FontName', 'Arial');
text(172, 0, '//', 'FontSize', 12, 'FontName', 'Arial');
[r, prob] = corr(DistancesMatrix(1:size(DistancesMatrix,1)*5)', BoutTypesNumBouts{1}(1:size(DistancesMatrix,1)*5)', 'rows', 'complete');
text(110, 90, {['r = ', num2str(r)]; ['p = ', num2str(prob)]}, 'FontSize', 8, 'FontName', 'Arial');

p.margintop = 5;
p.marginbottom = 10;
p.marginleft = 12;
p.marginright = 5;
set(gcf, 'PaperPositionMode', 'auto');

FileDateTimeStr = datestr(now, 'dd_mm_yyyy.hh_MM_ss');
print(fullfile(FinalFigureDir, ['ProportionDirSongs.', FileDateTimeStr '.png']), '-dpng', '-r300');

% In supp. figure
% Now plot number of undirected songs and DUN songs


figure;
set(gcf, 'Color', 'w');
set(gcf, 'Units', 'centimeters');
set(gcf, 'Position', [1 16 10.5 16]);
p = panel();
p.pack({1/2 1/2});

p(1).select();
hold on;
for i = 1:size(BoutTypesNumBouts{3}, 1),
    NotNaNValues = find(~isnan(BoutTypesNumBouts{3}(i,:)));
    plot(Distances(NotNaNValues), BoutTypesNumBouts{3}(i, NotNaNValues), 'ko-', 'Color', [0.75 0.75 0.75]);
end
for i = 1:6,
    errorbar(Distances(i), nanmean(BoutTypesNumBouts{3}(:,i)), nanstd(BoutTypesNumBouts{3}(:,i))/sqrt(length(find(~isnan(BoutTypesNumBouts{3}(:,i))))), 'ks', 'MarkerSize', 8, 'MarkerFaceColor', 'k');
end
ylabel('Proportion of undirected songs (%)');
axis([(Distances(1) - 5) (Distances(end)+5) 0 100]);
text(-28, 103, 'A', 'FontSize', 12, 'FontName', 'Arial');
text(172, 0, '//', 'FontSize', 12, 'FontName', 'Arial');
set(gca, 'XTick', Distances, 'XTickLabel', []);

p(2).select();
hold on;
for i = 1:size(BoutTypesNumBouts{2}, 1),
    NotNaNValues = find(~isnan(BoutTypesNumBouts{2}(i,:)));
    plot(Distances(NotNaNValues), BoutTypesNumBouts{2}(i, NotNaNValues), 'ko-', 'Color', [0.75 0.75 0.75]);
end
for i = 1:6,
    errorbar(Distances(i), nanmean(BoutTypesNumBouts{2}(:,i)), nanstd(BoutTypesNumBouts{2}(:,i))/sqrt(length(find(~isnan(BoutTypesNumBouts{2}(:,i))))), 'ks', 'MarkerSize', 8, 'MarkerFaceColor', 'k');
end
xlabel('Distance from female cage (cm)');
ylabel('Proportion of ambiguous (D/UN) songs (%)');
axis([(Distances(1) - 5) (Distances(end)+5) 0 100]);
set(gca, 'XTick', Distances, 'XTickLabel', DistanceString);

p.de.margin = 8;
p.fontsize = 8;
p.fontname = 'Arial';
p.margintop = 5;
p.marginbottom = 10;
p.marginleft = 12;
p.marginright = 5;
text(-28, 103, 'B', 'FontSize', 12, 'FontName', 'Arial');
text(172, 0, '//', 'FontSize', 12, 'FontName', 'Arial');

set(gcf, 'PaperPositionMode', 'auto');

FileDateTimeStr = datestr(now, 'dd_mm_yyyy.hh_MM_ss');
print(fullfile(FinalFigureDir, ['ProportionUnDir_DUN_Songs.', FileDateTimeStr '.png']), '-dpng', '-r300');

% close all;

% Next to plot the number of INs, number of motifs and first motif duration
% for all distances for D bouts alone
disp('Finished');
