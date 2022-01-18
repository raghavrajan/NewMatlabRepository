function [] = Convert_RavenData_To_NotMatFile(RawDataDir, SongFileList, FileType, RavenFileDir, RavenFile)
% ========================================================================================
% This is is a script to convert labels that have been put using Raven
% Software to labels in the .not.mat format so that it can also be opened
% with AutoSongSegmentLabel
% Usage:
%
% ========================================================================================
 
% Some common parameters
MinInt = 5; % in ms - minimum inter syllable duration
MinDur = 10; % in ms - minimum duration of a syllable

% First read in the raw files and get the duration of each file
Fid = fopen(fullfile(RawDataDir, SongFileList), 'r');
SongFiles = textscan(Fid, '%s', 'DeLimiter', '\n');
SongFiles = SongFiles{1};
fclose(Fid);

for i = 1:length(SongFiles),
    [RawData, Fs] = GetData(RawDataDir, SongFiles{i}, FileType, 0);
    FileLen(i) = length(RawData)/Fs;
end

% Now using the file durations, I will make two arrays that have the file
% onsets and file offsets
FileOffsets = cumsum(FileLen)*1000; % convert to ms
FileOffsets = FileOffsets(:);
FileOnsets = [1/Fs; FileOffsets(1:end-1)];

% Now read in the Raven file
Fid = fopen(fullfile(RavenFileDir, RavenFile), 'r');
RavenLabels = textscan(Fid, '%s', 'DeLimiter', '\n');
fclose(Fid);

RavenLabels = RavenLabels{1};
RavenLabels_Header = textscan(RavenLabels{1}, '%s', 'DeLimiter', '\t'); % First line is a header
RavenLabels_Header = RavenLabels_Header{1};

% Now locate the columns that I need to write to the .not.mat file -
% namely, "Begin Time (s)", "End Time (s)", "Tag".
% These correspond to onset, offset and label
OnsetColumn = strmatch('Begin Time (s)', RavenLabels_Header, 'exact');
OffsetColumn = strmatch('End Time (s)', RavenLabels_Header, 'exact');
TagColumn = strmatch('Tag', RavenLabels_Header, 'exact');

% Initialise onsets, offsets and labels
Onsets = [];
Offsets = [];
Labels = [];

% Now go through each Raven label line and extract the 
for i = 2:length(RavenLabels),
    TempLabels = textscan(RavenLabels{i}, '%s', 'DeLimiter', '\t');
    TempLabels = TempLabels{1};
    Onsets(end+1) = str2double(TempLabels{OnsetColumn})*1000; % convert to ms
    Offsets(end+1) = str2double(TempLabels{OffsetColumn})*1000; % convert to ms
    Labels{end+1} = TempLabels{TagColumn};
end

% Now, that all labels have been read in from the Raven File - 4 things
% need to be done
% 1) Split up the onsets, offsets and labels into the different file
% specific onsets and offsets. Some syllables are bound to be split between
% two files - these I have to check carefully and split them (see 2 below)
% 2) Some Raven labels are not single labels but include labels for
% multiple syllables - for example the entire motif is labelled with one
% label. This cannot be handled by AutoSongSegmentLabel as labels are
% always individual syllables. So, I have to split the underlying raw data
% into the appropriate number of syllables - this will be most convenient
% for further analysis
% 3) I should ignore all the 'bout#' labels as these are redundant labels
% for the entire bout
% 4) Then I write the data to a .not.mat file

% First remove all the labels with 'bout#' in them
BoutLabels = find(cellfun(@length, strfind(Labels, 'bout')));
Labels(BoutLabels) = [];
Onsets(BoutLabels) = [];
Offsets(BoutLabels) = [];

% First initialise a set of onsets, offsets and labels for each of the 
% First step 1 - split all onsets and offsets into the files. 
for i = 1:length(SongFiles),
    FileSpecific_Onsets{i} = [];
    FileSpecific_Offsets{i} = [];
    FileSpecific_Labels{i} = [];
end

% Will combine this with step 2 which can be in a nested loop Syllables
% with multiple labels that are spread across two files, I will take the
% raw data from the corresponding time points and split them accordingly
% right now in this step

for i = 1:length(Onsets),
    % First check if the onsets and offsets are spread across one file or
    % two
    % First for onsets
    OnsetTimeDiff = FileOnsets - Onsets(i);
    OffsetTimeDiff = FileOffsets - Onsets(i);
    OnsetOffsetTimeDiffs = OnsetTimeDiff.*OffsetTimeDiff;
    % If there is a time difference that is 0, then I will put it in the
    % file with the onset as this is onsets.
    % Otherwise, find the onset and offset pair that is negative and
    % positive, this would mean that the onset is between them.
    OnsetFile = find(OnsetOffsetTimeDiffs == 0);
    if (~isempty(OnsetFile))
        OnsetFile = OnsetFile(end);
    else
        OnsetFile = find(OnsetOffsetTimeDiffs < 0);
    end
    
    % Next for offsets
    OnsetTimeDiff = FileOnsets - Offsets(i);
    OffsetTimeDiff = FileOffsets - Offsets(i);
    OnsetOffsetTimeDiffs = OnsetTimeDiff.*OffsetTimeDiff;
    % If there is a time difference that is 0, then I will put it in the
    % file with the onset as this is onsets.
    % Otherwise, find the onset and offset pair that is negative and
    % positive, this would mean that the onset is between them.
    OffsetFile = find(OnsetOffsetTimeDiffs == 0);
    if (~isempty(OffsetFile))
        OffsetFile = OffsetFile(end);
    else
        OffsetFile = find(OnsetOffsetTimeDiffs < 0);
    end
    
    % Now check if there are multiple labels or just one label - if one
    % label and onset file and offset file are different then captialize
    % the label and split it between the two files
    
    if (length(Labels{i}) == 1)
        if (OnsetFile == OffsetFile)
            FileSpecific_Onsets{OnsetFile}(end+1) = Onsets(i) - FileOnsets(OnsetFile);
            FileSpecific_Offsets{OnsetFile}(end+1) = Offsets(i) - FileOnsets(OnsetFile);
            FileSpecific_Labels{OnsetFile}(end+1) = char(Labels{i});
        else
            FileSpecific_Onsets{OnsetFile}(end+1) = Onsets(i) - FileOnsets(OnsetFile);
            FileSpecific_Offsets{OnsetFile}(end+1) = FileOffsets(OnsetFile);
            FileSpecific_Labels{OnsetFile}(end+1) = upper(Labels(i));
            
            FileSpecific_Onsets{OnsetFile+1}(end+1) = FileOnsets(OnsetFile+1);
            FileSpecific_Offsets{OnsetFile+1}(end+1) = Offsets(i) - FileOnsets(OnsetFile+1);
            FileSpecific_Labels{OnsetFile+1}(end+1) = upper(Labels(i));
        end
    else
        % if there are multiple labels, then display the raw data the first
        % time and get the user to enter upper and lower thresholds to
        % split the syllables into as many labels
        if (OnsetFile == OffsetFile)
            [RawData, Fs] = GetData(RawDataDir, SongFiles{OnsetFile}, FileType, 0);
        else
            RawData = [];
            for j = OnsetFile:OffsetFile,
                [TempRawData, Fs] = GetData(RawDataDir, SongFiles{j}, FileType, 0);
                RawData = [RawData(:); TempRawData(:)];
            end
        end
        
        LabelSpecificData = RawData(round(Fs*(Onsets(i) - FileOnsets(OnsetFile))/1000):round(Fs*(Offsets(i) - FileOnsets(OnsetFile))/1000));
        % Now plot spectrogram and log amplitude below for user to decide
        % threshold
        ThresholdFig = figure;
        set(ThresholdFig, 'Position', [1 1 560 420]);
        p = panel();
        p.pack({1/2 1/2});dir 
        p(1).select();
        PlotSpectrogramInAxis_SongVar(LabelSpecificData, (1:1:length(LabelSpecificData))/Fs, Fs, gca);
        p(2).select();
        FFTWinSize = 5; % in ms
        [LogAmplitude] = ASSLCalculateLogAmplitudeAronovFee(LabelSpecificData, Fs, [], FFTWinSize, []);
        plot((1:1:length(LabelSpecificData))/Fs, LogAmplitude, 'b');
        
        Thresholds = inputdlg('Enter an upper and lower threshold separated by a comma', 'Upper and lower thresholds');
        Thresholds = textscan(Thresholds{1}, '%f', 'DeLimiter', ',');
        Thresholds = Thresholds{1};
        close(ThresholdFig);
        
        [TempSyllOnsets, TempSyllOffsets] = ASSLSegmentDataAronovFee(LogAmplitude, Fs, MinInt, MinDur, Thresholds);
        for j = 1:length(TempSyllOnsets),
            FileSpecific_Onsets{OnsetFile}(end+1) = Onsets(i) - FileOnsets(OnsetFile) + TempSyllOnsets(j);
            FileSpecific_Offsets{OnsetFile}(end+1) = Onsets(i) - FileOnsets(OnsetFile) + TempSyllOffsets(j);
            if (j <= length(Labels{i}))
                FileSpecific_Labels{OnsetFile}(end+1) = Labels{i}(j);
            else
                FileSpecific_Labels{OnsetFile}(end+1) = Labels{i}(end);
            end
        end            
    end
end

% Have to convert all labels to char instead of double
for i = 1:length(FileSpecific_Labels),
    FileSpecific_Labels{i} = char(FileSpecific_Labels{i});
end

% Write note files
for i = 1:length(FileSpecific_Labels),
    onsets = FileSpecific_Onsets{i};
    offsets = FileSpecific_Offsets{i};
    labels = FileSpecific_Labels{i};
    sm_win = 8;
    min_int = 5;
    min_dur = 10;
    threshold = [-40 -50];
    if (~exist(fullfile(RawDataDir, 'ASSLNoteFiles'), 'dir'))
        PresentDir = pwd;
        cd(RawDataDir);
        !mkdir ASSLNoteFiles
        cd(PresentDir);
    end
    save(fullfile(RawDataDir, 'ASSLNoteFiles', [SongFiles{i}, '.not.mat']), 'onsets', 'offsets', 'labels', 'threshold', 'min_int', 'min_dur', 'sm_win');
end
disp(['Finished converting Raven file : ', RavenFile, ' to .not.mat files']);