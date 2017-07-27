function [] = SplitSyllsAcrossFiles(FileList, FileType, NoteFileDir, OriginalSyll)

% This file is used to split syllables that have been merged during
% automatic segmenting
FFTWinSize = 5; % in ms
PrePostPadding = 10; % in ms

Fid = fopen(FileList, 'r');
TempNames = textscan(Fid, '%s', 'DeLimiter', '\n');
FileNames = TempNames{1};
fclose(Fid);

% First get the raw amplitude waveforms of all the syllables that have to
% be split and plot them along with the local minimas

PresentDir = pwd;
Index = 1;
for i = 1:length(FileNames),
    NoteInfo(i) = load(fullfile(PresentDir, NoteFileDir, [FileNames{i}, '.not.mat']));
    OriginalSyllIndices = find(NoteInfo(i).labels == OriginalSyll);
    
    for j = 1:length(OriginalSyllIndices),
        [RawSong, Fs] = GetData(PresentDir, FileNames{i}, FileType, 0);
        [LogAmplitude] = ASSLCalculateLogAmplitudeAronovFee(RawSong, Fs, [], FFTWinSize, []);
        
        if (round(NoteInfo(i).onsets(OriginalSyllIndices(j)) * Fs/1000) > 0)
            if (round(NoteInfo(i).offsets(OriginalSyllIndices(j)) * Fs/1000) <= length(LogAmplitude))
                OriginalSyllAmp{Index} = LogAmplitude(round(NoteInfo(i).onsets(OriginalSyllIndices(j)) * Fs/1000):round(NoteInfo(i).offsets(OriginalSyllIndices(j)) * Fs/1000));
            else
                OriginalSyllAmp{Index} = LogAmplitude(round(NoteInfo(i).onsets(OriginalSyllIndices(j)) * Fs/1000):end);
            end
        else
            OriginalSyllAmp{Index} = LogAmplitude(1:round(NoteInfo(i).offsets(OriginalSyllIndices(j)) * Fs/1000));
        end
        FileNo(Index) = i;
        SyllIndex(Index) = OriginalSyllIndices(j);
        if (round((NoteInfo(i).onsets(OriginalSyllIndices(j)) - PrePostPadding) * Fs/1000) > 0)
            if (round((NoteInfo(i).offsets(OriginalSyllIndices(j)) + PrePostPadding) * Fs/1000) <= length(RawSong))
                SyllData{Index} = RawSong(round((NoteInfo(i).onsets(OriginalSyllIndices(j)) - PrePostPadding) * Fs/1000):round((NoteInfo(i).offsets(OriginalSyllIndices(j)) + PrePostPadding) * Fs/1000));
            else
                SyllData{Index} = RawSong(round((NoteInfo(i).onsets(OriginalSyllIndices(j)) - PrePostPadding) * Fs/1000):end);
            end
        else
            SyllData{Index} = RawSong(1:round((NoteInfo(i).offsets(OriginalSyllIndices(j)) + PrePostPadding) * Fs/1000));
        end
        Index = Index + 1;
    end
end

% Plot all the amplitudes with one representative spectrogram and let the
% user choose where to split and how many parts to split into
figure;
subplot(2,1,1);
PlotSpectrogramInAxis_SongVar(SyllData{1}, (1:1:length(SyllData{1}))/Fs, Fs, gca);

subplot(2,1,2);
hold on;
IndicesToPlot = randperm(length(OriginalSyllAmp));
IndicesToPlot = IndicesToPlot(1:(min(length(IndicesToPlot), 50)));

for i = IndicesToPlot(:)',
    plot((1:1:length(OriginalSyllAmp{i}))/Fs, OriginalSyllAmp{i}, 'b');
end
axis tight;
xlabel('Time (sec)', 'FontSize', 12);
ylabel('Log amplitude (dB)', 'FontSize', 12);
set(gcf, 'Color', 'w');
set(gca, 'FontSize', 14);
set(gcf, 'Position', [680 89 800 900]);
AxisLimits = axis;

% Now ask the user to use left click to put the centers for where syllables
% should be segmented and right click when done.
title('Instructions: Left click to place a location wherever the syllable needs to be split; right click when done; q to quit', 'FontSize', 9);
Flag = 1;
PointsToSplit = [];
while(Flag)
    [x, y, button] = ginput(1);
    if (button == 3)
        Flag = 0;
        break;
    end
    if (button == 113)
        return;
    end
    PointsToSplit(end+1) = x(1)*1000;
    plot([x x], AxisLimits(3:4), 'r--', 'LineWidth', 2);
end

NewSyllLabels = inputdlg('Enter the labels for the new syllables separated by commas', 'New SyllLabels');

if (isempty(NewSyllLabels))
    while (isempty(NewSyllLabels))
        NewSyllLabels = inputdlg('Enter the labels for the new syllables separated by commas', 'New SyllLabels');
    end
end

TempNewSyllLabels = NewSyllLabels{1};
TempNewSyllLabels = TempNewSyllLabels(find(~isspace(TempNewSyllLabels)));
clear NewSyllLabels;
CommaIndices = find(TempNewSyllLabels == ',');


Index = 1;
if (isempty(CommaIndices))
    NewSyllLabels(Index) = TempNewSyllLabels;
else
    NewSyllLabels(Index) = TempNewSyllLabels(1:CommaIndices(1)-1);
    Index = Index + 1;
    for i = 1:length(CommaIndices)-1,
        NewSyllLabels(Index) = TempNewSyllLabels((CommaIndices(i)+1):CommaIndices(i+1) - 1);
        Index = Index + 1;
    end
    NewSyllLabels(Index) = TempNewSyllLabels((CommaIndices(end)+1):end);
end

% Now to go through each syllable, split it by finding the local minimum in
% the log amplitude waveform near the point where it needs to be split and
% then cut it about 2ms before and after this point. 2ms is arbitrary -
% since I will anyway later adjust boundaries, I'm guessing this should not
% matter
% After cutting, the note times should be added into the original note file
% and should be written to file.
% While finding local minima, I will look in a window +/-5ms around the
% point that I have chosen for splitting

PrePostPaddingForBoundary = 5; % in ms
PrePostPaddingForSplit = 2; % in ms
for j = 1:length(PointsToSplit),
    SplitBoundaries(j,:) = PointsToSplit(j) + [-PrePostPaddingForBoundary PrePostPaddingForBoundary];
end

ShortSyllFlag = 0;
UniqueFiles = unique(FileNo);
for Files = 1:length(UniqueFiles),
    SyllsToSplit = find(FileNo == UniqueFiles(Files));
    TempSyllOnsets = [];
    TempSyllLabels = [];
    for i = SyllsToSplit(:)',
        NewSyllOnsets = [];
        for j = 1:length(PointsToSplit),
            if (j == 1)
                NewSyllOnsets(end+1,1) = 0;  % onset of new syllable
            end
            if (length(OriginalSyllAmp{i}) < (round(SplitBoundaries(j,2) * Fs/1000)))
                ShortSyllFlag = 1;
                break;
            end
            [MinVal, MinValLocation] = min(OriginalSyllAmp{i}(round(SplitBoundaries(j,1) * Fs/1000):round(SplitBoundaries(j,2) * Fs/1000)));

            MinValLocation = MinValLocation + (round(SplitBoundaries(j,1) * Fs/1000) - 1);
            % Now new boundary is 2ms before and after this point
            NewSyllOnsets(end,2) = (MinValLocation * 1000/Fs) - PrePostPaddingForSplit; % new offset
            NewSyllOnsets(end+1,1) = (MinValLocation * 1000/Fs) + PrePostPaddingForSplit; %
            TempSyllLabels(end+1) = NewSyllLabels(j);
        end
        TempSyllLabels(end+1) = NewSyllLabels(end);
        
        if (ShortSyllFlag == 1)
            ShortSyllFlag = 0;
            disp('Skipped one syllable');
            continue;
        end

        NewSyllOnsets(end,2) = length(OriginalSyllAmp{i})*1000/Fs;

        NewSyllOnsets = NewSyllOnsets + NoteInfo(FileNo(i)).onsets(SyllIndex(i));
        TempSyllOnsets = [TempSyllOnsets; NewSyllOnsets];
    end
    % Delete original syllable
    NoteInfo(UniqueFiles(Files)).onsets(SyllIndex(SyllsToSplit)) = [];
    NoteInfo(UniqueFiles(Files)).offsets(SyllIndex(SyllsToSplit)) = [];
    NoteInfo(UniqueFiles(Files)).labels(SyllIndex(SyllsToSplit)) = [];
    
    % Now add new syllables
    NoteInfo(UniqueFiles(Files)).onsets = [NoteInfo(UniqueFiles(Files)).onsets(:); TempSyllOnsets(:,1)];
    NoteInfo(UniqueFiles(Files)).offsets = [NoteInfo(UniqueFiles(Files)).offsets(:); TempSyllOnsets(:,2)];
    NoteInfo(UniqueFiles(Files)).labels = [NoteInfo(UniqueFiles(Files)).labels(:)' char(TempSyllLabels(:)')];
    
    % Now sort onsets, offsets and labels
    [SortedVals, SortedIndices] = sort(NoteInfo(UniqueFiles(Files)).onsets);
    NoteInfo(UniqueFiles(Files)).onsets = NoteInfo(UniqueFiles(Files)).onsets(SortedIndices);
    NoteInfo(UniqueFiles(Files)).offsets = NoteInfo(UniqueFiles(Files)).offsets(SortedIndices);
    NoteInfo(UniqueFiles(Files)).labels = NoteInfo(UniqueFiles(Files)).labels(SortedIndices);
    
    % Now write data to note file. Save fields as individual variables
    % to be consistent with original note files
    TempNoteInfo = NoteInfo(UniqueFiles(Files));
    save(fullfile(PresentDir, NoteFileDir, [FileNames{UniqueFiles(Files)}, '.not.mat']), '-struct', 'TempNoteInfo');
end

OutputString = ['Finished splitting ', num2str(length(OriginalSyllAmp)), ' ', OriginalSyll, ' into '];
for i = 1:length(NewSyllLabels),
    OutputString = [OutputString, NewSyllLabels(i)];
    if (i < (length(NewSyllLabels) - 1))
        OutputString = [OutputString, ', '];
    else
        if (i == (length(NewSyllLabels) - 1))
            OutputString = [OutputString, ' and '];
        end
    end
end
disp(OutputString);