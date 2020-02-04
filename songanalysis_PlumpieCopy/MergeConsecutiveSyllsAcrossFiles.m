10function [] = MergeConsecutiveSyllsAcrossFiles(FileList, NoteFileDir, FirstSyll, SecondSyll, MergedSyllLabel)

% This file is used to merge consecutive syllables that have been separated
% during labelling. 
Fid = fopen(FileList, 'r');
TempNames = textscan(Fid, '%s', 'DeLimiter', '\n');
FileNames = TempNames{1};
fclose(Fid);

% First get the time interval between consecutive syllables to merge
PresentDir = pwd;
SyllableDetails = [];
Intervals = [];
for i = 1:length(FileNames),
    NoteInfo = load(fullfile(PresentDir, NoteFileDir, [FileNames{i}, '.not.mat']));
    FirstSyllIndices = find(NoteInfo.labels == FirstSyll);
    
    for j = 1:length(FirstSyllIndices),
        if (FirstSyllIndices(j) < length(NoteInfo.labels))
            if (NoteInfo.labels(FirstSyllIndices(j)+1) == SecondSyll)
                SyllableDetails(end+1,:) = [i FirstSyllIndices(j)];
                Intervals(end+1) = NoteInfo.onsets(FirstSyllIndices(j)+1) - NoteInfo.offsets(FirstSyllIndices(j));
            end
        end
    end
end

% Plot all the intervals and let the user choose if consecutive syllables
% below a certain interval are the only ones that need to be merged
figure;
hold on;
Edges = 0:10:200;
plot(Edges, histc(Intervals, Edges)/length(Intervals), 'ko-');
title(['Intervals between all consecutive ', FirstSyll, ' and ', SecondSyll], 'FontSize', 16);
xlabel('Time (ms)', 'FontSize', 14);
ylabel('Fraction', 'FontSize', 14);
set(gcf, 'Color', 'w');
set(gca, 'FontSize', 14);
set(gcf, 'Position', [839 32 500 400]);

% Now ask the user what is the time-interval below which the two
% consecutive syllables should be merged
MaxInterval = inputdlg('Enter the maximum interval for two syllables to be merged', 'Maximum Interval for Merging');
if (isempty(MaxInterval))
    while (isempty(MaxInterval))
        MaxInterval = inputdlg('Enter the maximum interval for two syllables to be merged', 'Maximum Interval for Merging');
    end
end
MaxInterval = str2double(MaxInterval{1});

% Now find all intervals below MaxInterval and merge them. The vector
% syllable details is of same length as Intervals and has details of the
% file # and the syllable # for the corresponding intervals. This can be
% used for merging in each file once the valid intervals have been
% identified

LongIntervals = find(Intervals > MaxInterval);
disp(['Merging ', num2str(length(Intervals) - length(LongIntervals)), ' out of ', num2str(length(Intervals)), ' consecutive syllables with interval <= ', num2str(MaxInterval)]);

Intervals(LongIntervals) = [];
SyllableDetails(LongIntervals,:) = [];
for i = 1:length(FileNames),
    FileIndices = find(SyllableDetails(:,1) == i);
    if (~isempty(FileIndices))
        FirstSyllIndices = SyllableDetails(FileIndices,2);
        NoteInfo = load(fullfile(PresentDir, NoteFileDir, [FileNames{i}, '.not.mat']));
        
        % Change labels of the first syllable to the new label and delete
        % next syllable
        NoteInfo.labels(FirstSyllIndices) = MergedSyllLabel;
        NoteInfo.labels(FirstSyllIndices + 1) = [];
        
        % Now keep first set of onsets and consecutive syllable offsets
        NoteInfo.offsets(FirstSyllIndices) = [];
        NoteInfo.onsets(FirstSyllIndices + 1) = [];
        
        % Now write data to note file. Save fields as individual variables
        % to be consistent with original note files
        save(fullfile(PresentDir, NoteFileDir, [FileNames{i}, '.not.mat']), '-struct', 'NoteInfo');
    end
end
disp(['Finished merging all ', FirstSyll, ' and ', SecondSyll]);