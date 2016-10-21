function [AllLabels, AllOnsets, AllOffsets, AllUnAdjustedOnsets, AllUnAdjustedOffsets, OnsetFileNos, OffsetFileNos] = CombineContinuousDataNoteFiles(DirName, Files, NoteFileDir, FileType)

AllLabels = [];
AllOnsets = [];
AllOffsets = [];

AllUnAdjustedOnsets = [];
AllUnAdjustedOffsets = [];

OnsetFileNos = [];
OffsetFileNos = [];

CumulativeFileTime = 0;

for i = 1:length(Files),
    [RawData, Fs] = GetData(DirName, Files{i}, FileType, 0);
    NoteInfo{i} = load(fullfile(NoteFileDir, [Files{i}, '.not.mat']));
    if (~isempty(NoteInfo{i}.labels))
        AllLabels = [AllLabels NoteInfo{i}.labels];
        
        AllOnsets = [AllOnsets(:); (NoteInfo{i}.onsets(:) + CumulativeFileTime)];
        AllUnAdjustedOnsets = [AllUnAdjustedOnsets(:); (NoteInfo{i}.onsets(:))];
        
        AllOffsets = [AllOffsets(:); (NoteInfo{i}.offsets(:) + CumulativeFileTime)];
        AllUnAdjustedOffsets = [AllUnAdjustedOffsets(:); (NoteInfo{i}.offsets(:))];
        
        OnsetFileNos = [OnsetFileNos(:); ones(size(NoteInfo{i}.offsets(:)))*i];
        OffsetFileNos = [OffsetFileNos(:); ones(size(NoteInfo{i}.offsets(:)))*i];
    end
    CumulativeFileTime = CumulativeFileTime + (length(RawData)*1000/Fs);
end

% Now delete all labels with '0' as they're just noise
IndicesToDelete = find(AllLabels == '0');
AllLabels(IndicesToDelete) = [];

AllOnsets(IndicesToDelete) = [];
AllUnAdjustedOnsets(IndicesToDelete) = [];

AllOffsets(IndicesToDelete) = [];
AllUnAdjustedOffsets(IndicesToDelete) = [];

OnsetFileNos(IndicesToDelete) = [];
OffsetFileNos(IndicesToDelete) = [];

% Now fix all consecutive capital letters as they're the same syllable
% split over two files
ConsecutiveCaps = regexp(AllLabels, '[A-Z]');
SyllsToMerge = [];
for i = 1:length(ConsecutiveCaps)-1,
    if (ConsecutiveCaps(i+1) == (ConsecutiveCaps(i) + 1))
        if (AllLabels(ConsecutiveCaps(i)) == AllLabels(ConsecutiveCaps(i+1)))
            SyllsToMerge(end+1) = ConsecutiveCaps(i);
        end
    end
end

AllLabels(SyllsToMerge) = lower(AllLabels(SyllsToMerge));
AllLabels(SyllsToMerge+1) = [];
AllOnsets(SyllsToMerge+1) = [];
AllUnAdjustedOnsets(SyllsToMerge+1) = [];
OnsetFileNos(SyllsToMerge+1) = [];

AllOffsets(SyllsToMerge) = [];
AllUnAdjustedOffsets(SyllsToMerge) = [];
OffsetFileNos(SyllsToMerge) = [];