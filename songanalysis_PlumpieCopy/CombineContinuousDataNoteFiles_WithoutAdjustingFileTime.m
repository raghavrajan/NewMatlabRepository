function [AllLabels, AllOnsets, AllOffsets, OnsetFileNos, OffsetFileNos, OnsetFileLabelNos, OffsetFileLabelNos] = CombineContinuousDataNoteFiles_WithoutAdjustingFileTime(DirName, Files, NoteFileDir, FileType)

AllLabels = [];
AllOnsets = [];
AllOffsets = [];
OnsetFileNos = [];
OffsetFileNos = [];
OnsetFileLabelNos = [];
OffsetFileLabelNos = [];

for i = 1:length(Files),
    [RawData, Fs] = GetData(DirName, Files{i}, FileType, 0);
    NoteInfo{i} = load(fullfile(NoteFileDir, [Files{i}, '.not.mat']));
    if (~isempty(NoteInfo{i}.labels))
        AllLabels = [AllLabels NoteInfo{i}.labels];
        AllOnsets = [AllOnsets(:); (NoteInfo{i}.onsets(:))];
        AllOffsets = [AllOffsets(:); (NoteInfo{i}.offsets(:))];
        OnsetFileNos = [OnsetFileNos(:); ones(size(NoteInfo{i}.offsets(:)))*i];
        OffsetFileNos = [OffsetFileNos(:); ones(size(NoteInfo{i}.offsets(:)))*i];
        OnsetFileLabelNos = [OnsetFileLabelNos(:); (1:1:length(NoteInfo{i}.labels))'];
        OffsetFileLabelNos = [OffsetFileLabelNos(:); (1:1:length(NoteInfo{i}.labels))'];
    end
end

% Now delete all labels with '0' as they're just noise
IndicesToDelete = find(AllLabels == '0');
AllLabels(IndicesToDelete) = [];
AllOnsets(IndicesToDelete) = [];
AllOffsets(IndicesToDelete) = [];
OnsetFileNos(IndicesToDelete) = [];
OffsetFileNos(IndicesToDelete) = [];
OnsetFileLabelNos(IndicesToDelete) = [];
OffsetFileLabelNos(IndicesToDelete) = [];

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
OnsetFileNos(SyllsToMerge+1) = [];
OnsetFileLabelNos(SyllsToMerge + 1) = [];

AllOffsets(SyllsToMerge) = [];
OffsetFileNos(SyllsToMerge) = [];
OffsetFileLabelNos(SyllsToMerge) = [];