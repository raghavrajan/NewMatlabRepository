function [NoteOnsets, NoteOffsets, NoteLabels] = GetNoteInformation(DirectoryName, FileNames)

cd(DirectoryName);
% Now load up all the information about syllable onsets and offsets from
% the song bout files

NoteOnsets = cell(length(FileNames), 1);
NoteOffsets = cell(length(FileNames), 1);
NoteLabels = cell(length(FileNames), 1);

for i = 1:length(FileNames),
    Notes = load([FileNames{i},'.not.mat']);
    
%   Note onsets and offsets from uisonganal are in ms, so that has to be
%   converted to seconds
    NoteOnsets{i} = Notes.onsets/1000;
    NoteOffsets{i} = Notes.offsets/1000;
    NoteLabels{i} = Notes.labels;
    clear Notes;
end