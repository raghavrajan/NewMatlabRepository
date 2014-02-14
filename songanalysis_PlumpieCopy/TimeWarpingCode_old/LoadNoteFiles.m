function [NoteOnsets, NoteOffsets, NoteLabels, NoofMotifs] = LoadNoteFiles(DirectoryName, FileNames, RecordLengths, Motif)

cd(DirectoryName);
% Now load up all the information about syllable onsets and offsets from
% the song bout files

RecTime = 0;
NoteOnsets = [];
NoteOffsets = [];
NoteLabels = [];
NoofMotifs = [];

for i = 1:length(FileNames),
    Notes = load([FileNames{i},'.not.mat']);
    
%   Note onsets and offsets from uisonganal are in ms, so that has to be
%   converted to seconds
    Notes.onsets = Notes.onsets/1000 + RecTime;
    Notes.offsets = Notes.offsets/1000 + RecTime;
    RecTime = RecTime + RecordLengths(i);
    
    NoteOnsets = [NoteOnsets; Notes.onsets]; 
    NoteOffsets = [NoteOffsets; Notes.offsets];
    NoteLabels = [NoteLabels [Notes.labels]];
    
    NoofMotifs(i) = length(strfind(Notes.labels,Motif));
    clear Notes;
    
end