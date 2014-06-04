function [BoutOnsets, BoutOffsets, BoutLens] = MA_LoadBoutNoteFiles(SongFileName, NoteFileDir)

PresentDir = pwd;
if (exist(NoteFileDir, 'dir'))
    cd(NoteFileDir);
else
    BoutOnsets = [];
    BoutOffsets = [];
    BoutLens = [];
    return;
end

if (exist([SongFileName, '.not.mat']))
    Temp = load([SongFileName, '.not.mat']);
    BoutIndices = find(Temp.labels == 'B');
    BoutOnsets = Temp.onsets(BoutIndices);
    BoutOffsets = Temp.offsets(BoutIndices);
    BoutLens = BoutOffsets - BoutOnsets;
else
    BoutOnsets = [];
    BoutOffsets = [];
    BoutLens = [];
end

cd(PresentDir);
