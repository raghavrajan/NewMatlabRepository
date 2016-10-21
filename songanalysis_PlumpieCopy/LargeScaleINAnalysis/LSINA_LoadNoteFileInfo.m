function [NoteInfo, FileLen] = LSINA_LoadNoteFileInfo(BirdParameters)

% Now for each of the filelists, load up the labels and onsets and offsets
for j = 1:length(BirdParameters.SongFileNames),
    NoteInfo{j} = load(fullfile(BirdParameters.DataDirectory, 'ASSLNoteFiles', [BirdParameters.SongFileNames{j}, '.not.mat']));
    % Remove all the '0' labels
    Indices = find(NoteInfo{j}.labels == '0');
    if (~isempty(Indices))
        NoteInfo{j}.labels(Indices) = [];
        NoteInfo{j}.onsets(Indices) = [];
        NoteInfo{j}.offsets(Indices) = [];
    end
    [RawData, Fs] = GetData(fullfile(BirdParameters.DataDirectory), BirdParameters.SongFileNames{j}, BirdParameters.FileType, 0);
    FileLen(j) = length(RawData)*1000/Fs;
end