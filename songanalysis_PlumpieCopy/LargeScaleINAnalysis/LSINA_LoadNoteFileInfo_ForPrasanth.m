function [BoutNoteInfo] = LSINA_LoadNoteFileInfo_ForPrasanth(BirdParameters)

% Now for each of the filelists, load up the labels and onsets and offsets
% from two different sets of ASSLNoteFiles directories

% The regular set is already loaded up, so now I will only load up the
% extra note files

% Now the set corresponding to bout labels
if (exist(fullfile(BirdParameters.DataDirectory, 'bout_label_ASSLNoteFiles'), 'dir'))
    for j = 1:length(BirdParameters.SongFileNames),
        BoutNoteInfo{j} = load(fullfile(BirdParameters.DataDirectory, 'bout_label_ASSLNoteFiles', [BirdParameters.SongFileNames{j}, '.not.mat']));
        % Remove all the '0' labels
        Indices = find(BoutNoteInfo{j}.labels == '0');
        if (~isempty(Indices))
            BoutNoteInfo{j}.labels(Indices) = [];
            BoutNoteInfo{j}.onsets(Indices) = [];
            BoutNoteInfo{j}.offsets(Indices) = [];
        end
    end
else
    BoutNoteInfo = [];
end