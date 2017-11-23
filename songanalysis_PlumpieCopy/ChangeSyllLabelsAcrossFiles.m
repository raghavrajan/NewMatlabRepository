function [] = ChangeSyllLabelsAcrossFiles(FileList, NoteFileDir, OriginalSyllLabel, NewSyllLabel)

PresentDir = pwd;
% This file is used to merge consecutive syllables that have been separated
% during labelling. 
Fid = fopen(FileList, 'r');
TempNames = textscan(Fid, '%s', 'DeLimiter', '\n');
FileNames = TempNames{1};
fclose(Fid);

for i = 1:length(FileNames),
    NoteInfo = load(fullfile(PresentDir, NoteFileDir, [FileNames{i}, '.not.mat']));
    % Change labels of the first syllable to the new label and delete
    % next syllable
    OldSylls = find(NoteInfo.labels == OriginalSyllLabel);
    if (~isempty(OldSylls))
        NoteInfo.labels(OldSylls) = NewSyllLabel;

        % Now write data to note file. Save fields as individual variables
        % to be consistent with original note files
        save(fullfile(PresentDir, NoteFileDir, [FileNames{i}, '.not.mat']), '-struct', 'NoteInfo');
    end
end
disp(['Finished changing all ', OriginalSyllLabel, ' into ', NewSyllLabel]);