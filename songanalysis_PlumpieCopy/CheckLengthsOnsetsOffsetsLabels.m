function [Flag] = CheckLengthsOnsetsOffsetsLabels(Directory, FileList)

% Script to check if the length of labels, onsets and offsets are the same
% in a given note file

Flag = 0;

PresentDir = pwd;

cd(Directory);
Fid = fopen(FileList, 'r');
Files = textscan(Fid, '%s', 'DeLimiter', '\n');
Files = Files{1};
fclose(Fid);

cd(fullfile(Directory, 'ASSLNoteFiles'));
for i = 1:length(Files),
    Temp = load([Files{i}, '.not.mat']);
    if ((length(Temp.labels) ~= length(Temp.onsets)) || (length(Temp.labels) ~= length(Temp.offsets)) || (length(Temp.offsets) ~= length(Temp.onsets)))
        disp(Files{i});
        Flag = 1;
    end
end

cd(PresentDir);
