function [SongFileNames] = MA_ExtractSongFileNames(SongFileList)

%==========================================================================
% Function for reading a text file and getting all the lines into a cell
% array
% Raghav Rajan - 29th November 2013
%==========================================================================

Fid = fopen(SongFileList, 'r');
Temp = textscan(Fid, '%s', 'DeLimiter', '\n');
SongFileNames = Temp{1};
fclose(Fid);

for i = 1:length(SongFileNames),
    SlashIndex = find((SongFileNames{i} == '/') | (SongFileNames{i} == '\'));
    if (~isempty(SlashIndex))
        SongFileNames{i} = SongFileNames{i}(SlashIndex(end)+1:end);
    end
end