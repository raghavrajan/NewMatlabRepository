function [RawData, Fs] = ASSLGetRawData(DirectoryName, FileName, FileType, SongChanNo)

% 
% SlashIndex = find((FileName == '/') | (FileName == '\'));
% if (~isempty(SlashIndex))
%     FileName = FileName(SlashIndex(end)+1:end);
% end
% 
% FileSep = filesep;
% if (DirectoryName(end) ~= FileSep)
%     DirectoryName(end+1) = FileSep;
% end

[RawData, Fs] = GetData(DirectoryName, FileName, FileType, SongChanNo);
