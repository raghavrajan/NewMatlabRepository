function [RawData, Fs] = ASSLGetRawData(DirectoryName, FileName, FileType, SongChanNo)


SlashIndex = find((FileName == '/') | (FileName == '\'));
if (~isempty(SlashIndex))
    FileName = FileName(SlashIndex(end)+1:end);
end

switch FileType
    case 'okrank'
        [RawData, Fs] = ReadOKrankData(DirectoryName, FileName, SongChanNo);
    
    case 'obs'
        [RawData, Fs] = soundin_copy(DirectoryName, FileName, 'obs');
        RawData = RawData * 5/32768;

    case 'wav'
        PresentDir = pwd;
        cd(DirectoryName);
        [RawData, Fs] = wavread(FileName);
        cd(PresentDir);
end
