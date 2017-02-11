function [RawData, Fs] = GetData(DirectoryName, FileName, FileType, SongChanNo)

%========== Help text for GetData.m =======================================
% Description: This is a function used to get data from one channel of any
% type of data file. Currently it supports 'wav', 'observer' or 'okrank'
% format. I will include support for Intan format too soon
% Usage: [RawData, Fs] = GetData(DirectoryName, FileName, FileType, SongChanNo)
% Inputs: 
%   1) DirectoryName - directory where the data file is present
%   2) FileName - name of raw data file from which data has to be read
%   3) FileType - currently one of 3 options
%        a) 'wav' - for wave file
%        b) 'okrank' - for OKrank file
%        c) 'obs' - for Observer file
%   4) SongChanNo - the channel # within the file for which the data has to
%   be read
% Outputs:
%   1) RawData - the data read in a column vector
%   2) Fs - sampling rate of data
% Eg of usage:
%   [RawData, Fs] = GetData('/home/raghav/y50g89', 'y50g89_082410101110',
%   'okrank', 1);
%
% =========================================================================

SlashIndex = find((FileName == '/') | (FileName == '\'));
if (~isempty(SlashIndex))
    FileName = FileName(SlashIndex(end)+1:end);
end

FileSep = filesep;
if (DirectoryName(end) ~= FileSep)
    DirectoryName(end+1) = FileSep;
end

switch FileType
    case 'okrank'
        [RawData, Fs] = ReadOKrankData(DirectoryName, FileName, SongChanNo);
    
    case 'obs'
        [RawData, Fs] = soundin_copy(DirectoryName, FileName, ['obs', num2str(SongChanNo), 'r']);
        RawData = RawData * 5/32768;

    case 'wav'
        PresentDir = pwd;
        cd(DirectoryName);
        [RawData, Fs] = audioread(FileName);
        cd(PresentDir);
end
