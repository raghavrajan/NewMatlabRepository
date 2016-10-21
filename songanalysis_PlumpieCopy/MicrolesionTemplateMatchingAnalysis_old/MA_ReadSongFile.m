function [Song, Fs] = MA_ReadSongFile(DataDir, SongFileName, FileType)

%==========================================================================
% Function for reading song data from a file
% Raghav Rajan - 29th November 2013
%==========================================================================
FileSep = filesep;

if (strfind(FileType, 'okrank'))
    [Song, Fs] = ReadOKrankData(DataDir, SongFileName, 1);
    Song = Song/10;
else
    if (strfind(FileType, 'wav'))
        PresentDir = pwd;
        cd(DataDir);
        [Song, Fs] = wavread(SongFileName);
        cd(PresentDir);
    else
        if (strfind(FileType, 'obs'))
            channel_string = strcat('obs',num2str(0),'r');
            [Song, Fs] = soundin_copy([DataDir, FileSep], SongFileName, channel_string);
            % Convert to uV - 5V on the data acquisition is 32768
            Song = Song * 5/32768;
            % Convert to +1 to -1V - same as wav files
            Song = Song/5;
        end
    end
end