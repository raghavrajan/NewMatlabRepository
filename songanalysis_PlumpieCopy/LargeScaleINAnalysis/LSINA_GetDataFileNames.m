function [SongFileNames] =LSINA_GetDataFileNames(BirdParameters)

PresentDir = pwd;

cd(BirdParameters.DataDirectory);

if (isfield(BirdParameters, 'UndirSongFileList'))
    Fid = fopen(BirdParameters.UndirSongFileList, 'r');
else
    Fid = fopen(BirdParameters.SongFileList, 'r');
end
TempNames = textscan(Fid, '%s', 'DeLimiter', '\n');
SongFileNames = TempNames{1};

fclose(Fid);

cd(PresentDir);
