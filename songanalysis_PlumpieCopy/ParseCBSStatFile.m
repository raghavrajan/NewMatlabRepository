function [] = ParseCBSStatFile(FileName, BirdName)

Fid = fopen(FileName, 'r');
tline = fgetl(Fid);
tline = fgetl(Fid);
tline = fgetl(Fid);
tline = fgetl(Fid);

SongFileIndex = 1;
while (strfind(tline, BirdName))
    if (strfind(tline, 'SAVE'))
        SpaceIndex = find(tline == ' ');
        SongFiles{SongFileIndex} = tline(1:(SpaceIndex-1));
        SongFileIndex = SongFileIndex + 1;
    end
    tline = fgetl(Fid);
end
fclose(Fid);

Fid = fopen([FileName, '.txt'], 'w');
for i = 1:length(SongFiles),
    fprintf(Fid, '%s\n', SongFiles{i});
end
fclose(Fid);