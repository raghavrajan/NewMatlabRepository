function [] = MakeDirFileList(FileList, UnDirFileList)

Files = fopen(FileList, 'r');
tline = fgetl(Files);
SongFileIndex = 1;
while (ischar(tline(1)))
    SongFiles{SongFileIndex} = tline;
    SongFileIndex = SongFileIndex + 1;
    tline = fgetl(Files);
end
fclose(Files);

Files = fopen(UnDirFileList, 'r');
tline = fgetl(Files);
SongFileIndex = 1;
while (ischar(tline(1)))
    UnDirSongFiles{SongFileIndex} = tline;
    SongFileIndex = SongFileIndex + 1;
    tline = fgetl(Files);
end
fclose(Files);

DirSongFiles = setdiff(SongFiles, UnDirSongFiles);

Fid = fopen([FileList, '.dirsongs.txt'], 'w');
for i = 1:length(DirSongFiles),
    fprintf(Fid, '%s\n', DirSongFiles{i});
end
fclose(Fid);