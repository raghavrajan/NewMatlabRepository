function [] = FindValidSongFiles(BatchFileName,Motif)

fid = fopen(BatchFileName,'r');
DirectoryName = pwd;
DirectoryName(end + 1) = '/';

Silence = wavread('Silence.wav');

ValidSongFiles = [];
ValidSongFileIndex = 1;

while (~(feof(fid)))
    tline = fgetl(fid);
    NoteFile = [tline,'.not.mat'];
    FiltFile = [tline,'.filt'];
    
    if ~(exist(NoteFile,'file'))
        continue;
    end
    
    Notes = load(NoteFile);
    
    for i = 1:length(Motif),
        MotifNotes = find(Notes.labels == Motif(i));
        if (length(MotifNotes) > 0)
            ValidSongFiles{ValidSongFileIndex} = tline;
            ValidSongFileIndex = ValidSongFileIndex + 1;
            break;
        end
    end
end

fclose(fid);

OutputFileName = [BatchFileName, '.ValidSongFiles.txt'];
fid = fopen(OutputFileName,'w');
for i = 1:length(ValidSongFiles),
    for j = 1:length(ValidSongFiles{i});
        fprintf(fid,'%c',ValidSongFiles{i}(j));
    end
    fprintf(fid,'\n');
end
fclose(fid);
        