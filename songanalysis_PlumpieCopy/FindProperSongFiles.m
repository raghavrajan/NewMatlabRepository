function [SongFilesList] = FindProperSongFiles(NotesBatchFile,motif)

ActualNotesFile = ['Actual',NotesBatchFile];

fid = fopen(NotesBatchFile);
NotesFilesList = textscan(fid,'%s');
fclose(fid);

SongFilesList = [];

for i = 1:size(NotesFilesList{1},1),
    Notes = load(cell2mat(NotesFilesList{1}(i)));
    flag = 0;
    for j = 1:length(motif),
        MotifLabels = find(Notes.labels == motif(j));
        if (length(MotifLabels) > 0)
            flag = 1;
            break;
        end
    end
    if (flag == 1)
        SongFilesList{end + 1} = NotesFilesList{1}(i);
    end
end

if (size(SongFilesList,1) > 0)
    fid = fopen(ActualNotesFile,'w');
    for i = 1:length(SongFilesList),
        if (i == length(SongFilesList))
            fprintf(fid,'%s',cell2mat(SongFilesList{i}));
        else
            fprintf(fid,'%s\n',cell2mat(SongFilesList{i}));
        end
    end
    fclose(fid);
else
    disp('None of the files in the batch file had the motif');
end

