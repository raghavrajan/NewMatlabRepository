function [] = MakeContextInfoFile(FileList, DirectedFiles, FemaleInOutFiles, MotifSyllables)

Index = strfind(FileList, '.txt');
ContextFileName = [FileList(1:Index), 'Context.txt'];
ContentFileName = [FileList(1:Index), 'Content.txt'];
DirFileName = [FileList(1:Index), 'Dir.txt'];
UnDirFileName = [FileList(1:Index), 'UnDir.txt'];

Fid = fopen(FileList, 'r');
Temp = textscan(Fid, '%s', 'DeLimiter', '\n');
Temp = Temp{1};
fclose(Fid);

Fid = fopen(DirFileName, 'w');
for i = 1:length(DirectedFiles),
    fprintf(Fid, '%s\n', DirectedFiles{i});
    DirFileIndices(i) = find(cellfun(@length, strfind(Temp, DirectedFiles{i})));
end
fclose(Fid);

UnDirectedFiles = Temp;
UnDirectedFiles(DirFileIndices) = [];
Fid = fopen(UnDirFileName, 'w');
for i = 1:length(UnDirectedFiles),
    fprintf(Fid, '%s\n', UnDirectedFiles{i});
    UnDirFileIndices(i) = find(cellfun(@length, strfind(Temp, UnDirectedFiles{i})));
end
fclose(Fid);

Fid = fopen(ContextFileName, 'w');
for i = 1:length(Temp),
    fprintf(Fid, '%s\t', Temp{i});
    if (isempty(find(DirFileIndices == i)))
        fprintf(Fid, 'U\n');
    else
        fprintf(Fid, 'D\n');
    end
end
fclose(Fid);

if (exist('AdjustedNoteFiles', 'dir'))
    Fid = fopen(ContentFileName, 'w');
    for i = 1:length(Temp),
        fprintf(Fid, '%s\t', Temp{i});
        if (ispc)
            Notes = load(['AdjustedNoteFiles\', Temp{i}, '.not.mat']);
        else
            Notes = load(['AdjustedNoteFiles/', Temp{i}, '.not.mat']);
        end
        Notes.labels(find(Notes.labels == '0')) = [];
        
        if (isempty(Notes.labels))
            fprintf(Fid, 'B\n');
        else
            MotifSyllableFlag = 0;
            for j = 1:length(MotifSyllables)
                if (~isempty(find(Notes.labels == MotifSyllables(j))))
                    MotifSyllableFlag = 1;
                    break;
                end
            end
            if (MotifSyllableFlag == 1)
                fprintf(Fid, 'S\n');
            else
                fprintf(Fid, 'C\n');
            end
        end
    end
    fclose(Fid);
end

disp('Finished writing text files for context and content information');
