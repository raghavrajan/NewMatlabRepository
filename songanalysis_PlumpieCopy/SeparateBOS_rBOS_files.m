function [SongFiles] = SeparateBOS_rBOS_files(FileName,FileExt,FileNos,FileTypes)

DirectoryName = pwd;
DirectoryName(end+1) = '/';

for SongTypes = 1:size(FileTypes,1),
    SongFiles(SongTypes).Files = [];
    SongFiles(SongTypes).Name = FileTypes{SongTypes};
end

for j = 1:length(FileNos),
    i = FileNos(j);
    if (i < 10);
        DataFileName = strcat(FileName,'.00',num2str(i),FileExt);
    else
        if (i < 100);
            DataFileName = strcat(FileName,'.0',num2str(i),FileExt);
        else
            DataFileName = strcat(FileName,'.',num2str(i),FileExt);
        end
    end
    
    dot_index = find(DataFileName == '.');
    dot_index = dot_index(end);
    RecFileName = strcat(DataFileName(1:dot_index),'rec');

    fid = fopen(RecFileName,'r');
    index = 0;
    while (index < 5)
        if (feof(fid))
            break;
        end
        tline = fgetl(fid);
        if (strfind(tline,'Stimulus') > 0)
            index = 5;
            for SongTypes = 1:size(FileTypes,1),
                if (strfind(tline, FileTypes{SongTypes}))
                    SongFiles(SongTypes).Files(end + 1) = FileNos(j);
                    break;
                end
            end
        end
    end

    if (index ~= 5)
        fclose(fid);
        disp('There appears to be no song stimulus in this file');
        return;
    end
    fclose(fid);
end
