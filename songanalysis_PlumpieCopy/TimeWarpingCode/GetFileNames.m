function [DataFileNames,RecTime] = GetFileNames(DirectoryName,FileName,FileExt,FileNos,FileType)

DataFileNames = [];
RecTime = [];

cd(DirectoryName);

for i = 1:length(FileNos),
    if (strfind(FileType,'obs'))
        if (FileNos(i) < 10)
            DataFile = [FileName,'.00',num2str(FileNos(i))];
        else
            if (FileNos(i) < 100)
                DataFile = [FileName,'.0',num2str(FileNos(i))];            
            else
                DataFile = [FileName,'.',num2str(FileNos(i))];
            end
        end
        if (FileExt(1) == '.')
            DataFileNames{i} = [DataFile,FileExt];
        else
            DataFileNames{i} = [DataFile,'.',FileExt];
        end
    else
        if (strfind(FileType, 'okrank'))
            if (FileNos(i) < 100000)
                DataFile = [FileName,'0', num2str(FileNos(i))];
            else
                DataFile = [FileName, num2str(FileNos(i))];
            end
            
            DataFileNames{i} = DataFile;
        end
    end
    
    RecFileName = strcat(DataFile,'.rec');
    disp(RecFileName);
    
    fid = fopen(RecFileName,'r');
    index = 0;
    while ~(feof(fid))
        tline = fgetl(fid);
        if (strfind(tline,'rec end') > 0)
            break;
        end
        if (strfind(tline, 'ai_freq'))
            StartIndex = find(tline ==':');
            Fs = str2double(tline((StartIndex + 1):end));
        end
        
        if (strfind(tline, 'n_samples'))
            break;
        end
    end

    if (strfind(FileType, 'obs'))
        StartIndex = find(tline == '=');
        EndIndex = strfind(tline,'ms');
        RecTime(i) = str2double(tline((StartIndex + 1):(EndIndex - 1)));
        RecTime(i) = RecTime(i)/1000;
    else
        if (strfind(FileType, 'okrank'))
            StartIndex = find(tline == ':');
            RecTime(i) = str2double(tline((StartIndex + 1):end));
            RecTime(i) = RecTime(i)/Fs;
        end
    end
    fclose(fid);
end