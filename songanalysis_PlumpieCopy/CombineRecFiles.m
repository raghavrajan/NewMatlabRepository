function [RecTime, TotalSamples] = CombineRecFiles(DirectoryName,FileName,FileNos,FileType)

TotalSamples = 0;
RecTime = 0;

cd(DirectoryName);

for i = 1:length(FileNos),  
    if (strfind(FileType, 'obs'))
        if (FileNos(i) < 10)
            DataFile = [FileName,'.00',num2str(FileNos(i))];
        else
            if (FileNos(i) < 100)
                DataFile = [FileName,'.0',num2str(FileNos(i))];            
            else
                DataFile = [FileName,'.',num2str(FileNos(i))];
            end
        end
    else
        if (strfind(FileType, 'okrank'))
            if (FileNos < 100000)
                DataFile = [FileName, '0', num2str(FileNos(i))];
            else
                DataFile = [FileName, num2str(FileNos(i))];
            end
        end
    end

    RecFileName = strcat(DataFile,'.rec');

    if (~exist(RecFileName,'file'))
        continue;
    end
    
    fid = fopen(RecFileName,'r');
    index = 0;
    while (index < 5)
        if (feof(fid))
            break;
        end
        tline = fgetl(fid);
        if (strfind(tline,'rec end') > 0)
            StartIndex = find(tline == '=');
            EndIndex = strfind(tline,'ms');
            RecTime = RecTime + str2double(tline((StartIndex + 1):(EndIndex - 1)));
        end
        
        if (strfind(FileType, 'obs'))
            if (strfind(tline, 'Samples') > 0)
                StartIndex = find(tline == '=');
                TotalSamples = TotalSamples + str2double(tline((StartIndex + 1):end));
            end
        else
            if (strfind(FileType, 'okrank'))
                if (strfind(tline, 'n_samples') > 0)
                    StartIndex = find(tline == ':');
                    TotalSamples = TotalSamples + str2double(tline((StartIndex + 1):end));
                end
            end
        end
    end

    fclose(fid);
end