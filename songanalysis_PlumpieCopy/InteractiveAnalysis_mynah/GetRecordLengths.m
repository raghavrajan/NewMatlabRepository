function [RecordLengths, NChans] = GetRecordLengths(DirectoryName, FileNames, FileType)

cd(DirectoryName);
RecordLengths = cell(length(FileNames), 1);
NChans = cell(length(FileNames), 1);

for i = 1:length(FileNames),
    if (strfind(FileType, 'observer'))
        DotIndex = find(FileNames{i} == '.');
        RecFileName = [FileNames{i}(1:(DotIndex(end) - 1)), '.rec'];
        RecFid = fopen(RecFileName, 'r');
        index = 0;
        while (index < 5)
            if (feof(RecFid))
                break;
            end
            tline = fgetl(RecFid);
            if ((strfind(tline, 'rec end')) > 0)
                StartIndex = find(tline == '=');
                EndIndex = strfind(tline,'ms');
                RecordLengths{i} = str2double(tline((StartIndex + 1):(EndIndex - 1)));
                RecordLengths{i} = RecordLengths{i}/1000;
            end
            if ((strfind(tline, 'Chans')) > 0)
                StartIndex = find(tline == '=');
                NChans{i} = str2double(tline((StartIndex + 1):end));
                index = 5;
            end
        end
        fclose(RecFid);
    else
        if (strfind(FileType, 'okrank'))
            RecFileName = [FileNames{i}, '.rec'];
        end
        RecFid = fopen(RecFileName, 'r');
        index = 0;
        while (index < 5)
            if (feof(RecFid))
                break;
            end
            tline = fgetl(RecFid);
            if ((strfind(tline, 'ai_freq')) > 0)
                ColonIndex = find(tline == ':');
                Fs = str2double(tline((ColonIndex + 1):end));
            end
            if ((strfind(tline, 'n_ai_chan')) > 0)
                ColonIndex = find(tline == ':');
                NChans{i} = str2double(tline((ColonIndex + 1):end));
            end
            if ((strfind(tline, 'n_samples')) > 0)
                index = 5;
                ColonIndex = find(tline == ':');
                RecordLengths{i} = str2double(tline((ColonIndex + 1):end));
                RecordLengths{i} = RecordLengths{i}/Fs;
            end
        end
        fclose(RecFid);
    end
end