function [] = SplitSpikeSortFiles(DirectoryName, SortedSpikeFileName, FileName, FileNos, FileExt)

PresentDirectory = pwd;

cd(DirectoryName);

SortedSpikeTimes = load(SortedSpikeFileName);

RecTime = 0;

for i = 1:length(FileNos),
    
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
        DataFileName = [DataFile,FileExt];
    else
        DataFileName = [DataFile,'.',FileExt];
    end
    
    RecFileName = strcat(DataFile,'.rec');

    fid = fopen(RecFileName,'r');
    index = 0;
    while (index < 5)
        if (feof(fid))
            break;
        end
        tline = fgetl(fid);
        if (strfind(tline,'rec end') > 0)
            index = 5;
        end
    end

    StartIndex = find(tline == '=');
    EndIndex = strfind(tline,'ms');
    TempRecTime = str2double(tline((StartIndex + 1):(EndIndex - 1)));
    TempRecTime = TempRecTime/1000;
    fclose(fid);
    
    OutputFileName = [DataFileName, SortedSpikeFileName((end - 7):end)];
    TempSpikes = SortedSpikeTimes(find((SortedSpikeTimes > RecTime) & (SortedSpikeTimes <= (RecTime + TempRecTime))));
    TempSpikes = TempSpikes - RecTime;
    save(OutputFileName, 'TempSpikes','-ASCII');
    RecTime = RecTime + TempRecTime;
end
