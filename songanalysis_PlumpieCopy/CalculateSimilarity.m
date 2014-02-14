function [DirSim, UnDirSim] = CalculateSimilarity(BirdName, Motif, FileType)

cd(['/media/FreeAgent GoFlex Drive/RaghavData/HVC_Microlesions/', BirdName]);

Files = dir(['*.', Motif, '.*AllFiles*.mat']);
for i = 1:length(Files),
    Data{i} = load(Files(i).name);
end

for i = 1:length(Data),
    Temp = cell2mat(Data{i}.DirBout.MaxBoutSeqMatchVal);
    DirSim(i,:) = [mean(Temp(2:2:end)) std(Temp(2:2:end))];
    
    UnDirTrials = [];
    clear UnDirFileNo;
    for j = 1:length(Data{i}.UnDirBout.FileName),
        DotIndex = find(Data{i}.UnDirBout.FileName{j} == '.');
        if (strfind(FileType, 'obs'))
            UnDirFileNo(j) = str2double(Data{i}.UnDirBout.FileName{j}(DotIndex(end-1)+1:DotIndex(end)-1));
        else
            if (strfind(FileType, 'wav'))
                UnDirFileNo(j) = str2double(Data{i}.UnDirBout.FileName{j}(DotIndex(end)+1:end));
            else
                if (strfind(FileType, 'okrank'))
                   UnDirFileNo(j) = str2double(Data{i}.UnDirBout.FileName{j}((end-5):(end-4))) + str2double(Data{i}.UnDirBout.FileName{j}((end-3):(end-2)))/60 + str2double(Data{i}.UnDirBout.FileName{j}((end-1):(end)))/3600;
                end
            end
        end
    end
    for j = 1:length(Data{i}.DirBout.FileName),
        DotIndex = find(Data{i}.DirBout.FileName{j} == '.');
        if (strfind(FileType, 'obs'))
            FileNo = str2double(Data{i}.DirBout.FileName{j}(DotIndex(end-1)+1:DotIndex(end)-1));
        else
            if (strfind(FileType, 'wav'))
                FileNo = str2double(Data{i}.DirBout.FileName{j}(DotIndex(end)+1:end));
            else
                if (strfind(FileType, 'okrank'))
                    FileNo = str2double(Data{i}.DirBout.FileName{j}((end-5):(end-4))) + str2double(Data{i}.DirBout.FileName{j}((end-3):(end-2)))/60 + str2double(Data{i}.DirBout.FileName{j}((end-1):(end)))/3600;
                end
            end
        end
        
        LastUnDirFileBeforeDir = find(UnDirFileNo < FileNo, 1, 'last');
        if (~isempty(LastUnDirFileBeforeDir))
            if (LastUnDirFileBeforeDir > 5)
                UnDirTrials = [UnDirTrials [LastUnDirFileBeforeDir-5:1:LastUnDirFileBeforeDir]];
            else
                UnDirTrials = [UnDirTrials [1:1:LastUnDirFileBeforeDir]];
            end
        end
        
        FirstUnDirFileAfterDir = find(UnDirFileNo > FileNo, 1, 'first');
        if (~isempty(FirstUnDirFileAfterDir))
            if (FirstUnDirFileAfterDir < (length(UnDirFileNo) - 5))
                UnDirTrials = [UnDirTrials [FirstUnDirFileAfterDir:1:FirstUnDirFileAfterDir+4]];
            else
                UnDirTrials = [UnDirTrials [FirstUnDirFileAfterDir:1:length(UnDirFileNo)]];
            end
        end
    end
    if (~isempty(UnDirTrials))
        UnDirTrials = unique(UnDirTrials);
    else
        UnDirTrials = 1:1:length(UnDirFileNo);
    end
    Temp = cell2mat(Data{i}.UnDirBout.MaxBoutSeqMatchVal(UnDirTrials));
    UnDirSim(i,:) = [mean(Temp(2:2:end)) std(Temp(2:2:end))];
end        