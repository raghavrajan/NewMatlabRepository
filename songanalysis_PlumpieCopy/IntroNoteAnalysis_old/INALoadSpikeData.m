function [SpikeData, ASpikeData] = INALoadSpikeData(FileNames, SpikeFileDir, ClusterNos, SpikeChanNo, FileLength, ContinuousOrNot)

PresentDir = pwd;
cd(SpikeFileDir);

if (strfind(ContinuousOrNot, 'continuous'))
    CSpikeData = [];
else
    SpikeData = [];
end

FTime = 0;
for i = 1:length(FileNames),
    Temp = load(['Chan', num2str(SpikeChanNo), '_', FileNames{i}, '.spk']);
    disp(['Loaded ', ['Chan', num2str(SpikeChanNo), '_', FileNames{i}, '.spk']]);
    
    if (Temp(end,2) > FileLength{i})
        Temp(:,2) = Temp(:,2)/1000;
    end
    
    TempSpikeIndices = [];
    for j = 1:length(ClusterNos),
        TempSpikeIndices = [TempSpikeIndices; find(Temp(:,1) == ClusterNos(j))];
    end
    SpikeData{i} = Temp(TempSpikeIndices, 2);
    if (strfind(ContinuousOrNot, 'continuous'))
        CSpikeData{1} = [CSpikeData{1}; (Temp(TempSpikeIndices,2) + FTime)];
        FTime = FTime + FileLength{i};
    end
end

ASpikeData = SpikeData;

if (strfind(ContinuousOrNot, 'continuous'))
    SpikeData = CSpikeData;
end

cd(PresentDir);