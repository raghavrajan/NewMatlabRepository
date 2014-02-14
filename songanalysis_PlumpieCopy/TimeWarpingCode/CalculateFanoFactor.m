function [FanoFactor_30, FanoFactor_100] = CalculateFanoFactor(FileInfo,MedianMotif,Latency)

PST_30 = [];

for i = 1:length(FileInfo.WSpikeTrain),
    TempSpikeTrain = FileInfo.WSpikeTrain{i};
    Index = 1;
    for j = -0.1:0.001:(MedianMotif.Length + Latency - 0.03),
        PST_30(i,Index) = length(find((TempSpikeTrain > j) & (TempSpikeTrain < (j + 0.03))));
        Index = Index + 1;
    end
end

FanoFactor_30 = var(PST_30)./mean(PST_30);

PST_100 = [];

for i = 1:length(FileInfo.WSpikeTrain),
    TempSpikeTrain = FileInfo.WSpikeTrain{i};
    Index = 1;
    for j = -0.1:0.001:(MedianMotif.Length + Latency - 0.1),
        PST_100(i,Index) = length(find((TempSpikeTrain >= j) & (TempSpikeTrain < (j + 0.1))));
        Index = Index + 1;
    end
end

FanoFactor_100 = var(PST_100)./mean(PST_100);
