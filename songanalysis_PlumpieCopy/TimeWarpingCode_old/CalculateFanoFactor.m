function [FanoFactor] = CalculateFanoFactor(FileInfo,MedianMotif,Latency)

PST = [];
for i = 1:length(FileInfo.SpikeTrain),
    TempSpikeTrain = FileInfo.SpikeTrain{i};
    Index = 1;
    for j = 0:0.001:(MedianMotif.Length + Latency),
        PST(i,Index) = length(find((TempSpikeTrain > j) & (TempSpikeTrain < (j + 0.03))));
        Index = Index + 1;
    end
end

FanoFactor = var(PST)./mean(PST);
