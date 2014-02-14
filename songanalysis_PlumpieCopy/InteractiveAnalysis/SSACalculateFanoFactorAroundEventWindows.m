function [FanoFactor] = SSACalculateFanoFactorAroundEventWindows(SpikeTrain, MedianMotif, WindowWidth, PreSongStartDuration, PreSongEndDuration, EventTimes, EventWindowWidth, ContextString)

PST = [];

for i = 1:length(SpikeTrain),
    TempSpikeTrain = SpikeTrain{i};
    for k = 1:length(EventTimes),
        Index = 1;
        for j = (EventTimes(k) - EventWindowWidth):0.001:(EventTimes(k) + EventWindowWidth - WindowWidth),
            PST{k}(i,Index) = length(find((TempSpikeTrain >= j) & (TempSpikeTrain < (j + WindowWidth))));
            Index = Index + 1;
        end
    end
end

for k = 1:length(EventTimes),
    VarPST{k} = var(PST{k});
    MeanPST{k} = mean(PST{k});
    FanoFactor(k) = mean(VarPST{k}(find(MeanPST{k}))./MeanPST{k}(find(MeanPST{k})));
end

disp([ContextString, ': The window width is ', num2str(WindowWidth * 1000), ' ms and the mean Fano factor is ', num2str(mean(FanoFactor))]);
FanoFactor = [WindowWidth mean(FanoFactor)];
