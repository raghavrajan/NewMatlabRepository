function [ActualCrossings] = FindThresholdCrossings(Data, Threshold)

Crossings = find(abs(Data - Threshold) <= 0.5);

PotentialCrossings = find(diff(Crossings) > 1);

if (isempty(PotentialCrossings))
    ActualCrossings = [];
    return;
end

[MinVal, Index] = min(abs(Data(Crossings(1:PotentialCrossings(1)))));

ActualCrossings(1) = Crossings(Index);

for i = 2:length(PotentialCrossings),
    [MinVal, Index] = min(abs(Data(Crossings((PotentialCrossings(i-1)+1):PotentialCrossings(i)))));
    ActualCrossings(i) = Crossings(Index + PotentialCrossings(i-1)); 
end

[MinVal, Index] = min(abs(Data(Crossings((PotentialCrossings(end)+1):end))));
ActualCrossings(end+1) = Crossings(Index + PotentialCrossings(end));

