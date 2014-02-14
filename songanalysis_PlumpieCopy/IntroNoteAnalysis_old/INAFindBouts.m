function [BoutIndices, BoutLengths] = INAFindBouts(Onsets, Offsets, InterBoutInterval)

BoutIndices = [];
BoutLengths = [];

InterSyllIntervals = Onsets(2:end) - Offsets(1:(end-1));
Bouts = find(InterSyllIntervals > InterBoutInterval);

for i = 0:length(Bouts),
    if (i ~= length(Bouts))
        if (i == 0)
            BoutIndices = [BoutIndices; 1; Bouts(i+1)];
            BoutLengths = [BoutLengths; (Offsets(Bouts(i+1)) - Onsets(1))];
        else
            BoutIndices = [BoutIndices; (Bouts(i) + 1); Bouts(i+1)];
            BoutLengths = [BoutLengths; (Offsets(Bouts(i+1)) - Onsets(Bouts(i) + 1))];
        end
    else
        if (i == 0)
            BoutIndices = [BoutIndices; 1; length(Onsets)];
            BoutLengths = [BoutLengths; (Offsets(end) - Onsets(1))];
        else
            BoutIndices = [BoutIndices; (Bouts(i) + 1); length(Onsets)];
            BoutLengths = [BoutLengths; (Offsets(end) - Onsets(Bouts(i) + 1))];
        end
    end
end