function [Bouts] = ASSLGetBoutEdges(SyllOnsets, SyllOffsets, InterBoutInterval, BoutPaddingTime, FileDur)

Bouts = [];

Intervals = SyllOnsets(2:end) - SyllOffsets(1:end-1);
LongIntervals = find(Intervals >= InterBoutInterval);
if (~isempty(LongIntervals))
    for i = 1:length(LongIntervals),
        if (i == 1)
            BoutOnset = (SyllOnsets(1) - BoutPaddingTime);
            BoutOffset = (SyllOffsets(LongIntervals(i)) + BoutPaddingTime);    
        else
            BoutOnset = (SyllOnsets(LongIntervals(i-1)+1) - BoutPaddingTime);
            BoutOffset = (SyllOffsets(LongIntervals(i)) + BoutPaddingTime);
        end
        if (BoutOnset > (InterBoutInterval - BoutPaddingTime))
            Bouts(end+1,:) = [BoutOnset BoutOffset];
        end
    end
    BoutOnset = (SyllOnsets(LongIntervals(end)+1) - BoutPaddingTime);
    BoutOffset = (SyllOffsets(end) + BoutPaddingTime);
    if (BoutOffset < (FileDur*1000 - (InterBoutInterval - BoutPaddingTime)))
        Bouts(end+1,:) = [BoutOnset BoutOffset];
    end
else
    BoutOnset = (SyllOnsets(1) - BoutPaddingTime);
    BoutOffset = (SyllOffsets(end) + BoutPaddingTime);
    if ((BoutOnset > (InterBoutInterval - BoutPaddingTime)) && (BoutOffset < (1000*FileDur - (InterBoutInterval - BoutPaddingTime))))
        Bouts(end+1,:) = [BoutOnset BoutOffset];
    end
end

BoutsToBeRemoved = [];
for i = 1:size(Bouts,1),
    NumSyllables = find((SyllOnsets >= Bouts(i,1)) & (SyllOnsets <= Bouts(i,2)));
    if (length(NumSyllables) < 3)
        BoutsToBeRemoved(end+1) = i;
    end
end

if (~isempty(BoutsToBeRemoved))
    Bouts(BoutsToBeRemoved,:) = [];
end