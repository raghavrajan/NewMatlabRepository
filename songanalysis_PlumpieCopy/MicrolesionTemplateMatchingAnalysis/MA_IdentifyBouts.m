function [Bouts] = MA_IdentifyBouts(GapDurs, Lens, Fs, Onsets, Offsets, InterBoutInterval, Index)

Onsets = Onsets(:);
Offsets = Offsets(:);

LongGaps = find(GapDurs >= InterBoutInterval);

if (~isempty(LongGaps))
    Bouts = [Onsets(1) Onsets(LongGaps + 1)'; Offsets(LongGaps)' Offsets(end)]';
    Bouts = [Bouts [1:1:size(Bouts,1)]' ones(size(Bouts,1),1)*Index];
else
    Bouts = [-100 -100 1 Index];
end