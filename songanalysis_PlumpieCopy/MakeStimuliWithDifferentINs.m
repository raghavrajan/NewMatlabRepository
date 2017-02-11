function [OutputData] = MakeStimuliWithDifferentINs(Motif, Gap, IN, NumINs, NumBouts, GapBetweenBouts, Side, Fs)

% GapBetweenBouts in sec.

SongBout = [];
IN = IN(:);
Motif = Motif(:);
Gap = Gap(:);

for i = NumINs:-1:1,
    SongBout = [SongBout; IN*(1 - (i-1)*0.1); Gap];
end

SongBout = [SongBout; Motif];

BoutGap = [];
for i = 1:round(GapBetweenBouts/(length(Gap)/Fs));
    BoutGap = [BoutGap; Gap];
end

Song = [];
for i = 1:NumBouts,
    Song = [Song; SongBout; BoutGap];
end

Song = Song/max(Song);

Silence = repmat(Gap, floor(size(Song,1)/size(Gap,1)), 1);
Silence = [Silence; Gap(1:(size(Song,1) - size(Silence,1)))];

if (strfind(Side, 'Left'))
    OutputData = [Song Silence];
else
    OutputData = [Silence Song];
end
    
