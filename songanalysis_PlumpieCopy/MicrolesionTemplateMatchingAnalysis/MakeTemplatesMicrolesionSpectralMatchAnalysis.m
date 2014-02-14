function [MotifTemplate] = MakeTemplatesMicrolesionSpectralMatchAnalysis(RawSong, Fs, XDur, Sequence, TimeStretch, FreqStretch)

Time = (0:1:length(RawSong)-1)/Fs;
sm_win = 8; % in ms
FFTWinSize = sm_win;
FFTWinOverlap = 4; % step size in ms
Motif = RawSong(find((Time >= XDur(1)) & (Time <= XDur(2))));

Index = 0;
for i = 1:length(TimeStretch),
    for j = 1:length(FreqStretch), 
        Index = Index + 1;
        MotifTemplate(Index).Label = Sequence;
        MotifTemplate(Index).FFTWinSize = FFTWinSize;
        MotifTemplate(Index).FFTWinOverlap = FFTWinOverlap;
        MotifTemplate(Index).TimeStretch = TimeStretch(i);
        
        % [MotifTemplate(Index).MotifTemplate] = WarpMotifTemplatePhaseVocoder(Motif, Fs, (1 - TimeStretch(i)/100), (1 - FreqStretch(j)/100));
        [MotifTemplate(Index).MotifTemplate] = StretchCompressTemplates(Motif, Fs, (1 + TimeStretch(i)/100), (1 - FreqStretch(j)/100));
        MotifTemplate(Index).FreqStretch = FreqStretch(j);
    end
end
