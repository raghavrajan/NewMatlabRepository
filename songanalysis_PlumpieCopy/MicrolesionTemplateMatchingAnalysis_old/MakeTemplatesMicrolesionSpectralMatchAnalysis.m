function [MotifTemplate] = MakeTemplatesMicrolesionSpectralMatchAnalysis(RawSong, Fs, XDur, Sequence, TimeStretch, FreqStretch, varargin)

Time = (0:1:length(RawSong)-1)/Fs;

if (nargin > 6)
    sm_win = varargin{1};
    FFTWinSize = sm_win;
    FFTWinOverlap = varargin{2}; % step size in ms
else
    sm_win = 8; % in ms
    FFTWinSize = sm_win;
    FFTWinOverlap = 4; % step size in ms
end

Motif = RawSong(find((Time >= XDur(1)) & (Time <= XDur(2))));

Index = 0;
for i = 1:length(TimeStretch),
    for j = 1:length(FreqStretch), 
        Index = Index + 1;
        MotifTemplate(Index).Label = Sequence;
        MotifTemplate(Index).FFTWinSize = FFTWinSize;
        MotifTemplate(Index).FFTWinOverlap = FFTWinOverlap;
        MotifTemplate(Index).TimeStretch = TimeStretch(i);
        
%        [MotifTemplate(Index).MotifTemplate] = WarpMotifTemplatePhaseVocoder(Motif, Fs, (1 - TimeStretch(i)/100), (1 - FreqStretch(j)/100));
        [PitchShiftedMotif] = PitchShiftSound(Motif, Fs, (1 - FreqStretch(j)/100));
        
        if ((i == 1) && (j == 1))
            MotifTemplate(Index).Fs = Fs;
            MotifTemplate(Index).RawSound = Motif;
        end
        
        MotifTemplate(Index).PitchShiftedSound = PitchShiftedMotif;
        [MotifTemplate(Index).MotifTemplate] = StretchCompressTemplates(PitchShiftedMotif, Fs, (1 + TimeStretch(i)/100), (1 - FreqStretch(j)/100), FFTWinSize, FFTWinOverlap);
        MotifTemplate(Index).FreqStretch = FreqStretch(j);
    end
end
