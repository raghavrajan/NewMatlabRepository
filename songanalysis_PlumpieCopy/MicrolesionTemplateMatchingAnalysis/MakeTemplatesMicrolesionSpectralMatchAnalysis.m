function [SyllableTemplates, TemplatePNGFig_Syllables, TemplatePNGFig_SyllLabelTimes] = MakeTemplatesMicrolesionSpectralMatchAnalysis(RawSong, Fs, OnsetsOffsets, labels, UniqueLabels, TimeStretch, FreqStretch, varargin)

Time = (0:1:length(RawSong)-1)/Fs;

if (nargin > 7)
    sm_win = varargin{1};
    FFTWinSize = sm_win;
    FFTWinOverlap = varargin{2}; % step size in ms
    Normalization = varargin{3};
else
    sm_win = 8; % in ms
    FFTWinSize = sm_win;
    FFTWinOverlap = 4; % step size in ms
    Normalization = 1;
end

if (nargin > 10)
    PitchShiftedSong = varargin{4};
    S = varargin{5};
    T = varargin{6};
else
    Index = 0;
    for j = 1:length(FreqStretch),
        fprintf('f%d ', j);
        PitchShiftedSong{j} = PitchShiftSound(RawSong, Fs, (1 - FreqStretch(j)/100));
        [P, F, S1, T{j}] = CalculateMultiTaperSpectrogram(PitchShiftedSong{j}, Fs, FFTWinSize, FFTWinOverlap, 1.5);
        Freq1 = find((F >= 860) & (F <= 8600));
        S{j} = log10(abs(S1(Freq1,:)));
    end
end

fprintf('\n');

PaddingTime = 0.1; % in ms
TemplatePNGFig_Syllables = [zeros(round(PaddingTime * Fs),1)];
TempSyllTime = PaddingTime;
TemplatePNGFig_SyllLabelTimes = [];

for k = 1:length(UniqueLabels),
    Indices = find(labels == UniqueLabels(k));
    fprintf('%c ', UniqueLabels(k));
    for SyllIndex = 1:min(3,length(Indices)),
        fprintf('>');
        XDur = [OnsetsOffsets(Indices(SyllIndex), 1) OnsetsOffsets(Indices(SyllIndex), 2)];
        Index = 0;
        for j = 1:length(FreqStretch),
            for i = 1:length(TimeStretch),
                Index = Index + 1;
                
                SyllableTemplates{k}{SyllIndex}.MotifTemplate(Index).Label = UniqueLabels(k);
                SyllableTemplates{k}{SyllIndex}.MotifTemplate(Index).FFTWinSize = FFTWinSize;
                SyllableTemplates{k}{SyllIndex}.MotifTemplate(Index).FFTWinOverlap = FFTWinOverlap;
                SyllableTemplates{k}{SyllIndex}.MotifTemplate(Index).TimeStretch = TimeStretch(i);

        %        [MotifTemplate(Index).MotifTemplate] = WarpMotifTemplatePhaseVocoder(Motif, Fs, (1 - TimeStretch(i)/100), (1 - FreqStretch(j)/100));

        %        [PitchShiftedMotif] = PitchShiftSound(Motif, Fs, (1 - FreqStretch(j)/100));

                % Shortened the offset by 12.5% since it seemed to be consistently off by a
                % about 12.5% of the total duration
                Motif = RawSong(find((Time >= XDur(1)) & (Time <= (XDur(2) - (0.125*(XDur(2) - XDur(1)))))));
                PitchShiftedMotif = PitchShiftedSong{j}(find((Time >= XDur(1)) & (Time <= (XDur(2) - (0.125*(XDur(2) - XDur(1)))))));

                % Add motif to templates if time stretch = 0
                if ((SyllIndex == 1) && (i == 1))
                    TemplatePNGFig_Syllables = [TemplatePNGFig_Syllables; Motif; zeros(round(PaddingTime*Fs),1)];
                    TemplatePNGFig_SyllLabelTimes(end+1) = TempSyllTime;
                    TempSyllTime = TempSyllTime + length(Motif)/Fs + PaddingTime;
                end
                
                if ((i == 1) && (j == 1))
                    SyllableTemplates{k}{SyllIndex}.MotifTemplate(Index).Fs = Fs;
                    SyllableTemplates{k}{SyllIndex}.MotifTemplate(Index).RawSound = Motif;
                end

                SyllableTemplates{k}{SyllIndex}.MotifTemplate(Index).PitchShiftedSound = PitchShiftedMotif;

                PitchShiftedMotifSpectrogram = S{j}(:,find((T{j} >= XDur(1)) & (T{j} <= (XDur(2) - (0.125*(XDur(2) - XDur(1)))))));
                PitchShiftedMotifTime = T{j}(find((T{j} >= XDur(1)) & (T{j} <= (XDur(2) - (0.125*(XDur(2) - XDur(1)))))));

                [SyllableTemplates{k}{SyllIndex}.MotifTemplate(Index).MotifTemplate] = StretchCompressTemplates(PitchShiftedMotifSpectrogram, PitchShiftedMotifTime, (1 + TimeStretch(i)/100), Normalization);
                SyllableTemplates{k}{SyllIndex}.MotifTemplate(Index).Normalization = Normalization;
                SyllableTemplates{k}{SyllIndex}.MotifTemplate(Index).FreqStretch = FreqStretch(j);
            end
        end
    end
    fprintf('\n');
end

