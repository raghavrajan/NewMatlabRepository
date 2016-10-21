function [WarpedMotifTemplate] = WarpMotifTemplatePhaseVocoder(RawData, Fs, TimeStretch, FreqStretch)

if (FreqStretch ~= 1)
    % WarpedMotifRawData = pvoc(RawData, FreqStretch); 
    % Using another algorithm to shift pitch - seems to work pretty well
    % It expresses changes in pitch in terms of cents
    % To calculate cents I used the formula
    % cent change = 1200 * log(f2/f1) / log(2)
    parameter.fsAudio = Fs;
    WarpedMotifRawData = pitchShiftViaTSM(RawData, 1200 * log(FreqStretch) / log(2), parameter); 
    % WarpedMotifRawData = resample(WarpedMotifRawData, length(RawData), length(WarpedMotifRawData));
else
    WarpedMotifRawData = RawData;
end

if (TimeStretch ~= 1)
    WarpedMotifRawData = pvoc(WarpedMotifRawData, TimeStretch);
end

WinSize = round(8/1000 * Fs);
WinOverlap = round(0.5 * WinSize);
% [S1, F, T, P] = spectrogram(WarpedMotifRawData, hamming(WinSize), WinOverlap, WinSize, Fs);
[P, F, S1, T] = CalculateMultiTaperSpectrogram(WarpedMotifRawData(:)', Fs, 8, 4, 1.5);
Freq1 = find((F >= 860) & (F <= 8600));
S = log10(abs(S1(Freq1,:)));
S = (S - mean(reshape(S, 1, size(S,1)*size(S,2))))/std(reshape(S, 1, size(S,1)*size(S,2)));
WarpedMotifTemplate = S;
