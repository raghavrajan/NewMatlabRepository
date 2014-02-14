function [WarpedMotifTemplate] = WarpMotifTemplatePhaseVocoder(RawData, Fs, TimeStretch, FreqStretch)

if (FreqStretch ~= 1)
    WarpedMotifRawData = pvoc(RawData, FreqStretch); 
    WarpedMotifRawData = resample(WarpedMotifRawData, length(RawData), length(WarpedMotifRawData));
else
    WarpedMotifRawData = RawData;
end

if (TimeStretch ~= 1)
    WarpedMotifRawData = pvoc(WarpedMotifRawData, TimeStretch);
end

WinSize = round(8/1000 * Fs);
WinOverlap = round(0.5 * WinSize);
[S1, F, T, P] = spectrogram(WarpedMotifRawData, hamming(WinSize), WinOverlap, WinSize, Fs);
Freq1 = find((F >= 860) & (F <= 8600));
S = log10(abs(S1(Freq1,:)));
S = (S - mean(reshape(S, 1, size(S,1)*size(S,2))))/std(reshape(S, 1, size(S,1)*size(S,2)));
WarpedMotifTemplate = S;
