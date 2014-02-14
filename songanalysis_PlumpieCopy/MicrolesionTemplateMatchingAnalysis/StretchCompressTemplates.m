function [WarpedMotifTemplate] = StretchCompressTemplates(RawData, Fs, TimeStretch, FreqStretch)

WinSize = 8; % in ms
WinOverlap = 4; % step size in ms
%[S1, F, T, P] = spectrogram(RawData, hamming(WinSize), WinOverlap, WinSize, Fs);
[P, F, S1, T] = CalculateMultiTaperSpectrogram(RawData, Fs, WinSize, WinOverlap, 1.5);
Freq1 = find((F >= 860) & (F <= 8600));
S = log10(abs(S1(Freq1,:)));

x = T;
xx = linspace(T(1), T(end), round(TimeStretch*length(T)));
for i = 1:size(S, 1),
    WarpedS(i,:) = spline(x, S(i,:), xx);
%    WarpedS(i,:) = interp1(x, S(i,:), xx);
end
WarpedS = (WarpedS - mean(WarpedS(:)))/std(WarpedS(:));
WarpedMotifTemplate = WarpedS;

