function [LogAmplitude] = ASSLCalculateLogAmplitudeKao(RawData, Fs, Time, FFTWinSize, FFTWinOverlap, varargin)

% WinSize = round(FFTWinSize * Fs/1000);
% WinOverlap = round(WinSize * FFTWinOverlap);
% [S, F, T, P] = spectrogram(RawData, hamming(WinSize), WinOverlap, WinSize, Fs);
% Freq = find((F >= 300) & (F <= 8000));
% LogAmplitude = 10*log10(sum((S(Freq,:)).*conj((S(Freq,:)))));
% LogAmplitude = spline(T, LogAmplitude, Time);

if (nargin <= 5)
    F_High = 300; % high pass filter cut-off
    F_Low = 8000; % low pass filter cut-off
else
    F_High = varargin{1}; % high pass filter cut-off
    F_Low = varargin{2}; % low pass filter cut-off
end
SmoothWinSize = 2; % in ms
SmoothWin = ones(round(SmoothWinSize * Fs/1000),1);
SmoothWin = SmoothWin/sum(SmoothWin);

LogAmplitude = conv(abs(RawData), SmoothWin, 'same');
