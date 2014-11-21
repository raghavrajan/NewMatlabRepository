function [TempBoutFFT] = MA_CalcRhythmSpectrogram(SongBout, Fs, DiscreteBout_Fs, Freq)

[LogAmplitude] = ASSLCalculateLogAmplitudeAronovFee(SongBout, Fs, (1:1:length(SongBout))/Fs, 8, 0.9);
%LogAmplitude = spline((1:1:length(LogAmplitude))/Fs, LogAmplitude, (1:1:length(LogAmplitude)*DiscreteBout_Fs/Fs)/DiscreteBout_Fs);

DiscreteBout_Fs = Fs;

if (mod(length(LogAmplitude),2) == 1)
    LogAmplitude(end) = [];
end

LogAmplitude = LogAmplitude - mean(LogAmplitude(:));

NFFT = 8192*4;
%NFFT = length(LogAmplitude);

HanningWin = hanning(NFFT);
%LogAmplitude = HanningWin(:).*LogAmplitude(:);

TempBoutFFT = fft(LogAmplitude, NFFT)/NFFT;
%TempBoutFFT = TempBoutFFT/RMS_SongBout;
%TempBoutFFT = TempBoutFFT/(NFFT/DiscreteBout_Fs);
%TempBoutFFT = (TempBoutFFT).*conj(TempBoutFFT);
TempBoutFFT = 2*abs(TempBoutFFT(1:NFFT/2 + 1));
%TempBoutFFT = TempBoutFFT * sqrt((NFFT/sum(HanningWin)).^2);

%TempBoutFFT = TempBoutFFT/sum(TempBoutFFT);
%TempBoutFFT = TempBoutFFT/(NFFT/DiscreteBout_Fs);
%TempBoutFFT = TempBoutFFT/NFFT;
%TempBoutFFT = abs(TempBoutFFT).^2/NFFT;
%TempBoutFFT(2:end-1) = 2 * TempBoutFFT(2:end-1);

TempFreq = linspace(0, 1, NFFT/2+1) * DiscreteBout_Fs/2;

% Normalization routine

% disp(['RMS of log amplitude is ', num2str(RMS_SongBout)]);
% TempBoutFFT = TempBoutFFT/(sum(TempBoutFFT) * (DiscreteBout_Fs/NFFT)); % Normalization by (sum of all magnitudes times the frequency resolution)
%TempBoutFFT = TempBoutFFT/(NFFT);
%TempBoutFFT = TempBoutFFT/(sum(TempBoutFFT) * (DiscreteBout_Fs/NFFT)); % Normalization by (sum of all magnitudes times the frequency resolution)
% TempBoutFFT = TempBoutFFT/(RMS_SongBout);

TempBoutFFT = interp1(TempFreq, TempBoutFFT, Freq);