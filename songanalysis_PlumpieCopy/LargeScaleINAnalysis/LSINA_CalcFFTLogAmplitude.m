function [LogAmplitude] = LSINA_CalcFFTLogAmplitude(RawSong, Fs, FFTWindowSize, FFTWindowStep)

FiltSong = bandpass(RawSong,Fs,300,10000);

nfft=round(Fs*FFTWindowSize);
spect_win = hanning(nfft);
noverlap = nfft - round(FFTWindowStep * Fs); %number of overlapping points       

%now calculate power spectrum
[spect, freq, time_song, Power] = spectrogram(FiltSong,spect_win,noverlap,nfft,Fs,'yaxis');

NeccFreqs = find((freq >= 860) & (freq <= 8600)); % am using the power between 860 and 8600 Hz only

LogAmplitude = 10*log10(sum(Power(NeccFreqs,:)));
       