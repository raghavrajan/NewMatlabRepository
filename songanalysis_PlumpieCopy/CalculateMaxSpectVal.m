function [MaxSpectVal] = CalculateMaxSpectVal(rawsong, Fs)

% Time axis
time = 0:1/Fs:(length(rawsong)/Fs);
time(end) = [];

filtsong = bandpass(rawsong,Fs,300,8000);

nfft=round(Fs*8/1000);
nfft = 2^nextpow2(nfft);
spect_win = hanning(nfft);
noverlap = round(0.9*length(spect_win)); %number of overlapping points       

[spect, freq, time_song] = spectrogram(filtsong,spect_win,noverlap,nfft,Fs,'yaxis');
MaxSpectVal = max(max(abs(spect)))/2;