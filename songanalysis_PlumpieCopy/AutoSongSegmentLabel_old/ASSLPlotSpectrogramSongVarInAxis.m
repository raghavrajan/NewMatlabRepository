function [] = ASSLPlotSpectrogramSongVarInAxis(rawsong, Fs, SpecAxis)

% Time axis
time = 0:1/Fs:(length(rawsong)/Fs);
time(end) = [];

start_time_index = 1;
end_time_index = length(time);

% Plot song spectrogram

time = time(start_time_index:end_time_index);
rawsong = rawsong(start_time_index:end_time_index);

filtsong = bandpass(rawsong,Fs,300,10000);
%filtsong = bandpass_fft_filter(rawsong,8000,1500,Fs);
%filtnoise = bandpass_fft_filter(rawsong,1200,100,Fs);
Len = round(Fs*2/1000);
h = ones(1, Len)/Len;

nfft=round(Fs*8/1000);
nfft = 2^nextpow2(nfft);
spect_win = hanning(nfft);
noverlap = round(0.9*length(spect_win)); %number of overlapping points       

axes(SpecAxis);
%now calculate spectrogram
%     [spect, freq, time_song] = specgram(filtsong, nfft, Fs, spect_win, noverlap);
[spect, freq, time_song] = spectrogram(filtsong,spect_win,noverlap,nfft,Fs,'yaxis');
idx_spect=scale_spect(spect);  %calculate index array for spectrogram
f_min = freq(1);
f_max = freq(length(freq));
freq_spect = [f_min, f_max];

t_min = time(1); %convert to ms
t_max = time(end); %convert to ms
%adjust time axis for spectrogram offset (1/2 window duration in ms)
%t_min = t_min + 0.5*nfft/Fs;  
%t_max = t_max + 0.5*nfft/Fs;  
time_spect = [t_min, t_max];   
cm = disp_idx_spect(idx_spect, time_spect, freq_spect, -50, ...
        10, 1.2, 'hot', 'classic');
axis([t_min t_max 300 8000]);
set(gca, 'FontSize', 8);
xlabel('Time (sec)', 'FontSize', 10);
ylabel('Frequency (Hz)', 'FontSize', 10);
zoom xon;
hold on;
%plot(time, (filtsong * 1000) + 12000);
%plot(time, (filtsong2 * 1000) + 14000,'r');
%NoteTimes = SmoothSong > 10*RMS;
%plot(time,NoteTimes * 14000,'k');
disp('Finished');