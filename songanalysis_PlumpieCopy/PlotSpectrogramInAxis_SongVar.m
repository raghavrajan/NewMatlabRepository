function [] = PlotSpectrogramInAxis_SongVar(SongVar, Time, Fs, Axes)

% Plot song spectrogram

filtsong = bandpass(SongVar,Fs,300,10000);
Len = round(Fs*2/1000);
h = ones(1, Len)/Len;

nfft=round(Fs*8/1000);
nfft = 2^nextpow2(nfft);
spect_win = hanning(nfft);
noverlap = round(0.95*length(spect_win)); %number of overlapping points       

%now calculate spectrogram
%     [spect, freq, time_song] = specgram(filtsong, nfft, Fs, spect_win, noverlap);
[spect, freq, time_song] = spectrogram(filtsong,spect_win,noverlap,nfft,Fs,'yaxis');
idx_spect=scale_spect(spect);  %calculate index array for spectrogram
f_min = freq(1);
f_max = freq(length(freq));
freq_spect = [f_min, f_max];

t_min = Time(1); %convert to ms
t_max = Time(end); %convert to ms
%adjust time axis for spectrogram offset (1/2 window duration in ms)
%t_min = t_min + 0.5*nfft/Fs;  
%t_max = t_max + 0.5*nfft/Fs;  
Time_spect = [t_min, t_max];   
axes(Axes);
hold off;
cm = disp_idx_spect(idx_spect, Time_spect, freq_spect, -50, ...
        10, 1.2, 'hot', 'classic');
axis([t_min t_max 300 8000]);
set(gca, 'FontSize', 10);
set(gca, 'XTick', []);
ylabel('Frequency (Hz)', 'FontSize', 12);
zoom xon;
hold on;
%plot(time, (filtsong * 1000) + 12000);
%plot(time, (filtsong2 * 1000) + 14000,'r');
%NoteTimes = SmoothSong > 10*RMS;
%plot(time,NoteTimes * 14000,'k');
disp('Finished');