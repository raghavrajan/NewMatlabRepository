function [] = PlotSingleChannelAndSongData(pathname,filename,ChannelNo);

% Get Spike Data from observer file
channel_string = strcat('obs',num2str(ChannelNo),'r');
[rawsong,Fs] = soundin_copy(pathname,filename,channel_string);

% Convert to uV - 5V on the data acquisition is 32768
rawsong = rawsong * 5/32768;
time = 0:1/Fs:length(rawsong)/Fs;
time(end) = [];


figure;
subplot(2,1,1);
plot(time,rawsong,'b');
axis tight;

ChannelNo = 0;

% Get Song Data from observer file
channel_string = strcat('obs',num2str(0),'r');
[rawsong,Fs] = soundin_copy(pathname,filename,channel_string);

% Convert to uV - 5V on the data acquisition is 32768
rawsong = rawsong * 1/32768;

% Time axis
time = 0:1/32000:(length(rawsong)/32000);
time(end) = [];

subplot(2,1,2);

filtsong = bandpass(rawsong,Fs,300,10000);
    
nfft=round(Fs*8/1000);
nfft = 2^nextpow2(nfft);
spect_win = hanning(nfft);
noverlap = round(0.9*length(spect_win)); %number of overlapping points       
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
%t_min = t_min + 0.5*nfft*1000/Fs;  
%t_max = t_max + 0.5*nfft*1000/Fs;  
t_min = t_min;  
t_max = t_max;  
time_spect = [t_min, t_max];                

cm = disp_idx_spect(idx_spect, time_spect, freq_spect, -70, ...
    5, 1.5, 'gray', 'classic');
axis([time(1) time(end) 300 10000]);
set(gca,'ytick',[]);
set(gca,'Visible','off');
set(gca,'xtick',[]);
% ylabel('Frequency (Hz)','FontSize',16);
disp('Finished');