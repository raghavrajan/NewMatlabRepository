function [spike_times] = threshold_observer_spike_data(pathname,filename,channel_no,upper_threshold,lower_threshold, varargin)

if (nargin > 5)
	Range(1) = varargin{1};
    Range(2) = varargin{2};
end

channel_string = ['obs',num2str(channel_no),'r'];

[rawdata,Fs] = soundin_copy(pathname,filename,channel_string);
rawdata = rawdata * 500/32768;

% upper_thresh_cross = find(rawdata > upper_threshold);
% lower_thresh_cross = find(rawdata > lower_threshold);
% 
% real_spikes = [];
% spike_index = 1;
% for i = 1:length(upper_thresh_cross),
%     index = find(lower_thresh_cross > upper_thresh_cross(i),1,'first');
%     if (index < 24)
%         real_spikes(spike_index) = upper_thresh_cross(i);
%         spike_index = spike_index + 1;
%     end
% end


[spike_times, waveforms] = FindSpikes(rawdata,lower_threshold,upper_threshold,32,16,0,Fs);
waveforms = waveforms';

time = (1:1:(length(rawdata)))/Fs;
if (~exist('Range', 'var'))
    Range(1) = time(1);
    Range(2) = time(end);
end

RangeIndex = round(Range * Fs);

figure;
set(gcf,'Color','w');

axes('position',[0.15 0.565 0.8 0.2]);
plot(time,rawdata,'b');
axis tight;
hold on;

if (length(spike_times) > 0)
    for i = 1:length(spike_times),
        index = find(time <= spike_times(i), 1, 'last');
        plot(time((index-8):(index + 23)), rawdata((index - 8):(index + 23)),'r');
    end
end
axis tight;
ylabel('Amplitude (\muV)','FontSize',16,'FontWeight','bold');
set(gca,'FontSize',14,'FontWeight','bold');
set(gca,'XTickLabel',[]);
% set(gca,'ytick',[]);

axes('position',[0.15 0.34 0.8 0.2]);
plot(time,rawdata);
axis tight;
set(gca,'FontSize',14,'FontWeight','bold');

axes('position',[0.15 0.065 0.8 0.2]);
set(gca,'FontSize',18);
channel_string = strcat('obs0r');
% pathname = pwd;
% pathname = strcat(pathname,'/');
[RawData,Fs] = soundin_copy(pathname,filename,channel_string);
Times = 0:1/Fs:length(RawData)/Fs;
Times(end) = [];
% plot(Times,RawData);
% axis([-0.2 t_max min(RawData) max(RawData)]);

filtsong = bandpass(RawData,Fs,300,10000);

nfft=round(Fs*8/1000);
nfft = 2^nextpow2(nfft);
spect_win = hanning(nfft);
noverlap = round(0.9*length(spect_win)); %number of overlapping points       
%now calculate spectrogram
% [spect, freq, time] = specgram(filtsong, nfft, Fs, spect_win, noverlap);
[spect, freq, time] = spectrogram(filtsong,spect_win,noverlap,nfft,Fs,'yaxis');

idx_spect=scale_spect(spect);  %calculate index array for spectrogram
f_min = freq(1);
f_max = freq(length(freq));
freq_spect = [f_min, f_max];
time = time;
t_min = time(1); 
t_max = time(length(time)); 
%adjust time axis for spectrogram offset (1/2 window duration in ms)
t_min = t_min + 0.5*nfft*1/Fs;  
t_max = t_max + 0.5*nfft*1/Fs;  
t_min = t_min;  
t_max = t_max;  
time_spect = [t_min, t_max];                

disp_idx_spect(idx_spect, time_spect, freq_spect, -55, ...
	0, 1.5, 'gray', 'classic');
axis([0 Times(end) 300 10000]);
xlabel('Time (sec)','FontSize',16,'FontWeight','bold');
% set(gca,'ytick',[1000 5000 10000]);
% set(gca,'xtick',[]);
% [rawsong,Fs] = soundin_copy(pathname,filename,'obs0r');
% plot(time,rawsong);
% hold on;
% y_axis = ones(length(spike_times),1) * (max(rawsong) + 100);
% axis tight;
hold on;
set(gca,'Visible','off');

axes('position',[0.15 0.79 0.35 0.15]);
if (length(spike_times) > 0)
    isi = diff(spike_times);
    edges = 0:0.001:0.1;
    isi_hist = histc(isi,edges);
    bar(edges,isi_hist,'histc');
    axis tight;
    set(gca,'FontSize',14,'FontWeight','bold');
end

axes('position',[0.6 0.79 0.3 0.15]);
hold on;
for i = 1:size(waveforms,1);
    plot((1:1:32)/32,waveforms(i,:));
end
axis tight;
set(gca,'FontSize',18);

disp([num2str(length(spike_times)) ' spikes found']);

RangeSpikeTimes = find((spike_times >= Range(1)) & (spike_times <= Range(2)));
disp(['Firing rate is ', num2str(length(RangeSpikeTimes)/(Range(2) - Range(1))), ' Hz']);