function [] = threshold_observer_spike_data(pathname,filename,channel_string,upper_threshold,lower_threshold)

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


[spike_times] = find_spikes(rawdata,lower_threshold,upper_threshold,16,0);

T = spike_times/3.2;

waveforms = load('test');
if (length(waveforms) > 0)
    waveforms = waveforms(:,2);
    
    temp = reshape(waveforms,32,length(waveforms)/32);
    temp = reshape(temp',length(waveforms)/32,1,32);
    wv = temp;
    wv(:,2,:) = wv(:,1,:);
    wv(:,3,:) = wv(:,1,:);
    wv(:,4,:) = wv(:,1,:);
else
    wv = [];
end

file_index = strfind(filename,'.');

if (filename(file_index(end) + 1) == 'b')
    output_filename = strcat(pathname,filename(1:file_index(end)),'BBin.mat');
else
    output_filename = strcat(pathname,filename(1:file_index(end)),'CBin.mat');
end

save(output_filename,'T','wv');

output_filename = strcat(pathname,filename,'.channo',channel_string,'.spiketimes');

spike_times = spike_times/(Fs);
save(output_filename,'spike_times','-ASCII');

time = 0:1/32000:(length(rawdata)/32000);time(end) = [];

waveforms = load('test');

figure;
set(gcf,'Color','w');

axes('position',[0.15 0.565 0.8 0.2]);
plot(time,rawdata,'b');
axis tight;
hold on;
if (length(waveforms) > 0)
    plot(waveforms(:,1),waveforms(:,2),'r');
end
axis tight;
ylabel('Amplitude (\muV)','FontSize',18);
set(gca,'FontSize',18);
set(gca,'xtick',[]);
set(gca,'ytick',[]);

axes('position',[0.15 0.34 0.8 0.2]);
plot(time,rawdata);
axis tight;
set(gca,'FontSize',18);
xlabel('Time (sec)','FontSize',18);

axes('position',[0.15 0.79 0.8 0.2]);
set(gca,'FontSize',18);
channel_string = strcat('obs0r');
pathname = pwd;
pathname = strcat(pathname,'/');
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
[spect, freq, time] = specgram(filtsong, nfft, Fs, spect_win, noverlap);
idx_spect=scale_spect(spect);  %calculate index array for spectrogram
f_min = freq(1);
f_max = freq(length(freq));
freq_spect = [f_min, f_max];
time = time - 0.2;
t_min = time(1); %convert to ms
t_max = time(length(time)); %convert to ms
%adjust time axis for spectrogram offset (1/2 window duration in ms)
%t_min = t_min + 0.5*nfft*1000/Fs;  
%t_max = t_max + 0.5*nfft*1000/Fs;  
t_min = t_min;  
t_max = t_max;  
time_spect = [t_min, t_max];                

disp_idx_spect(idx_spect, time_spect, freq_spect, -55, ...
	0, 2, 'gray', 'classic');
axis([-0.2 t_max 300 10000]);
% set(gca,'ytick',[1000 5000 10000]);
% set(gca,'xtick',[]);
% [rawsong,Fs] = soundin_copy(pathname,filename,'obs0r');
% plot(time,rawsong);
% hold on;
% y_axis = ones(length(spike_times),1) * (max(rawsong) + 100);
% axis tight;
hold on;
set(gca,'Visible','off');

axes('position',[0.15 0.065 0.35 0.2]);
if (length(spike_times) > 0)
    isi = diff(spike_times);
    edges = 0:0.001:0.1;
    isi_hist = histc(isi,edges);
    bar(edges,isi_hist,'histc');
    axis tight;
    set(gca,'FontSize',18);
end

axes('position',[0.6 0.065 0.3 0.2]);
hold on;
for i = 1:size(temp,1);
    plot(linspace(0,0.1,32),temp(i,:));
end
axis tight;
set(gca,'FontSize',18);

disp([num2str(length(spike_times)) ' spikes found']);
