function [] = make_raster_for_song_playback(SpikeFile,RawDataFile,ChannelNo,bin_size,DataFile)

SpikeTimes = load(SpikeFile);

channel_string = strcat('obs',num2str(ChannelNo),'r');
pathname = pwd;
pathname = strcat(pathname,'/');
[RawData,Fs] = soundin_copy(pathname,RawDataFile,channel_string);
Times = 0:1/Fs:length(RawData)/Fs;
Times(end) = [];
clear RawData;

dot_index = find(RawDataFile == '.');
dot_index = dot_index(end);
RecFileName = strcat(RawDataFile(1:dot_index),'rec');

fid = fopen(RecFileName,'r');
index = 0;
while (index < 5)
    if (feof(fid))
        break;
    end
    tline = fgetl(fid);
    if (strfind(tline,'stim') > 0)
        index = 5;
    end
end

if (index ~= 5)
    fclose(fid);
    disp('There appears to be no song stimulus in this file');
    return;
end

StimTime =  str2num(tline((end - 7):(end - 3)))/1000;
tline = fgetl(fid);
RecEndTime = str2num(tline((end - 7):(end - 3)))/1000;
fclose(fid);

disp(['Stimulation start time in each file is at ',num2str(StimTime),' seconds']);
disp(['File end time in each file is at ',num2str(RecEndTime),' seconds']);

NoOfTrials = length(Times)/(RecEndTime * Fs);

for i = 1:NoOfTrials,
    TrialOnsets(i) = ((i-1) * (RecEndTime)) + StimTime;
end

figure(3);
RasterSpikeTimes = [];
for i = 1:NoOfTrials,
    figure(3);
    TrialSpikeTimes = SpikeTimes(find((SpikeTimes > (TrialOnsets(i) - 0.5)) & (SpikeTimes < (TrialOnsets(i) + 1.5)))) - TrialOnsets(i);
    y_value = ones(length(TrialSpikeTimes),1) * i;
    RasterSpikeTimes = [RasterSpikeTimes; [TrialSpikeTimes y_value]];
end

axes('position',[0.1 0.25 0.8 0.7]);
edges = -0.5:bin_size:1.5;
pst = histc(RasterSpikeTimes(:,1),edges);
pst = pst/(bin_size * NoOfTrials);
%pst = pst/10;
bar(edges,pst,'histc');
RasterSpikeTimes(:,2) = RasterSpikeTimes(:,2) + max(pst) + 5;
hold on;
plot(RasterSpikeTimes(:,1),RasterSpikeTimes(:,2),'w+','MarkerSize',2);
marker_string = repmat('|',size(RasterSpikeTimes,1),1);
text(RasterSpikeTimes(:,1),RasterSpikeTimes(:,2),marker_string,'HorizontalAlignment','center','VerticalAlignment','middle','FontSize',4);
axis([-0.5 1.5 0 (max(RasterSpikeTimes(:,2)) + 0.2)]);
%set(gca,'ytick',[]);
ylabel('Trial No','FontSize',14,'FontWeight','bold');

channel_string = strcat('obs0r');
pathname = pwd;
pathname = strcat(pathname,'/');
[RawData,Fs] = soundin_copy(pathname,DataFile,channel_string);
Times = 0:1/Fs:length(RawData)/Fs;
Times(end) = [];
Indices = (find((Times >= 1.5) & (Times <= 3.5)));

axes('position',[0.1 0.05 0.8 0.15]);
%plot((Times(Indices) - Notes.onsets(syllable(median_motif,1))),RawData(Indices),'b');
filtsong = bandpass(RawData(Indices),Fs,300,10000);

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
time = time - 0.5;
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
axis([t_min t_max 300 10000]);
set(gca,'ytick',[1000 5000 10000]);
set(gca,'xtick',[]);
xlabel('Time (sec)','FontSize',14,'FontWeight','bold');