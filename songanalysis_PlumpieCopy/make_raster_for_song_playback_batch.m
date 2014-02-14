function [] = make_raster_for_song_playback_batch(FileName,FileExt,FileNos,ChannelNo,bin_size,upper_threshold,lower_threshold)

% The thresholded spikes figure and waveform characteristics
%===============================================================
DirectoryName = pwd;
DirectoryName(end+1) = '/';

figure(1);
RawTrace = axes('position',[0.1 0.75 0.8 0.2]);
ThresholdedSpikes = axes('position',[0.1 0.525 0.8 0.2]);
SoundAmplitude = axes('position',[0.1 0.3 0.8 0.2]);
ISI = axes('position',[0.1 0.025 0.35 0.2]);
Waveforms = axes('position',[0.55 0.025 0.35 0.2]);

RecStartTime = 0;
RecEndTime = 0;
RasterSpikeTimes = [];
  
StartNo = FileNos(1);
EndNo = FileNos(end);

for j = 1:length(FileNos),
    i = FileNos(j);
    if (i < 10);
        DataFileName = strcat(FileName,'.00',num2str(i),FileExt);
    else
        if (i < 100);
            DataFileName = strcat(FileName,'.0',num2str(i),FileExt);
        else
            DataFileName = strcat(FileName,'.',num2str(i),FileExt);
        end
    end
    
    dot_index = find(DataFileName == '.');
    dot_index = dot_index(end);
    RecFileName = strcat(DataFileName(1:dot_index),'rec');

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
    RecEndTime = RecStartTime + str2num(tline((end - 7):(end - 3)))/1000;
    RecordLength = str2num(tline((end - 7):(end - 3)))/1000;
    fclose(fid);

    disp(['Stimulation start time in each file is at ',num2str(StimTime),' seconds']);
    disp(['File end time in each file is at ',num2str(RecEndTime),' seconds']);
    
    channel_string = strcat('obs',num2str(ChannelNo),'r');
    pathname = pwd;
    pathname = strcat(pathname,'/');

    [rawdata,Fs] = soundin_copy(pathname,DataFileName,channel_string);
    rawdata = rawdata * 500/32768;
    
    Time = RecStartTime:1/Fs:RecEndTime;
    Time(end) = [];
    
    threshold_observer_spike_data(DirectoryName,DataFileName,ChannelNo,upper_threshold,lower_threshold);
    SpikeFileName = [DataFileName,'.channo',num2str(ChannelNo),'.spiketimes'];
    TempSpikeTimes = load(SpikeFileName);
    [TempSpikeAmplitudes,TempSpikeWaveforms] = GetSpikeAmplitudes(DirectoryName,DataFileName,['obs',num2str(ChannelNo),'r'],TempSpikeTimes);
    SpikeAmplitudes = [TempSpikeAmplitudes];
    SpikeWaveforms = [TempSpikeWaveforms];
    
    TempSpikeTimes = TempSpikeTimes + RecStartTime;
    SpikeTimes = [TempSpikeTimes];
    SpikeWaveformsTimes = zeros(size(SpikeWaveforms));
    SpikeWaveformsTimes(:,9) = TempSpikeTimes;
    
    for ColumnNo = 1:8,
        SpikeWaveformsTimes(:,ColumnNo) = SpikeWaveformsTimes(:,9) - (9 - ColumnNo)/32000;
    end

    for ColumnNo = 10:32,
        SpikeWaveformsTimes(:,ColumnNo) = SpikeWaveformsTimes(:,9) + ColumnNo/32000;
    end
    
    waveforms = load('test');
    
%     figure(1);
%     
%     axes(RawTrace);
%     plot(Time,rawdata);
%     axis tight;
%     hold on;
%     ylabel('Amplitude','FontSize',14,'FontWeight','bold');
%     set(gca,'xtick',[]);
%     
%     axes(ThresholdedSpikes);
%     plot(Time,rawdata,'b');
%     hold on;
% %     waveforms(:,1) = waveforms(:,1) + RecStartTime;
% %     plot(waveforms(:,1),waveforms(:,2),'r');
%     for PlotNo = 1:size(SpikeWaveforms,1),
%         plot(SpikeWaveformsTimes(PlotNo,:),SpikeWaveforms(PlotNo,:),'r');
%     end
%     axis tight;
%     ylabel('Amplitude','FontSize',14,'FontWeight','bold');
%     set(gca,'xtick',[]);
% 
%     axes(SoundAmplitude);
%     [rawsong,Fs] = soundin_copy(pathname,DataFileName,'obs0r');
%     plot(Time,rawsong);
%     hold on;
%     axis tight;
%     xlabel('Time (sec)','FontSize',14,'FontWeight','bold');
%     ylabel('Amplitude','FontSize',14,'FontWeight','bold');
% 
%     axes(ISI);
%     isi = diff(SpikeTimes);
%     edges = 0:0.001:0.1;
%     isi_hist = histc(isi,edges);
%     bar(edges,isi_hist,'histc');
%     axis tight;
% 
%     axes(Waveforms);
%     hold on;
%     if (length(waveforms) > 0)
%         waveforms = waveforms(:,2);
%         temp = reshape(waveforms,32,length(waveforms)/32);
%         temp = reshape(temp',length(waveforms)/32,1,32);
%     
%         for waveform_index = 1:size(temp,1);
%             plot(linspace(0,0.1,32),temp(waveform_index,:));
%         end
%         axis tight;
%     end

    disp([num2str(length(SpikeTimes)) ' spikes found']);

    TrialOnset = RecStartTime + StimTime;
    RecStartTime = RecEndTime;
    
    TrialSpikeTimes = SpikeTimes - TrialOnset;
    y_value = ones(length(TrialSpikeTimes),1) * (j/5);
    RasterSpikeTimes = [RasterSpikeTimes; [TrialSpikeTimes y_value]];
end

figure;
set(gcf,'Color','w');
axes('position',[0.1 0.4 0.8 0.525]);
edges = -(StimTime):bin_size:(RecordLength - StimTime);
pst = histc(RasterSpikeTimes(:,1),edges);

NoOfTrials = EndNo - StartNo + 1;

pst = pst/(bin_size * NoOfTrials);
% pst = pst/10;
bar(edges,pst,'histc');
if (max(pst) > 10)
    RasterSpikeTimes(:,2) = RasterSpikeTimes(:,2) * 5;
end

RasterSpikeTimes(:,2) = RasterSpikeTimes(:,2) + max(pst) + 5;
hold on;
plot(RasterSpikeTimes(:,1),RasterSpikeTimes(:,2),'w+','MarkerSize',2);
marker_string = repmat('|',size(RasterSpikeTimes,1),1);
text(RasterSpikeTimes(:,1),RasterSpikeTimes(:,2),marker_string,'HorizontalAlignment','center','VerticalAlignment','middle','FontSize',4);
axis([-(StimTime) (RecordLength - StimTime) 0 (max(RasterSpikeTimes(:,2)) + 0.2)]);
%set(gca,'ytick',[]);
ylabel('Firing Rate (Hz)','FontSize',16,'FontWeight','bold');
set(gca,'FontSize',14,'FontWeight','bold')
channel_string = strcat('obs0r');
pathname = pwd;
pathname = strcat(pathname,'/');
[RawData,Fs] = soundin_copy(pathname,DataFileName,channel_string);

axes('position',[0.1 0.05 0.8 0.25]);
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
t_min = time(1); %convert to ms
t_max = time(length(time)); %convert to ms
%adjust time axis for spectrogram offset (1/2 window duration in ms)
t_min = t_min + 0.5*nfft*1/Fs;  
t_max = t_max + 0.5*nfft*1/Fs;  
t_min = t_min;  
t_max = t_max;  
time_spect = [t_min, t_max];                
time_spect = time_spect - StimTime;

disp_idx_spect(idx_spect, time_spect, freq_spect, -55, ...
	20, 1.5, 'gray', 'classic');
axis([t_min t_max 300 10000]);
% set(gca,'ytick',[]);
% set(gca,'xtick',[]);
set(gca,'Visible','off');
xlabel('Time (sec)','FontSize',14,'FontWeight','bold');