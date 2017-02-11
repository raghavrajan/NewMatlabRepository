function [] = make_raster_for_song_playback_batch_for_windows(FileName,FileExt,FileNos,ChannelNo,bin_size,upper_threshold,lower_threshold,site,string,varargin)

% The thresholded spikes figure and waveform characteristics
%===============================================================

if (length(FileNos) == 0)
   StartNum = cell2mat(varargin(1));
   EndNum = cell2mat(varargin(2));
   [FileTypes] = filesort(FileName,'.rec',StartNum, EndNum,string);
   TrialTypes = fieldnames(FileTypes);
else
   TrialTypes = {'Playback'};
end

ISITimes = [];
for TrialIndices = 1:length(TrialTypes),
    if (exist('FileTypes'))
        FileNos = FileTypes.(cell2mat(TrialTypes(TrialIndices)));
    end
    
    figure;
    set(gcf,'Visible','off');
    RawTrace = axes('position',[0.1 0.75 0.8 0.2]);
    ThresholdedSpikes = axes('position',[0.1 0.525 0.8 0.2]);
    SoundAmplitude = axes('position',[0.1 0.3 0.8 0.2]);
    ISI = axes('position',[0.1 0.025 0.35 0.2]);
    Waveforms = axes('position',[0.55 0.025 0.35 0.2]);

    RecStartTime = 0;
    RecEndTime = 0;
    RasterSpikeTimes = [];

    for i = 1:length(FileNos),
        if (FileNos(i) < 10);
            DataFileName = strcat(FileName,'.00',num2str(FileNos(i)),FileExt);
        else
            if (FileNos(i) < 100);
                DataFileName = strcat(FileName,'.0',num2str(FileNos(i)),FileExt);
            else
                DataFileName = strcat(FileName,'.',num2str(FileNos(i)),FileExt);
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
        
        disp(RecFileName);
        disp(['Stimulation start time in each file is at ',num2str(StimTime),' seconds']);
        disp(['File end time in each file is at ',num2str(RecEndTime),' seconds']);

        channel_string = strcat('obs',num2str(ChannelNo),'r');
        pathname = pwd;
        pathname = strcat(pathname,'/');

        [rawdata,Fs] = soundin_copy(pathname,DataFileName,channel_string);
        rawdata = rawdata * 500/32768;
        
        Time = RecStartTime:1/Fs:RecEndTime;
        Time(end) = [];

        [SpikeTimes] = find_spikes(rawdata,lower_threshold,upper_threshold,16,0);
        if (length(SpikeTimes) > 0)
            waveforms = load('test');
            SpikeTimes = SpikeTimes/(Fs);
            SpikeTimes = SpikeTimes + RecStartTime;
        else
            continue;
        end

        axes(RawTrace);
        plot(Time,rawdata);
        axis tight;
        hold on;
        ylabel('Amplitude','FontSize',14,'FontWeight','bold');
        set(gca,'xtick',[]);

        axes(ThresholdedSpikes);
        plot(Time,rawdata,'b');
        hold on;
        waveforms(:,1) = waveforms(:,1) + RecStartTime;
        plot(waveforms(:,1),waveforms(:,2),'r');
        axis tight;
        ylabel('Amplitude','FontSize',14,'FontWeight','bold');
        set(gca,'xtick',[]);

        axes(SoundAmplitude);
        [rawsong,Fs] = soundin_copy(pathname,DataFileName,'obs0r');
        rawsong = rawsong * 1/32768;
        plot(Time,rawsong);
        hold on;
        axis tight;
        xlabel('Time (sec)','FontSize',14,'FontWeight','bold');
        ylabel('Amplitude','FontSize',14,'FontWeight','bold');

        ISITimes = [ISITimes; diff(SpikeTimes)];
       
        axes(Waveforms);
        hold on;
        if (length(waveforms) > 0)
            waveforms = waveforms(:,2);
            temp = reshape(waveforms,32,length(waveforms)/32);
            temp = reshape(temp',length(waveforms)/32,1,32);

            for waveform_index = 1:size(temp,1);
                plot(linspace(0,0.1,32),temp(waveform_index,:));
            end
            axis tight;
        end

        disp([num2str(length(SpikeTimes)) ' spikes found']);

        TrialOnset = RecStartTime + StimTime;
        RecStartTime = RecEndTime;

        TrialSpikeTimes = SpikeTimes - TrialOnset;
        y_value = ones(length(TrialSpikeTimes),1) * i;
        RasterSpikeTimes = [RasterSpikeTimes; [TrialSpikeTimes y_value]];
    end

    set(gcf,'Visible','on');
    axes(ISI);
    edges = 0:0.001:0.1;
    if (length(ISITimes) > 0)
        isi_hist = histc(ISITimes,edges);
        bar(edges,isi_hist,'histc');
        axis tight;
    end

    axes(RawTrace);
    title(strcat(site,'-','obs',num2str(ChannelNo),'-',num2str(StartNum),'-',num2str(EndNum),'-', cell2mat(TrialTypes(TrialIndices))),'FontWeight','bold','FontSize',14);
    figure;
    set(gcf,'Color','w');
    axes('position',[0.15 0.15 0.8 0.6]);
    set(gca,'FontSize',18);
    edges = -(StimTime):bin_size:(RecordLength - StimTime);
    if (length(RasterSpikeTimes) > 0)
        pst = histc(RasterSpikeTimes(:,1),edges);

        NoOfTrials = length(FileNos);

        pst = pst/(bin_size * NoOfTrials);
        if (max(pst) > 50)
            RasterSpikeTimes(:,2) = RasterSpikeTimes(:,2) * 10;
        else
            if (max(pst) > 10)
                RasterSpikeTimes(:,2) = RasterSpikeTimes(:,2) * 5;
            end
        end
        tempbar = bar(edges,pst,'histc');
        set(tempbar,'EdgeColor','k');
        set(tempbar,'FaceColor','w');
        RasterSpikeTimes(:,2) = RasterSpikeTimes(:,2) + max(pst) + 5;
        hold on;
        plot(RasterSpikeTimes(:,1),RasterSpikeTimes(:,2),'w+','MarkerSize',2);
        marker_string = repmat('|',size(RasterSpikeTimes,1),1);
        text(RasterSpikeTimes(:,1),RasterSpikeTimes(:,2),marker_string,'HorizontalAlignment','center','VerticalAlignment','middle','FontSize',5,'FontName','FixedWidth','FontUnits','pixels');
        axis([-(StimTime) (RecordLength - StimTime) 0 (max(RasterSpikeTimes(:,2)) + 0.2)]);
        %set(gca,'ytick',[]);
    end
    ylabel('Firing rate (Hz)','FontSize',18);
    xlabel('Time (sec)','FontSize',18);
    
%     title(strcat(site,'-','obs',num2str(ChannelNo),'-',num2str(StartNum),'-',num2str(EndNum),'-',cell2mat(TrialTypes(TrialIndices))),'FontWeight','bold','FontSize',14);
    
    channel_string = strcat('obs0r');
    pathname = pwd;
    pathname = strcat(pathname,'/');
    [RawData,Fs] = soundin_copy(pathname,DataFileName,channel_string);

    axes('position',[0.15 0.8 0.8 0.15]);
    filtsong = bandpass(RawData,Fs,300,8000);

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
        5, 1.5, 'gray', 'classic');
    axis([t_min t_max 300 8000]);
    set(gca,'ytick',[]);
    set(gca,'xtick',[]);
    set(gca,'Visible','off');
%     xlabel(strcat('Time (sec)','___','Threshold_',num2str(upper_threshold),'/',num2str(lower_threshold)),'FontSize',12,'FontWeight','bold');
end

figs = findobj('Type','figure');
figs = sort(figs);

for FigNo = 1:length(figs),
    if (length(varargin) ~= 0)
        if (FigNo < 3)
            OutputFileName = strcat(site,'_',num2str(StartNum),'_',num2str(EndNum),'_',cell2mat(TrialTypes(1)));
        else
            OutputFileName = strcat(site,'_',num2str(StartNum),'_',num2str(EndNum),'_',cell2mat(TrialTypes(2)));            
        end
    else
        OutputFileName = strcat(site,'_',num2str(FileNos(1)),'_',num2str(FileNos(end)));
    end
    saveas(figure(figs(FigNo)),strcat(OutputFileName,'fig',num2str(FigNo),'.fig'));
end
