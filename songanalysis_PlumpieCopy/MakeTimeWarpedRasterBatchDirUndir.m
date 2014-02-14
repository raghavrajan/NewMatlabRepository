function [SpikeTrain,median_syllable_starting,median_gap_starting] = MakeTimeWarpedRasterBatchDirUndir(DirectoryName,FileName,FileExt,FileNos,ChannelString,BinSize,Motif,latency,SpikeSortOrThreshold)

DirSpikeTimes = [];
DirSpikeAmplitudes = [];
DirSpikeWaveforms = [];
DirActualSpikeWaveforms = [];
DirNoteOnsets = [];
DirNoteOffsets = [];
DirNoteLabels = [];
DirRecEndTime = 0;
DirRawData = [];
DirRawSongData = [];
DirUnwarpedRaster = [];
DirMotifRelatedFiringRate = [];

UndirSpikeTimes = [];
UndirSpikeAmplitudes = [];
UndirSpikeWaveforms = [];
UndirActualSpikeWaveforms = [];
UndirNoteOnsets = [];
UndirNoteOffsets = [];
UndirNoteLabels = [];
UndirRecEndTime = 0;
UndirRawData = [];
UndirRawSongData = [];
UndirUnwarpedRaster = [];
UndirMotifRelatedFiringRate = [];

if (strfind(SpikeSortOrThreshold,'spikesort'))
    MaxClusters = input('Type the maximum number of clusters in the sorted data excluding the outliers: ');
    ClusterNos = input('Type the cluster nos in the sorted data that you want analysed: ');
    IncludeOutliers = input('Type yes if you want the outliers to be included, else type no: ');
end

if (strfind(SpikeSortOrThreshold,'threshold'))
    HighThresh = input('Type the higher threshold in uV: ');
    LowThresh = input('Type the lower threshold in uV: ');
end

cd(DirectoryName);

% First get all the data from the Directed bouts

for i = 1:length(FileNos.Directed),
    if (FileNos.Directed(i) < 10)
        DataFile = [FileName,'.00',num2str(FileNos.Directed(i))];
    else
        if (FileNos.Directed(i) < 100)
            DataFile = [FileName,'.0',num2str(FileNos.Directed(i))];            
        else
            DataFile = [FileName,'.',num2str(FileNos.Directed(i))];
        end
    end
    if (FileExt(1) == '.')
        DataFileName = [DataFile,FileExt];
    else
        DataFileName = [DataFile,'.',FileExt];
    end
    
    if (ispc)
        if (DirectoryName(end) ~= '\')
            DirectoryName = [DirectoryName,'\'];
        end
    else
        if (isunix)
            if (DirectoryName(end) ~= '/')
                DirectoryName = [DirectoryName,'/'];
            end
        end 
    end
    
    RecFileName = strcat(DataFile,'.rec');

    fid = fopen(RecFileName,'r');
    index = 0;
    while (index < 5)
        if (feof(fid))
            break;
        end
        tline = fgetl(fid);
        if (strfind(tline,'rec end') > 0)
            index = 5;
        end
    end

    StartIndex = find(tline == '=');
    EndIndex = strfind(tline,'ms');
    RecTime = str2double(tline((StartIndex + 1):(EndIndex - 1)));
    RecTime = RecTime/1000;
    
    fclose(fid);
    
    if (strfind(SpikeSortOrThreshold,'spikesort'))
        SpikeFileName = [DataFile,'.spk'];
        TempSpikeTimes = load(SpikeFileName);
        for ClusterIndex = 1:length(ClusterNos),
            TempSpikes = TempSpikeTimes(find((TempSpikeTimes(:,1) == ClusterNos(ClusterIndex))),2)/1000;
            [TempSpikeAmplitudes,TempSpikeWaveforms] = GetSpikeAmplitudes(DirectoryName,DataFileName,ChannelString,TempSpikes);
            DirSpikeAmplitudes = [DirSpikeAmplitudes; TempSpikeAmplitudes];
            DirSpikeWaveforms = [DirSpikeWaveforms; TempSpikeWaveforms];
            TempSpikes = TempSpikes + DirRecEndTime;
            DirSpikeTimes = [SpikeTimes; TempSpikes];
            if (strfind(IncludeOutliers,'yes'))
                TempSpikes = TempSpikeTimes(find((TempSpikeTimes(:,1) == (MaxClusters + ClusterNos(ClusterIndex) + 1))),2)/1000;
                [TempSpikeAmplitudes, TempSpikeWaveforms] = GetSpikeAmplitudes(DirectoryName,DataFileName,ChannelString,TempSpikes);
                DirSpikeAmplitudes = [DirSpikeAmplitudes; TempSpikeAmplitudes];
                DirSpikeWaveforms = [DirSpikeWaveforms; TempSpikeWaveforms];
                TempSpikes = TempSpikes + DirRecEndTime;
                DirSpikeTimes = [DirSpikeTimes; TempSpikes];
            end
        end
    else
        threshold_observer_spike_data(DirectoryName,DataFileName,ChannelString,HighThresh,LowThresh);
        SpikeFileName = [DataFileName,'.channo',ChannelString,'.spiketimes'];
        TempSpikeTimes = load(SpikeFileName);
        [TempSpikeAmplitudes,TempSpikeWaveforms] = GetSpikeAmplitudes(DirectoryName,DataFileName,ChannelString,TempSpikeTimes);
        DirSpikeAmplitudes = [DirSpikeAmplitudes; TempSpikeAmplitudes];
        DirSpikeWaveforms = [DirSpikeWaveforms; TempSpikeWaveforms];
        TempSpikeTimes = TempSpikeTimes + DirRecEndTime;
        DirSpikeTimes = [DirSpikeTimes; TempSpikeTimes];
        
    end
    [DirSpikeTimes,DirSpikeTimeIndices] = sort(DirSpikeTimes);
    if ((length(DirSpikeAmplitudes) == length(DirSpikeTimes)))
        DirSpikeAmplitudes = DirSpikeAmplitudes(DirSpikeTimeIndices);
        DirSpikeWaveforms = DirSpikeWaveforms(DirSpikeTimeIndices,:);
    end
    
    NoteFileName = [DataFileName,'.not.mat'];
    Notes = load(NoteFileName);

    Notes.onsets = Notes.onsets/1000;
    Notes.offsets = Notes.offsets/1000;
    
    Notes.onsets = Notes.onsets + DirRecEndTime;
    Notes.offsets = Notes.offsets + DirRecEndTime;
    
    DirNoteOnsets = [DirNoteOnsets; [Notes.onsets]];
    DirNoteOffsets = [DirNoteOffsets; [Notes.offsets]];
    DirNoteLabels = [DirNoteLabels [Notes.labels]];

    pathname = pwd;
    pathname = strcat(pathname,'/');
    [TempRawData,Fs] = soundin_copy(pathname,DataFileName,'obs0r');
    DirRawSongData = [DirRawSongData; TempRawData];
    
    pathname = pwd;
    pathname = strcat(pathname,'/');
    [TempRawData,Fs] = soundin_copy(pathname,DataFileName,ChannelString);
    DirRawData = [DirRawData; TempRawData];
    
    DirRecEndTime = DirRecEndTime + RecTime;
end

% Now, get all the data from the Undirected bouts

for i = 1:length(FileNos.Undirected),
    if (FileNos.Undirected(i) < 10)
        DataFile = [FileName,'.00',num2str(FileNos.Undirected(i))];
    else
        if (FileNos.Undirected(i) < 100)
            DataFile = [FileName,'.0',num2str(FileNos.Undirected(i))];            
        else
            DataFile = [FileName,'.',num2str(FileNos.Undirected(i))];
        end
    end
    if (FileExt(1) == '.')
        DataFileName = [DataFile,FileExt];
    else
        DataFileName = [DataFile,'.',FileExt];
    end
    
    if (ispc)
        if (DirectoryName(end) ~= '\')
            DirectoryName = [DirectoryName,'\'];
        end
    else
        if (isunix)
            if (DirectoryName(end) ~= '/')
                DirectoryName = [DirectoryName,'/'];
            end
        end 
    end
    
    RecFileName = strcat(DataFile,'.rec');

    fid = fopen(RecFileName,'r');
    index = 0;
    while (index < 5)
        if (feof(fid))
            break;
        end
        tline = fgetl(fid);
        if (strfind(tline,'rec end') > 0)
            index = 5;
        end
    end

    StartIndex = find(tline == '=');
    EndIndex = strfind(tline,'ms');
    RecTime = str2double(tline((StartIndex + 1):(EndIndex - 1)));
    RecTime = RecTime/1000;
    
    fclose(fid);
    
    if (strfind(SpikeSortOrThreshold,'spikesort'))
        SpikeFileName = [DataFile,'.spk'];
        TempSpikeTimes = load(SpikeFileName);
        for ClusterIndex = 1:length(ClusterNos),
            TempSpikes = TempSpikeTimes(find((TempSpikeTimes(:,1) == ClusterNos(ClusterIndex))),2)/1000;
            [TempSpikeAmplitudes,TempSpikeWaveforms] = GetSpikeAmplitudes(DirectoryName,DataFileName,ChannelString,TempSpikes);
            UndirSpikeAmplitudes = [UndirSpikeAmplitudes; TempSpikeAmplitudes];
            UndirSpikeWaveforms = [UndirSpikeWaveforms; TempSpikeWaveforms];
            TempSpikes = TempSpikes + UndirRecEndTime;
            UndirSpikeTimes = [SpikeTimes; TempSpikes];
            if (strfind(IncludeOutliers,'yes'))
                TempSpikes = TempSpikeTimes(find((TempSpikeTimes(:,1) == (MaxClusters + ClusterNos(ClusterIndex) + 1))),2)/1000;
                [TempSpikeAmplitudes, TempSpikeWaveforms] = GetSpikeAmplitudes(DirectoryName,DataFileName,ChannelString,TempSpikes);
                UndirSpikeAmplitudes = [UndirSpikeAmplitudes; TempSpikeAmplitudes];
                UndirSpikeWaveforms = [UndirSpikeWaveforms; TempSpikeWaveforms];
                TempSpikes = TempSpikes + UndirRecEndTime;
                UndirSpikeTimes = [UndirSpikeTimes; TempSpikes];
            end
        end
    else
        threshold_observer_spike_data(DirectoryName,DataFileName,ChannelString,HighThresh,LowThresh);
        SpikeFileName = [DataFileName,'.channo',ChannelString,'.spiketimes'];
        TempSpikeTimes = load(SpikeFileName);
        [TempSpikeAmplitudes,TempSpikeWaveforms] = GetSpikeAmplitudes(DirectoryName,DataFileName,ChannelString,TempSpikeTimes);
        UndirSpikeAmplitudes = [UndirSpikeAmplitudes; TempSpikeAmplitudes];
        UndirSpikeWaveforms = [UndirSpikeWaveforms; TempSpikeWaveforms];
        TempSpikeTimes = TempSpikeTimes + UndirRecEndTime;
        UndirSpikeTimes = [UndirSpikeTimes; TempSpikeTimes];
        
    end
    [UndirSpikeTimes,UndirSpikeTimeIndices] = sort(UndirSpikeTimes);
    if ((length(UndirSpikeAmplitudes) == length(UndirSpikeTimes)))
        UndirSpikeAmplitudes = UndirSpikeAmplitudes(UndirSpikeTimeIndices);
        UndirSpikeWaveforms = UndirSpikeWaveforms(UndirSpikeTimeIndices,:);
    end
    
    NoteFileName = [DataFileName,'.not.mat'];
    Notes = load(NoteFileName);

    Notes.onsets = Notes.onsets/1000;
    Notes.offsets = Notes.offsets/1000;
    
    Notes.onsets = Notes.onsets + UndirRecEndTime;
    Notes.offsets = Notes.offsets + UndirRecEndTime;
    
    UndirNoteOnsets = [UndirNoteOnsets; [Notes.onsets]];
    UndirNoteOffsets = [UndirNoteOffsets; [Notes.offsets]];
    UndirNoteLabels = [UndirNoteLabels [Notes.labels]];

    pathname = pwd;
    pathname = strcat(pathname,'/');
    [TempRawData,Fs] = soundin_copy(pathname,DataFileName,'obs0r');
    UndirRawSongData = [UndirRawSongData; TempRawData];
    
    pathname = pwd;
    pathname = strcat(pathname,'/');
    [TempRawData,Fs] = soundin_copy(pathname,DataFileName,ChannelString);
    UndirRawData = [UndirRawData; TempRawData];
    
    UndirRecEndTime = UndirRecEndTime + RecTime;
end

% =========================================================================
% Finished collecting data
disp('Finished collecting data');
% =========================================================================

DirRawData = DirRawData * 5000/32768;
DirRawSongData = DirRawSongData * 1/32768;

DirTimes = 0:1/Fs:length(DirRawSongData)/Fs;
DirTimes(end) = [];

DirSpikeAmplitudes = DirSpikeAmplitudes * 5/32768 * 1000;
DirSpikeWaveforms = DirSpikeWaveforms * 5/32768 * 1000;

UndirRawData = UndirRawData * 5000/32768;
UndirRawSongData = UndirRawSongData * 1/32768;

UndirTimes = 0:1/Fs:length(UndirRawSongData)/Fs;
UndirTimes(end) = [];

UndirSpikeAmplitudes = UndirSpikeAmplitudes * 5/32768 * 1000;
UndirSpikeWaveforms = UndirSpikeWaveforms * 5/32768 * 1000;

% IntroNotes = find(NoteLabels == 'i');
% Temp = find(diff(IntroNotes) > 1);
% Temp = Temp + 1;
% BoutBeginnings = [1 IntroNotes(Temp)];

clear Notes;

if (exist('latency','var'))
    DirSpikeTimes = DirSpikeTimes + latency;
    UndirSpikeTimes = UndirSpikeTimes + latency;
end

% SpikeWaveformFigure = figure;
% hold on;

RasterFigure = figure;
set(gcf,'Color','w');
RasterPlot = axes('position',[0.1 0.15 0.8 0.6]);
set(gca,'FontSize',18);
MotifPlot =  axes('position',[0.1 0.775 0.8 0.2]);
set(gca,'FontSize',18);
MaxRaster = 0;

RawDataFigure = figure;
set(gcf,'Color','w');
RawDataPlot = axes('position',[0.15 0.1 0.8 0.8]);
hold on;
set(gca,'FontSize',18);

SpikeAmplitudeFigure = figure;
set(gcf,'Color','w');
AmplitudePlot = axes('position',[0.1 0.15 0.8 0.6]);
set(gca,'FontSize',18);
SAMotifPlot =  axes('position',[0.1 0.775 0.8 0.2]);
set(gca,'FontSize',18);
motif_count = 1;
MaxRawData = 0;

BoutOnsetOffsetFigure = figure;
set(gcf,'Color','w');
BoutOnsetPlot = axes('position',[0.1 0.15 0.8 0.3]);
BoutOffsetPlot =  axes('position',[0.1 0.65 0.8 0.3]);

CVFigure = figure;
set(gcf,'Color','w');
CVPlot = axes('position',[0.1 0.15 0.8 0.2]);
CVRasterPlot = axes('position',[0.1 0.375 0.8 0.425]);
CVMotifPlot = axes('position',[0.1 0.825 0.8 0.15]);

for MotifNo = 1:length(Motif),
    if (MotifNo > 1)
        break;
    end
    
    syllable = [];
    syllable_lengths = [];
    gap_lengths = [];
    motif = Motif(1:(end - (MotifNo - 1)));

    Dirmotifs = strfind(DirNoteLabels,motif);
    
    if (length(Dirmotifs) == 0)
        continue;
    end
    
    for i = 1:length(dirmotif),
        Dirsyllable(:,i) = Dirmotifs + (i - 1);
        Dirsyllable_lengths(:,i) = DirNoteOffsets(Dirsyllable(:,i)) - DirNoteOnsets(syllable(:,i));
    end

    for i = 1:(length(motif)-1),
        Dirgap_lengths(:,i) = DirNoteOnsets(syllable(:,i+1)) - DirNoteOffsets(syllable(:,i));
    end

    if (MotifNo == 1)
        motif_durations = NoteOffsets(syllable(:,end)) - NoteOnsets(syllable(:,1));
        [sorted_motif_durations,sorted_indices] = sort(motif_durations);

        if (mod(length(motif_durations),2) == 0)
            median_motif = sorted_indices(length(motif_durations)/2);
        else
            median_motif = sorted_indices((length(motif_durations) + 1)/2);
        end
        
        if (exist('syllable_lengths','var'))
            median_syllable_lengths = syllable_lengths(median_motif,:);
        end
        if (length(gap_lengths) > 0)
            median_gap_lengths = gap_lengths(median_motif,:);
        end
        if (exist('syllable_lengths','var'));
            for i = 1:length(median_syllable_lengths),
                if (i == 1)
                    median_syllable_starting(i) = 0;
                else
                    median_syllable_starting(i) = sum(median_syllable_lengths(1:(i-1))) + sum(median_gap_lengths(1:(i-1)));
                end
            end
        end
        if (length(gap_lengths) > 0)
            for i = 1:length(median_gap_lengths),
                if (i == 1)
                    median_gap_starting(i) = median_syllable_lengths(i);
                else
                    median_gap_starting(i) = sum(median_syllable_lengths(1:i)) + sum(median_gap_lengths(1:(i-1)));
                end
            end
        end
    end
    
    bin_size = BinSize;

%     figure;
%     UnwarpedRaster = [];
    Raster = [];
    Amplitudes = [];
    trial_pst = [];
    for i = 1:size(syllable_lengths,1),
        figure(RawDataFigure);
        axes(RawDataPlot);
        if ((i < 7) && (MotifNo == 1))
            TempRawData = RawData(find((Times > (NoteOnsets(syllable(i,1)) - 0.2)) & (Times < NoteOffsets(syllable(i,end))))) - NoteOnsets(syllable(i,1));
            TempRawData = TempRawData + MaxRawData;
            MaxRawData = MaxRawData + max(TempRawData) - min(TempRawData);
            TempTimes = -0.2:1/32000:(NoteOffsets(syllable(i,end)) - NoteOnsets(syllable(i,1)));
            plot(TempTimes(1:length(TempRawData)),TempRawData,'k');
            axis tight;
            xlabel('Time (sec)','FontSize',18);
            ylabel('Amplitude (\muV)','FontSize',18);
            title(FigTitle);
        end
        
               
        TrialSpikeTimes = SpikeTimes(find((SpikeTimes > (NoteOnsets(syllable(i,1)) - 0.2)) & (SpikeTimes < NoteOnsets(syllable(i,1)))));
        TrialSpikeAmplitudes = SpikeAmplitudes(find((SpikeTimes > (NoteOnsets(syllable(i,1)) - 0.2)) & (SpikeTimes < NoteOnsets(syllable(i,1)))));
        TrialSpikeWaveforms = SpikeWaveforms(find((SpikeTimes > (NoteOnsets(syllable(i,1)) - 0.2)) & (SpikeTimes < NoteOnsets(syllable(i,1)))),:);
        TrialSpikeTimes = TrialSpikeTimes - NoteOnsets(syllable(i,1));
        y_value = ones(length(TrialSpikeTimes),1) * i/5;
        Raster = [Raster; [TrialSpikeTimes y_value]];
%         UnwarpedRaster = [UnwarpedRaster; [TrialSpikeTimes y_value]];
        Amplitudes = [Amplitudes; [TrialSpikeTimes TrialSpikeAmplitudes]];
        ActualSpikeWaveforms = [ActualSpikeWaveforms; TrialSpikeWaveforms];

        syll_count = 0;
        gap_count = 0;
        for j = 1:(2*length(motif) - 1),
            if (mod(j,2) == 1)
                syll_count = syll_count + 1;
                TrialSpikeTimes = SpikeTimes(find((SpikeTimes > NoteOnsets(syllable(i,syll_count))) & (SpikeTimes < NoteOffsets(syllable(i,syll_count)))));
                y_value = ones(length(TrialSpikeTimes),1) * i/5;
                TrialSpikeAmplitudes = SpikeAmplitudes(find((SpikeTimes > NoteOnsets(syllable(i,syll_count))) & (SpikeTimes < NoteOffsets(syllable(i,syll_count)))));
                TrialSpikeWaveforms = SpikeWaveforms(find((SpikeTimes > NoteOnsets(syllable(i,syll_count))) & (SpikeTimes < NoteOffsets(syllable(i,syll_count)))),:);
                TrialSpikeTimes = (TrialSpikeTimes - NoteOnsets(syllable(i,syll_count)))/syllable_lengths(i,syll_count);
                TrialSpikeTimes = (TrialSpikeTimes * median_syllable_lengths(1,syll_count)) + median_syllable_starting(1,syll_count);
                UnwarpedRaster = [UnwarpedRaster; TrialSpikeTimes];
                Raster = [Raster; [TrialSpikeTimes y_value]];
                Amplitudes = [Amplitudes; [TrialSpikeTimes TrialSpikeAmplitudes]];
                ActualSpikeWaveforms = [ActualSpikeWaveforms; TrialSpikeWaveforms];
            else
                gap_count = gap_count + 1;
                TrialSpikeTimes = SpikeTimes(find((SpikeTimes > NoteOffsets(syllable(i,gap_count))) & (SpikeTimes < NoteOnsets(syllable(i,(gap_count + 1))))));
                y_value = ones(length(TrialSpikeTimes),1) * i/5;
                TrialSpikeAmplitudes = SpikeAmplitudes(find((SpikeTimes > NoteOffsets(syllable(i,gap_count))) & (SpikeTimes < NoteOnsets(syllable(i,(gap_count + 1))))));                
                TrialSpikeWaveforms = SpikeWaveforms(find((SpikeTimes > NoteOffsets(syllable(i,gap_count))) & (SpikeTimes < NoteOnsets(syllable(i,(gap_count + 1))))),:);                
                TrialSpikeTimes = (TrialSpikeTimes - NoteOffsets(syllable(i,gap_count)))/gap_lengths(i,gap_count);
                TrialSpikeTimes = (TrialSpikeTimes * median_gap_lengths(1,gap_count)) + median_gap_starting(1,gap_count);
                UnwarpedRaster = [UnwarpedRaster; TrialSpikeTimes];
                Raster = [Raster; [TrialSpikeTimes y_value]];
                Amplitudes = [Amplitudes; [TrialSpikeTimes TrialSpikeAmplitudes]];                
                ActualSpikeWaveforms = [ActualSpikeWaveforms; TrialSpikeWaveforms];                
            end
        end
%         edges = -0.2:bin_size:motif_durations(median_motif);
%         temp_trial_pst = histc(Raster(:,1),edges);
%         temp_trial_pst(end) = [];
%         trial_pst = [trial_pst; temp_trial_pst'];
        MotifRelatedFiringRate(motif_count) = length(UnwarpedRaster)/motif_durations(motif_count);
        UnwarpedRaster = UnwarpedRaster * 1000;
        SpikeTrain{motif_count} = num2cell(UnwarpedRaster');
        motif_count = motif_count + 1;
        UnwarpedRaster = [];
    end
    figure(RasterFigure);
    axes(RasterPlot);
%     edges = -0.2:bin_size:motif_durations(median_motif);
%     pst = histc(Raster(:,1),edges);
%     pst(end) = [];
%     edges(end) = [];
%     pst = pst/size(motif_durations,1);
%     pst = pst/bin_size;
%     pst = pst/100;
    %bar(edges,pst,'histc');
    hold on;
    %Raster(:,2) = Raster(:,2) + max(pst) + 1;
    Raster(:,2) = Raster(:,2) + MaxRaster;
    plot(Raster(:,1),Raster(:,2),'w+');
    marker_string = repmat('|',size(Raster,1),1);
    text(Raster(:,1),Raster(:,2),marker_string,'HorizontalAlignment','center','VerticalAlignment','middle','FontUnits','points','FontSize',2);
%     line(Raster(:,1)',Raster(:,2)','Color','k','LineWidth',1);
%     if (length(Raster) > 0)
%         plot([0 (motif_durations(median_motif) - 0.2)], [(max(Raster(:,2)) + 0.1) (max(Raster(:,2)) + 0.1)],'r');
    axis([-0.2 motif_durations(median_motif) 0 (max(Raster(:,2)) + 0.5)]);
%         MaxRaster = max(Raster(:,2)) + 0.1;
%     else
%         plot([0 (motif_durations(median_motif) - 0.2)], [(MaxRaster + 0.1) (MaxRaster + 0.1)],'r');
%         axis([-0.2 motif_durations(median_motif) 0 (MaxRaster + 0.5)]);
%     end
%     set(gca,'ytick',[0,max(pst)/2,max(pst)]);
    set(gca,'ytick',[]);
    set(gca,'FontSize',16);
    ylabel('Trials','FontSize',18);
    xlabel('Time (sec)','FontSize',18);
    title(FigTitle);
    
    figure(SpikeAmplitudeFigure);
    axes(AmplitudePlot);
    hold on;
    plot(Amplitudes(:,1),Amplitudes(:,2),'k+','MarkerSize',2);
    if (length(Amplitudes) > 1)
        axis([-0.2 motif_durations(median_motif) 0 (max(Amplitudes(:,2)))]);
    end
    set(gca,'FontSize',16);
    ylabel('Amplitude(\muV)','FontSize',18);
    xlabel('Time (sec)','FontSize',18);
    title(FigTitle);
    
    if (MotifNo == 1)
        channel_string = strcat('obs0r');
        pathname = pwd;
        pathname = strcat(pathname,'/');
    %     [RawData,Fs] = soundin_copy(pathname,DataFileName,channel_string);
        Indices = (find((Times >= (NoteOnsets(syllable(median_motif,1)))) & (Times <= NoteOffsets(syllable(median_motif,length(motif))))));

        %plot((Times(Indices) - Notes.onsets(syllable(median_motif,1))),RawData(Indices),'b');
        filtsong = bandpass(RawSongData(Indices),Fs,300,10000);

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
        time = time;
        t_min = time(1); %convert to ms
        t_max = time(length(time)); %convert to ms
        %adjust time axis for spectrogram offset (1/2 window duration in ms)
        %t_min = t_min + 0.5*nfft*1000/Fs;  
        %t_max = t_max + 0.5*nfft*1000/Fs;  
        t_min = t_min;  
        t_max = t_max;  
        time_spect = [t_min, t_max];                
        
        figure(RasterFigure);
        axes(MotifPlot);

        disp_idx_spect(idx_spect, time_spect, freq_spect, -55, ...
            5, 1.5, 'gray', 'classic');
        axis([-0.2 t_max 300 10000]);
        set(gca,'FontSize',16);
        set(gca,'ytick',[]);
        set(gca,'xtick',[]);
        set(gca,'Visible','off');
        
        figure(SpikeAmplitudeFigure);
        axes(SAMotifPlot);

        disp_idx_spect(idx_spect, time_spect, freq_spect, -55, ...
            5, 1.5, 'gray', 'classic');
        axis([-0.2 t_max 300 10000]);
        set(gca,'FontSize',16);        
        set(gca,'ytick',[]);
        set(gca,'xtick',[]);
        set(gca,'Visible','off');        
        
        figure(CVFigure);
        axes(CVMotifPlot);

        disp_idx_spect(idx_spect, time_spect, freq_spect, -55, ...
            5, 1.5, 'gray', 'classic');
        axis([0 t_max 300 10000]);
        set(gca,'FontSize',16);        
        set(gca,'ytick',[]);
        set(gca,'xtick',[]);
        set(gca,'Visible','off');        
        title(FigTitle);
    end
    
    % Do rasters for Bout onsets and Bout offsets - interval of 3s or more
    % between two consecutive syllables is considered as two different
    % bouts
    
    bout_no = 1;
    BoutOnsetRaster = [];
    BoutOffsetRaster = [];
    if (MotifNo == 1)
        IntroNotes = find(NoteLabels == 'i');
        BoutOnsetSpikes = SpikeTimes(find((SpikeTimes > (NoteOnsets(1) - 3)) & (SpikeTimes < (NoteOnsets(1) + 0.2))));
        BoutOnsetSpikes = BoutOnsetSpikes - NoteOnsets(1);
        BoutOnsetRaster = [BoutOnsetRaster; [BoutOnsetSpikes (ones(length(BoutOnsetSpikes),1) * bout_no/5)]];

        BoutDifferences = NoteOnsets(2:end) - NoteOffsets(1:(end - 1));
        BoutOffsets = find(BoutDifferences > 3);
        BoutOnsets = BoutOffsets + 1;
        [TempOnsets, TempOnsetIndices] = intersect(BoutOnsets,IntroNotes);
        BoutOnsets = BoutOnsets(TempOnsetIndices);
        BoutOffsets = BoutOffsets(TempOnsetIndices);
        for i = 1:length(BoutOnsets),
            bout_no = bout_no + 1;
            BoutOnsetSpikes = SpikeTimes(find((SpikeTimes > (NoteOnsets(BoutOnsets(i)) - 3)) & (SpikeTimes < (NoteOnsets(BoutOnsets(i)) + 0.2))));
            BoutOnsetSpikes = BoutOnsetSpikes - NoteOnsets(BoutOnsets(i));
            BoutOnsetRaster = [BoutOnsetRaster; [BoutOnsetSpikes (ones(length(BoutOnsetSpikes),1) * bout_no/5)]];
        end
        
        figure(BoutOnsetOffsetFigure);
        axes(BoutOnsetPlot);
        hold on;
        set(gca,'FontSize',18);
        xlabel('Time (sec)');
        edges = -3:0.1:0.2;
        pst = histc(BoutOnsetRaster(:,1),edges);
        pst = pst/bout_no;
        pst = pst/0.1;
        
        if (max(pst) > 100)
            pst = pst/100;
            ylabel('Firing rate (/100) (Hz)');
        else
            if (max(pst) > 10)
                pst = pst/10;
                ylabel('Firing rate (/10) (Hz)');
            else
                ylabel('Firing rate (Hz)');                
            end
        end
        
        BoutOnsetPST = bar(edges,pst,'histc');
        set(BoutOnsetPST,'EdgeColor','k');
        set(BoutOnsetPST,'FaceColor','w');
        BoutOnsetRaster(:,2) = BoutOnsetRaster(:,2) + max(pst) + 0.5;
        plot(BoutOnsetRaster(:,1),BoutOnsetRaster(:,2),'w+');
        marker_string = repmat('|',size(BoutOnsetRaster,1),1);
        text(BoutOnsetRaster(:,1),BoutOnsetRaster(:,2),marker_string,'HorizontalAlignment','center','VerticalAlignment','middle','FontSize',1.5,'FontName','FixedWidth','FontUnits','pixels');
        axis([-3 0.2 0 (BoutOnsetRaster(end,2) + 0.2)]);
        
        bout_no = 0;
        for i = 1:length(BoutOffsets),
            bout_no = bout_no + 1;
            BoutOffsetSpikes = SpikeTimes(find((SpikeTimes > (NoteOffsets(BoutOffsets(i)) - 0.2)) & (SpikeTimes < (NoteOnsets(BoutOffsets(i)) + 3))));
            BoutOffsetSpikes = BoutOffsetSpikes - NoteOffsets(BoutOffsets(i));
            BoutOffsetRaster = [BoutOffsetRaster; [BoutOffsetSpikes (ones(length(BoutOffsetSpikes),1) * bout_no/5)]];
        end
        
        bout_no = bout_no + 1;
        BoutOffsetSpikes = SpikeTimes(find((SpikeTimes > (NoteOffsets(end) - 0.2)) & (SpikeTimes < (NoteOnsets(end) + 3))));
        BoutOffsetSpikes = BoutOffsetSpikes - NoteOffsets(end);
        BoutOffsetRaster = [BoutOffsetRaster; [BoutOffsetSpikes (ones(length(BoutOffsetSpikes),1) * bout_no/5)]];
        
        figure(BoutOnsetOffsetFigure);
        axes(BoutOffsetPlot);
        hold on;
        set(gca,'FontSize',18);
        xlabel('Time (sec)');
        edges = -0.2:0.1:3;
        pst = histc(BoutOffsetRaster(:,1),edges);
        pst = pst/bout_no;
        pst = pst/0.1;
        
        if (max(pst) > 100)
            pst = pst/100;
            ylabel('Firing rate (/100) (Hz)');
        else
            if (max(pst) > 10)
                pst = pst/10;
                ylabel('Firing rate (/10) (Hz)');
            else
                ylabel('Firing rate (Hz)');
            end
        end
        
        BoutOffsetPST = bar(edges,pst,'histc');
        set(BoutOffsetPST,'EdgeColor','k');
        set(BoutOffsetPST,'FaceColor','w');
        BoutOffsetRaster(:,2) = BoutOffsetRaster(:,2) + max(pst) + 0.5;
        plot(BoutOffsetRaster(:,1),BoutOffsetRaster(:,2),'w+');
        marker_string = repmat('|',size(BoutOffsetRaster,1),1);
        text(BoutOffsetRaster(:,1),BoutOffsetRaster(:,2),marker_string,'HorizontalAlignment','center','VerticalAlignment','middle','FontSize',1.5,'FontName','FixedWidth','FontUnits','pixels');
        axis([-0.2 3 0 (BoutOffsetRaster(end,2) + 0.2)]);
        title(FigTitle);
    end
        
    % do the CV plot - Fano Factor = variance/mean of spike count in a
    % sliding window of size 30ms - slid across the song in 1ms intervals
    
    if (MotifNo == 1)
        CVBinPST = [];
        for SpikeTrainNo = 1:length(SpikeTrain),
            TempSpikeTrain = cell2mat(SpikeTrain{SpikeTrainNo});
            for i = 0:1:round(motif_durations(median_motif)*1000),
                CVBinPST(SpikeTrainNo,(i+1)) = length(find((TempSpikeTrain > i) & (TempSpikeTrain < (i + 30))));
            end
        end
        figure(CVFigure);
        axes(CVPlot);
        set(gca,'FontSize',18);
        CV = var(CVBinPST)./mean(CVBinPST);
        CV(isnan(CV)) = 0;
        x_axis = (0:1:round(motif_durations(median_motif) * 1000))/1000;
        plot(x_axis,CV,'k');
        axis([0 motif_durations(median_motif) 0 max(CV) + 0.2]);
        xlabel('Time (sec)','FontSize',18)
        ylabel('Fano factor','FontSize',18)
        
        axes(CVRasterPlot);
        plot(Raster((find(Raster(:,1) >= 0)),1),Raster((find(Raster(:,1) >= 0)),2),'w+');
        marker_string = repmat('|',size(Raster(find(Raster(:,1) >= 0)),1),1);
        text(Raster((find(Raster(:,1) >= 0)),1),Raster((find(Raster(:,1) >= 0)),2),marker_string,'HorizontalAlignment','center','VerticalAlignment','middle','FontSize',3,'FontName','FixedWidth','FontUnits','pixels');
        axis([0 motif_durations(median_motif) 0 max(Raster(:,2)) + 0.2]);
        set(gca,'FontSize',18);
        set(gca,'xtick',[]);
        ylabel('Trials','FontSize',18)
        title(FigTitle);
    end
    
        
    temp_motifs = [];
    for i = 1:length(motif),
        temp_motifs = [temp_motifs; (motifs + (i - 1))'];
    end
%     NoteLabels(temp_motifs) = [];
%     NoteOnsets(temp_motifs) = [];
%     NoteOffsets(temp_motifs) = [];
end

disp(['Total number of motifs is ',num2str(length(motif_durations))]);
disp(['Mean firing rate during song is ',num2str(mean(MotifRelatedFiringRate)),' Hz']);
disp('Finished');