function [] = PlotSpikeSortData(DirectoryName,FileName,FileExt,FileNo,ChannelNo,ClusterNos,MaxClusters, varargin)

if (nargin > 7)
    Range(1) = varargin{1};
    Range(2) = varargin{2};
end

Colours = ['rgkmcy'];

for FileNos = 1:length(FileNo),

    if (FileNo(FileNos) < 10)
        DataFile = [FileName,'.00',num2str(FileNo(FileNos))];
    else
        if (FileNo < 100)
            DataFile = [FileName,'.0',num2str(FileNo(FileNos))];            
        else
            DataFile = [FileName,'.',num2str(FileNo(FileNos))];
        end
    end
    if (FileExt(1) == '.')
        DataFileName = [DataFile,FileExt];
    else
        DataFileName = [DataFile,'.',FileExt];
    end

    if (DirectoryName(end) ~= '/')
        DirectoryName = [DirectoryName,'/'];
    end
    
    cd(DirectoryName);

    ChannelString = ['obs',num2str(ChannelNo),'r'];

    MainFigure = figure;
    [RawData,Fs] = soundin_copy(DirectoryName,DataFileName,'obs0r');
    plot_motif_spectrogram(RawData, Fs, MainFigure,subplot(2,1,2))
    axes(subplot(2,1,2));
    SpecAxis = axis;
    
    [RawData,Fs] = soundin_copy(DirectoryName,DataFileName,ChannelString);
    RawData = RawData * 500/32768;
    
    Times = 0:1/Fs:length(RawData)/Fs;
    Times(end) = [];

    if (~exist('Range', 'var'))
        Range(1) = Times(1);
        Range(2) = Times(end);
    end
    
    TempSpikeTimes = load([DataFile,'.spk']);
    
    if (length(ClusterNos) < 1)
        ClusterNos = 0:1:MaxClusters;
    end
    
    SpikeTimes = struct('Times',{});
%     SpikeTimes(1).Times = TempSpikeTimes;
    
    for i = 1:length(ClusterNos),
        SpikeTimes(i).Times = [];
        TempSpikes = TempSpikeTimes(find((TempSpikeTimes(:,1) == ClusterNos(i))),2);
        if (~isempty(TempSpikes))
            if (TempSpikes(end) > Times(end))
                TempSpikes = TempSpikes/1000;
            end
        end
        SpikeTimes(i).Times = [SpikeTimes(i).Times; TempSpikes];
        
        %TempSpikes = TempSpikeTimes(find((TempSpikeTimes(:,1) == (MaxClusters + ClusterNos(i) + 1))),2)/1000;
        %SpikeTimes(i).Times = [SpikeTimes(i).Times; TempSpikes];
    end

    subplot(2,1,1);
    plot(Times,RawData,'b');
    hold on;
    for i = 1:length(SpikeTimes),
        for j = 1:length(SpikeTimes(i).Times),
            Temp = SpikeTimes(i).Times;
            SpikeIndex = find(Times < Temp(j),1,'last');
            SpikeWaveformIndices = (SpikeIndex - 8):(SpikeIndex + 23);
            plot(Times(SpikeWaveformIndices),RawData(SpikeWaveformIndices),Colours(i));
        end
    end
    axis tight;
    RawDataAxis = axis;
    RawDataAxis(1:2) = SpecAxis(1:2);
    axis(RawDataAxis);
end
RangeSpikeTimes = [];
for i = 1:length(SpikeTimes),
    RangeSpikeTimes = [RangeSpikeTimes; find((SpikeTimes(i).Times >= Range(1)) & (SpikeTimes(i).Times <= Range(2)))];
end
    
disp(['Firing rate is ', num2str(length(RangeSpikeTimes)/(Range(2) - Range(1))), ' Hz']);