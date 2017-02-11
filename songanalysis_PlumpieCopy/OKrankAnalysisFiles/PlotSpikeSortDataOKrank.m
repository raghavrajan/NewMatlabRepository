function [FiringRate] = PlotSpikeSortDataOKrank(DirectoryName, FileName, SpikeFile, ChannelNo, ClusterNos, MaxClusters, FileType, varargin)

if (nargin > 7)
    Range(1) = varargin{1};
    Range(2) = varargin{2};
end

Colours = ['rgkmcy'];

cd(DirectoryName);

MainFigure = figure;
set(gcf, 'Color', 'w');
if ~isempty(strfind(FileType, 'okrank'))
    [RawData, Fs] = ReadOKrankData(DirectoryName, FileName, 1);
else
    if (~isempty(strfind(handles.SS.FileType, 'obs')))
        [RawData, Fs] = soundin_copy([DirectoryName, '/'], FileName, 'obs0r');
        RawData = RawData * 5/32768;
    end
end
plot_motif_spectrogram(RawData, Fs, MainFigure,subplot(1,1,1))
SpecAxis = axis;
hold on;

if ~isempty(strfind(FileType, 'okrank'))
    [RawData, Fs] = ReadOKrankData(DirectoryName, FileName, ChannelNo);
else
    if (~isempty(strfind(handles.SS.FileType, 'obs')))
        [RawData, Fs] = soundin_copy([DirectoryName, '/'], FileName, ['obs', num2str(5-ChannelNo), 'r']);
        RawData = RawData * 5/32768;
    end
end
RawData = RawData * 2000;
RawData = RawData + 11000 + abs(min(RawData));
    
Times = 0:1/Fs:length(RawData)/Fs;
Times(end) = [];

if (~exist('Range', 'var'))
    Range(1) = Times(1);
    Range(2) = Times(end);
end

%TempSpikeTimes = load(['Chan', num2str(ChannelNo), '_', FileName,'.spk']);
TempSpikeTimes = load(SpikeFile);

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
end

plot(Times,RawData,'b');
hold on;
for i = 1:length(SpikeTimes),
    Temp = SpikeTimes(i).Times;
    for j = 1:length(SpikeTimes(i).Times),
        SpikeIndex = find(Times <= Temp(j),1,'last');
        if ((SpikeIndex > 8) && (SpikeIndex < (length(Times) - 23)))
            SpikeWaveformIndices = (SpikeIndex - 8):(SpikeIndex + 23);
        else
            if (SpikeIndex <= 8)
                SpikeWaveformIndices = 1:(SpikeIndex + 23);
            else
                SpikeWaveformIndices = (SpikeIndex - 8):length(Times);
            end
        end
        plot(Times(SpikeWaveformIndices),RawData(SpikeWaveformIndices),Colours(i));
    end
end
axis tight;
RawDataAxis = axis;
RawDataAxis(1:2) = SpecAxis(1:2);
axis(RawDataAxis);
RangeSpikeTimes = [];
for i = 1:length(SpikeTimes),
    RangeSpikeTimes = [RangeSpikeTimes; find((SpikeTimes(i).Times >= Range(1)) & (SpikeTimes(i).Times <= Range(2)))];
end

FiringRate = (length(RangeSpikeTimes)/(Range(2) - Range(1)));
disp(['Firing rate is ', num2str(length(RangeSpikeTimes)/(Range(2) - Range(1))), ' Hz']);