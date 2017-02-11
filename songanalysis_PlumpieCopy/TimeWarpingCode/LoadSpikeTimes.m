function [SpikeTimes, SpikeAmplitudes, SpikeWaveforms, ClusterParameters] = LoadSpikeTimes(DirectoryName,FileNames,ChannelNo,RecordLengths, FileType)

cd(DirectoryName);
% Now load up all the information about syllable onsets and offsets from
% the song bout files

SpikeTimes = [];
SpikeWaveforms = [];
SpikeAmplitudes = [];

Prompt = {'Enter the cluster nos. that have to be loaded','Enter the total number of clusters','Type yes if you want to include outliers, otherwise type no'};
DialogTitle = 'Input for cluster parameters';
DefaultParameters = {'0','1','no'};

ClusterParameters = inputdlg(Prompt, DialogTitle, 1, DefaultParameters);

ClusterNos = [];

for i = 1:length(ClusterParameters{1})
    if ~(isnan(str2double(ClusterParameters{1}(i))))
        ClusterNos = [ClusterNos str2double(ClusterParameters{1}(i))];
    end
end

MaxClusters = str2double(ClusterParameters{2});
OutlierInclude = [ClusterParameters{3}];

RecTime = 0;

for i = 1:length(FileNames),
    disp([FileNames{i}]);
    DotIndex = find(FileNames{i} == '.');
    if (length(DotIndex) < 1)
        DotIndex = length(FileNames{i}) + 1;
    end
    SpikeTimeFile = [FileNames{i}(1:(DotIndex(end) - 1)),'.spk'];
    TempSpikeTimes = load(SpikeTimeFile);
 
    for ClusterIndex = 1:length(ClusterNos),
            TempSpikes = TempSpikeTimes(find((TempSpikeTimes(:,1) == ClusterNos(ClusterIndex))),2)/1000;
%            TempSpikes = TempSpikeTimes(find((TempSpikeTimes(:,1) == ClusterNos(ClusterIndex))),2);
            [TempSpikeAmplitudes,TempSpikeWaveforms] = GetSpikeAmplitudes(DirectoryName,FileNames{i},ChannelNo,TempSpikes, FileType);
            TempSpikes = TempSpikes + RecTime;
            SpikeAmplitudes = [SpikeAmplitudes; TempSpikeAmplitudes];
            SpikeWaveforms = [SpikeWaveforms; TempSpikeWaveforms];
            SpikeTimes = [SpikeTimes; TempSpikes];
%           The following section is to include outliers
            if (strfind(OutlierInclude,'yes'))
                TempSpikes = TempSpikeTimes(find((TempSpikeTimes(:,1) == (MaxClusters + ClusterNos(ClusterIndex) + 1))),2)/1000;
 %               TempSpikes = TempSpikeTimes(find((TempSpikeTimes(:,1) == (MaxClusters + ClusterNos(ClusterIndex) + 1))),2);
                [TempSpikeAmplitudes, TempSpikeWaveforms] = GetSpikeAmplitudes(DirectoryName,FileNames{i},ChannelNo,TempSpikes, FileType);
                TempSpikes = TempSpikes + RecTime;
                SpikeAmplitudes = [SpikeAmplitudes; TempSpikeAmplitudes];
                SpikeWaveforms = [SpikeWaveforms; TempSpikeWaveforms];
                SpikeTimes = [SpikeTimes; TempSpikes];
            end
    end
    RecTime = RecTime + RecordLengths(i);    
end
[SortedSpikeTimes, SortedIndices] = sort(SpikeTimes);
SpikeTimes = SpikeTimes(SortedIndices);
SpikeAmplitudes = SpikeAmplitudes(SortedIndices);
SpikeWaveforms = SpikeWaveforms(SortedIndices,:);
