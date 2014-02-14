function [SpikeTimes, SpikeAmplitudes, SpikeWaveforms, ClusterParameters] = LoadMClustSpikeTimes(DirectoryName,FileNames,ChannelNo,RecordLengths)

cd(DirectoryName);
% Now load up all the information about syllable onsets and offsets from
% the song bout files

SpikeTimes = [];
SpikeWaveforms = [];
SpikeAmplitudes = [];

Prompt = {'Enter the cluster nos. that have to be loaded'};
DialogTitle = 'Input for cluster parameters';
DefaultParameters = {'0'};

ClusterParameters = inputdlg(Prompt, DialogTitle, 1, DefaultParameters);

ClusterNos = [];

for i = 1:length(ClusterParameters{1})
    if ~(isnan(str2double(ClusterParameters{1}(i))))
        ClusterNos = [ClusterNos str2double(ClusterParameters{1}(i))];
    end
end

RecTime = 0;

for i = 1:length(FileNames),
    disp([FileNames{i}]);
    SpikeTimeFile = [FileNames{i},'_',num2str(ClusterNos),'.times'];
    TempSpikes = load(SpikeTimeFile);
    [TempSpikeAmplitudes,TempSpikeWaveforms] = GetSpikeAmplitudes(DirectoryName,FileNames{i},ChannelNo,TempSpikes);
    TempSpikes = TempSpikes + RecTime;
    SpikeAmplitudes = [SpikeAmplitudes; TempSpikeAmplitudes];
    SpikeWaveforms = [SpikeWaveforms; TempSpikeWaveforms];
    SpikeTimes = [SpikeTimes; TempSpikes];
    RecTime = RecTime + RecordLengths(i);    
end

[SortedSpikeTimes, SortedIndices] = sort(SpikeTimes);
SpikeTimes = SpikeTimes(SortedIndices);
SpikeAmplitudes = SpikeAmplitudes(SortedIndices);
SpikeWaveforms = SpikeWaveforms(SortedIndices,:);
