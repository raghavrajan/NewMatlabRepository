function [SpikeAmplitudes, SpikeWaveforms] = GetSpikeAmplitudes(DirectoryName,DataFileName,ChannelString,SpikeTimes)

if (DirectoryName(end) ~= '/')
    DirectoryName = [DirectoryName,'/'];
end

[RawData,Fs] = soundin_copy(DirectoryName,DataFileName,ChannelString);
Times = 0:1/Fs:length(RawData)/Fs;
Times(end) = [];

SpikeAmplitudes = [];
SpikeWaveforms = [];
for i = 1:length(SpikeTimes),
    Temp = SpikeTimes(i);
    SpikeIndex = find(Times < Temp,1,'last');
    SpikeAmplitudes(i,1) = max(RawData((SpikeIndex - 8):(SpikeIndex + 23))) - min(RawData((SpikeIndex - 8):(SpikeIndex + 23)));
    SpikeWaveforms(i,:) = (RawData((SpikeIndex - 8):(SpikeIndex + 23)));
end