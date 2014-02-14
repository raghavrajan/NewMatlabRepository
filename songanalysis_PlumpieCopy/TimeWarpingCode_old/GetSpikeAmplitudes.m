function [SpikeAmplitudes, SpikeWaveforms] = GetSpikeAmplitudes(DirectoryName,DataFileName,ChannelNo,SpikeTimes)

if (DirectoryName(end) ~= '/')
    DirectoryName = [DirectoryName,'/'];
end

if (ChannelNo(1) ~= 'o')
    ChannelString = ['obs',num2str(ChannelNo),'r'];
else
    ChannelString = ChannelNo;
end

[RawData,Fs] = soundin_copy(DirectoryName,DataFileName,ChannelString);
RawData = RawData * 500/32768;
Times = 0:1/Fs:length(RawData)/Fs;
Times(end) = [];

SpikeAmplitudes = [];
SpikeWaveforms = [];
for i = 1:length(SpikeTimes),
    Temp = SpikeTimes(i);
    SpikeIndex = find(Times < Temp,1,'last');
    if (((SpikeIndex - 8) < 1) || ((SpikeIndex + 47) > length(RawData)))
        SpikeAmplitudes(i,1) = max(RawData((SpikeIndex - 8):(SpikeIndex + 23))) - min(RawData((SpikeIndex - 8):(SpikeIndex + 23)));        
        SpikeWaveforms(i,:) = zeros(1,56);        
    else
        SpikeAmplitudes(i,1) = max(RawData((SpikeIndex - 8):(SpikeIndex + 23))) - min(RawData((SpikeIndex - 8):(SpikeIndex + 23)));
        SpikeWaveforms(i,:) = (RawData((SpikeIndex - 8):(SpikeIndex + 47)));
    end
end