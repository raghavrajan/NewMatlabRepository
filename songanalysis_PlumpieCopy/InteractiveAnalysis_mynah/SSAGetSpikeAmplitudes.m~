function [SpikeAmplitudes, SpikeWaveforms] = SSAGetSpikeAmplitudes(DirectoryName, RecFileDirectory, DataFileName,ChannelNo,SpikeTimes, FileType)

if (DirectoryName(end) ~= '/')
    DirectoryName = [DirectoryName,'/'];
end

if (strfind(FileType, 'obs'))
    if (ChannelNo(1) ~= 'o')
        ChannelString = ['obs',num2str(ChannelNo),'r'];
    else
        ChannelString = ChannelNo;
    end
end

if (strfind(FileType, 'obs'))
    [RawData,Fs] = SSASoundIn(DirectoryName, RecFileDirectory, DataFileName, ChannelString);
    RawData = RawData * 500/32768;
else
    if (strfind(FileType, 'okrank'))
        PresentDirectory = pwd;
        [RawData,Fs] = SSAReadOKrankData(DirectoryName, RecFileDirectory, DataFileName, ChannelNo);
        [b, a] = butter(3, [300/16000 6000/16000]);
        RawData = filtfilt(b, a, RawData);
        RawData = RawData * 100;
    end
end

Times = 0:1/Fs:length(RawData)/Fs;
Times(end) = [];

SpikeAmplitudes = [];
SpikeWaveforms = [];
for i = 1:length(SpikeTimes),
    Temp = SpikeTimes(i);
    SpikeIndex = find(Times < Temp,1,'last');
    if (((SpikeIndex - 16) < 1) || ((SpikeIndex + 47) > length(RawData)))
        SpikeAmplitudes(i,1) = max(RawData((SpikeIndex - 16):(SpikeIndex + 23))) - min(RawData((SpikeIndex - 8):(SpikeIndex + 23)));        
        SpikeWaveforms(i,:) = zeros(1,64);        
    else
        SpikeAmplitudes(i,1) = max(RawData((SpikeIndex - 16):(SpikeIndex + 23))) - min(RawData((SpikeIndex - 8):(SpikeIndex + 23)));
        SpikeWaveforms(i,:) = (RawData((SpikeIndex - 16):(SpikeIndex + 47)));
    end
end