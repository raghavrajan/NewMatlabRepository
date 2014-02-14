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
        RawData = RawData * 100;
    end
end

Wp = [700 8000] * 2 / Fs;
Ws = [500 10000] * 2 / Fs;
[N, Wn] = buttord(Wp, Ws, 3, 20);
[b, a] = butter(N, Wn);
%   [b, a] = butter(4, [600/16000 6000/16000]);
RawData = filtfilt(b, a, RawData);

Times = [1:1:length(RawData)]/Fs;

SpikeAmplitudes = [];
SpikeWaveforms = [];
for i = 1:length(SpikeTimes),
    Temp = SpikeTimes(i);
    SpikeIndex = find(Times < Temp,1,'last');
    if (((SpikeIndex - 8) < 1) || ((SpikeIndex + 47) > length(RawData)))
        if ((SpikeIndex - 8) < 1)
            SpikeAmplitudes(i,1) = max(RawData(1:(SpikeIndex + 23))) - min(RawData(1:(SpikeIndex + 23)));
        else
            SpikeAmplitudes(i,1) = max(RawData((SpikeIndex - 8):end)) - min(RawData((SpikeIndex - 8):end));
        end
        SpikeWaveforms(i,:) = zeros(1,56);        
    else
        SpikeAmplitudes(i,1) = max(RawData((SpikeIndex - 8):(SpikeIndex + 23))) - min(RawData((SpikeIndex - 8):(SpikeIndex + 23)));
        SpikeWaveforms(i,:) = (RawData((SpikeIndex - 8):(SpikeIndex + 47)));
    end
end