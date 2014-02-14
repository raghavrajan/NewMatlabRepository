function [] = SpikeDetect(DirectoryName, FileName, FileType, ChannelNo, LowerThreshold, UpperThreshold, WindowSize, SkipWindowSize)

T = [];
wv = [];

if (strfind(FileType, 'obs'))
    ChannelString = ['obs',num2str(ChannelNo),'r'];
    [RawData, Fs] = soundin_copy(DirectoryName, FileName, ChannelString);
    RawData = RawData * 500/32768;
end

[SpikeTimes, SpikeWaveforms] = DetectSpikes(RawData,LowerThreshold,UpperThreshold,WindowSize,SkipWindowSize,0,Fs);

SpikeWaveforms = SpikeWaveforms';

T = SpikeTimes * 10000;

wv = [wv; zeros(length(SpikeTimes),4,32)];

wv(:,1,:) = SpikeWaveforms;
wv(:,2,:) = SpikeWaveforms;
wv(:,3,:) = SpikeWaveforms;
wv(:,4,:) = SpikeWaveforms;

OutputFileName = [FileName,'.mat'];
    
save(OutputFileName,'T','wv');

disp('Finished');



