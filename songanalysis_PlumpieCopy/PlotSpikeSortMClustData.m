function [] = PlotSpikeSortMClustData(DirectoryName, DataFile, FileExt, FileType, ChannelNo, ClusterNo)

Colours = ['rgkmcy'];

if (strfind(FileType, 'obs'))
    if (FileExt(1) == '.')
        DataFileName = [DataFile,FileExt];
    else
        DataFileName = [DataFile,'.',FileExt];
    end
else
    DataFileName = [DataFile, FileExt];
end

if (DirectoryName(end) ~= '/')
    DirectoryName = [DirectoryName,'/'];
end
    
cd(DirectoryName);

if (strfind(FileType, 'obs'))
    ChannelString = ['obs',num2str(ChannelNo),'r'];
    [RawData,Fs] = soundin_copy(DirectoryName,DataFileName,ChannelString);
else
    if (strfind(FileType, 'okrank'))
        [RawData, Fs] = ReadOKrankData(DirectoryName, DataFileName, ChannelNo);
    end
end

Times = 0:1/Fs:length(RawData)/Fs;
Times(end) = [];

    
SpikeTimes = struct('Times',{});
%     SpikeTimes(1).Times = TempSpikeTimes;
    
for i = 1:length(ClusterNo),
    SpikeTimes(i).Times = load([DataFile,'_',num2str(ClusterNo),'.times']);
%    SpikeTimes(i).Times = SpikeTimes(i).Times/10000;
end

figure;
plot(Times,RawData,'b');
hold on;
for i = 1:length(SpikeTimes),
    Temp = SpikeTimes(i).Times;
    for j = 1:length(SpikeTimes(i).Times),
        SpikeIndex = find(Times < Temp(j),1,'last');
        SpikeWaveformIndices = (SpikeIndex - 8):(SpikeIndex + 23);
        plot(Times(SpikeWaveformIndices),RawData(SpikeWaveformIndices),Colours(i));
    end
end