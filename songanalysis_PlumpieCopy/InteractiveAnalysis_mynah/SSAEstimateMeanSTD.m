function [Mean, STD] = SSAEstimateMeanSTD(DirectoryName, RecFileDirectory, FileInfo, ChannelNo, FileType)

Mean = [];
STD = [];

if (strfind(FileType, 'obs'))
    if (ChannelNo(1) ~= 'o')
        ChannelString = ['obs',num2str(ChannelNo),'r'];
    else
        ChannelString = ChannelNo;
    end
end

for i = 1:length(FileInfo.FileNames),
    DataFileName = FileInfo.FileNames{i};
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
    RawData2 = RawData;
    Times = 0:1/Fs:length(RawData)/Fs;
    Times(end) = [];

    SpikeTimes = FileInfo.SpikeData.Times{i};
    
    for i = 1:length(SpikeTimes),
        Temp = SpikeTimes(i);
        SpikeIndex = find(Times < Temp,1,'last');
        RawData((SpikeIndex - 8):(SpikeIndex + 23)) = -1000;
    end
    RawData(find(RawData < -999)) = [];
    Mean(i,1) = mean(RawData2);
    Mean(i,2) = mean(RawData);
    STD(i,1) = std(RawData2);
    STD(i,2) = std(RawData);
    disp(DataFileName);
    disp(['Mean and std of complete data is ', num2str(Mean(i,1)), ' +/- ', num2str(STD(i,1))]);
    disp(['Mean and std of data with spikes removed is ', num2str(Mean(i,2)), ' +/- ', num2str(STD(i,2))]);
end
