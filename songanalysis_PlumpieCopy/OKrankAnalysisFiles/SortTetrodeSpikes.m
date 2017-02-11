function [] = SortTetrodeSpikes(DirectoryName, FileList, ChanNo, ClusterNos, PreNSamples, PostNSamples, ArtefactV, ThreshMultiplier, LoPassFilt, HiPassFilt)

SpikeCount = 0;
TotalSpikeCount = 0;

T = [];
wv = [];

fid = fopen(FileList, 'r');
FileName = fgetl(fid);
RootFileName = FileName;
FileIndex = 0;

while (ischar(FileName(1)))
    FileIndex = FileIndex + 1;
    for i = 1:length(ChanNo),
        [RawData(:,i), Fs] = ReadOKrankData(DirectoryName, FileName, ChanNo(i));
        [b, a] = butter(3, [HiPassFilt/16000 LoPassFilt/16000]);
        RawData(:,i) = filtfilt(b, a, RawData(:,i));
    end
    RawData = sum(((RawData*(inv(cov(RawData)))).*RawData)');
    Threshold(FileIndex) = ThreshMultiplier * median(abs(RawData))/0.6745;
    FileName = fgetl(fid);
    clear RawData;
end
fclose(fid);

Threshold = median(median(Threshold));
disp(['Threshold is ', num2str(Threshold)]);

fid = fopen(FileList, 'r');
FileName = fgetl(fid);
RootFileName = FileName;
StartTime = 0;
FileIndex = 0;

while (ischar(FileName(1)))
    clear RawData;
    FileIndex = FileIndex + 1;
    for i = 1:length(ChanNo),
        [RawData(:,i), Fs] = ReadOKrankData(DirectoryName, FileName, ChanNo(i));
        RawData(:,i) = RawData(:,i) * 100;
    end
    [ST, Wv1, Wv2, Wv3, Wv4] = MedianThreshold_TetrodeFindSpikes([RawData/100], StartTime, Fs, Threshold, PreNSamples, PostNSamples, ArtefactV, LoPassFilt, HiPassFilt);
    FalseSpikes = find(ST < StartTime);
    ST(FalseSpikes) = [];
    Wv1(:,FalseSpikes) = [];
    Wv2(:,FalseSpikes) = [];
    Wv3(:,FalseSpikes) = [];
    Wv4(:,FalseSpikes) = [];
    
%     ShortISIs = find(diff(ST) < 0.0003);
%     ST(ShortISIs) = [];
%     Wv1(:,ShortISIs) = [];
%     Wv2(:,ShortISIs) = [];
%     Wv3(:,ShortISIs) = [];
%     Wv4(:,ShortISIs) = [];
%     disp(['Dropped ', num2str(length(ShortISIs)), ' spikes that had ISIs less than 0.3ms']);
    
    StartTime = StartTime + (length(RawData))/Fs;
    
    SpikeCount = length(ST);
    
    if (length(ST) > 0)
        T = [T; ST];
        wv(((TotalSpikeCount + 1):(TotalSpikeCount + SpikeCount)),1,:) = Wv1';
        wv(((TotalSpikeCount + 1):(TotalSpikeCount + SpikeCount)),2,:) = Wv2';
        wv(((TotalSpikeCount + 1):(TotalSpikeCount + SpikeCount)),3,:) = Wv3';
        wv(((TotalSpikeCount + 1):(TotalSpikeCount + SpikeCount)),4,:) = Wv4';
    end
    
  
    TotalSpikeCount = TotalSpikeCount + SpikeCount;
    
    %ContinueYN = menu(['Continue or not'],'Continue', 'Quit');
    %if (ContinueYN == 2)
    %    return;
    %end
    
    fprintf('Loaded %i spikes from datafile %s\n', length(ST), FileName)
%    disp(['Max is ', num2str(max(max(max(RawData1), max(RawData2)), max(max(RawData3), max(RawData4))))]);
    FileName = fgetl(fid);
end

PeakIndex = PreNSamples + 1;
[ClusterIds] = PlotTetrodeFeatures(T, wv, FileName, ClusterNos, PeakIndex, Fs);
if (size(ClusterIds,1) < size(ClusterIds,2))
    ClusterIds = ClusterIds';
end
if (size(T,1) < size(T,2))
    T = T';
end

for i = 1:max(ClusterIds),
    ClusterIndices = find(ClusterIds == i);
    OutputTimes = [T(ClusterIndices)];
    outputfile = [RootFileName,'_chan',int2str(ChanNo),'_',num2str(i),'.times'];
    save(outputfile, 'OutputTimes', '-ASCII');
end

disp(['No of ISIs that are smaller than 0 is ', num2str(length(find(diff(T) <= 0)))]);
outputfile = [RootFileName,'_chan',int2str(ChanNo),'.mat'];
T = T * 10000;
T = T';

save(outputfile,'T','wv');