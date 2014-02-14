function [] = ConvertOKrankData(DirectoryName, FileList, ChanNo, UpperThreshold)

SpikeCount = 0;
TotalSpikeCount = 0;

T = [];
wv = [];

fid = fopen(FileList, 'r');
FileName = fgetl(fid);
RootFileName = FileName;
StartTime = 0;
FileIndex = 0;

TEOPeakValues = [];
while (ischar(FileName(1)))
    FileIndex = FileIndex + 1;
    [RawData, Fs] = ReadOKrankData(DirectoryName, FileName, ChanNo);
    RawData = RawData * 100;
   
    WindowToSkip = 32;
    
%    [ST, Wv1, Wv2, Wv3, Wv4] = tetrode_find_spikes(RawData1, RawData2, RawData3, RawData4,UpperThreshold, 32, StartTime, WindowToSkip, Fs);

    if (~exist('Threshold', 'var'))
        Threshold = 0;
    end
    [ST, Wv1, Wv2, Wv3, Wv4, Threshold, TEOPeaks] = TEO_FindSpikes([RawData/100], UpperThreshold, 32, StartTime, WindowToSkip, Fs, Threshold);
    if (size(TEOPeaks,1) < size(TEOPeaks,2))
        TEOPeakValues = [TEOPeakValues; TEOPeaks'];
    else
        TEOPeakValues = [TEOPeakValues; TEOPeaks];
    end
    
    FalseSpikes = find(ST < StartTime);
    ST(FalseSpikes) = [];
    Wv1(:,FalseSpikes) = [];
    Wv2(:,FalseSpikes) = [];
    Wv3(:,FalseSpikes) = [];
    Wv4(:,FalseSpikes) = [];
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
    
    fprintf('Loaded %i spikes from datafile %s\n', length(ST), FileName)
%    disp(['Max is ', num2str(max(max(max(RawData1), max(RawData2)), max(max(RawData3), max(RawData4))))]);
    FileName = fgetl(fid);
end

disp(['No of ISIs that are smaller than 0 is ', num2str(length(find(diff(T) <= 0)))]);
outputfile = [RootFileName,'_chan',int2str(ChanNo),'.mat'];
T = T * 10000;
T = T';

save(outputfile,'T','wv');