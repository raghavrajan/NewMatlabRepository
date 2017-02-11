function [] = SortSpikes(DirectoryName, FileList, ChanNo, ClusterNos, PreNSamples, PostNSamples, ArtefactV, ThreshMultiplier, LoPassFilt, HiPassFilt, PosNegFirst)

Time = clock;
LogFileName = [DirectoryName, '/SortSpikes_', num2str(Time(1)), num2str(Time(2)), num2str(Time(3)), num2str(Time(4)), num2str(Time(5)), num2str(round(Time(6))), '.log'];
LogFid = fopen(LogFileName, 'w');
fprintf(LogFid, 'Sort Spikes log file with parameters\n');
fprintf(LogFid, 'Directory name : %s\n', DirectoryName);    
fprintf(LogFid, 'File list file : %s\n', FileList);
fprintf(LogFid, 'Spike channel no : %i\n', ChanNo);
fprintf(LogFid, 'No of kmeans clusters : %i\n', ClusterNos);
fprintf(LogFid, 'No of pre samples : %i\n', PreNSamples);
fprintf(LogFid, 'No of post samples : %i\n', PostNSamples);
fprintf(LogFid, 'Artefact voltages : %i and %i\n', ArtefactV(1), ArtefactV(2));
fprintf(LogFid, 'Threshold multiplier : %i\n', ThreshMultiplier);
fprintf(LogFid, 'Low pass filter (Hz) : %i\n', LoPassFilt);
fprintf(LogFid, 'High pass filter (Hz) : %i\n', HiPassFilt);
fprintf(LogFid, 'Positive first or negative first waveform : %s\n', PosNegFirst);

UpSampleFactor = 4;
SpikeCount = 0;
TotalSpikeCount = 0;

T = [];
wv = [];

fid = fopen(FileList, 'r');
FileName = fgetl(fid);
RootFileName = FileName;
FileIndex = 0;

while (ischar(FileName(1)))
    fprintf(LogFid, '%s\n', FileName);
    FileIndex = FileIndex + 1;
    [RawData, Fs] = ReadOKrankData(DirectoryName, FileName, ChanNo);
    [b, a] = butter(3, [HiPassFilt/16000 LoPassFilt/16000]);
    RawData = filtfilt(b, a, RawData);

    if (strfind(PosNegFirst, 'negative'))
        RawData = -RawData;
    end
    
    Threshold(FileIndex) = ThreshMultiplier * median(abs(RawData))/0.6745;
%    FiltData = sqrt(filter(ones(1, 6)/6, 1, RawData.*RawData));
%    Threshold(FileIndex) = ThreshMultiplier * median(FiltData)/0.6745;
%     clear TEO ATEO;
%     for j =  3:6,
%         TEO_K = j;
%         for i = 1:size(RawData, 2),
%             TEO{j-2}(:,i) = (RawData((1 + TEO_K):(end-TEO_K),i).*RawData((1 + TEO_K):(end-TEO_K),i)) - (RawData((1:(end - 2*TEO_K)),i).*(RawData(((2*TEO_K + 1):end),i)));
%             TEO{j-2}(:,i) = conv(TEO{j-2}(:,i), hamming(4*TEO_K + 1), 'same');
%         end
%     end
% 
%     for i = 1:size(RawData,2),
%         TempTEO = [];
%         for j = 1:length(TEO),
%             if (size(TEO{j},1) < size(TEO{1},1))
%                 SizeDifference = size(TEO{1},1) - size(TEO{j},1);
%                 TempTEO = [TempTEO [zeros(SizeDifference/2,1); TEO{j}(:,i); zeros(SizeDifference/2,1)]];
%             else
%                 TempTEO = [TempTEO TEO{j}(:,i)];
%             end
%         end
%         ATEO(:,i) = max(TempTEO, [], 2);
%     end
%     Threshold(FileIndex) = ThreshMultiplier*mean(ATEO);
    FileName = fgetl(fid);
end
fclose(fid);

Threshold = median(Threshold);
disp(['Threshold is ', num2str(Threshold*100), ' uV']);
fprintf(LogFid, 'Threshold is %guV\n', Threshold*100);

fid = fopen(FileList, 'r');
FileName = fgetl(fid);
RootFileName = FileName;
StartTime = 0;
FileIndex = 0;

while (ischar(FileName(1)))
    FileIndex = FileIndex + 1;
    [RawData, Fs] = ReadOKrankData(DirectoryName, FileName, ChanNo);
    RawData = RawData * 100;
   
    [ST, Wv1, Wv2, Wv3, Wv4] = MedianThreshold_FindSpikes([RawData/100], StartTime, Fs, Threshold, PreNSamples, PostNSamples, ArtefactV, LoPassFilt, HiPassFilt, PosNegFirst, UpSampleFactor, LogFid);
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
    
    fprintf('Loaded %i spikes from datafile %s\n', length(ST), FileName);
    fprintf(LogFid, 'Loaded %i spikes from datafile %s\n', length(ST), FileName);
%    disp(['Max is ', num2str(max(max(max(RawData1), max(RawData2)), max(max(RawData3), max(RawData4))))]);
    FileName = fgetl(fid);
end

UpSampleFactor = 1;

PeakIndex = PreNSamples*UpSampleFactor + 1;
[ClusterIds] = PlotFeatures(T, wv, FileName, ClusterNos, PeakIndex, Fs, UpSampleFactor, LogFid);
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
fprintf(LogFid, 'No of ISIs that are smaller than 0 is %i\n', length(find(diff(T) <= 0)));
outputfile = [RootFileName,'_chan',int2str(ChanNo),'.mat'];
fprintf(LogFid, 'Saved data to file : %s\n', outputfile);
fclose(LogFid);
T = T * 10000;
T = T';

save(outputfile,'T','wv');