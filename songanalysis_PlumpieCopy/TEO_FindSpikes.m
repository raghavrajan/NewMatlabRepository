function [ST, Wv1, Wv2, Wv3, Wv4, Threshold, TEOPeaks] = TEO_FindSpikes(RawData,UpperThreshold, WindowSize, StartTime, WindowToSkip, Fs, Threshold)

[b, a] = butter(3, [300/16000 6000/16000]);
RawData = filtfilt(b, a, RawData);

NoofTEOs = 5;

ST = [];
Wv1 = [];
Wv2 = [];
Wv3 = [];
Wv4 = [];

Time(:,1) = (1:1:size(RawData,1))/Fs;

for j =  1:NoofTEOs,
    TEO_K = j;
    for i = 1:size(RawData, 2),
        TEO{j}(:,i) = (RawData((1 + TEO_K):(end-TEO_K),i).*RawData((1 + TEO_K):(end-TEO_K),i)) - (RawData((1:(end - 2*TEO_K)),i).*(RawData(((2*TEO_K + 1):end),i)));
        TEO{j}(:,i) = conv(TEO{j}(:,i), hamming(4*TEO_K + 1), 'same');
    end
end

for i = 1:size(RawData,2),
    TempTEO = [];
    for j = 1:NoofTEOs,
        if (size(TEO{j},1) < size(TEO{1},1))
            SizeDifference = size(TEO{1},1) - size(TEO{j},1);
            TempTEO = [TempTEO [zeros(SizeDifference/2,1); TEO{j}(:,i); zeros(SizeDifference/2,1)]];
        else
            TempTEO = [TempTEO TEO{j}(:,i)];
        end
    end
    ATEO(:,i) = max(TempTEO, [], 2);
end

% for i = 1:size(TEO,2),
%     SortedTEO(:,i) = sort(TEO(:,i));
% end

TEO = ATEO;
if (Threshold == 0)
    Threshold = 10;
end

%Threshold = SortedTEO(round(0.99*size(SortedTEO,1)),:);

SpikeOnsets = [];
for i = 1:size(RawData,2),
    [Pks, Locs] = findpeaks(ATEO(:,i));
    ActualLocs = find((Pks > Threshold) & (Pks < 200));
    SpikeOnsets = [SpikeOnsets; Time(Locs(ActualLocs))];
end

TEOPeaks = Pks(ActualLocs);

% h=[1 -1];
% SpikeOnsets = [];
% for i = 1:size(RawData,2),
%     Temp = zeros(size(ThresholdCrossing,1), 1);
%     Temp(find(ThresholdCrossing(:,i) > 0)) = 1;
%     Transitions = conv(Temp, h, 'same');
%     SpikeOnsets = [SpikeOnsets; Time(find(Transitions > 0))];
% end

SpikeOnsets = sort(SpikeOnsets);
SpikeOnsets = unique(SpikeOnsets);
DoubleSpikeIndices = find(diff(SpikeOnsets) <= 0.0004);
SpikeOnsets(DoubleSpikeIndices + 1) = [];

RegSpikeWaveforms = [];
InterpSpikeWaveforms = [];

for i = 1:length(SpikeOnsets),
    SpikeIndex = round(SpikeOnsets(i) * Fs);
    if (SpikeIndex < 17)
        continue;
    end
    
    if (SpikeIndex > (size(RawData,1) - 32))
        continue;
    end
    SpikeWaveforms = RawData((SpikeIndex - 8):(SpikeIndex + 8),:);
    WFIndex = [0:1:47];
    RegSpikeWaveforms(i,:) = RawData((SpikeIndex - 16):(SpikeIndex + 31),:);
    InterpWFIndex = [0:0.1:47];
    InterpSpikeWaveforms(i,:) = interp1(WFIndex, RegSpikeWaveforms(i,:), InterpWFIndex);
    
    [Max, MaxIndices] = max(SpikeWaveforms);
    [Max, MaxChanNo] = max(max(SpikeWaveforms));
    ActualSpikeIndex = SpikeIndex - 8 + MaxIndices(MaxChanNo) - 1;
    if ((ActualSpikeIndex < 9) || (ActualSpikeIndex > (size(RawData,1) - 23)))
        continue;
    end
    ST(i,1) = Time(ActualSpikeIndex);
    Wv1(:,i) = RawData((ActualSpikeIndex - 8):(ActualSpikeIndex + 23), 1);
    Wv2(:,i) = RawData((ActualSpikeIndex - 8):(ActualSpikeIndex + 23), 1);
    Wv3(:,i) = RawData((ActualSpikeIndex - 8):(ActualSpikeIndex + 23), 1);
    Wv4(:,i) = RawData((ActualSpikeIndex - 8):(ActualSpikeIndex + 23), 1);    
end

[SpikeTimes, UniqueIndices] = unique(ST);
ST = ST(UniqueIndices);
Wv1 = Wv1(:, UniqueIndices);
Wv2 = Wv2(:, UniqueIndices);
Wv3 = Wv3(:, UniqueIndices);
Wv4 = Wv4(:, UniqueIndices);

ST = ST + StartTime;

disp('Finished');


