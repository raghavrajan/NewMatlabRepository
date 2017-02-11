function [ST, Wv1, Wv2, Wv3, Wv4, Threshold] = MedianThreshold_TetrodeFindSpikes(RawData, StartTime, Fs, Threshold, PreNSamples, PostNSamples, ArtefactV, LoPassFilt, HiPassFilt)

[b, a] = butter(3, [HiPassFilt/16000 LoPassFilt/16000]);
for i = 1:size(RawData,2),
    RawData(:,i) = filtfilt(b, a, RawData(:,i));
end

ST = [];
Wv1 = [];
Wv2 = [];
Wv3 = [];
Wv4 = [];

Time(:,1) = (1:1:size(RawData,1))/Fs;
CovRawData = sum(((RawData*(inv(cov(RawData)))).*RawData)');

if (size(CovRawData,1) < size(CovRawData,2))
    CovRawData = CovRawData';
end

ThresholdCrossing = CovRawData > Threshold;

h=[1 -1];
SpikeOnsets = [];
Temp = zeros(size(ThresholdCrossing));
Temp(ThresholdCrossing > 0) = 1;

Transitions = conv(Temp, h, 'same');
SpikeOnsets = [SpikeOnsets; Time(Transitions > 0)];

SpikeOnsets = sort(SpikeOnsets);
SpikeOnsets = unique(SpikeOnsets);
DoubleSpikeIndices = find(diff(SpikeOnsets) <= 0.0001);
SpikeOnsets(DoubleSpikeIndices + 1) = [];

RegSpikeWaveforms = [];
InterpSpikeWaveforms = [];

SpikeCount = 0;
OmittedSpikes = 0;

for i = 1:length(SpikeOnsets),
    SpikeIndex = round(SpikeOnsets(i) * Fs);
    
    if (SpikeIndex > (size(RawData,1) - 11))
        continue;
    end
    
    SpikeWaveforms = RawData((SpikeIndex):(SpikeIndex + 10),:);
    [Max, MaxIndices] = max(SpikeWaveforms);
    [Max, MaxChanNo] = max(max(SpikeWaveforms));
    ActualSpikeIndex = SpikeIndex + MaxIndices(MaxChanNo) - 1;
        
    if ((ActualSpikeIndex < 2*PreNSamples) || (ActualSpikeIndex > (size(RawData,1) - 2*PostNSamples)))
        OmittedSpikes = OmittedSpikes + 1;
        continue;
    end
    
    SpikeIndex = ActualSpikeIndex + 1;
    
    SpikeWaveforms = RawData((SpikeIndex - PreNSamples*2):(SpikeIndex + PostNSamples*2),:);
    SWIndex = linspace(0, length(SpikeWaveforms), length(SpikeWaveforms));
    ISWIndex = linspace(0, length(SpikeWaveforms), length(SpikeWaveforms)*10);
    ISpikeWaveforms = interp1(SWIndex, SpikeWaveforms, ISWIndex, 'spline');
    StartI = find(ISpikeWaveforms(1:(PreNSamples*2*10),MaxChanNo) < Threshold/100, 1, 'last') + 1;
    EndI = PreNSamples*2*10 + find(ISpikeWaveforms(PreNSamples*2*10:end, MaxChanNo) < Threshold/100, 1, 'first') - 1;
    if (isempty(StartI) || isempty(EndI))
        OmittedSpikes = OmittedSpikes + 1;
        continue;
    end
    MaxIndex = sum(ISWIndex(StartI:EndI).*(ISpikeWaveforms(StartI:EndI, MaxChanNo)' - Threshold/100))/sum(ISpikeWaveforms(StartI:EndI, MaxChanNo)' - Threshold/100);
    [MinValue, MinIndex] = min(abs(ISWIndex - MaxIndex));
    SpikeIndex = MinIndex;
    
    if ((SpikeIndex + PostNSamples*10 - 10) > length(ISpikeWaveforms) || ((SpikeIndex - PreNSamples*10 - 9) < 1))
        OmittedSpikes = OmittedSpikes + 1;
        continue;
    end
    
    SpikeCount = SpikeCount + 1;
    ST(SpikeCount,1) = Time(ActualSpikeIndex);

    Wv1(:,SpikeCount) = decimate(ISpikeWaveforms((SpikeIndex - PreNSamples*10 - 9):(SpikeIndex + PostNSamples*10 - 10),1), 10);
    Wv2(:,SpikeCount) = decimate(ISpikeWaveforms((SpikeIndex - PreNSamples*10 - 9):(SpikeIndex + PostNSamples*10 - 10),2), 10);
    Wv3(:,SpikeCount) = decimate(ISpikeWaveforms((SpikeIndex - PreNSamples*10 - 9):(SpikeIndex + PostNSamples*10 - 10),3), 10);
    Wv4(:,SpikeCount) = decimate(ISpikeWaveforms((SpikeIndex - PreNSamples*10 - 9):(SpikeIndex + PostNSamples*10 - 10),4), 10);    
end

disp(['Omitted ', num2str(OmittedSpikes), ' spikes']);
ArtefactIndices = [(find(max(Wv1) > ArtefactV(1))) (find(min(Wv1) < ArtefactV(2)))];
ArtefactIndices = unique(ArtefactIndices);
ST(ArtefactIndices) = [];
Wv1(:, ArtefactIndices) = [];
Wv2(:, ArtefactIndices) = [];
Wv3(:, ArtefactIndices) = [];
Wv4(:, ArtefactIndices) = [];

[SpikeTimes, UniqueIndices] = unique(ST);
ST = ST(UniqueIndices);
Wv1 = Wv1(:, UniqueIndices);
Wv2 = Wv2(:, UniqueIndices);
Wv3 = Wv3(:, UniqueIndices);
Wv4 = Wv4(:, UniqueIndices);

ST = ST + StartTime;

