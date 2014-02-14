function [ST, Wv1, Wv2, Wv3, Wv4, Threshold] = MedianThreshold_FindSpikes(RawData, StartTime, Fs, Threshold, PreNSamples, PostNSamples, ArtefactV, LoPassFilt, HiPassFilt, PosNegFirst, UpSampleFactor, LogFid)

[b, a] = butter(3, [HiPassFilt/16000 LoPassFilt/16000]);
RawData = filtfilt(b, a, RawData);

if (strfind(PosNegFirst, 'negative'))
    RawData = -RawData;
end

ST = [];
Wv1 = [];
Wv2 = [];
Wv3 = [];
Wv4 = [];

Time(:,1) = (1:1:size(RawData,1))/Fs;

ITime = linspace(Time(1), Time(end), length(Time)*UpSampleFactor);
SmoothData = spline(Time, RawData, ITime);

[Pks, Locs] = findpeaks(RawData, 'MinPeakHeight', Threshold);

SpikeOnsets = Time(Locs);

SpikeOnsets = sort(SpikeOnsets);
SpikeOnsets = unique(SpikeOnsets);

RegSpikeWaveforms = [];
InterpSpikeWaveforms = [];

SpikeCount = 0;
for i = 1:length(SpikeOnsets),
    SpikeIndex = round(SpikeOnsets(i) * Fs);
    
    if ((SpikeIndex <= PreNSamples) || (SpikeIndex > (size(RawData,1) - PostNSamples))) 
        continue;
    end
    
    SpikeWaveforms = RawData((SpikeIndex - PreNSamples):(SpikeIndex + PostNSamples - 1),:);
    SWIndex = Time(SpikeIndex - PreNSamples):1/Fs:Time(SpikeIndex + PostNSamples - 1);
    ISWIndex = linspace(SWIndex(1), SWIndex(end), length(SWIndex) * UpSampleFactor);
    ISpikeWaveforms = spline(SWIndex, SpikeWaveforms, ISWIndex);
    [MaxValue, MaxIndex] = max(ISpikeWaveforms);
    
    ISpikeIndex = find(ITime <= ISWIndex(MaxIndex), 1, 'last');
    
    if (((ISpikeIndex) <= PreNSamples*UpSampleFactor) || ((ISpikeIndex) >= (length(SmoothData) - PostNSamples*UpSampleFactor + 1)))
        continue;
    end
    SpikeCount = SpikeCount + 1;
    ST(SpikeCount,1) = ISWIndex(MaxIndex);
    
    ISpikeIndex = find(ITime <= ST(SpikeCount, 1), 1, 'last');
        
    Wv1(:,SpikeCount) = decimate(SmoothData((ISpikeIndex - PreNSamples*UpSampleFactor):(ISpikeIndex + PostNSamples*UpSampleFactor - 1)), UpSampleFactor);
    Wv2(:,SpikeCount) = decimate(SmoothData((ISpikeIndex - PreNSamples*UpSampleFactor):(ISpikeIndex + PostNSamples*UpSampleFactor - 1)), UpSampleFactor);
    Wv3(:,SpikeCount) = decimate(SmoothData((ISpikeIndex - PreNSamples*UpSampleFactor):(ISpikeIndex + PostNSamples*UpSampleFactor - 1)), UpSampleFactor);
    Wv4(:,SpikeCount) = decimate(SmoothData((ISpikeIndex - PreNSamples*UpSampleFactor):(ISpikeIndex + PostNSamples*UpSampleFactor - 1)), UpSampleFactor);
end

ISpikeIndices = round(ST * Fs * UpSampleFactor);

ShortISIs = find(diff(ST) <= 1/Fs);
for i = 1:length(ShortISIs),
    [MinValue, RemovalIndices(i)] = min([SmoothData(ISpikeIndices(ShortISIs(i))) SmoothData(ISpikeIndices(ShortISIs(i) + 1))]);
    if (RemovalIndices(i) == 1)
        RemovalIndices(i) = ShortISIs(i);
    else
        RemovalIndices(i) = ShortISIs(i) + 1;
    end
end
if (exist('RemovalIndices', 'var'))
    ST(RemovalIndices) = [];
    Wv1(:, RemovalIndices) = [];
    Wv2(:, RemovalIndices) = [];
    Wv3(:, RemovalIndices) = [];
    Wv4(:, RemovalIndices) = [];
    disp(['Threw out ', num2str(length(RemovalIndices)), ' spikes with too short ISIs out of a total of ', num2str(length(ST)), ' spikes']);
    fprintf(LogFid, 'Threw out %i spikes with too short ISIs out of a total of %i spikes\n', length(RemovalIndices), (length(ST) + length(RemovalIndices)));
else
    disp(['Threw out 0 spikes with too short ISIs out of a total of ', num2str(length(ST)), ' spikes']);
    fprintf(LogFid, 'Threw out 0 spikes with too short ISIs out of a total of %i spikes\n', (length(ST)));
end

ArtefactIndices = [(find(max(Wv1) > ArtefactV(1))) (find(min(Wv1) < ArtefactV(2)))];
ArtefactIndices = unique(ArtefactIndices);
ST(ArtefactIndices) = [];
Wv1(:, ArtefactIndices) = [];
Wv2(:, ArtefactIndices) = [];
Wv3(:, ArtefactIndices) = [];
Wv4(:, ArtefactIndices) = [];
disp(['Threw out ', num2str(length(ArtefactIndices)), ' artefact spikes out of a total of ', num2str(length(ST)), ' spikes']);
fprintf(LogFid, 'Threw out %i artefact spikes out of a total of %i spikes\n', length(ArtefactIndices), (length(ST) + length(ArtefactIndices)));

[SpikeTimes, UniqueIndices] = unique(ST);
ST = ST(UniqueIndices);
Wv1 = Wv1(:, UniqueIndices);
Wv2 = Wv2(:, UniqueIndices);
Wv3 = Wv3(:, UniqueIndices);
Wv4 = Wv4(:, UniqueIndices);

ST = ST + StartTime;
