function [AlignedWaveforms] = AlignSpikeWaveforms(TempWaveforms, MaxShift, Threshold, PreNSamples)

UpSamplingFactor = 4;
x = 1:1:size(TempWaveforms,2);
xx = linspace(x(1), x(end), length(x)*UpSamplingFactor);

UnAlignedIndex = 0;

for i = 1:size(TempWaveforms,1),
    Wv_xx = spline(x, TempWaveforms(i,:), xx);
    [Pk_xx, PkIndex_xx] = max(Wv_xx(14*UpSamplingFactor:18*UpSamplingFactor));
    PkIndex_xx = PkIndex_xx + (14*UpSamplingFactor - 1);
    StartIndex = find(Wv_xx(1:PkIndex_xx) >= Threshold, 1, 'first');
    EndIndex = PkIndex_xx + find(Wv_xx((PkIndex_xx + 1):end) <= Threshold, 1, 'first');
    EndIndex = EndIndex - 1;
    AlignIndex = sum([StartIndex:1:EndIndex].*(Wv_xx(StartIndex:EndIndex) - Threshold))/sum(Wv_xx(StartIndex:EndIndex) - Threshold);
    if ((abs(AlignIndex/UpSamplingFactor - PreNSamples)) <= MaxShift)
        AlignedWaveforms(i,:) = spline((xx - (AlignIndex/UpSamplingFactor - PreNSamples)), Wv_xx, xx);
        AlIndex(i) = AlignIndex - PreNSamples;
        AlignmentStatus(i) = 1;
    else
        AlignedWaveforms(i,:) = Wv_xx;
        UnAlignedIndex = UnAlignedIndex + 1;
        AlIndex(i) = 0;
        AlignmentStatus(i) = 0;
    end
end
AlignedWaveforms(:,1:MaxShift) = [];
AlignedWaveforms(:,(end - MaxShift + 1):end) = [];
AlignedWaveforms = AlignedWaveforms;
disp(['Finished aligning waveforms: ', num2str(UnAlignedIndex), ' of ', num2str(size(TempWaveforms,1)), ' had max shift greater or less than ', num2str(MaxShift)]);
disp(['Max shift is ', num2str(max(AlIndex)), ', min shift is ', num2str(min(AlIndex)), ' and median shift is ', num2str(median(AlIndex))]);