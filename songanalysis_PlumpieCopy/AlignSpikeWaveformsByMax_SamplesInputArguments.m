function [AlignedWaveforms] = AlignSpikeWaveformsByMax_SamplesInputArguments(TempWaveforms, MaxShift, Threshold, PreNSamples, FinalPreNSamples, FinalPostNSamples)

UpSamplingFactor = 4;
x = 1:1:size(TempWaveforms,2);
xx = linspace(x(1), x(end), length(x)*UpSamplingFactor);

UnAlignedIndex = 0;

Wv_xx = spline(x, TempWaveforms, xx);
[Pk_xx, PkIndex_xx] = max(Wv_xx(:, 15*UpSamplingFactor:25*UpSamplingFactor), [], 2);
PkIndex_xx = PkIndex_xx + 15*UpSamplingFactor - 1;
AlignIndex = PkIndex_xx;
for i = 1:size(Wv_xx,1),
    if (((AlignIndex(i) - (FinalPreNSamples*UpSamplingFactor)) < 1) || ((AlignIndex(i) + (FinalPostNSamples*UpSamplingFactor - 1)) > size(Wv_xx,2)))
        UnAlignedIndex = UnAlignedIndex + 1;
    else
        AlignedWaveforms(i,:) = Wv_xx(i, (AlignIndex(i) - (FinalPreNSamples*UpSamplingFactor)):(AlignIndex(i) + (FinalPostNSamples*UpSamplingFactor - 1)));
    end
end
AlignedWaveforms = AlignedWaveforms - repmat(mean(AlignedWaveforms,2), 1, size(AlignedWaveforms, 2));
AlIndex = AlignIndex/UpSamplingFactor - PreNSamples;
AlignmentStatus = ones(size(AlignIndex))*1;

disp(['Finished aligning waveforms: ', num2str(UnAlignedIndex), '(', num2str(100 * UnAlignedIndex/size(Wv_xx,1)), '%) could not be aligned as there was not enough as there was not enough data']);