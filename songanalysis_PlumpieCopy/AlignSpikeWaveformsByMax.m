function [AlignedWaveforms] = AlignSpikeWaveformsByMax(TempWaveforms, MaxShift, Threshold, PreNSamples)

UpSamplingFactor = 4;
x = 1:1:size(TempWaveforms,2);
xx = linspace(x(1), x(end), length(x)*UpSamplingFactor);

UnAlignedIndex = 0;

FinalPreNSamples = 16;
FinalPostNSamples = 32;

Wv_xx = spline(x, TempWaveforms, xx);
[Pk_xx, PkIndex_xx] = max(Wv_xx(:, 20*UpSamplingFactor:28*UpSamplingFactor), [], 2);
PkIndex_xx = PkIndex_xx + 20*UpSamplingFactor - 1;
AlignIndex = PkIndex_xx;
AlignedWaveforms = Wv_xx(:, (AlignIndex - (FinalPreNSamples*UpSamplingFactor)):(AlignIndex + (FinalPostNSamples*UpSamplingFactor - 1)));
AlignedWaveforms = AlignedWaveforms - repmat(mean(AlignedWaveforms,2), 1, size(AlignedWaveforms, 2));
AlIndex = AlignIndex/UpSamplingFactor - PreNSamples;
AlignmentStatus = ones(size(AlignIndex))*1;

disp(['Finished aligning waveforms: ', num2str(UnAlignedIndex), ' of ', num2str(size(TempWaveforms,1)), ' had max shift greater or less than ', num2str(MaxShift)]);
disp(['Max shift is ', num2str(max(AlIndex)), ', min shift is ', num2str(min(AlIndex)), ' and median shift is ', num2str(median(AlIndex))]);