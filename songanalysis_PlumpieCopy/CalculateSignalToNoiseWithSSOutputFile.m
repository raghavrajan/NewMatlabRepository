function [SignalToNoise] = CalculateSignalToNoiseWithSSOutputFile(SSOutputFile, ClusterNos)

load(SSOutputFile);

MeanWF = mean(SSOutput.AlignedWaveforms);
Residuals = SSOutput.AlignedWaveforms - repmat(MeanWF, size(SSOutput.AlignedWaveforms,1), 1);

SignalToNoise = (max(MeanWF) - min(MeanWF))/std(Residuals(:));
