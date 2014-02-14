function [PeakVals] = CalculateMaxMatches(Bout)

Width = 0.005;
GaussianLen = 4;
Fs = Bout.BoutAmplitudeFs{1};

XGauss = 1:1:(1 + round(2 * GaussianLen * Width * (Fs)));
XGauss = XGauss - (length(XGauss) + 1)/2;
GaussWin = (1/((Width * Fs) * sqrt(2 * pi))) * exp(-(XGauss.*XGauss)/(2 * (Width * Fs) * (Width * Fs)));

PeakVals = [];

for i = 1:length(Bout.MaxBoutSeqMatch),
    Temp = conv(Bout.MaxBoutSeqMatch{i}, GaussWin, 'same');
    PeakVals = [PeakVals; max(Temp)];
end