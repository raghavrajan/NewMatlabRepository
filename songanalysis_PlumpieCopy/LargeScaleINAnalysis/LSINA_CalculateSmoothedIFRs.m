function [ClusterIFR, OptionalClusterIFR] = LSINA_CalculateSmoothedIFRs(SURecordingDetails, Width, GaussianLen, Fs)

% Gaussian window for correlation
XGauss = 1:1:(1 + round(2 * GaussianLen * Width * (Fs)));
XGauss = XGauss - (length(XGauss) + 1)/2;
GaussWin = (1/((Width * Fs) * sqrt(2 * pi))) * exp(-(XGauss.*XGauss)/(2 * (Width * Fs) * (Width * Fs)));

% Initializing all the output variables to empty arrays
for i = 1:length(SURecordingDetails.ClusterSpikeTimes),
    ClusterIFR{i} = [];
    OptionalClusterIFR{i} = [];
end

for i = 1:length(SURecordingDetails.ClusterSpikeTimes),
    if (~isempty(SURecordingDetails.ClusterSpikeTimes{i}))
        Time = 0:1/Fs:SURecordingDetails.FileLen(i)/1000;
        ClusterIFR{i} = CalculateIFR(SURecordingDetails.ClusterSpikeTimes{i}/1000, Time);
        ClusterIFR{i} = conv(ClusterIFR{i}, GaussWin, 'same');
    end
    if ((~isempty(SURecordingDetails.ClusterSpikeTimes{i})) && (~isempty(SURecordingDetails.OptionalClusterSpikeTimes{i})))
        Time = 0:1/Fs:SURecordingDetails.FileLen(i)/1000;
        OptionalClusterIFR{i} = CalculateIFR([SURecordingDetails.ClusterSpikeTimes{i}/1000; SURecordingDetails.OptionalClusterSpikeTimes{i}/1000], Time);
        OptionalClusterIFR{i} = conv(OptionalClusterIFR{i}, GaussWin, 'same');
    end
end

