function [] = FindMatches(Bout, Threshold, OutputFileName, Color)

Fid = fopen(OutputFileName, 'w');

BoutPeaks = [];

Width = 0.005;
GaussianLen = 4;
Fs = Bout.BoutAmplitudeFs{1};

XGauss = 1:1:(1 + round(2 * GaussianLen * Width * (Fs)));
XGauss = XGauss - (length(XGauss) + 1)/2;
GaussWin = (1/((Width * Fs) * sqrt(2 * pi))) * exp(-(XGauss.*XGauss)/(2 * (Width * Fs) * (Width * Fs)));


for i = 1:length(Bout.MaxBoutSeqMatch),
    Matches = Bout.MaxBoutSeqMatch{i} > Threshold;
    h=[1 -1];
    temp = zeros(size(Matches,1), size(Matches,2));
    temp(find(Matches > 0)) = 1;
    trans=conv(h, temp, 'same');
    onsets=find(trans > 0);
    offsets=find(trans < 0);
    if (~isempty(onsets))
        for j = 1:length(onsets),
            [maxval, maxind] = max(Bout.MaxBoutSeqMatch{i}(onsets(j):offsets(j)));
            if (isfield(Bout, 'FileName'))
                fprintf(Fid, '%s\t%g\t%g\n', Bout.FileName{i}, Bout.T{i}(onsets(j) + maxind - 1), maxval);
            end
        end
    else
        if (isfield(Bout, 'FileName'))
            disp(['There were no values above the threshold ', num2str(Threshold), ' in file ', Bout.FileName{i}]);
        end
    end
    Temp = conv(Bout.MaxBoutSeqMatch{i}, GaussWin, 'same');
    
    [Peaks, Locs] = findpeaks(Temp);
    BoutPeaks = [BoutPeaks Peaks];
end

Edges = linspace(0, max(BoutPeaks)*1.1, 100);
if (exist('Color', 'var'))
    plot(Edges, histc(BoutPeaks, Edges)/sum(histc(BoutPeaks, Edges)), Color, 'LineWidth', 2);
else
    figure;
    BoutPeakBar = bar(Edges, histc(BoutPeaks, Edges));
end
fclose(Fid);
