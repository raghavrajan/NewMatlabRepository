function [PeakVals, PeakIntervals, TotalTime, SongTimes, FileNames] = MicrolesionCalculateNoofMatches(Bout, Threshold, varargin)

if (nargin > 2)
    SongTime = varargin{1};
end

disp(['Threshold is ', num2str(Threshold)]);

Width = 0.005;
GaussianLen = 4;
Fs = Bout.BoutAmplitudeFs{1};

XGauss = 1:1:(1 + round(2 * GaussianLen * Width * (Fs)));
XGauss = XGauss - (length(XGauss) + 1)/2;
GaussWin = (1/((Width * Fs) * sqrt(2 * pi))) * exp(-(XGauss.*XGauss)/(2 * (Width * Fs) * (Width * Fs)));

PeakVals = [];
PeakIntervals = [];
TotalTime = 0;
SongTimes = [];
FileNames = [];

for i = 1:length(Bout.MaxBoutSeqMatch),
    Temp = conv(Bout.MaxBoutSeqMatch{i}, GaussWin, 'same');
    [Peaks, Locs] = findpeaks(Temp, 'MINPEAKHEIGHT', Threshold);
    
    if (size(Peaks,1) < size(Peaks,2))
        Peaks = Peaks';
    end
    
    if (size(Locs,1) < size(Locs,2))
        Locs = Locs';
    end
    
    PeakVals = [PeakVals; Peaks];
    PeakIntervals = [PeakIntervals; diff(Bout.T{i}(Locs)')];
    TotalTime = TotalTime + Bout.T{i}(end);
    if (nargin > 2)
        SongTimes = [SongTimes; ones(size(Peaks))*SongTime(i)];
    end
    FileNames = [FileNames; repmat({Bout.FileName{i}}, size(Peaks,1), size(Peaks,2))];
end