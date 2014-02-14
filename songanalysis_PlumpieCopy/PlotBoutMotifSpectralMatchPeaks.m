function [Threshold] = PlotBoutMotifSpectralMatchPeaks(FileName, Colour, DirOrUndir, Percentile)

load(FileName);

UnDirPeaks = [];
if (strfind(DirOrUndir, 'undir'))
    for i = 1:length(UnDirBout.Peaks),
        UnDirPeaks = [UnDirPeaks UnDirBout.Peaks{i}];
    end
else
    for i = 1:length(DirBout.Peaks),
        UnDirPeaks = [UnDirPeaks DirBout.Peaks{i}];
    end
end

Edges = linspace(0,1.1*max(UnDirPeaks), 100);
plot(Edges, histc(UnDirPeaks, Edges)/sum(histc(UnDirPeaks, Edges)), Colour);
disp([FileName, ': mean is ', num2str(mean(UnDirPeaks)), ' and std is ', num2str(std(UnDirPeaks))]);
disp([FileName, ': median is ', num2str(median(UnDirPeaks)), ' and mad is ', num2str(mad(UnDirPeaks))]);

Temp = sort(UnDirPeaks);
Threshold = Temp(round(Percentile*length(UnDirPeaks)/100));
disp([FileName, ': ', num2str(Percentile), '% percentile is ', num2str(Temp(round(Percentile*length(UnDirPeaks)/100)))]);