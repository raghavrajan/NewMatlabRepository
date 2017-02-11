function [MeanBoutPeaks, StdBoutPeaks, PeakVals, TotalFileTime] = CalculateMeanRandomMatches(FileName, Colour, BirdName, PlotOption, ThresholdMultiplier)

load(FileName);

BoutPeaks = [];
MaxBoutPeaks = [];
TotalFileTime = 0;
PeakVals = [];

for i = 1:length(WavBout.Peaks),
    for j = 1:length(BirdName),
        if (isempty((strfind(WavBout.FileName{i}, BirdName{j}))))
            IncludeFile = 1;
        else
            IncludeFile = 0;
            break;
        end
    end
            
    if (IncludeFile == 1)
        BoutPeaks = [BoutPeaks WavBout.Peaks{i}];
        MaxBoutPeaks = [MaxBoutPeaks; max(WavBout.Peaks{i})];
        TotalFileTime = TotalFileTime + WavBout.FileLength{i};
    end
end

for i = 1:length(ObsBout.Peaks),
    for j = 1:length(BirdName),
        if (isempty((strfind(ObsBout.FileName{i}, BirdName{j}))))
            IncludeFile = 1;
        else
            IncludeFile = 0;
            break;
        end
    end
            
    if (IncludeFile == 1)
        BoutPeaks = [BoutPeaks ObsBout.Peaks{i}];
        MaxBoutPeaks = [MaxBoutPeaks; max(ObsBout.Peaks{i})];
        TotalFileTime = TotalFileTime + ObsBout.FileLength{i};
    end
end

for i = 1:length(OKrankBout.Peaks),
    for j = 1:length(BirdName),
        if (isempty((strfind(OKrankBout.FileName{i}, BirdName{j}))))
            IncludeFile = 1;
        else
            IncludeFile = 0;
            break;
        end
    end
    if (IncludeFile == 1)
        BoutPeaks = [BoutPeaks OKrankBout.Peaks{i}];
        MaxBoutPeaks = [MaxBoutPeaks; max(OKrankBout.Peaks{i})];
        TotalFileTime = TotalFileTime + OKrankBout.FileLength{i};
    end
end

MeanBoutPeaks = mean(BoutPeaks);
StdBoutPeaks = std(BoutPeaks);
PeakVals = BoutPeaks(find(BoutPeaks > (mean(BoutPeaks) + ThresholdMultiplier*std(BoutPeaks))));

if (strfind(PlotOption, 'on'))
    Edges = linspace(0, max(BoutPeaks)*1.1, 100);
    plot(Edges, histc(BoutPeaks, Edges)/sum(histc(BoutPeaks, Edges)), Colour, 'LineWidth', 2);
    hold on;
    axis tight;
    temp = axis;
    plot([mean(BoutPeaks) mean(BoutPeaks)], [0 temp(4)], [Colour '--'], 'LineWidth', 3);
    plot([(mean(BoutPeaks) + ThresholdMultiplier*std(BoutPeaks)) (mean(BoutPeaks) + ThresholdMultiplier*std(BoutPeaks))], [0 temp(4)], [Colour ':'], 'LineWidth', 3);
end
disp(['Mean of random bout peaks is ', num2str(mean(BoutPeaks)), ' and std is ', num2str(std(BoutPeaks))]);
disp(['Mean of random max bout peaks is ', num2str(mean(MaxBoutPeaks)), ' and std is ', num2str(std(MaxBoutPeaks))]);
disp(['# of bout peaks > mean + ', num2str(ThresholdMultiplier), ' * std of bout peaks is ', num2str(length(find(BoutPeaks > (mean(BoutPeaks) + ThresholdMultiplier*std(BoutPeaks)))))]);
disp(['# of max bout peaks > mean + ', num2str(ThresholdMultiplier), ' * std of bout peaks is ', num2str(length(find(MaxBoutPeaks > (mean(BoutPeaks) + ThresholdMultiplier*std(BoutPeaks)))))]);
disp(['# of max bout peaks is ', num2str(length(MaxBoutPeaks))]);
disp(['Total file time is ', num2str(TotalFileTime)]);
