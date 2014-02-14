function [] = PlotExampleSyllables(TemplateFile, FileName, RawDataDir, FileType, Threshold, PreDur, PostDur)

load(TemplateFile);
load(FileName);

UnDirPeaks = [];
for i = 1:length(UnDirBout.Peaks),
    UnDirPeaks = [UnDirPeaks UnDirBout.Peaks{i}];
end
figure;
subplot(2,1,1);
Edges = linspace(0,1.1*max(UnDirPeaks), 100);
plot(Edges, histc(UnDirPeaks, Edges)/sum(histc(UnDirPeaks, Edges)));
hold on;
for i = 1:length(Threshold),
    plot([Threshold(i) Threshold(i)], [0 0.0015], 'k', 'LineWidth', 2);
end
axis tight;
temp = axis;
temp(3:4) = [0 0.0015];
axis(temp);

subplot(2,1,2);
contourf(MotifTemplate);

MatchIndex = 0;
Index = 0;
for i = 1:length(UnDirBout.Peaks),
    if (length(Threshold) == 1)
        Matches = find(UnDirBout.Peaks{i} > Threshold);
    else
        Matches = find((UnDirBout.Peaks{i} > Threshold(1)) & (UnDirBout.Peaks{i} < Threshold(2)));
    end
    if (~isempty(Matches))
        Diff = diff(UnDirBout.PeakTimes{i}(Matches));
        Matches(find(Diff < 0.02)) = [];
        if (~isempty(find(Diff < 0.02)))
            disp(['Dropped ', num2str(find(Diff < 0.02)), ' syllables']);
        end
        
        for j = 1:length(Matches),
            MatchIndex = MatchIndex + 1;
            Peaks(MatchIndex) = UnDirBout.Peaks{i}(Matches(j));
            PeakTimes(MatchIndex) = UnDirBout.PeakTimes{i}(Matches(j));
            FileNames{MatchIndex} = UnDirBout.FileName{i};
        end
    end
end

[SortedPeaks, SortedIndices] = sort(Peaks);

for i = 1:length(SortedIndices),
    if (mod(Index, 48) == 0)
        figure;
        set(gcf, 'Color', 'w');
        set(gcf, 'Position', [242 211 918 711]);
        Index = 0;
    end
    Index = Index + 1;
    subplot(6, 8, Index);
    PlotSpectrogram2([RawDataDir, '/'], FileNames{SortedIndices(i)}, FileType, 0, (PeakTimes(SortedIndices(i)) - PreDur), (PeakTimes(SortedIndices(i)) + PostDur));
    set(gca, 'XTick', []);
    set(gca, 'YTick', []);
    xlabel([]);
    ylabel([]);
    hold on;
    plot([PeakTimes(SortedIndices(i)) PeakTimes(SortedIndices(i))], [0 10000], 'b--', 'LineWidth', 2);
    title(num2str(Peaks(SortedIndices(i))), 'FontSize', 12);
end
