function [] = PlotSequentialMotifLengths(Dir, Motif, BirdName, Percentile)

cd(Dir);
for i = 1:length(Motif),
    clear Files;
    clear Data;
    Files = dir(['*Syll', Motif(1:i), '.mat']);
    if (~isempty(Files))
        
        for j = 1:length(Files),
            Data{j} = load(Files(j).name);
        end

        figure;
        Threshold = PlotBoutMotifSpectralMatchPeaks(Files(1).name, 'b', 'undir', Percentile);
        axis([0 6 0 0.0015]);
        hold on;
        plot([Threshold Threshold], [0 0.0015], 'k--', 'LineWidth', 2);
        title(['Sequence - ', Motif(1:i)]);
        zoom xon;
        uiwait(gcf);
        Threshold = input(['What threshold do you want to use for sequence ', Motif(1:i), ': ']);
        for j = 1:length(Data),
            UnDirPeaks = [];
            UnDirPeakTimes = [];
            FileTime = 0;
            for k = 1:length(Data{j}.UnDirBout.Peaks),
                UnDirPeaks = [UnDirPeaks Data{j}.UnDirBout.Peaks{k}];
                UnDirPeakTimes = [UnDirPeakTimes (Data{j}.UnDirBout.PeakTimes{k} + FileTime)];
                FileTime = FileTime + Data{j}.UnDirBout.FileLength{k};
            end
            Matches = find(UnDirPeaks > Threshold);
            Diff = diff(UnDirPeakTimes(Matches));
            UnDirMatches(j,i) = length(Matches) - length(find(Diff < 0.05));
            

            DirPeaks = [];
            DirPeakTimes = [];
            FileTime = 0;
            for k = 1:length(Data{j}.DirBout.Peaks),
                DirPeaks = [DirPeaks Data{j}.DirBout.Peaks{k}];
                DirPeakTimes = [DirPeakTimes (Data{j}.DirBout.PeakTimes{k} + FileTime)];
                FileTime = FileTime + Data{j}.DirBout.FileLength{k};
            end
            Matches = find(DirPeaks > Threshold);
            Diff = diff(DirPeakTimes(Matches));
            DirMatches(j,i) = length(Matches) - length(find(Diff < 0.05));
        end
    end
end

for i = 1:size(UnDirMatches, 1),
    UnDirMatches(i,:) = UnDirMatches(i,:)/UnDirMatches(i,1);
    DirMatches(i,:) = DirMatches(i,:)/DirMatches(i,1);
end

figure;
set(gcf, 'Color', 'w');
subplot(2,1,1);
plot(UnDirMatches', 's-');
XLabels = [];
for i = 1:length(Motif),
    XLabels{i} = Motif(1:i);
end
set(gca, 'FontSize', 10);
set(gca, 'XTick', [1:1:length(Motif)], 'XTickLabel', XLabels);
ylabel('# of sequences/# of a', 'FontSize', 12, 'FontName', 'Arial');
title([BirdName, ' - UnDirected'], 'FontSize', 12, 'FontName', 'Arial');
LegendString{1} = 'Pre';
for i = 2:size(UnDirMatches,1),
    LegendString{i} = ['Post ', num2str(i-1)];
end
legend(LegendString);

subplot(2,1,2);
plot(DirMatches', 's-');
set(gca, 'FontSize', 10);
set(gca, 'XTick', [1:1:length(Motif)], 'XTickLabel', XLabels);
ylabel('# of sequences/# of a', 'FontSize', 12, 'FontName', 'Arial');
title([BirdName, ' - Directed'], 'FontSize', 12, 'FontName', 'Arial');
legend(LegendString);

disp('Finished');