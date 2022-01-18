function [] = PlotSyllOnsetRasterPSTs(SyllOnsetCombinedData, TitleString, SyllLabels)

figure;
p = panel();
p.pack(1,length(SyllOnsetCombinedData.SmoothedPST));

for Sylls = 1:length(SyllOnsetCombinedData.SmoothedPST),
    p(1, Sylls).pack({4/5 1/5});
    p(1, Sylls, 1).select();
    hold on;
    TrialNum = 0;
    RandRasterTrials = randperm(length(SyllOnsetCombinedData.Raster{Sylls}));
    RandRasterTrials = RandRasterTrials(1:min(100, length(RandRasterTrials)));
    RandRasterTrials = sort(RandRasterTrials);
    for i = RandRasterTrials(:)',
        TrialNum = TrialNum + 1;
        if (~isempty(SyllOnsetCombinedData.Raster{Sylls}{i}))
            if (SyllOnsetCombinedData.DirUnDir{Sylls}(i) == 0)
                PlotRaster([SyllOnsetCombinedData.Raster{Sylls}{i}(:) ones(size(SyllOnsetCombinedData.Raster{Sylls}{i}(:)))], 'b', 0.5, TrialNum);
            else
                PlotRaster([SyllOnsetCombinedData.Raster{Sylls}{i}(:) ones(size(SyllOnsetCombinedData.Raster{Sylls}{i}(:)))], 'r', 0.5, TrialNum);
            end
        end
    end
    axis tight;
    Temp = axis;
    axis([SyllOnsetCombinedData.Time(1)*1000 SyllOnsetCombinedData.Time(end)*1000 0 TrialNum+2]);
    ylabel('Bout #');
    if (Sylls == round(length(SyllLabels)/2))
        title({TitleString; ['Syllable ', SyllLabels(Sylls)]});
    else
        title(['Syllable ', SyllLabels(Sylls)]);
    end
    text(-55, 1.03*(TrialNum+2), char(65), 'FontSize', 14, 'HorizontalAlignment', 'center', 'FontWeight', 'bold');
    p(1, Sylls, 2).select();
    hold on;
    DirTrials = find(SyllOnsetCombinedData.DirUnDir{Sylls} == 1);
    UnDirTrials = find(SyllOnsetCombinedData.DirUnDir{Sylls} == 0);

    patch([SyllOnsetCombinedData.Time*1000 fliplr(SyllOnsetCombinedData.Time*1000)], [(mean(SyllOnsetCombinedData.SmoothedPST{Sylls}(DirTrials,:)) - std(SyllOnsetCombinedData.SmoothedPST{Sylls}(DirTrials,:))/sqrt(length(DirTrials))) fliplr((mean(SyllOnsetCombinedData.SmoothedPST{Sylls}(DirTrials,:)) + std(SyllOnsetCombinedData.SmoothedPST{Sylls}(DirTrials,:))/sqrt(length(DirTrials))))], 'r', 'EdgeColor', 'none', 'FaceColor', 'r', 'FaceAlpha', 0.2);
    patch([SyllOnsetCombinedData.Time*1000 fliplr(SyllOnsetCombinedData.Time)*1000], [(mean(SyllOnsetCombinedData.SmoothedPST{Sylls}(UnDirTrials,:)) - std(SyllOnsetCombinedData.SmoothedPST{Sylls}(UnDirTrials,:))/sqrt(length(UnDirTrials))) fliplr((mean(SyllOnsetCombinedData.SmoothedPST{Sylls}(UnDirTrials,:)) + std(SyllOnsetCombinedData.SmoothedPST{Sylls}(UnDirTrials,:))/sqrt(length(UnDirTrials))))], 'b', 'EdgeColor', 'none', 'FaceColor', 'b', 'FaceAlpha', 0.2);
    
    DirPST = plot(1000*SyllOnsetCombinedData.Time, mean(SyllOnsetCombinedData.SmoothedPST{Sylls}(DirTrials,:)), 'r');
    UnDirPST = plot(1000*SyllOnsetCombinedData.Time, mean(SyllOnsetCombinedData.SmoothedPST{Sylls}(UnDirTrials,:)), 'b');

    xlabel('Time relative to start of motif (ms)');
    ylabel('Firing rate (Hz)');
    axis tight;
    Temp = axis;
    axis([SyllOnsetCombinedData.Time(1)*1000 SyllOnsetCombinedData.Time(end)*1000 0 1.35*Temp(4)]);
    text(-55, 1.35*(Temp(4)), char(66), 'FontSize', 14, 'HorizontalAlignment', 'center', 'FontWeight', 'bold');
    legend([DirPST UnDirPST], ['Directed (', num2str(length(DirTrials)), ')'], ['Undirected (', num2str(length(UnDirTrials)), ')'], 'FontSize', 6);
end
p.margintop = 15;
p.de.margin = 15;
p.fontsize = 8;
set(gcf, 'Position', [1 1 1000 500]);
set(gcf, 'PaperPositionMode', 'auto');