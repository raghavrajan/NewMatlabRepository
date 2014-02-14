function [RasterPSTFigure] = SSAPlotRasterPST(DirectoryName, RecFileDirectoryName, FileType, DirFileInfo, UnDirFileInfo, MedianMotif, Edges)

RasterPSTFigure = figure;
set(gcf, 'Color', 'w');
set(gcf, 'Position', [360 260 575 650]);
PSTAxis = axes('Position', [0.1 0.1 0.65 0.15]);
UnDirRasterAxis = axes('Position', [0.1 0.3 0.65 0.2]);
DirRasterAxis = axes('Position', [0.1 0.55 0.65 0.2]);
SpectrogramAxis = axes('Position', [0.1 0.8 0.65 0.15]);

DirISIAxis = axes('Position', [0.825 0.65 0.15 0.1]);
UnDirISIAxis = axes('Position', [0.825 0.4 0.15 0.1]);
WaveformAxis = axes('Position', [0.825 0.1 0.15 0.15]);

MaxPST = 0;
MaxWaveform = 0;
MaxISI = 0;

if (~isempty(DirFileInfo))
    axes(PSTAxis);
    plot(Edges, mean(DirFileInfo.PST),'r-', 'LineWidth', 2);
    hold on;
    MaxPST = max(mean(DirFileInfo.PST)) * 1.1;
    if (Edges(end) > MedianMotif.Length)
        axis([Edges(1) Edges(end) 0 MaxPST]);
    else
        axis([Edges(1) MedianMotif.Length 0 MaxPST]);
    end
    
    axes(DirRasterAxis);
    
    if (~isempty(UnDirFileInfo))
        MinRaster = min(DirFileInfo.SpikeRaster(end,2), UnDirFileInfo.SpikeRaster(end,2));
    else
        MinRaster = DirFileInfo.SpikeRaster(end,2);
    end
    
    PlotRaster(DirFileInfo.SpikeRaster, 'r', 1, 0, MinRaster);
    if (Edges(end) > MedianMotif.Length)
        axis([Edges(1) Edges(end) 0 (MinRaster + 1)]);
    else
        axis([Edges(1) MedianMotif.Length 0 (MinRaster + 1)]);
    end
    axes(DirISIAxis);
    ISIEdges = 0:0.1:1000;
    SpikeTimes = [];
    for i = 1:length(DirFileInfo.SpikeData.Times),
        SpikeTimes = [SpikeTimes; DirFileInfo.SpikeData.Times{i}];
    end
    SpikeTimes = SpikeTimes * 1000;
    DirISIs = histc(diff(SpikeTimes), ISIEdges);
    DirISIs = DirISIs/sum(DirISIs) * 100;
    DirISIBar = bar(ISIEdges, DirISIs, 'histc');
    MaxISI = max(DirISIs);
    
    axes(WaveformAxis);
    if (size(DirFileInfo.SpikeData.SongSpikeWaveforms,1) > 100)
        plot((1:1:(size(DirFileInfo.SpikeData.SongSpikeWaveforms,2)))/32, DirFileInfo.SpikeData.SongSpikeWaveforms(round(rand(100,1)*99 + 1),:), 'r');
    else
        plot((1:1:(size(DirFileInfo.SpikeData.SongSpikeWaveforms,2)))/32, DirFileInfo.SpikeData.SongSpikeWaveforms, 'r');
    end
    hold on;
    MaxWaveform = max(max(DirFileInfo.SpikeData.SongSpikeWaveforms)) * 1.1;
    axis tight;
end

if (~isempty(UnDirFileInfo))
    axes(PSTAxis);
    plot(Edges, mean(UnDirFileInfo.PST),'b-', 'LineWidth', 2);
    hold on;
    if (MaxPST < (max(mean(UnDirFileInfo.PST)) * 1.1))
        MaxPST = max(mean(UnDirFileInfo.PST)) * 1.1;
    end
    if (Edges(end) > MedianMotif.Length)
        axis([Edges(1) Edges(end) 0 MaxPST]);
    else
        axis([Edges(1) MedianMotif.Length 0 MaxPST]);
    end
    
    axes(UnDirRasterAxis);
    
    if (~isempty(DirFileInfo))
        MinRaster = min(DirFileInfo.SpikeRaster(end,2), UnDirFileInfo.SpikeRaster(end,2));
    else
        MinRaster = UnDirFileInfo.SpikeRaster(end,2);
    end
    
    PlotRaster(UnDirFileInfo.SpikeRaster, 'b', 1, 0, MinRaster);
    
	if (Edges(end) > MedianMotif.Length)
        axis([Edges(1) Edges(end) 0 (MinRaster + 1)]);
    else
        axis([Edges(1) MedianMotif.Length 0 (MinRaster + 1)]);
    end
    
    axes(UnDirISIAxis);
    ISIEdges = 0:0.1:1000;
    SpikeTimes = [];
    for i = 1:length(UnDirFileInfo.SpikeData.Times),
        SpikeTimes = [SpikeTimes; UnDirFileInfo.SpikeData.Times{i}];
    end
    SpikeTimes = SpikeTimes * 1000;    
    UnDirISIs = histc(diff(SpikeTimes), ISIEdges);
    UnDirISIs = UnDirISIs/sum(UnDirISIs) * 100;
    UnDirISIBar = bar(ISIEdges, UnDirISIs, 'histc');
    if (MaxISI < max(UnDirISIs))
        MaxISI = max(UnDirISIs);
    end
    
    axes(WaveformAxis);
    %hold off;
    if (size(UnDirFileInfo.SpikeData.SongSpikeWaveforms,1) > 100)
        plot(((1:1:(size(UnDirFileInfo.SpikeData.SongSpikeWaveforms,2)))/32 + 2), UnDirFileInfo.SpikeData.SongSpikeWaveforms(round(rand(100,1)*99 + 1),:), 'b');
    else
        plot(((1:1:(size(UnDirFileInfo.SpikeData.SongSpikeWaveforms,2)))/32 + 2), UnDirFileInfo.SpikeData.SongSpikeWaveforms, 'b');
    end
    if (MaxWaveform < max(max(UnDirFileInfo.SpikeData.SongSpikeWaveforms)));
        MaxWaveform = max(max(UnDirFileInfo.SpikeData.SongSpikeWaveforms)) * 1.1;
    end
    axis tight;
end

SSAPlotMotifSpectrogram(DirectoryName, RecFileDirectoryName, MedianMotif, RasterPSTFigure, SpectrogramAxis, FileType);

axes(PSTAxis);
set(gca, 'FontSize', 12, 'FontWeight', 'bold');
xlabel('Time (sec)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Firing Rate (Hz)', 'FontSize', 12, 'FontWeight', 'bold');
set(gca, 'Box', 'off');
if ((~isempty(DirFileInfo)) && (~isempty(UnDirFileInfo)))
    PSTLegend = legend('DIR', 'UNDIR');
else
    if (~isempty(DirFileInfo))
        PSTLegend = legend('DIR');
    else
        PSTLegend = legend('UNDIR');
    end
end
set(PSTLegend, 'Box', 'off');
set(gca, 'TickDir', 'out');

axes(DirRasterAxis);
set(gca, 'XColor', 'w');
set(gca, 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Directed Trial #', 'FontSize', 12, 'FontWeight', 'bold', 'Color', 'b');
set(gca, 'TickDir', 'out');

axes(UnDirRasterAxis);
set(gca, 'XColor', 'w');
set(gca, 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Undirected Trial #', 'FontSize', 12, 'FontWeight', 'bold', 'Color', 'k');
set(gca, 'TickDir', 'out');

axes(SpectrogramAxis);
if (Edges(end) > MedianMotif.Length)
    axis([Edges(1) Edges(end) 300 10000]);
else
    axis([Edges(1) MedianMotif.Length 300 10000]);
end
set(gca, 'Visible', 'on');
set(gca, 'XColor', 'w');
set(gca, 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Frequency (KHz)', 'FontSize', 12, 'FontWeight', 'bold', 'Color', 'k');
set(gca, 'YTick', [300 5000 10000], 'YTickLabel', [0.3 5 10]);
set(gca, 'TickDir', 'out');

axes(DirISIAxis);
axis auto;
temp = axis;
temp(1:2) = [0 1000];
temp(3) = 0;
temp(4) = MaxISI * 1.2;
axis(temp);
set(gca, 'Box', 'off');
set(gca, 'FontSize', 10, 'FontWeight', 'bold');
xlabel('ISI (ms)', 'FontSize', 10, 'FontWeight', 'bold', 'Color', 'k');
ylabel('Percent (%)', 'FontSize', 10, 'FontWeight', 'bold', 'Color', 'k');
set(gca, 'XTick', [0.1 1 10 1000], 'XTickLabel', [0.1 1 10 1000]);
title('Directed', 'FontSize', 10, 'FontWeight', 'bold', 'Color', 'b');
set(gca, 'TickDir', 'out');
set(gca, 'XScale', 'log');

axes(UnDirISIAxis);
axis auto;
temp = axis;
temp(1:2) = [0 1000];
temp(3) = 0;
temp(4) = MaxISI * 1.2;
axis(temp);
set(gca, 'Box', 'off');
set(gca, 'FontSize', 10, 'FontWeight', 'bold');
xlabel('ISI (ms)', 'FontSize', 10, 'FontWeight', 'bold', 'Color', 'k');
ylabel('Percent (%)', 'FontSize', 10, 'FontWeight', 'bold', 'Color', 'k');
set(gca, 'XTick', [0.1 1 10 1000], 'XTickLabel', [0.1 1 10 1000]);
title('Undirected', 'FontSize', 10, 'FontWeight', 'bold', 'Color', 'k');
set(gca, 'TickDir', 'out');
set(gca, 'XScale', 'log');

axes(WaveformAxis);
set(gca, 'FontSize', 10, 'FontWeight', 'bold');
xlabel('Time (ms)', 'FontSize', 10, 'FontWeight', 'bold', 'Color', 'k');
ylabel('Amplitude (uV)', 'FontSize', 10, 'FontWeight', 'bold', 'Color', 'k');
MaxWaveform = floor(MaxWaveform/100) * 100;
set(gca, 'YTick', [-MaxWaveform 0 MaxWaveform]);
set(gca, 'Box', 'off');
set(gca, 'TickDir', 'out');

if (~isempty(DirFileInfo))
    if (strfind(FileType, 'observer'))
        DotIndex = find(DirFileInfo.FileNames{1} == '.');
        TitleName = DirFileInfo.FileNames{1}(1:(DotIndex(1) - 1));
    else
        if (strfind(FileType, 'okrank'))
            TitleName = DirFileInfo.FileNames{1}(1:(end - 6));
        end
    end
else
    if (~isempty(UnDirFileInfo))
        if (strfind(FileType, 'observer'))
            DotIndex = find(UnDirFileInfo.FileNames{1} == '.');
            TitleName = UnDirFileInfo.FileNames{1}(1:(DotIndex(1) - 1));
        else
            if (strfind(FileType, 'okrank'))
                TitleName = UnDirFileInfo.FileNames{1}(1:(end - 6));
            end
        end
    end
end
axes(SpectrogramAxis);
title(TitleName, 'FontSize', 14, 'FontWeight', 'bold');