function [] = SSAAnalyseConvSeq(Syll, BinSize, FileInfo, TitleString)

Cols = ['rgbcm'];
ColNames{1} = 'red';
ColNames{2} = 'green';
ColNames{3} = 'blue';
ColNames{4} = 'cyan';
ColNames{5} = 'magenta';
Colors = [1 0 0; 0 1 0; 0 0 1; 0 1 1; 1 0 1];

for i = 1:length(FileInfo.NoteLabels),
    Matches{i} = strfind(FileInfo.NoteLabels{i}, Syll);
end

SyllMatchNo = 1;

for i = 1:length(Matches),
    for j = 1:length(Matches{i}),
        SyllFileIndex(SyllMatchNo) = i;
        SyllOnsetTime(SyllMatchNo) = FileInfo.NoteOnsets{i}(Matches{i}(j));
        SyllOffsetTime(SyllMatchNo) = FileInfo.NoteOffsets{i}(Matches{i}(j)+length(Syll)-1);
        PrevSyll = find((FileInfo.NoteOnsets{i} < SyllOnsetTime(SyllMatchNo)) & (FileInfo.NoteOnsets{i} >= (SyllOnsetTime(SyllMatchNo) - 0.5)), 1, 'last');
        
        if (~isempty(PrevSyll))
            SyllLabels(SyllMatchNo,:) = [FileInfo.NoteLabels{i}(PrevSyll) Syll];
            PrevSyllOnsetTime(SyllMatchNo) = FileInfo.NoteOnsets{i}(PrevSyll);
            PrevSyllOffsetTime(SyllMatchNo) = FileInfo.NoteOffsets{i}(PrevSyll);
        else
            SyllLabels(SyllMatchNo,:) = [' ' Syll];
            PrevSyllOnsetTime(SyllMatchNo) = -100;
            PrevSyllOffsetTime(SyllMatchNo) = -100;
        end
        SyllMatchNo = SyllMatchNo + 1;
    end
end
Edges = -1:BinSize:mean(SyllOffsetTime - SyllOnsetTime) + 0.1;
UniquePrevSylls = unique(SyllLabels(:,1));
RasterIncrement = 0;
RasterFig = figure;
RasterPlotAxes = axes('Position', [0.1 0.4 0.85 0.55]);
PSTPlotAxes = axes('Position', [0.1 0.1 0.85 0.25]);
for i = 1:length(UniquePrevSylls),
    SyllMatches = find(SyllLabels(:,1) == UniquePrevSylls(i));
    if (length(SyllMatches) < 5)
        continue;
    end
    [SortedLengths, SortedIndices] = sort(PrevSyllOnsetTime(SyllMatches) - SyllOnsetTime(SyllMatches));
    SyllMatches = SyllMatches(SortedIndices);
    Raster = [];
    PrevSyllOnsetTimeRaster = [];
    PrevSyllOffsetTimeRaster = [];    
    PST = [];
    for j = 1:length(SyllMatches),
        if (PrevSyllOffsetTime(SyllMatches(j)) > 0)
            SpikeTimeIndices = find((FileInfo.SpikeData.Times{SyllFileIndex(SyllMatches(j))} <= (SyllOffsetTime(SyllMatches(j)) + 0.1)) & (FileInfo.SpikeData.Times{SyllFileIndex(SyllMatches(j))} >= (PrevSyllOnsetTime(SyllMatches(j)) - 0.1)));
            SpikeTimes = FileInfo.SpikeData.Times{SyllFileIndex(SyllMatches(j))}(SpikeTimeIndices) - SyllOnsetTime(SyllMatches(j));
            PrevSyllOnsetTimeRaster = [PrevSyllOnsetTimeRaster; [(PrevSyllOnsetTime(SyllMatches(j)) - SyllOnsetTime(SyllMatches(j))) j]];
            PrevSyllOffsetTimeRaster = [PrevSyllOffsetTimeRaster; [(PrevSyllOffsetTime(SyllMatches(j)) - SyllOnsetTime(SyllMatches(j))) j]];
        else
            SpikeTimeIndices = find((FileInfo.SpikeData.Times{SyllFileIndex(SyllMatches(j))} <= (SyllOffsetTime(SyllMatches(j)) + 0.1)) & (FileInfo.SpikeData.Times{SyllFileIndex(SyllMatches(j))} >= (SyllOnsetTime(SyllMatches(j)) - 0.2)));
            SpikeTimes = FileInfo.SpikeData.Times{SyllFileIndex(SyllMatches(j))}(SpikeTimeIndices) - SyllOnsetTime(SyllMatches(j));            
        end
        if (size(SpikeTimes, 1) < size(SpikeTimes,2))
            SpikeTimes = SpikeTimes';
        end
        Raster = [Raster; [SpikeTimes ones(size(SpikeTimes))*j]];
        PST(j,:) = histc(SpikeTimes, Edges);
    end
    PST = PST/BinSize;
    disp([SyllLabels(SyllMatches(end),:), ' - ', ColNames{(mod(i,length(Cols)) + 1)}]);
    figure(RasterFig);
    axes(RasterPlotAxes); 
    hold on;
    if (~isempty(Raster))
        PlotRaster(Raster, Cols(mod(i,length(Cols)) + 1), 1, RasterIncrement);
    end
    PlotRaster([(SyllOffsetTime(SyllMatches) - SyllOnsetTime(SyllMatches))' [1:1:length(SyllMatches)]'], 'k', 2, RasterIncrement);
    if (~isempty(PrevSyllOnsetTimeRaster))
        PlotRaster(PrevSyllOnsetTimeRaster, 'k', 2, RasterIncrement);
        PlotRaster(PrevSyllOffsetTimeRaster, 'k', 2, RasterIncrement);
    end
    text(Edges(end)-0.02, (j/2 + RasterIncrement), SyllLabels(SyllMatches(end),:), 'FontSize', 12, 'FontWeight', 'bold', 'Color', Colors(mod(i, length(Cols)) + 1,:));
    RasterIncrement = RasterIncrement + j + 2;
    axis tight;
    axes(PSTPlotAxes);
    plot(Edges, mean(PST), Cols(mod(i,length(Cols)) + 1));
    hold on;
end
figure(RasterFig);
set(gcf, 'Position', [360 100 650 600]);
set(gcf, 'Color', 'w');
axes(RasterPlotAxes); 
set(gca, 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Trial #', 'FontSize', 14, 'FontWeight', 'bold');
set(gca, 'XColor', 'w');
axis tight;
temp = axis;
temp(2) = Edges(end);
if (temp(1) > 0)
    temp(1) = -0.1;
end
axis(temp);
plot([0 0], [temp(3) temp(4)], 'k--', 'LineWidth', 2);
axes(PSTPlotAxes); 
set(gca, 'FontSize', 12, 'FontWeight', 'bold');
set(gca, 'Box', 'off');
xlabel('Time (sec)', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('Firing Rate (Hz)', 'FontSize', 14, 'FontWeight', 'bold');
axis auto;
temp2 = axis;
temp2(1:2) = temp(1:2);
axis(temp2);
plot([0 0], [temp(3) temp(4)], 'k--', 'LineWidth', 2);
axes(RasterPlotAxes);
title([TitleString, '- Aligned to syllable onset'], 'FontSize', 14, 'FontWeight', 'bold');


Edges = -1:BinSize:(mean(SyllOffsetTime - SyllOnsetTime)+0.1);
UniquePrevSylls = unique(SyllLabels(:,1));
RasterIncrement = 0;
RasterFig2 = figure;
RasterPlotAxes2 = axes('Position', [0.1 0.4 0.85 0.55]);
PSTPlotAxes2 = axes('Position', [0.1 0.1 0.85 0.25]);
for i = 1:length(UniquePrevSylls),
    SyllMatches = find(SyllLabels(:,1) == UniquePrevSylls(i));
    if (length(SyllMatches) < 5)
        continue;
    end
    [SortedLengths, SortedIndices] = sort(PrevSyllOnsetTime(SyllMatches) - SyllOffsetTime(SyllMatches));
    SyllMatches = SyllMatches(SortedIndices);
    Raster = [];
    PrevSyllOnsetTimeRaster = [];
    PrevSyllOffsetTimeRaster = [];    
    PST = [];
    for j = 1:length(SyllMatches),
        if (PrevSyllOffsetTime(SyllMatches(j)) > 0)
            SpikeTimeIndices = find((FileInfo.SpikeData.Times{SyllFileIndex(SyllMatches(j))} <= (SyllOffsetTime(SyllMatches(j)) + 0.1)) & (FileInfo.SpikeData.Times{SyllFileIndex(SyllMatches(j))} >= (PrevSyllOnsetTime(SyllMatches(j)) - 0.1)));
            SpikeTimes = FileInfo.SpikeData.Times{SyllFileIndex(SyllMatches(j))}(SpikeTimeIndices) - SyllOffsetTime(SyllMatches(j));
            PrevSyllOnsetTimeRaster = [PrevSyllOnsetTimeRaster; [(PrevSyllOnsetTime(SyllMatches(j)) - SyllOffsetTime(SyllMatches(j))) j]];
            PrevSyllOffsetTimeRaster = [PrevSyllOffsetTimeRaster; [(PrevSyllOffsetTime(SyllMatches(j)) - SyllOffsetTime(SyllMatches(j))) j]];
        else
            SpikeTimeIndices = find((FileInfo.SpikeData.Times{SyllFileIndex(SyllMatches(j))} <= (SyllOffsetTime(SyllMatches(j)) + 0.1)) & (FileInfo.SpikeData.Times{SyllFileIndex(SyllMatches(j))} >= (SyllOnsetTime(SyllMatches(j)) - 0.2)));
            SpikeTimes = FileInfo.SpikeData.Times{SyllFileIndex(SyllMatches(j))}(SpikeTimeIndices) - SyllOffsetTime(SyllMatches(j));            
        end
        if (size(SpikeTimes, 1) < size(SpikeTimes,2))
            SpikeTimes = SpikeTimes';
        end
        Raster = [Raster; [SpikeTimes ones(size(SpikeTimes))*j]];
        PST(j,:) = histc(SpikeTimes, Edges);
    end
    PST = PST/BinSize;
    disp([SyllLabels(SyllMatches(end),:), ' - ', ColNames{(mod(i,length(Cols)) + 1)}]);
    figure(RasterFig2);
    axes(RasterPlotAxes2); 
    hold on;
    if (~isempty(Raster))
        PlotRaster(Raster, Cols(mod(i,length(Cols)) + 1), 1, RasterIncrement);
    end
    
    PlotRaster([(SyllOnsetTime(SyllMatches) - SyllOffsetTime(SyllMatches))' [1:1:length(SyllMatches)]'], 'k', 2, RasterIncrement);
    if (~isempty(PrevSyllOnsetTimeRaster))
        PlotRaster(PrevSyllOnsetTimeRaster, 'k', 2, RasterIncrement);
        PlotRaster(PrevSyllOffsetTimeRaster, 'k', 2, RasterIncrement);
    end
    text(Edges(end)-0.02, (j/2 + RasterIncrement), SyllLabels(SyllMatches(end),:), 'FontSize', 12, 'FontWeight', 'bold', 'Color', Colors(mod(i, length(Cols)) + 1,:));
    RasterIncrement = RasterIncrement + j + 2;
    axis tight;
    axes(PSTPlotAxes2);
    plot(Edges, mean(PST), Cols(mod(i,length(Cols)) + 1));
    hold on;
end
figure(RasterFig2);
set(gcf, 'Position', [360 100 650 600]);
set(gcf, 'Color', 'w');
axes(RasterPlotAxes2); 
set(gca, 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Trial #', 'FontSize', 14, 'FontWeight', 'bold');
set(gca, 'XColor', 'w');
axis tight;
temp = axis;
temp(2) = Edges(end);
if (temp(1) > 0)
    temp(1) = -0.1;
end
axis(temp);
plot([0 0], [temp(3) temp(4)], 'k--', 'LineWidth', 2);
axes(PSTPlotAxes2); 
set(gca, 'FontSize', 12, 'FontWeight', 'bold');
set(gca, 'Box', 'off');
xlabel('Time (sec)', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('Firing Rate (Hz)', 'FontSize', 14, 'FontWeight', 'bold');
axis auto;
temp2 = axis;
temp2(1:2) = temp(1:2);
axis(temp2);
plot([0 0], [temp(3) temp(4)], 'k--', 'LineWidth', 2);
axes(RasterPlotAxes2);
title([TitleString, '- Aligned to syllable offset'], 'FontSize', 14, 'FontWeight', 'bold');


disp('Finished analysing sequence');