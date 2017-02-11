function [] = SSAAnalyseSeq(Syll, BinSize, FileInfo, TitleString)

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
        NextSyll = find((FileInfo.NoteOnsets{i} >= SyllOffsetTime(SyllMatchNo)) & (FileInfo.NoteOnsets{i} <= (SyllOffsetTime(SyllMatchNo) + 0.5)), 1, 'first');
        
        if (~isempty(NextSyll))
            SyllLabels(SyllMatchNo,:) = [Syll FileInfo.NoteLabels{i}(NextSyll)];
            NextSyllOnsetTime(SyllMatchNo) = FileInfo.NoteOnsets{i}(NextSyll);
            NextSyllOffsetTime(SyllMatchNo) = FileInfo.NoteOffsets{i}(NextSyll);
        else
            SyllLabels(SyllMatchNo,:) = [Syll ' '];
            NextSyllOnsetTime(SyllMatchNo) = -100;
            NextSyllOffsetTime(SyllMatchNo) = -100;
        end
        SyllMatchNo = SyllMatchNo + 1;
    end
end
Edges = -0.6:BinSize:1;
UniqueNextSylls = unique(SyllLabels(:,end));
RasterIncrement = 0;
RasterFig = figure;
RasterPlotAxes = axes('Position', [0.1 0.4 0.85 0.55]);
PSTPlotAxes = axes('Position', [0.1 0.1 0.85 0.25]);
for i = 1:length(UniqueNextSylls),
    SyllMatches = find(SyllLabels(:,end) == UniqueNextSylls(i));
    if (length(SyllMatches) < 5)
        continue;
    end
    [SortedLengths, SortedIndices] = sort(NextSyllOnsetTime(SyllMatches) - SyllOnsetTime(SyllMatches));
    SyllMatches = SyllMatches(SortedIndices);
    Raster = [];
    NextSyllOnsetTimeRaster = [];
    NextSyllOffsetTimeRaster = [];    
    PST = [];
    for j = 1:length(SyllMatches),
        if (NextSyllOffsetTime(SyllMatches(j)) > 0)
            SpikeTimeIndices = find((FileInfo.SpikeData.Times{SyllFileIndex(SyllMatches(j))} >= (SyllOnsetTime(SyllMatches(j)) - 0.1)) & (FileInfo.SpikeData.Times{SyllFileIndex(SyllMatches(j))} <= (NextSyllOffsetTime(SyllMatches(j)) + 0.1)));
            SpikeTimes = FileInfo.SpikeData.Times{SyllFileIndex(SyllMatches(j))}(SpikeTimeIndices) - SyllOnsetTime(SyllMatches(j));
            NextSyllOnsetTimeRaster = [NextSyllOnsetTimeRaster; [(NextSyllOnsetTime(SyllMatches(j)) - SyllOnsetTime(SyllMatches(j))) j]];
            NextSyllOffsetTimeRaster = [NextSyllOffsetTimeRaster; [(NextSyllOffsetTime(SyllMatches(j)) - SyllOnsetTime(SyllMatches(j))) j]];
        else
            SpikeTimeIndices = find((FileInfo.SpikeData.Times{SyllFileIndex(SyllMatches(j))} >= (SyllOnsetTime(SyllMatches(j)) - 0.1)) & (FileInfo.SpikeData.Times{SyllFileIndex(SyllMatches(j))} <= (SyllOffsetTime(SyllMatches(j)) + 0.2)));
            SpikeTimes = FileInfo.SpikeData.Times{SyllFileIndex(SyllMatches(j))}(SpikeTimeIndices) - SyllOnsetTime(SyllMatches(j));            
        end
        if (size(SpikeTimes, 1) < size(SpikeTimes,2))
            SpikeTimes = SpikeTimes';
        end
        Raster = [Raster; [SpikeTimes ones(size(SpikeTimes))*j]];
        if (isempty(SpikeTimes))
            PST(j,:) = zeros(1, length(Edges));
        else
            PST(j,:) = histc(SpikeTimes, Edges);
        end
    end
    PST = PST/BinSize;
    disp([SyllLabels(SyllMatches(end),:), ' - ', ColNames{(mod(i,length(Cols)) + 1)}]);
    figure(RasterFig);
    axes(RasterPlotAxes); 
    hold on;
    if (~isempty(Raster))
        PlotRaster(Raster, Cols(mod(i,length(Cols)) + 1), 1, RasterIncrement);
        PlotRaster([(SyllOffsetTime(SyllMatches) - SyllOnsetTime(SyllMatches))' [1:1:length(SyllMatches)]'], 'k', 2, RasterIncrement);
        if (~isempty(NextSyllOnsetTimeRaster))
            PlotRaster(NextSyllOnsetTimeRaster, 'k', 2, RasterIncrement);
            PlotRaster(NextSyllOffsetTimeRaster, 'k', 2, RasterIncrement);
        end
        text(Edges(1)+0.02, (j/2 + RasterIncrement), SyllLabels(SyllMatches(end),:), 'FontSize', 12, 'FontWeight', 'bold', 'Color', Colors(mod(i, length(Cols)) + 1,:));
        RasterIncrement = RasterIncrement + j + 2;
        axis tight;
    end
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
temp(1) = Edges(1);
if (temp(2) < 0)
    temp(2) = 0.1;
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


Edges = -(mean(SyllOffsetTime - SyllOnsetTime)+0.6):BinSize:1;
UniqueNextSylls = unique(SyllLabels(:,end));
RasterIncrement = 0;
RasterFig2 = figure;
RasterPlotAxes2 = axes('Position', [0.1 0.4 0.85 0.55]);
PSTPlotAxes2 = axes('Position', [0.1 0.1 0.85 0.25]);
for i = 1:length(UniqueNextSylls),
    SyllMatches = find(SyllLabels(:,end) == UniqueNextSylls(i));
    if (length(SyllMatches) < 5)
        continue;
    end
    [SortedLengths, SortedIndices] = sort(NextSyllOnsetTime(SyllMatches) - SyllOffsetTime(SyllMatches));
    SyllMatches = SyllMatches(SortedIndices);
    Raster = [];
    NextSyllOnsetTimeRaster = [];
    NextSyllOffsetTimeRaster = [];    
    PST = [];
    for j = 1:length(SyllMatches),
        if (NextSyllOffsetTime(SyllMatches(j)) > 0)
            SpikeTimeIndices = find((FileInfo.SpikeData.Times{SyllFileIndex(SyllMatches(j))} >= (SyllOnsetTime(SyllMatches(j)) - 0.1)) & (FileInfo.SpikeData.Times{SyllFileIndex(SyllMatches(j))} <= (NextSyllOffsetTime(SyllMatches(j)) + 0.1)));
            SpikeTimes = FileInfo.SpikeData.Times{SyllFileIndex(SyllMatches(j))}(SpikeTimeIndices) - SyllOffsetTime(SyllMatches(j));
            NextSyllOnsetTimeRaster = [NextSyllOnsetTimeRaster; [(NextSyllOnsetTime(SyllMatches(j)) - SyllOffsetTime(SyllMatches(j))) j]];
            NextSyllOffsetTimeRaster = [NextSyllOffsetTimeRaster; [(NextSyllOffsetTime(SyllMatches(j)) - SyllOffsetTime(SyllMatches(j))) j]];
        else
            SpikeTimeIndices = find((FileInfo.SpikeData.Times{SyllFileIndex(SyllMatches(j))} >= (SyllOnsetTime(SyllMatches(j)) - 0.1)) & (FileInfo.SpikeData.Times{SyllFileIndex(SyllMatches(j))} <= (SyllOffsetTime(SyllMatches(j)) + 0.2)));
            SpikeTimes = FileInfo.SpikeData.Times{SyllFileIndex(SyllMatches(j))}(SpikeTimeIndices) - SyllOffsetTime(SyllMatches(j));            
        end
        if (size(SpikeTimes, 1) < size(SpikeTimes,2))
            SpikeTimes = SpikeTimes';
        end
        Raster = [Raster; [SpikeTimes ones(size(SpikeTimes))*j]];
        if (isempty(SpikeTimes))
            PST(j,:) = zeros(1, length(Edges));
        else
            PST(j,:) = histc(SpikeTimes, Edges);
        end
    end
    PST = PST/BinSize;
    disp([SyllLabels(SyllMatches(end),:), ' - ', ColNames{(mod(i,length(Cols)) + 1)}]);
    figure(RasterFig2);
    axes(RasterPlotAxes2); 
    hold on;
    if (~isempty(Raster))
        PlotRaster(Raster, Cols(mod(i,length(Cols)) + 1), 1, RasterIncrement);
        PlotRaster([(SyllOnsetTime(SyllMatches) - SyllOffsetTime(SyllMatches))' [1:1:length(SyllMatches)]'], 'k', 2, RasterIncrement);
        if (~isempty(NextSyllOnsetTimeRaster))
            PlotRaster(NextSyllOnsetTimeRaster, 'k', 2, RasterIncrement);
            PlotRaster(NextSyllOffsetTimeRaster, 'k', 2, RasterIncrement);
        end
        text(Edges(1)+0.02, (j/2 + RasterIncrement), SyllLabels(SyllMatches(end),:), 'FontSize', 12, 'FontWeight', 'bold', 'Color', Colors(mod(i, length(Cols)) + 1,:));
        RasterIncrement = RasterIncrement + j + 2;
        axis tight;
    end
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
temp(1) = Edges(1);
if (temp(2) < 0)
    temp(2) = 0.1;
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