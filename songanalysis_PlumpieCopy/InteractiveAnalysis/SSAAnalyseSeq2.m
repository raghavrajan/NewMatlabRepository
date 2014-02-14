function [] = SSAAnalyseSeq2(Syll, BinSize, FileInfo, TitleString, Motif)

InterBoutInterval = 0.5; % seconds
Cols = ['rgbcm'];
ColNames{1} = 'red';
ColNames{2} = 'green';
ColNames{3} = 'blue';
ColNames{4} = 'cyan';
ColNames{5} = 'magenta';
Colors = [1 0 0; 0 1 0; 0 0 1; 0 1 1; 1 0 1];

BoutNo = 1;
for i = 1:length(FileInfo.NoteLabels),
    Notes.onsets = FileInfo.NoteOnsets{i};
    Notes.offsets = FileInfo.NoteOffsets{i};
    Notes.labels = FileInfo.NoteLabels{i};
    
    BoutIndices = find((Notes.onsets(2:end) - Notes.offsets(1:end-1)) > InterBoutInterval);
    Bouts = [];
    if (~isempty(BoutIndices))
        if (Notes.onsets(1) > InterBoutInterval)
            Bouts(1,:) = [1 BoutIndices(1)];
            for i = 1:(length(BoutIndices)-1),
                Bouts(i+1,:) = [(BoutIndices(i)+1) (BoutIndices(i+1))];
            end
        else
            for i = 1:(length(BoutIndices)-1),
                Bouts(i,:) = [(BoutIndices(i)+1) (BoutIndices(i+1))];
            end
        end
        Bouts(end+1,:) = [(BoutIndices(end) + 1) length(Notes.labels)];
    else
        if (Notes.onsets(1) > InterBoutInterval)
            Bouts = [1 length(Notes.labels)];
        end
    end
    for i = 1:size(Bouts,1),
        Motifs = (strfind(Notes.labels(Bouts(i,1):Bouts(i,2)), Motif));
        if (isempty(Motifs))
            BoutStats.Motifs.NoofMotifs(BoutNo) = 0;
            continue;
        else
            Motifs = Motifs + Bouts(i,1) - 1;
            IntroNotes = find(Notes.labels(Bouts(i,1):Motifs(1)) == 'i');
            IntroNotes = IntroNotes + Bouts(i,1) - 1;
            if (isempty(IntroNotes))
                BoutStats.IntroNotes.NoofIntroNotes(BoutNo) = 0;
            else
                BoutStats.IntroNotes.NoofIntroNotes(BoutNo) = length(IntroNotes);
            end
            BoutStats.Motifs.NoofMotifs(BoutNo) = length(Motifs);
        end
        BoutNo = BoutNo + 1;
    end
end

Edges = -1:BinSize:1;
MaxIntroNotes = max(BoutStats.IntroNotes.NoofIntroNotes);
for i = 1:MaxIntroNotes,
    BoutNo(i) = 1;
    OnsetRaster{i} = [];
    OffsetRaster{i} = [];
    SyllOnsetRaster{i} = [];
    SyllOffsetRaster{i} = [];
    OnsetNextSyllOnsetRaster{i} = [];
    OffsetNextSyllOnsetRaster{i} = [];
end

for FileNo = 1:length(FileInfo.NoteLabels),
    Notes.onsets = FileInfo.NoteOnsets{FileNo};
    Notes.offsets = FileInfo.NoteOffsets{FileNo};
    Notes.labels = FileInfo.NoteLabels{FileNo};
    
    BoutIndices = find((Notes.onsets(2:end) - Notes.offsets(1:end-1)) > InterBoutInterval);
    Bouts = [];
    if (~isempty(BoutIndices))
        if (Notes.onsets(1) > InterBoutInterval)
            Bouts(1,:) = [1 BoutIndices(1)];
            for i = 1:(length(BoutIndices)-1),
                Bouts(i+1,:) = [(BoutIndices(i)+1) (BoutIndices(i+1))];
            end
        else
            for i = 1:(length(BoutIndices)-1),
                Bouts(i,:) = [(BoutIndices(i)+1) (BoutIndices(i+1))];
            end
        end
        Bouts(end+1,:) = [(BoutIndices(end) + 1) length(Notes.labels)];
    else
        if (Notes.onsets(1) > InterBoutInterval)
            Bouts = [1 length(Notes.labels)];
        end
    end
    for i = 1:size(Bouts,1),
        Motifs = (strfind(Notes.labels(Bouts(i,1):Bouts(i,2)), Motif));
        if (isempty(Motifs))
            BoutStats.Motifs.NoofMotifs(BoutNo) = 0;
            continue;
        else
            Motifs = Motifs + Bouts(i,1) - 1;
            if (Bouts(i,1) > (Motifs(1) - 8))
                IntroNotes = find(Notes.labels(Bouts(i,1):Motifs(1)) == 'i');
                IntroNotes = IntroNotes + Bouts(i,1) - 1;
            else
                IntroNotes = find(Notes.labels((Motifs(1) - 8):Motifs(1)) == 'i');
                IntroNotes = IntroNotes + Motifs(1) - 8 - 1;
            end
            if (~isempty(IntroNotes))
                IntroNotes = sort(IntroNotes, 'descend');
                for j = 1:length(IntroNotes),
                    SyllOnsetRaster{j} = [SyllOnsetRaster{j}; [(Notes.onsets(IntroNotes(j)) - Notes.offsets(IntroNotes(j))) BoutNo(j)]];
                    SyllOffsetRaster{j} = [SyllOffsetRaster{j}; [(Notes.offsets(IntroNotes(j)) - Notes.onsets(IntroNotes(j))) BoutNo(j)]];
                    OnsetNextSyllOnsetRaster{j} = [OnsetNextSyllOnsetRaster{j}; [(Notes.onsets(IntroNotes(j)+1) - Notes.onsets(IntroNotes(j))) BoutNo(j)]];
                    OffsetNextSyllOnsetRaster{j} = [OffsetNextSyllOnsetRaster{j}; [(Notes.onsets(IntroNotes(j)+1) - Notes.offsets(IntroNotes(j))) BoutNo(j)]];
                    SpikeTimeIndices = find((FileInfo.SpikeData.Times{FileNo} >= (Notes.onsets(IntroNotes(j)) - 0.1)) & (FileInfo.SpikeData.Times{FileNo} <= (Notes.offsets(IntroNotes(j)+1) + 0.2)));
                    SpikeTimes = FileInfo.SpikeData.Times{FileNo}(SpikeTimeIndices);
                    if (size(SpikeTimes, 1) < size(SpikeTimes,2))
                        SpikeTimes = SpikeTimes';
                    end
                    OnsetRaster{j} = [OnsetRaster{j}; [(SpikeTimes - Notes.onsets(IntroNotes(j))) ones(size(SpikeTimes))*BoutNo(j)]];
                    OffsetRaster{j} = [OffsetRaster{j}; [(SpikeTimes - Notes.offsets(IntroNotes(j))) ones(size(SpikeTimes))*BoutNo(j)]];
                    OnsetPST{j}(BoutNo(j),:) = histc((SpikeTimes - Notes.onsets(IntroNotes(j))), Edges);
                    OffsetPST{j}(BoutNo(j),:) = histc((SpikeTimes - Notes.offsets(IntroNotes(j))), Edges);
                    BoutNo(j) = BoutNo(j) + 1;
                end
            end
        end
    end
end
RasterFig = figure;
RasterPlotAxes = axes('Position', [0.1 0.4 0.85 0.55]);
PSTPlotAxes = axes('Position', [0.1 0.1 0.85 0.25]);

figure(RasterFig);
hold on;
RasterIncrement = 0;
for i = 1:length(OnsetRaster),
    if (~isempty(OnsetRaster{i}))
        axes(RasterPlotAxes);
        [SortedRaster, SortedIndices] = sort(OnsetNextSyllOnsetRaster{i}(:,1));
        for j = 1:length(SortedIndices),
            Indices = find(OnsetRaster{i}(:,2) == SortedIndices(j));
            OnsetRaster{i}(Indices,2) = j-1000;
            Indices = find(SyllOffsetRaster{i}(:,2) == SortedIndices(j));
            SyllOffsetRaster{i}(Indices,2) = j-1000;
            Indices = find(OnsetNextSyllOnsetRaster{i}(:,2) == SortedIndices(j));
            OnsetNextSyllOnsetRaster{i}(Indices,2) = j-1000;            
        end
        OnsetRaster{i}(:,2) = OnsetRaster{i}(:,2) + 1000;
        SyllOffsetRaster{i}(:,2) = SyllOffsetRaster{i}(:,2) + 1000;
        OnsetNextSyllOnsetRaster{i}(:,2) = OnsetNextSyllOnsetRaster{i}(:,2) + 1000;
        PlotRaster(OnsetRaster{i}, Cols(mod(i,length(Cols)) + 1), 1, RasterIncrement);
        PlotRaster(SyllOffsetRaster{i}, 'k', 2, RasterIncrement);
        PlotRaster(OnsetNextSyllOnsetRaster{i}, 'k', 2, RasterIncrement);
        RasterIncrement = RasterIncrement + max(OnsetRaster{i}(:,2)) + 2;
        axis tight;
        axes(PSTPlotAxes);
        OnsetPST{i} = OnsetPST{i}/BinSize;
        plot(Edges, mean(OnsetPST{i}), Cols(mod(i,length(Cols)) + 1));
        hold on;
    end
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
temp(1) = -0.1;
axis(temp);
axes(PSTPlotAxes); 
set(gca, 'FontSize', 12, 'FontWeight', 'bold');
set(gca, 'Box', 'off');
xlabel('Time (sec)', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('Firing Rate (Hz)', 'FontSize', 14, 'FontWeight', 'bold');
axis auto;
temp2 = axis;
temp2(1:2) = temp(1:2);
axis(temp2);
axes(RasterPlotAxes);
title([TitleString, '- Aligned to syllable onset'], 'FontSize', 14, 'FontWeight', 'bold');


RasterIncrement = 0;
RasterFig2 = figure;
RasterPlotAxes2 = axes('Position', [0.1 0.4 0.85 0.55]);
PSTPlotAxes2 = axes('Position', [0.1 0.1 0.85 0.25]);
figure(RasterFig2);
for i = 1:length(OffsetRaster),
    if (~isempty(OffsetRaster{i}))
        axes(RasterPlotAxes2); 
        [SortedRaster, SortedIndices] = sort(OffsetNextSyllOnsetRaster{i}(:,1));
        for j = 1:length(SortedIndices),
            Indices = find(OffsetRaster{i}(:,2) == SortedIndices(j));
            OffsetRaster{i}(Indices,2) = j - 1000;
            Indices = find(SyllOnsetRaster{i}(:,2) == SortedIndices(j));
            SyllOnsetRaster{i}(Indices,2) = j - 1000;
            Indices = find(OffsetNextSyllOnsetRaster{i}(:,2) == SortedIndices(j));
            OffsetNextSyllOnsetRaster{i}(Indices,2) = j - 1000;            
        end
        OffsetRaster{i}(:,2) = OffsetRaster{i}(:,2) + 1000;
        SyllOnsetRaster{i}(:,2) = SyllOnsetRaster{i}(:,2) + 1000;
        OffsetNextSyllOnsetRaster{i}(:,2) = OffsetNextSyllOnsetRaster{i}(:,2) + 1000;
        PlotRaster(OffsetRaster{i}, Cols(mod(i,length(Cols)) + 1), 1, RasterIncrement);
        PlotRaster(SyllOnsetRaster{i}, 'k', 2, RasterIncrement);
        PlotRaster(OffsetNextSyllOnsetRaster{i}, 'k', 2, RasterIncrement);
        RasterIncrement = RasterIncrement + max(OffsetRaster{i}(:,2)) + 2;
        axis tight;
        axes(PSTPlotAxes2);
        OffsetPST{i} = OffsetPST{i}/BinSize;
        plot(Edges, mean(OffsetPST{i}), Cols(mod(i,length(Cols)) + 1));
        hold on;
    end
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
temp(1) = -0.2;
if (temp(2) < 0)
    temp(2) = 0.1;
end
axis(temp);
axes(PSTPlotAxes2); 
set(gca, 'FontSize', 12, 'FontWeight', 'bold');
set(gca, 'Box', 'off');
xlabel('Time (sec)', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('Firing Rate (Hz)', 'FontSize', 14, 'FontWeight', 'bold');
axis auto;
temp2 = axis;
temp2(1:2) = temp(1:2);
axis(temp2);
axes(RasterPlotAxes2);
title([TitleString, '- Aligned to syllable offset'], 'FontSize', 14, 'FontWeight', 'bold');


disp('Finished analysing sequence');