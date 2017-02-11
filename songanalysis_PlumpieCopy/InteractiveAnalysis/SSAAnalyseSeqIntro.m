function [] = SSAAnalyseSeqIntro(Syll, BinSize, FileInfo, TitleString, Motif)

PreOnsetTime = 0.1; % seconds
PostOnsetTime = 0.1; % seconds
GaussianLen = 2; 
Width = 0.005; % seconds

InterBoutInterval = 0.5; % seconds
Cols = ['rgbc'];
ColNames{1} = 'red';
ColNames{2} = 'green';
ColNames{3} = 'blue';
ColNames{4} = 'cyan';
Colors = [1 0 0; 0 1 0; 0 0 1; 0 1 1];

BoutNo = 1;
for NoteFileNo = 1:length(FileInfo.NoteLabels),
    Notes.onsets = FileInfo.NoteOnsets{NoteFileNo};
    Notes.offsets = FileInfo.NoteOffsets{NoteFileNo};
    Notes.labels = FileInfo.NoteLabels{NoteFileNo};
    
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
for i = 1:MaxIntroNotes+1,
    BoutNo(i) = 1;
    OnsetRaster{i} = [];
    OffsetRaster{i} = [];
    SyllOnsetRaster{i} = [];
    SyllOffsetRaster{i} = [];
    OnsetNextSyllOnsetRaster{i} = [];
    OffsetNextSyllOnsetRaster{i} = [];
    OnsetPrevSyllOffsetRaster{i} = [];
    OffsetPrevSyllOffsetRaster{i} = [];
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
                if (size(IntroNotes,1) < size(IntroNotes,2))
                    IntroNotes = IntroNotes';
                end
                IntroNotes = [(IntroNotes(1)+1); IntroNotes];
                for j = 1:length(IntroNotes),
                    SyllOnsetRaster{j} = [SyllOnsetRaster{j}; [(Notes.onsets(IntroNotes(j)) - Notes.offsets(IntroNotes(j))) BoutNo(j)]];
                    SyllOffsetRaster{j} = [SyllOffsetRaster{j}; [(Notes.offsets(IntroNotes(j)) - Notes.onsets(IntroNotes(j))) BoutNo(j)]];
                    OnsetNextSyllOnsetRaster{j} = [OnsetNextSyllOnsetRaster{j}; [(Notes.onsets(IntroNotes(j)+1) - Notes.onsets(IntroNotes(j))) BoutNo(j)]];
                    OffsetNextSyllOnsetRaster{j} = [OffsetNextSyllOnsetRaster{j}; [(Notes.onsets(IntroNotes(j)+1) - Notes.offsets(IntroNotes(j))) BoutNo(j)]];
                    SpikeTimeIndices = find((FileInfo.SpikeData.Times{FileNo} >= (Notes.onsets(IntroNotes(j)) - PreOnsetTime)) & (FileInfo.SpikeData.Times{FileNo} <= (Notes.onsets(IntroNotes(j)+1) + PostOnsetTime)));
                    SpikeTimes = FileInfo.SpikeData.Times{FileNo}(SpikeTimeIndices);
                    if (size(SpikeTimes, 1) < size(SpikeTimes,2))
                        SpikeTimes = SpikeTimes';
                    end
                    OnsetRaster{j} = [OnsetRaster{j}; [(SpikeTimes - Notes.onsets(IntroNotes(j))) ones(size(SpikeTimes))*BoutNo(j)]];
                    OffsetRaster{j} = [OffsetRaster{j}; [(SpikeTimes - Notes.offsets(IntroNotes(j))) ones(size(SpikeTimes))*BoutNo(j)]];
                    OnsetPST{j}(BoutNo(j),:) = histc((SpikeTimes - Notes.onsets(IntroNotes(j))), Edges);
                    OffsetPST{j}(BoutNo(j),:) = histc((SpikeTimes - Notes.offsets(IntroNotes(j))), Edges);
                    OnsetSpikeTrain{j}{BoutNo(j)} = SpikeTimes - Notes.onsets(IntroNotes(j));
                    OffsetSpikeTrain{j}{BoutNo(j)} = SpikeTimes - Notes.offsets(IntroNotes(j));
                    if (IntroNotes(j) ~= Bouts(i,1))
                        OnsetPrevSyllOffsetRaster{j} = [OnsetPrevSyllOffsetRaster{j}; [(Notes.offsets(IntroNotes(j) - 1) - Notes.onsets(IntroNotes(j))) BoutNo(j)]];
                        OffsetPrevSyllOffsetRaster{j} = [OffsetPrevSyllOffsetRaster{j}; [(Notes.offsets(IntroNotes(j) - 1) - Notes.offsets(IntroNotes(j))) BoutNo(j)]];
                    end
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
        PlotRaster(OnsetNextSyllOnsetRaster{i}, 'm', 2, RasterIncrement);
        if (~isempty(OnsetPrevSyllOffsetRaster{i}))
            PlotRaster(OnsetPrevSyllOffsetRaster{i}, 'k', 2, RasterIncrement);
        end
        if (i == 1)
            text(-PreOnsetTime/2, (max(OnsetRaster{i}(:,2))/2 + RasterIncrement), Motif(1), 'FontSize', 12, 'FontWeight', 'bold', 'Color', Colors(mod(i, length(Cols)) + 1,:));
        else
            text(-PreOnsetTime/2, (max(OnsetRaster{i}(:,2))/2 + RasterIncrement), 'i', 'FontSize', 12, 'FontWeight', 'bold', 'Color', Colors(mod(i, length(Cols)) + 1,:));
        end
        RasterIncrement = RasterIncrement + max(OnsetRaster{i}(:,2)) + 2;
        axis tight;
        
        MaxRaster = -100;
        for j = 1:length(OnsetNextSyllOnsetRaster),
            if (~isempty(OnsetNextSyllOnsetRaster{j}))
                MaxRaster = max(MaxRaster, max(OnsetNextSyllOnsetRaster{j}(:,1)));
            end
        end
        Time = -PreOnsetTime:1/Fs:(MaxRaster + PostOnsetTime);
        
        for j = 1:length(OnsetSpikeTrain{i}),
            Indices = round(([OnsetSpikeTrain{i}{j}] + PreOnsetTime) * Fs);
            Indices(find(Indices == 0)) = 1;
            FR(j,Indices) = 1;   
        end

        XGauss = 1:1:(1 + round(2 * GaussianLen * Width * (Fs)));
        XGauss = XGauss - (length(XGauss) + 1)/2;
        GaussWin = (1/((Width * Fs) * sqrt(2 * pi))) * exp(-(XGauss.*XGauss)/(2 * (Width * Fs) * (Width * Fs)));

        % The 4 in the previous line refers to the fact that the gaussian window is
        % terminated at 4 times the standard deviation

        for j = 1:size(FR,1),
            ST(j,:) = conv(FR(j,:),GaussWin, 'same');
        end
        
        axes(PSTPlotAxes);
        OnsetPST{i} = OnsetPST{i}/BinSize;
        %plot(Edges, mean(OnsetPST{i}), Cols(mod(i,length(Cols)) + 1));
        plot(Time, mean(ST), Cols(mod(i,length(Cols)) + 1));
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
temp(1) = -PreOnsetTime;
MaxRaster = -100;
for i = 1:length(OnsetNextSyllOnsetRaster),
    if (~isempty(OnsetNextSyllOnsetRaster{i}))
        MaxRaster = max(MaxRaster, max(OnsetNextSyllOnsetRaster{i}(:,1)));
    end
end
temp(2) = MaxRaster + PostOnsetTime;
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
        if (~isempty(OffsetPrevSyllOffsetRaster{i}))
            PlotRaster(OffsetPrevSyllOffsetRaster{i}, 'm', 2, RasterIncrement);
        end
        MinRaster = 100;
        for j = 1:length(SyllOnsetRaster),
            if (~isempty(SyllOnsetRaster{j}))
                MinRaster = min(MinRaster, min(SyllOnsetRaster{j}(:,1)));
            end
        end
        if (i == 1)
            text(-(PreOnsetTime + MinRaster)/2, (max(OffsetRaster{i}(:,2))/2 + RasterIncrement), Motif(1), 'FontSize', 12, 'FontWeight', 'bold', 'Color', Colors(mod(i, length(Cols)) + 1,:));
        else
            text(-(PreOnsetTime + MinRaster)/2, (max(OffsetRaster{i}(:,2))/2 + RasterIncrement), 'i', 'FontSize', 12, 'FontWeight', 'bold', 'Color', Colors(mod(i, length(Cols)) + 1,:));
        end
        RasterIncrement = RasterIncrement + max(OffsetRaster{i}(:,2)) + 2;
        axis tight;
   
        MaxRaster = -100;
        for j = 1:length(OffsetNextSyllOnsetRaster),
            if (~isempty(OffsetNextSyllOnsetRaster{j}))
                MaxRaster = max(MaxRaster, max(OffsetNextSyllOnsetRaster{j}(:,1)));
            end
        end
        Fs = 10000;
        Time = -PreOnsetTime+MinRaster:1/Fs:(MaxRaster + PostOnsetTime);
        FR = zeros(length(OffsetSpikeTrain{i}),length(Time));
        ST = [];
        
        for j = 1:length(OffsetSpikeTrain{i}),
            Indices = round(([OffsetSpikeTrain{i}{j}] + PreOnsetTime + abs(MinRaster)) * Fs);
            Indices(find(Indices == 0)) = 1;
            FR(j,Indices) = 1;   
        end

        XGauss = 1:1:(1 + round(2 * GaussianLen * Width * (Fs)));
        XGauss = XGauss - (length(XGauss) + 1)/2;
        GaussWin = (1/((Width * Fs) * sqrt(2 * pi))) * exp(-(XGauss.*XGauss)/(2 * (Width * Fs) * (Width * Fs)));

        % The 4 in the previous line refers to the fact that the gaussian window is
        % terminated at 4 times the standard deviation

        for j = 1:size(FR,1),
            ST(j,:) = conv(FR(j,:),GaussWin, 'same');
        end
        
        axes(PSTPlotAxes2);
        OffsetPST{i} = OffsetPST{i}/BinSize;
%        plot(Edges, mean(OffsetPST{i}), Cols(mod(i,length(Cols)) + 1));
        plot(Time, mean(ST), Cols(mod(i,length(Cols)) + 1));
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
MinRaster = 100;

for i = 1:length(SyllOnsetRaster),
    if (~isempty(SyllOnsetRaster{i}))
        MinRaster = min(MinRaster, min(SyllOnsetRaster{i}(:,1)));
    end
end

temp(1) = -PreOnsetTime + MinRaster;

MaxRaster = -100;
for i = 1:length(OnsetNextSyllOnsetRaster),
    if (~isempty(OnsetNextSyllOnsetRaster{i}))
        MaxRaster = max(MaxRaster, max(OnsetNextSyllOnsetRaster{i}(:,1)));
    end
end
temp(2) = MaxRaster + PostOnsetTime;

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