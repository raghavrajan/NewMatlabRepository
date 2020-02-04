function [] = SSAAnalyseIntroNotes(RawDataDir, RecFileDir, PreBoutDuration, PostBoutDuration, FileInfo, TitleString, Motif, Motif2, FileType)

PostMotifOnsetDuration = 3;
InterBoutInterval = 0.5; % in seconds
GaussianLen = 2; 
Width = 0.0025; % seconds

Cols = ['rgbc'];
ColNames{1} = 'red';
ColNames{2} = 'green';
ColNames{3} = 'blue';
ColNames{4} = 'cyan';
Colors = [1 0 0; 0 1 0; 0 0 1; 0 1 1];

Syms{1} = '-';
Syms{2} = ':';
Syms{3} = '-.';
Syms{4} = '--';

AllFiles = dir('*.rec');
Edges = -PreBoutDuration:0.001:3;

for i = 1:length(AllFiles),
    FileTime{i} = AllFiles(i).name;
end

MotifNo = 1;
MotifRaster = [];
MotifPST = [];
MotifSpikeTrain = [];
OnsetsRaster = [];
OffsetsRaster = [];

for NoteFileNo = 1:length(FileInfo.NoteLabels),
    Notes.onsets = FileInfo.NoteOnsets{NoteFileNo};
    Notes.offsets = FileInfo.NoteOffsets{NoteFileNo};
    Notes.labels = FileInfo.NoteLabels{NoteFileNo};
    
    SongFile = FileInfo.FileNames{NoteFileNo};
    
    if (strfind(FileType, 'okrank'))
        [Song, Fs] = SSAReadOKrankData(RawDataDir, RecFileDir, SongFile, 1);
    else
        if (strfind(FileType, 'wav'))
            PresentDir = pwd;
            cd(RawDataDir);
            [Song, Fs] = wavread(SongFile);
            cd(PresentDir);
        else
            if (strfind(FileType, 'obs'));
                [Song, Fs] = SSASoundIn(RawDataDir, RecFileDir, SongFile, 'obs0r');
                Song = Song * 5/32768;
            end
        end
    end
    
    Time = (0:1:(length(Song)-1))/Fs;
%    [m_spec_deriv , m_AM, m_FM ,m_Entropy , m_amplitude , m_Freq, m_PitchGoodness , m_Pitch , Pitch_chose , Pitch_weight] = deriv(Song, Fs);
%    T = linspace(Time(1), Time(end), length(m_Entropy));

    if (strfind(FileType, 'okrank'))
        Index = find(cellfun(@length, strfind(FileTime, SongFile)));
    else
        if (strfind(FileType, 'obs'))
            DotIndex = union(strfind(SongFile, '.cbin') ,strfind(SongFile, '.bbin'));
            Index = find(cellfun(@length, strfind(FileTime, SongFile(1:DotIndex-1))));
        end
    end
    
    if (strfind(FileType, 'okrank'))
        FTime = FileTime{Index+1};
        ExtIndex = strfind(FTime, '.rec');
        FTime = str2double(FTime((ExtIndex-6):(ExtIndex-5)))*3600 + str2double(FTime((ExtIndex-4):(ExtIndex-3)))*60 + str2double(FTime((ExtIndex-2):(ExtIndex-1)));
        FTime = FTime - length(Song)/Fs;
    else
        if (strfind(FileType, 'obs'))
            FTime = FileTime{Index};
            PresentDir = pwd;
            cd(RecFileDir);
            RecFileFid = fopen(FTime, 'r');
            Tline = fgetl(RecFileFid);
            while (ischar(Tline(1)))
                if (strfind(Tline, 'created'))
                    ColonIndex = find(Tline == ':');
                    FTime = str2double(Tline((ColonIndex(2)-2):(ColonIndex(2)-1)))*3600 + str2double(Tline((ColonIndex(3)-2):(ColonIndex(3)-1)))*60 + str2double(Tline((ColonIndex(3)+1):(ColonIndex(3)+2)));
                    break;
                end
                Tline = fgetl(RecFileFid);
            end
            fclose(RecFileFid);
            cd(PresentDir);
        end
    end
    
    Motifs = union(strfind(Notes.labels, Motif), strfind(Notes.labels, Motif2));
    if (~isempty(Motifs))
        if (Motifs(1) == 1)
            MotifOnsetTime = Notes.onsets(Motifs(1));
            if (Notes.onsets(Motifs(1)) > 0.1)
                NoofINs(MotifNo) = 0;
                PrevSyllable(MotifNo) = 'Q';
                if (NoteFileNo == 1)
                    BoutInterval(MotifNo) = Notes.onsets(1) + 10000;
                else
                    BoutInterval(MotifNo) = FTime + Notes.onsets(1) - PrevOffsetTime;
                end
                MotifNo = MotifNo + 1;
            end
        else
            IntroNoteIndices = find(Notes.labels(1:Motifs(1)) == 'z');
            if (~isempty(IntroNoteIndices))
                if (length(IntroNoteIndices) > 1)
                    Intervals = Notes.onsets(IntroNoteIndices(2:end)) - Notes.offsets(IntroNoteIndices(1:end-1));
                    LongIntervals = find(Intervals > InterBoutInterval);
                    if (~isempty(LongIntervals))
                        IntroNoteIndices = IntroNoteIndices((LongIntervals(end) + 1):end);
                    end
                end
                NoofINs(MotifNo) = length(IntroNoteIndices);
                MotifOnsetTime = Notes.onsets(Motifs(1));
                if (IntroNoteIndices(1) == 1)
                    if (NoteFileNo == 1)
                        BoutInterval(MotifNo) = Notes.onsets(1);
                    end
                    PrevSyllable(MotifNo) = 'Q';
                else
                    BoutInterval(MotifNo) = Notes.onsets(IntroNoteIndices(1)) - Notes.offsets(IntroNoteIndices(1) - 1);
                    PrevSyllable(MotifNo) = Notes.labels(IntroNoteIndices(1) - 1);
                end

                IntroNoteIndices = sort(IntroNoteIndices, 'descend');
                for j = 1:length(IntroNoteIndices),
                    SpikeTimeIndices = find((FileInfo.SpikeData.Times{NoteFileNo} >= (Notes.onsets(IntroNoteIndices(j)) - 0.1)) & (FileInfo.SpikeData.Times{NoteFileNo} <= (Notes.offsets(IntroNoteIndices(j)) + 0.1)));
                    SpikeTimes = FileInfo.SpikeData.Times{NoteFileNo}(SpikeTimeIndices) - Notes.onsets(IntroNoteIndices(j));
                    if (size(SpikeTimes, 1) < size(SpikeTimes,2))
                        SpikeTimes = SpikeTimes';
                    end
                    MotifSpikeTrain{MotifNo}{j} = SpikeTimes;
                    if (IntroNoteIndices(j) ~= 1)
                        PrevSyllOffsetTime{MotifNo}{j} = Notes.offsets(IntroNoteIndices(j) - 1) - Notes.onsets(IntroNoteIndices(j));
                    else
                        PrevSyllOffsetTime{MotifNo}{j} = -1;
                    end
                    SyllOffsetTime{MotifNo}{j} = Notes.offsets(IntroNoteIndices(j)) - Notes.onsets(IntroNoteIndices(j));
                    NextSyllOnsetTime{MotifNo}{j} = Notes.onsets(IntroNoteIndices(j) + 1) - Notes.onsets(IntroNoteIndices(j));
                end
                MotifNo = MotifNo + 1;
            else
                NoofINs(MotifNo) = 0;
                MotifOnsetTime = Notes.onsets(Motifs(1));
                BoutInterval(MotifNo) = Notes.onsets(Motifs(1)) - Notes.offsets(Motifs(1) - 1);
                PrevSyllable(MotifNo) = Notes.labels(Motifs(1) - 1);
                MotifNo = MotifNo + 1;
            end
        end
        
        for i = 2:length(Motifs),
            IntroNoteIndices = find(Notes.labels(Motifs(i-1):Motifs(i)) == 'i');
            if (~isempty(IntroNoteIndices))
                IntroNoteIndices = IntroNoteIndices + Motifs(i-1) - 1;
                if (length(IntroNoteIndices) > 1)
                    Intervals = Notes.onsets(IntroNoteIndices(2:end)) - Notes.offsets(IntroNoteIndices(1:end-1));
                    LongIntervals = find(Intervals > InterBoutInterval);
                    if (~isempty(LongIntervals))
                        IntroNoteIndices = IntroNoteIndices((LongIntervals(end) + 1):end);
                    end
                end
                MotifOnsetTime = Notes.onsets(Motifs(i));
                NoofINs(MotifNo) = length(IntroNoteIndices);
                BoutInterval(MotifNo) = Notes.onsets(IntroNoteIndices(1)) - Notes.offsets(IntroNoteIndices(1)-1);
                PrevSyllable(MotifNo) = Notes.labels(IntroNoteIndices(1) - 1);

                IntroNoteIndices = sort(IntroNoteIndices, 'descend');
                for j = 1:length(IntroNoteIndices),
                    SpikeTimeIndices = find((FileInfo.SpikeData.Times{NoteFileNo} >= (Notes.onsets(IntroNoteIndices(j)) - 0.1)) & (FileInfo.SpikeData.Times{NoteFileNo} <= (Notes.offsets(IntroNoteIndices(j)) + 0.1)));
                    SpikeTimes = FileInfo.SpikeData.Times{NoteFileNo}(SpikeTimeIndices) - Notes.onsets(IntroNoteIndices(j));
                    if (size(SpikeTimes, 1) < size(SpikeTimes,2))
                        SpikeTimes = SpikeTimes';
                    end
                    MotifSpikeTrain{MotifNo}{j} = SpikeTimes;
                    if (IntroNoteIndices(j) ~= 1)
                        PrevSyllOffsetTime{MotifNo}{j} = Notes.offsets(IntroNoteIndices(j) - 1) - Notes.onsets(IntroNoteIndices(j));
                    else
                        PrevSyllOffsetTime{MotifNo}{j} = -1;
                    end
                    SyllOffsetTime{MotifNo}{j} = Notes.offsets(IntroNoteIndices(j)) - Notes.onsets(IntroNoteIndices(j));
                    NextSyllOnsetTime{MotifNo}{j} = Notes.onsets(IntroNoteIndices(j) + 1) - Notes.onsets(IntroNoteIndices(j));
                end
                
                MotifNo = MotifNo + 1;
            else
                NoofINs(MotifNo) = 0;
                MotifOnsetTime = Notes.onsets(Motifs(i));
                BoutInterval(MotifNo) = Notes.onsets(Motifs(i)) - Notes.offsets(Motifs(i)-1);
                PrevSyllable(MotifNo) = Notes.labels(Motifs(i) - 1);                
                MotifNo = MotifNo + 1;
            end
        end
    end
    PrevOffsetTime = Notes.offsets(end) + FTime;
end

MaxINotes = max(NoofINs);
%MaxINotes = 5;

RasterFig = figure;
PlotGap = 0.04;
PlotWidth = (0.85 - (MaxINotes - 1)*PlotGap)/MaxINotes;

SeqSummFig = figure;

figure(RasterFig);
for i = 1:MaxINotes,
    RasterPlotAxes(i) = axes('Position', [(0.1 + (i-1)*PlotWidth + (i-1)*PlotGap) 0.4 PlotWidth 0.55]);
    PSTPlotAxes(i) = axes('Position', [(0.1 + (i-1)*PlotWidth + (i-1)*PlotGap) 0.1 PlotWidth 0.25]);
end
RasterIncrement = 0;

Fs = 1000;
XGauss = 1:1:(1 + round(2 * GaussianLen * Width * (Fs)));
XGauss = XGauss - (length(XGauss) + 1)/2;
GaussWin = (1/((Width * Fs) * sqrt(2 * pi))) * exp(-(XGauss.*XGauss)/(2 * (Width * Fs) * (Width * Fs)));

Edges = -0.1:0.001:0.3;

figure(RasterFig);
for i = 1:MaxINotes,
    for j = 1:i,
        Raster{j} = [];
        PrevOFFRaster{j} = [];
        OFFRaster{j} = [];
        NextONRaster{j} = [];
        PST{j} = [];
        PreMotifSpikeTrain{j} = [];
        DurMotifSpikeTrain{j} = [];
        FullMotifSpikeTrain{j} = [];
        PreCorr{j} = [];
        DurCorr{j} = [];
        FullCorr{j} = [];
    end
    PreFiringRate = [];
    DurFiringRate = [];
    FullFiringRate = [];
    NoteLens = [];
    MedianNoteLen = [];
    
    Indices = find(NoofINs == i);
    if (~isempty(Indices))
        for j = 1:length(Indices),
            for k = 1:i,
                NoteLens(k,j) = SyllOffsetTime{Indices(j)}{k};
            end
        end
        
        for k = 1:i,
            MedianNoteLen(k) = median(NoteLens(k,:));
            MedianNote{k}.Length = median(NoteLens(k,:));
        end
        
        for j = 1:length(Indices),
            for k = 1:i,
                Raster{k} = [Raster{k}; [MotifSpikeTrain{Indices(j)}{k} ones(size(MotifSpikeTrain{Indices(j)}{k}))*j]];
                PrevOFFRaster{k} = [PrevOFFRaster{k}; [PrevSyllOffsetTime{Indices(j)}{k} j]];
                OFFRaster{k} = [OFFRaster{k}; [SyllOffsetTime{Indices(j)}{k} j]];
                NextONRaster{k} = [NextONRaster{k}; [NextSyllOnsetTime{Indices(j)}{k} j]];
                PST{k}(j,:) = conv(histc(MotifSpikeTrain{Indices(j)}{k}, Edges)/0.001, GaussWin, 'same');
                PreFiringRate(j,k) = length(find((MotifSpikeTrain{Indices(j)}{k} >= -0.05) & (MotifSpikeTrain{Indices(j)}{k} < 0)))/0.05;
                DurFiringRate(j,k) = length(find((MotifSpikeTrain{Indices(j)}{k} >= 0) & (MotifSpikeTrain{Indices(j)}{k} < SyllOffsetTime{Indices(j)}{k})))/(SyllOffsetTime{Indices(j)}{k});
                FullFiringRate(j,k) = length(find((MotifSpikeTrain{Indices(j)}{k} >= -0.05) & (MotifSpikeTrain{Indices(j)}{k} < (SyllOffsetTime{Indices(j)}{k}))))/(SyllOffsetTime{Indices(j)}{k} + 0.05);
                
                SpikeIndices = find((MotifSpikeTrain{Indices(j)}{k} >= -0.05) & (MotifSpikeTrain{Indices(j)}{k} < 0));
                PreMotifSpikeTrain{k}{j} = MotifSpikeTrain{Indices(j)}{k}(SpikeIndices);
        
                SpikeIndices = find((MotifSpikeTrain{Indices(j)}{k} >= 0) & (MotifSpikeTrain{Indices(j)}{k} < SyllOffsetTime{Indices(j)}{k}));
                DurMotifSpikeTrain{k}{j} = (MotifSpikeTrain{Indices(j)}{k}(SpikeIndices) * MedianNoteLen(k))/SyllOffsetTime{Indices(j)}{k};

                FullMotifSpikeTrain{k}{j} = [PreMotifSpikeTrain{k}{j}; DurMotifSpikeTrain{k}{j}];
            end
        end
        for j = 1:i,
            axes(RasterPlotAxes(j));
            if (~isempty(Raster{j}))
                PlotRaster(Raster{j}, Cols(mod(i,length(Cols)) + 1), 2, RasterIncrement);
            end
            %PlotRaster(PrevOFFRaster{j}, 'k', 4, RasterIncrement);
            PlotRaster(OFFRaster{j}, 'k', 3, RasterIncrement);
            PlotRaster(NextONRaster{j}, 'k', 3, RasterIncrement);
            if (j == 1)
                text(-0.2, (max(NextONRaster{j}(:,2))/2 + RasterIncrement), num2str(i), 'FontSize', 12, 'FontWeight', 'bold', 'Color', Colors(mod(i, length(Cols)) + 1,:));
            end
            axes(PSTPlotAxes(j));
            %if (size(PST{j},1) >= 5)
                plot(Edges, mean(PST{j}), [Cols(mod(i,length(Cols)) + 1), Syms{mod(floor(i/length(Cols)), length(Syms)) + 1}], 'LineWidth', 2);
                hold on;
            %end
            [PreCorr{j}] = SSACalculateCorrGaussSameSize(PreMotifSpikeTrain{j}, MedianNote{j}, 0.001, 0.05, MedianNoteLen(j), 'Undirected - pre correlation', 4);
            [DurCorr{j}] = SSACalculateCorrGaussSameSize(DurMotifSpikeTrain{j}, MedianNote{j}, 0.001, 0, 0, 'Undirected - dur correlation', 4);
            [FullCorr{j}] = SSACalculateCorrGaussSameSize(FullMotifSpikeTrain{j}, MedianNote{j}, 0.001, 0.05, 0, 'Undirected - full correlation', 4);
        end
        figure(SeqSummFig);
        subplot(3,2,1);
        hold on;
        errorbar([1:1:i], mean(PreFiringRate, 1), std(PreFiringRate/sqrt(size(PreFiringRate,1)), [], 1), [Cols(mod(i,length(Cols)) + 1), Syms{mod(floor(i/length(Cols)), length(Syms)) + 1}, 's']);
        
        subplot(3,2,3);
        hold on;
        errorbar([1:1:i], mean(DurFiringRate, 1), std(DurFiringRate, [], 1)/sqrt(size(DurFiringRate,1)), [Cols(mod(i,length(Cols)) + 1), Syms{mod(floor(i/length(Cols)), length(Syms)) + 1}, 's']);
        
        subplot(3,2,5);
        hold on;
        errorbar([1:1:i], mean(FullFiringRate, 1), std(FullFiringRate, [], 1)/sqrt(size(FullFiringRate,1)), [Cols(mod(i,length(Cols)) + 1), Syms{mod(floor(i/length(Cols)), length(Syms)) + 1}, 's']);
        
        subplot(3,2,2);
        hold on;
        Corr = [];
        StdCorr = [];
        for j = 1:i,
            Corr(j) = PreCorr{j}(2);
            StdCorr(j) = PreCorr{j}(3);
        end
        errorbar([1:1:i], Corr, StdCorr, [Cols(mod(i,length(Cols)) + 1), Syms{mod(floor(i/length(Cols)), length(Syms)) + 1}, 's']);
        
        subplot(3,2,4);
        hold on;
        Corr = [];
        StdCorr = [];
        for j = 1:i,
            Corr(j) = DurCorr{j}(2);
            StdCorr(j) = DurCorr{j}(3);
        end
        errorbar([1:1:i], Corr, StdCorr, [Cols(mod(i,length(Cols)) + 1), Syms{mod(floor(i/length(Cols)), length(Syms)) + 1}, 's']);
        
        subplot(3,2,6);
        hold on;
        Corr = [];
        StdCorr = [];
        for j = 1:i,
            Corr(j) = FullCorr{j}(2);
            StdCorr(j) = FullCorr{j}(3);
        end
        errorbar([1:1:i], Corr, StdCorr, [Cols(mod(i,length(Cols)) + 1), Syms{mod(floor(i/length(Cols)), length(Syms)) + 1}, 's']);
                
        RasterIncrement = RasterIncrement + max(NextONRaster{j}(:,2)) + 2;
    end
end

INSummFig = figure;

for i = 1:MaxINotes,
    PreFiringRate = [];
    DurFiringRate = [];
    FullFiringRate = [];
    PreMotifSpikeTrain = [];
    DurMotifSpikeTrain = [];
    FullMotifSpikeTrain = [];
    
    Indices = find(NoofINs >= i);
    NoteLens = [];
    for j = 1:length(Indices),
        NoteLens(j) = SyllOffsetTime{Indices(j)}{i};
    end
    MedianNoteLen = median(NoteLens);
    
    for j = 1:length(Indices),
        PreFiringRate(j,1) = length(find((MotifSpikeTrain{Indices(j)}{i} >= -0.05) & (MotifSpikeTrain{Indices(j)}{i} < 0)))/0.05;
        DurFiringRate(j,1) = length(find((MotifSpikeTrain{Indices(j)}{i} >= 0) & (MotifSpikeTrain{Indices(j)}{i} < SyllOffsetTime{Indices(j)}{i})))/(SyllOffsetTime{Indices(j)}{i});
        FullFiringRate(j,1) = length(find((MotifSpikeTrain{Indices(j)}{i} >= -0.05) & (MotifSpikeTrain{Indices(j)}{i} < (SyllOffsetTime{Indices(j)}{i}))))/(SyllOffsetTime{Indices(j)}{i} + 0.05);      
        
        SpikeIndices = find((MotifSpikeTrain{Indices(j)}{i} >= -0.05) & (MotifSpikeTrain{Indices(j)}{i} < 0));
        PreMotifSpikeTrain{j} = MotifSpikeTrain{Indices(j)}{i}(SpikeIndices);
        
        SpikeIndices = find((MotifSpikeTrain{Indices(j)}{i} >= 0) & (MotifSpikeTrain{Indices(j)}{i} < SyllOffsetTime{Indices(j)}{i}));
        DurMotifSpikeTrain{j} = (MotifSpikeTrain{Indices(j)}{i}(SpikeIndices) * MedianNoteLen)/SyllOffsetTime{Indices(j)}{i};
        
        FullMotifSpikeTrain{j} = [PreMotifSpikeTrain{j}; DurMotifSpikeTrain{j}];
    end
    MedianNote.Length = median(NoteLens);
    [PreCorr] = SSACalculateCorrGaussSameSize(PreMotifSpikeTrain, MedianNote, 0.001, 0.05, MedianNoteLen, 'Undirected - pre correlation', 4);
    [DurCorr] = SSACalculateCorrGaussSameSize(DurMotifSpikeTrain, MedianNote, 0.001, 0, 0, 'Undirected - dur correlation', 4);
    [FullCorr] = SSACalculateCorrGaussSameSize(FullMotifSpikeTrain, MedianNote, 0.001, 0.05, 0, 'Undirected - full correlation', 4);
    
    figure(INSummFig);
    subplot(3,2,1);
    hold on;
    errorbar(i, mean(PreFiringRate), std(PreFiringRate)/sqrt(length(Indices)), 'bs');
    
    subplot(3,2,3);
    hold on;
    errorbar(i, mean(DurFiringRate), std(DurFiringRate)/sqrt(length(Indices)), 'bs');
    
    subplot(3,2,5);
    hold on;
    errorbar(i, mean(FullFiringRate), std(FullFiringRate)/sqrt(length(Indices)), 'bs');
    
    subplot(3,2,2);
    hold on;
    errorbar(i, PreCorr(2), PreCorr(3), 'bs');
    
    subplot(3,2,4);
    hold on;
    errorbar(i, DurCorr(2), DurCorr(3), 'bs');
    
    subplot(3,2,6);
    hold on;
    errorbar(i, FullCorr(2), FullCorr(3), 'bs');
    
end

figure(RasterFig);
set(gcf, 'Position', [360 100 650 600]);
set(gcf, 'Color', 'w');

MaxPSTAxis = -100;

for i = 1:MaxINotes,
    axes(RasterPlotAxes(i)); 
    set(gca, 'FontSize', 12, 'FontWeight', 'bold');
    %if (i == 1)
    %    ylabel('Trials', 'FontSize', 14, 'FontWeight', 'bold');
    %end
    set(gca, 'YColor', 'w');
    set(gca, 'YTick', []);
    set(gca, 'XColor', 'w');
    axis tight;
    Temp = axis;
    Temp(1:3) = [Edges(1) 0.15 0];
    plot([0 0], [0 Temp(4)], 'k--', 'LineWidth', 1.5);
    axis(Temp);
    %set(gca, 'Box', 'on');
    title(['i-', num2str(i)], 'FontSize', 14, 'FontWeight', 'bold');
    
    axes(PSTPlotAxes(i));
    if (i == 3)
        xlabel('Time (sec)', 'FontSize', 14, 'FontWeight', 'bold');
    else
        set(gca, 'XTick', []);
    end
    if (i == 1)
        ylabel('Firing Rate (Hz)', 'FontSize', 14, 'FontWeight', 'bold');
    else
        set(gca, 'YColor', 'w');
    end
    
    axis tight;
    Temp = axis;
    Temp(1:3) = [Edges(1) 0.15 0];
    axis(Temp);
    MaxPSTAxis = max(MaxPSTAxis, Temp(4));
    plot([0 0], [0 Temp(4)], 'k--', 'LineWidth', 1.5);
    set(gca, 'FontSize', 12, 'FontWeight', 'bold');
    set(gca, 'Box', 'off');
end

for i = 1:MaxINotes,
    axes(PSTPlotAxes(i));
    plot([0 0], [0 MaxPSTAxis], 'k--', 'LineWidth', 1.5);
    axis([Edges(1) 0.15 0 MaxPSTAxis]);
end
%save([FileInfo.FileNames{1}, '.motif_transitions.mat'], 'MotifRaster', 'OnsetsRaster', 'MotifPST', 'NoofINs', 'BoutInterval');
disp('Finished analysing notes');