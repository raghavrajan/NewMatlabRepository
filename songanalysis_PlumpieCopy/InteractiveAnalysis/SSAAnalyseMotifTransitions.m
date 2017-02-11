function [] = SSAAnalyseMotifTransitions(RawDataDir, RecFileDir, PreBoutDuration, PostBoutDuration, FileInfo, TitleString, Motif, Motif2, FileType)

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
            DotIndex = strfind(SongFile, '.cbin');
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
                SpikeTimeIndices = find((FileInfo.SpikeData.Times{NoteFileNo} >= (MotifOnsetTime - PreBoutDuration)) & (FileInfo.SpikeData.Times{NoteFileNo} <= (MotifOnsetTime + PostMotifOnsetDuration)));
                SpikeTimes = FileInfo.SpikeData.Times{NoteFileNo}(SpikeTimeIndices) - MotifOnsetTime;
                if (size(SpikeTimes, 1) < size(SpikeTimes,2))
                    SpikeTimes = SpikeTimes';
                end
                MotifRaster = [MotifRaster; [SpikeTimes ones(size(SpikeTimes))*MotifNo]];
                MotifSpikeTrain{MotifNo} = SpikeTimes;
                MotifPST(MotifNo,:) = histc(SpikeTimes, Edges);

                OnsetsRaster = [OnsetsRaster; [(Time(1) - MotifOnsetTime) MotifNo]];

                MotifNo = MotifNo + 1;
            end
        else
            IntroNoteIndices = find(Notes.labels(1:Motifs(1)) == 'i');
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
                        BoutInterval(MotifNo) = Notes.onsets(1) + 5;
                    else
                        BoutInterval(MotifNo) = FTime + Notes.onsets(1) - PrevOffsetTime;
                    end
                    PrevSyllable(MotifNo) = 'Q';
                else
                    BoutInterval(MotifNo) = Notes.onsets(IntroNoteIndices(1)) - Notes.offsets(IntroNoteIndices(1) - 1);
                    PrevSyllable(MotifNo) = Notes.labels(IntroNoteIndices(1) - 1);
                end

                if (IntroNoteIndices(1) == 1)
                    SpikeTimeIndices = find((FileInfo.SpikeData.Times{NoteFileNo} >= (MotifOnsetTime - PreBoutDuration)) & (FileInfo.SpikeData.Times{NoteFileNo} <= (MotifOnsetTime + PostMotifOnsetDuration)));
                else
                    if ((MotifOnsetTime - PreBoutDuration) <= Notes.offsets(IntroNoteIndices(1) - 1))
                        SpikeTimeIndices = find((FileInfo.SpikeData.Times{NoteFileNo} >= (Notes.offsets(IntroNoteIndices(1) - 1))) & (FileInfo.SpikeData.Times{NoteFileNo} <= (MotifOnsetTime + PostMotifOnsetDuration)));
                    else
                        SpikeTimeIndices = find((FileInfo.SpikeData.Times{NoteFileNo} >= (MotifOnsetTime - PreBoutDuration)) & (FileInfo.SpikeData.Times{NoteFileNo} <= (MotifOnsetTime + PostMotifOnsetDuration)));
                    end
                end
                SpikeTimes = FileInfo.SpikeData.Times{NoteFileNo}(SpikeTimeIndices) - MotifOnsetTime;
                if (size(SpikeTimes, 1) < size(SpikeTimes,2))
                    SpikeTimes = SpikeTimes';
                end
                MotifRaster = [MotifRaster; [SpikeTimes ones(size(SpikeTimes))*MotifNo]];
                MotifSpikeTrain{MotifNo} = SpikeTimes;
                MotifPST(MotifNo,:) = histc(SpikeTimes, Edges);

                if (IntroNoteIndices(1) == 1)
                    OnsetsRaster = [OnsetsRaster; [(Time(1) - MotifOnsetTime) MotifNo]];
                else
                    if ((MotifOnsetTime - PreBoutDuration) <= Notes.offsets(IntroNoteIndices(1) - 1))
                        OnsetsRaster = [OnsetsRaster; [(Notes.offsets(IntroNoteIndices(1) - 1) - MotifOnsetTime) MotifNo]];
                    else
                        OnsetsRaster = [OnsetsRaster; [(-PreBoutDuration) MotifNo]];
                    end
                end
                MotifNo = MotifNo + 1;
            else
                NoofINs(MotifNo) = 0;
                MotifOnsetTime = Notes.onsets(Motifs(1));
                BoutInterval(MotifNo) = Notes.onsets(Motifs(1)) - Notes.offsets(Motifs(1) - 1);
                PrevSyllable(MotifNo) = Notes.labels(Motifs(1) - 1);
                if ((MotifOnsetTime - PreBoutDuration) <= Notes.offsets(Motifs(1) - 1))
                    SpikeTimeIndices = find((FileInfo.SpikeData.Times{NoteFileNo} >= (Notes.offsets(Motifs(1) - 1))) & (FileInfo.SpikeData.Times{NoteFileNo} <= (MotifOnsetTime + PostMotifOnsetDuration)));
                else
                    SpikeTimeIndices = find((FileInfo.SpikeData.Times{NoteFileNo} >= (MotifOnsetTime - PreBoutDuration)) & (FileInfo.SpikeData.Times{NoteFileNo} <= (MotifOnsetTime + PostMotifOnsetDuration)));
                end
                SpikeTimes = FileInfo.SpikeData.Times{NoteFileNo}(SpikeTimeIndices) - MotifOnsetTime;
                if (size(SpikeTimes, 1) < size(SpikeTimes,2))
                    SpikeTimes = SpikeTimes';
                end
                MotifRaster = [MotifRaster; [SpikeTimes ones(size(SpikeTimes))*MotifNo]];
                MotifSpikeTrain{MotifNo} = SpikeTimes;
                MotifPST(MotifNo,:) = histc(SpikeTimes, Edges);

                if ((MotifOnsetTime - PreBoutDuration) <= Notes.offsets(Motifs(1) - 1))
                    OnsetsRaster = [OnsetsRaster; [(Notes.offsets(Motifs(1)) - MotifOnsetTime) MotifNo]];
                else
                    OnsetsRaster = [OnsetsRaster; [(-PreBoutDuration) MotifNo]];
                end
 
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
                if ((MotifOnsetTime - PreBoutDuration) <= Notes.offsets(IntroNoteIndices(1) - 1))
                    SpikeTimeIndices = find((FileInfo.SpikeData.Times{NoteFileNo} >= (Notes.offsets(IntroNoteIndices(1) - 1))) & (FileInfo.SpikeData.Times{NoteFileNo} <= (MotifOnsetTime + PostMotifOnsetDuration)));
                else
                    SpikeTimeIndices = find((FileInfo.SpikeData.Times{NoteFileNo} >= (MotifOnsetTime - PreBoutDuration)) & (FileInfo.SpikeData.Times{NoteFileNo} <= (MotifOnsetTime + PostMotifOnsetDuration)));
                end
                SpikeTimes = FileInfo.SpikeData.Times{NoteFileNo}(SpikeTimeIndices) - MotifOnsetTime;
                if (size(SpikeTimes, 1) < size(SpikeTimes,2))
                    SpikeTimes = SpikeTimes';
                end
                MotifRaster = [MotifRaster; [SpikeTimes ones(size(SpikeTimes))*MotifNo]];
                MotifSpikeTrain{MotifNo} = SpikeTimes;
                MotifPST(MotifNo,:) = histc(SpikeTimes, Edges);

                if ((MotifOnsetTime - PreBoutDuration) <= Notes.offsets(IntroNoteIndices(1) - 1))
                    OnsetsRaster = [OnsetsRaster; [(Notes.offsets(IntroNoteIndices(1) - 1) - MotifOnsetTime) MotifNo]];
                else
                    OnsetsRaster = [OnsetsRaster; [(-PreBoutDuration) MotifNo]];
                end
                MotifNo = MotifNo + 1;
            else
                NoofINs(MotifNo) = 0;
                MotifOnsetTime = Notes.onsets(Motifs(i));
                BoutInterval(MotifNo) = Notes.onsets(Motifs(i)) - Notes.offsets(Motifs(i)-1);
                PrevSyllable(MotifNo) = Notes.labels(Motifs(i) - 1);                
                if ((MotifOnsetTime - PreBoutDuration) <= Notes.offsets(Motifs(i) - 1))
                    SpikeTimeIndices = find((FileInfo.SpikeData.Times{NoteFileNo} >= (Notes.offsets(Motifs(i) - 1))) & (FileInfo.SpikeData.Times{NoteFileNo} <= (MotifOnsetTime + PostMotifOnsetDuration)));
                else
                    SpikeTimeIndices = find((FileInfo.SpikeData.Times{NoteFileNo} >= (MotifOnsetTime - PreBoutDuration)) & (FileInfo.SpikeData.Times{NoteFileNo} <= (MotifOnsetTime + PostMotifOnsetDuration)));
                end
                SpikeTimes = FileInfo.SpikeData.Times{NoteFileNo}(SpikeTimeIndices) - MotifOnsetTime;
                if (size(SpikeTimes, 1) < size(SpikeTimes,2))
                    SpikeTimes = SpikeTimes';
                end
                MotifRaster = [MotifRaster; [SpikeTimes ones(size(SpikeTimes))*MotifNo]];
                MotifSpikeTrain{MotifNo} = SpikeTimes;
                MotifPST(MotifNo,:) = histc(SpikeTimes, Edges);

                if ((MotifOnsetTime - PreBoutDuration) <= Notes.offsets(Motifs(i) - 1))
                    OnsetsRaster = [OnsetsRaster; [(Notes.offsets(Motifs(i) - 1) - MotifOnsetTime) MotifNo]];
                else
                    OnsetsRaster = [OnsetsRaster; [(-PreBoutDuration) MotifNo]];
                end
                MotifNo = MotifNo + 1;
            end
        end
    end
    PrevOffsetTime = Notes.offsets(end) + FTime;
end

MotifPST = MotifPST/0.001;
RasterFig = figure;
RasterPlotAxes = axes('Position', [0.1 0.4 0.85 0.55]);
PSTPlotAxes = axes('Position', [0.1 0.1 0.85 0.25]);

figure(RasterFig);
hold on;
RasterIncrement = 0;

Fs = 1000;
XGauss = 1:1:(1 + round(2 * GaussianLen * Width * (Fs)));
XGauss = XGauss - (length(XGauss) + 1)/2;
GaussWin = (1/((Width * Fs) * sqrt(2 * pi))) * exp(-(XGauss.*XGauss)/(2 * (Width * Fs) * (Width * Fs)));

for i = 0:max(NoofINs),
    INRaster = [];
    ONRaster = [];
    Indices = find(NoofINs == i);
    if (~isempty(Indices))
        axes(RasterPlotAxes);
        SpikeTrain = [];
        PST = [];
        for j = 1:length(Indices),
            RasterIndices = find(MotifRaster(:,2) == Indices(j));
            INRaster = [INRaster; MotifRaster(RasterIndices,:)];
            RasterIndices = find(OnsetsRaster(:,2) == Indices(j));
            ONRaster = [ONRaster; OnsetsRaster(RasterIndices,:)];
            SpikeTrain{j} = MotifSpikeTrain{Indices(j)};
            PST(j,:) = conv(MotifPST(Indices(j),:), GaussWin, 'same');
        end
        if (~isempty(INRaster))
            Intervals = BoutInterval(Indices);
            [SortedIntervals, SortedIndices] = sort(Intervals);
            for j = 1:length(SortedIndices),
                RasterIndices = find(INRaster(:,2) == Indices(SortedIndices(j)));
                INRaster(RasterIndices,2) = j - 1000;
                RasterIndices = find(ONRaster(:,2) == Indices(SortedIndices(j)));
                ONRaster(RasterIndices,2) = j - 1000;
            end
            INRaster(:,2) = INRaster(:,2) + 1000;
            ONRaster(:,2) = ONRaster(:,2) + 1000;
            PlotRaster(INRaster, Cols(mod(i,length(Cols)) + 1), 2, RasterIncrement);
            PlotRaster(ONRaster, 'k', 4, RasterIncrement);
            text(-PreBoutDuration/2, (max(INRaster(:,2))/2 + RasterIncrement), num2str(i), 'FontSize', 12, 'FontWeight', 'bold', 'Color', Colors(mod(i, length(Cols)) + 1,:));
            RasterIncrement = RasterIncrement + max(INRaster(:,2)) + 2;
        end
        
        axes(PSTPlotAxes);
        if (size(PST,1) >= 5)
            plot(Edges, mean(PST), [Cols(mod(i,length(Cols)) + 1), Syms{mod(floor(i/length(Cols)), length(Syms)) + 1}]);
            hold on;
        end
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
Temp = axis;
Temp(1:2) = [-PreBoutDuration 0];
axis(Temp);
title([TitleString, '- Ordered by # of intro notes - Motif Onset at time 0'], 'FontSize', 14, 'FontWeight', 'bold');
plot([0 0], [0 Temp(4)], 'k--', 'LineWidth', 1.5);

axes(PSTPlotAxes);
set(gca, 'FontSize', 12, 'FontWeight', 'bold');
set(gca, 'Box', 'off');
xlabel('Time (sec)', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('Firing Rate (Hz)', 'FontSize', 14, 'FontWeight', 'bold');
axis tight;
Temp = axis;
Temp(1:2) = [-PreBoutDuration 0];
axis(Temp);
plot([0 0], [0 Temp(4)], 'k--', 'LineWidth', 1.5);

RasterFig2 = figure;
RasterPlotAxes2 = axes('Position', [0.1 0.4 0.85 0.55]);
PSTPlotAxes2 = axes('Position', [0.1 0.1 0.85 0.25]);

figure(RasterFig2);
hold on;
RasterIncrement = 0;

Fs = 1000;
XGauss = 1:1:(1 + round(2 * GaussianLen * Width * (Fs)));
XGauss = XGauss - (length(XGauss) + 1)/2;
GaussWin = (1/((Width * Fs) * sqrt(2 * pi))) * exp(-(XGauss.*XGauss)/(2 * (Width * Fs) * (Width * Fs)));

% The 4 in the previous line refers to the fact that the gaussian window is
% terminated at 4 times the standard deviation

INRaster = MotifRaster;
ONRaster = OnsetsRaster;
PST = [];

if (~isempty(INRaster))
    [SortedIntervals, SortedIndices] = sort(BoutInterval);
    for j = 1:length(SortedIndices),
        RasterIndices = find(INRaster(:,2) == SortedIndices(j));
        INRaster(RasterIndices,2) = j - 1000;
        RasterIndices = find(ONRaster(:,2) == SortedIndices(j));
        ONRaster(RasterIndices,2) = j - 1000;
        PST(j,:) = conv(MotifPST(SortedIndices(j),:), GaussWin, 'same');
    end
    INRaster(:,2) = INRaster(:,2) + 1000;
    ONRaster(:,2) = ONRaster(:,2) + 1000;
    
    Units = [0 0.2 1 10 100 1000];
    for i = 1:length(Units);
        if (i ~= length(Units))
            Indices = find((SortedIntervals >= Units(i) & (SortedIntervals < Units(i+1))));
        else
            Indices = find(SortedIntervals >= Units(end));
        end
        BINRaster = [];
        BONRaster = [];
        if (~isempty(Indices))
            for j = 1:length(Indices),
                NewIndices = find(INRaster(:,2) == Indices(j));
                BINRaster = [BINRaster; INRaster(NewIndices,:)];
                NewIndices = find(ONRaster(:,2) == Indices(j));
                BONRaster = [BONRaster; ONRaster(NewIndices,:)];
            end
            BINRaster(:,2) = BINRaster(:,2) - min(BINRaster(:,2)) + 1;
            BONRaster(:,2) = BONRaster(:,2) - min(BONRaster(:,2)) + 1;
            
            axes(RasterPlotAxes2);
            PlotRaster(BINRaster, Cols(mod(i,length(Cols)) + 1), 2, RasterIncrement);
            PlotRaster(BONRaster, 'k', 4, RasterIncrement);
            
            if (i ~= length(Units))
                text(-PreBoutDuration/2, (max(BINRaster(:,2))/2 + RasterIncrement), [num2str(Units(i)), '-', num2str(Units(i+1))], 'FontSize', 12, 'FontWeight', 'bold', 'Color', Colors(mod(i, length(Cols)) + 1,:));
            else
                text(-PreBoutDuration/2, (max(BINRaster(:,2))/2 + RasterIncrement), ['>=', num2str(Units(i))], 'FontSize', 12, 'FontWeight', 'bold', 'Color', Colors(mod(i, length(Cols)) + 1,:));
            end

            RasterIncrement = RasterIncrement + max(BINRaster(:,2)) + 2;
            
            axes(PSTPlotAxes2);
            if (length(Indices) >= 5)
                plot(Edges, mean(PST(Indices,:)), [Cols(mod(i,length(Cols)) + 1), Syms{mod(floor(1/length(Cols)), length(Syms)) + 1}]);
                hold on;
            end
        end
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
Temp = axis;
Temp(1:2) = [-PreBoutDuration 0];
axis(Temp);
title([TitleString, '- Ordered by bout interval - onset of 1st intro note at t = 0'], 'FontSize', 14, 'FontWeight', 'bold');
plot([0 0], [0 Temp(4)], 'k--', 'LineWidth', 1.5);

axes(PSTPlotAxes2);
set(gca, 'FontSize', 12, 'FontWeight', 'bold');
set(gca, 'Box', 'off');
xlabel('Time (sec)', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('Firing Rate (Hz)', 'FontSize', 14, 'FontWeight', 'bold');
axis tight;
Temp = axis;
Temp(1:2) = [-PreBoutDuration 0];
axis(Temp);
plot([0 0], [0 Temp(4)], 'k--', 'LineWidth', 1.5);

save([FileInfo.FileNames{1}, '.motif_transitions.mat'], 'MotifRaster', 'OnsetsRaster', 'MotifPST', 'NoofINs', 'BoutInterval');
disp('Finished analysing notes');