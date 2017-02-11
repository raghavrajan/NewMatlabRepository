function [] = SSAAnalyseBoutTransitions(RawDataDir, RecFileDir, PreBoutDuration, PostBoutDuration, FileInfo, TitleString, Motif, Motif2, FileType, ContinuousOrNot)

InterBoutInterval = 2; % in seconds
GaussianLen = 2; 
Width = 0.01; % seconds

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
OnsetEdges = -PreBoutDuration:0.001:3;
OffsetEdges = -3:0.001:PostBoutDuration;

for i = 1:length(AllFiles),
    FileTime{i} = AllFiles(i).name;
end

BoutNo = 1;
BoutOnsetRaster = [];
BoutOffsetRaster = [];
BoutLength = [];
BoutOnsetPST = [];
BoutOffsetPST = [];
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
    Bouts = [];
    Bouts = find((Notes.onsets(2:end) - Notes.offsets(1:(end-1))) > InterBoutInterval);
    if (isempty(Bouts))
        continue;
    end
    for BoutIndex = 0:length(Bouts),
        if (BoutIndex == 0)
            StartIndex = 1;
        else
            StartIndex = Bouts(BoutIndex) + 1;
        end
        
        if (BoutIndex == length(Bouts))
            EndIndex = length(Notes.onsets);
        else
            EndIndex = Bouts(BoutIndex + 1);
        end
        if ~((Notes.onsets(StartIndex) > InterBoutInterval) && (Notes.offsets(EndIndex) < (Time(end) - InterBoutInterval)))
            continue;
        end
        Motifs = union(strfind(Notes.labels(StartIndex:EndIndex), Motif), strfind(Notes.labels(StartIndex:EndIndex), Motif2));
        if (isempty(Motifs))
            continue;
        else
            BoutLength(BoutNo) = Notes.offsets(EndIndex) - Notes.onsets(StartIndex);
            if (StartIndex == 1)
                SpikeTimeIndices = find((FileInfo.SpikeData.Times{NoteFileNo} >= Time(1)) & (FileInfo.SpikeData.Times{NoteFileNo} <= (Notes.onsets(StartIndex) + 3)));
                SpikeTimes = FileInfo.SpikeData.Times{NoteFileNo}(SpikeTimeIndices) - Notes.onsets(StartIndex);
                OffsetsRaster = [OffsetsRaster; [(Time(1) - Notes.onsets(StartIndex)) BoutNo]];
            else
                SpikeTimeIndices = find((FileInfo.SpikeData.Times{NoteFileNo} >= (Notes.offsets(StartIndex - 1))) & (FileInfo.SpikeData.Times{NoteFileNo} <= (Notes.onsets(StartIndex) + 3)));
                SpikeTimes = FileInfo.SpikeData.Times{NoteFileNo}(SpikeTimeIndices) - Notes.onsets(StartIndex);
                OffsetsRaster = [OffsetsRaster; [(Notes.offsets(StartIndex - 1) - Notes.onsets(StartIndex)) BoutNo]];
            end
            if (size(SpikeTimes, 1) < size(SpikeTimes,2))
                SpikeTimes = SpikeTimes';
            end        
            BoutOnsetRaster = [BoutOnsetRaster; [SpikeTimes ones(size(SpikeTimes))*BoutNo]];
            BoutOnsetPST(BoutNo,:) = histc(SpikeTimes, OnsetEdges);
 
            if (EndIndex == length(Notes.onsets))
                SpikeTimeIndices = find((FileInfo.SpikeData.Times{NoteFileNo} >= (Notes.offsets(EndIndex) - 3)) & (FileInfo.SpikeData.Times{NoteFileNo} <= Time(end)));
                SpikeTimes = FileInfo.SpikeData.Times{NoteFileNo}(SpikeTimeIndices) - Notes.offsets(EndIndex);            
                OnsetsRaster = [OnsetsRaster; [(Time(end) - Notes.offsets(EndIndex)) BoutNo]];
            else
                SpikeTimeIndices = find((FileInfo.SpikeData.Times{NoteFileNo} >= (Notes.offsets(EndIndex) - 3)) & (FileInfo.SpikeData.Times{NoteFileNo} <= (Notes.onsets(EndIndex + 1))));
                SpikeTimes = FileInfo.SpikeData.Times{NoteFileNo}(SpikeTimeIndices) - Notes.offsets(EndIndex);            
                OnsetsRaster = [OnsetsRaster; [(Notes.onsets(EndIndex + 1) - Notes.offsets(EndIndex)) BoutNo]];
            end
            if (size(SpikeTimes, 1) < size(SpikeTimes,2))
                SpikeTimes = SpikeTimes';
            end
            BoutOffsetRaster = [BoutOffsetRaster; [SpikeTimes ones(size(SpikeTimes))*BoutNo]];
            BoutOffsetPST(BoutNo,:) = histc(SpikeTimes, OffsetEdges);
            BoutNo = BoutNo + 1;
        end
    end
end

BoutOnsetPST = BoutOnsetPST/0.001;
BoutOffsetPST = BoutOffsetPST/0.001;

[SortedBoutLengths, SortedIndices] = sort(BoutLength);
for i = 1:length(SortedIndices),
    Indices = find(BoutOnsetRaster(:,2) == SortedIndices(i));
    BoutOnsetRaster(Indices,2) = SortedIndices(i) - 1000;
    Indices = find(OffsetsRaster(:,2) == SortedIndices(i));
    OffsetsRaster(Indices,2) = SortedIndices(i) - 1000;
    
    Indices = find(BoutOffsetRaster(:,2) == SortedIndices(i));
    BoutOffsetRaster(Indices,2) = SortedIndices(i) - 1000;
    Indices = find(OnsetsRaster(:,2) == SortedIndices(i));
    OnsetsRaster(Indices,2) = SortedIndices(i) - 1000;
end
BoutOnsetRaster(:,2) = BoutOnsetRaster(:,2) + 1000;
OffsetsRaster(:,2) = OffsetsRaster(:,2) + 1000;

BoutOffsetRaster(:,2) = BoutOffsetRaster(:,2) + 1000;
OnsetsRaster(:,2) = OnsetsRaster(:,2) + 1000;

BoutOnsetFig = figure;
RasterPlotAxes = axes('Position', [0.1 0.4 0.85 0.55]);
PSTPlotAxes = axes('Position', [0.1 0.1 0.85 0.25]);

figure(BoutOnsetFig);
hold on;
RasterIncrement = 0;

Fs = 1000;
XGauss = 1:1:(1 + round(2 * GaussianLen * Width * (Fs)));
XGauss = XGauss - (length(XGauss) + 1)/2;
GaussWin = (1/((Width * Fs) * sqrt(2 * pi))) * exp(-(XGauss.*XGauss)/(2 * (Width * Fs) * (Width * Fs)));


axes(RasterPlotAxes);
PlotRaster(BoutOnsetRaster, 'b', 2, RasterIncrement);
PlotRaster(OffsetsRaster, 'k', 4, RasterIncrement);
PST = [];
for i = 1:size(BoutOnsetPST,1),
    PST(i,:) = conv(BoutOnsetPST(i,:), GaussWin, 'same');
end
axes(PSTPlotAxes);
plot(OnsetEdges, mean(PST), 'b');
   
figure(BoutOnsetFig);
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
title([TitleString, 'Bout onset at time 0'], 'FontSize', 14, 'FontWeight', 'bold');
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

BoutOffsetFig = figure;
RasterPlotAxes2 = axes('Position', [0.1 0.4 0.85 0.55]);
PSTPlotAxes2 = axes('Position', [0.1 0.1 0.85 0.25]);

figure(BoutOffsetFig);
hold on;
RasterIncrement = 0;

Fs = 1000;
XGauss = 1:1:(1 + round(2 * GaussianLen * Width * (Fs)));
XGauss = XGauss - (length(XGauss) + 1)/2;
GaussWin = (1/((Width * Fs) * sqrt(2 * pi))) * exp(-(XGauss.*XGauss)/(2 * (Width * Fs) * (Width * Fs)));

% The 4 in the previous line refers to the fact that the gaussian window is
% terminated at 4 times the standard deviation

axes(RasterPlotAxes2);
PlotRaster(BoutOffsetRaster, 'b', 2, RasterIncrement);
PlotRaster(OnsetsRaster, 'k', 4, RasterIncrement);
PST = [];
for i = 1:size(BoutOffsetPST,1),
    PST(i,:) = conv(BoutOffsetPST(i,:), GaussWin, 'same');
end
axes(PSTPlotAxes2);
plot(OffsetEdges, mean(PST), 'b');

figure(BoutOffsetFig);
set(gcf, 'Position', [360 100 650 600]);
set(gcf, 'Color', 'w');
axes(RasterPlotAxes2); 
set(gca, 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Trial #', 'FontSize', 14, 'FontWeight', 'bold');
set(gca, 'XColor', 'w');
axis tight;
Temp = axis;
Temp(1:2) = [0 PostBoutDuration];
axis(Temp);
title([TitleString, 'Bout offset at t = 0'], 'FontSize', 14, 'FontWeight', 'bold');
plot([0 0], [0 Temp(4)], 'k--', 'LineWidth', 1.5);

axes(PSTPlotAxes2);
set(gca, 'FontSize', 12, 'FontWeight', 'bold');
set(gca, 'Box', 'off');
xlabel('Time (sec)', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('Firing Rate (Hz)', 'FontSize', 14, 'FontWeight', 'bold');
axis tight;
Temp = axis;
Temp(1:2) = [0 PostBoutDuration];
axis(Temp);
plot([0 0], [0 Temp(4)], 'k--', 'LineWidth', 1.5);

disp('Finished analysing notes');