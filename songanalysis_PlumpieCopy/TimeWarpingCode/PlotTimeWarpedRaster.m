function [DirFileInfo, UnDirFileInfo, MedianMotif] = PlotTimeWarpedRaster(DirectoryName, FileName, FileExt, FileNos, SpikeSortOrThreshold, FileType, ChannelNo, Motif, BinSize, Latency)

% Code to plot singing related activity aligned and time warped to one of
% the motifs of the different song bouts for directed and undirected song

PresentDirectory = pwd;

if (DirectoryName(end) ~= '/')
    DirectoryName = [DirectoryName,'/'];
end

DirFileInfo.BinSize = BinSize;
DirFileInfo.Latency = Latency;
DirFileInfo.ChannelNo = ChannelNo;
DirFileInfo.Motif = Motif;

UnDirFileInfo.BinSize = BinSize;
UnDirFileInfo.Latency = Latency;
UnDirFileInfo.ChannelNo = ChannelNo;
UnDirFileInfo.Motif = Motif;

% =========================================================================
% Variables and initialisation of the variables
% =========================================================================
% Variables related to the directed song bouts
% DirFileInfo is a structure that will contain all the data

% Variables related to the undirected song bouts

% =========================================================================
% Get File name and record length
% =========================================================================

[FileNames, RecordLengths] = GetFileNames(DirectoryName,FileName,FileExt,FileNos.Directed,FileType);
DirFileInfo.FileNames = FileNames;
DirFileInfo.RecordLengths = RecordLengths;

clear FileNames;
clear RecordLengths;

[FileNames, RecordLengths] = GetFileNames(DirectoryName,FileName,FileExt,FileNos.Undirected, FileType);
UnDirFileInfo.FileNames = FileNames;
UnDirFileInfo.RecordLengths = RecordLengths;

clear FileNames RecordLengths;

% =========================================================================
% Load up information about syllable onsets and offsets and spike times
% =========================================================================

% Directed song bouts
[NoteOnsets, NoteOffsets, NoteLabels] = LoadNoteFiles(DirectoryName, DirFileInfo.FileNames,DirFileInfo.RecordLengths);

DirFileInfo.Notes.NoteOnsets = NoteOnsets;
DirFileInfo.Notes.NoteOffsets = NoteOffsets;
DirFileInfo.Notes.NoteLabels = NoteLabels;

clear NoteOnsets NoteOffsets NoteLabels;

[Syllables, Gaps] = ProcessNoteFiles(DirFileInfo, Motif);

DirFileInfo.Syllables = Syllables;
DirFileInfo.Gaps = Gaps;

clear Syllables Gaps;

if (strfind(SpikeSortOrThreshold,'spikesort'))
    [SpikeTimes,SpikeAmplitudes,SpikeWaveforms,DirFileInfo.ClusterParameters] = LoadSpikeTimes(DirectoryName,DirFileInfo.FileNames,ChannelNo,DirFileInfo.RecordLengths, FileType);
    DirFileInfo.SpikeData.Times = SpikeTimes;
    DirFileInfo.SpikeData.Waveforms = SpikeWaveforms;
    DirFileInfo.SpikeData.Amplitudes = SpikeAmplitudes;
end

if (strfind(SpikeSortOrThreshold,'threshold'))
    [SpikeTimes,SpikeAmplitudes,SpikeWaveforms,DirFileInfo.ThresholdingParameters] = ThresholdSpikeTimes(DirectoryName,DirFileInfo.FileNames,ChannelNo,DirFileInfo.RecordLengths, FileType);
    DirFileInfo.SpikeData.Times = SpikeTimes;
    DirFileInfo.SpikeData.Waveforms = SpikeWaveforms;
    DirFileInfo.SpikeData.Amplitudes = SpikeAmplitudes;
end

clear SpikeTimes SpikeWaveforms SpikeAmplitudes;

% Undirected song bouts
[NoteOnsets, NoteOffsets, NoteLabels] = LoadNoteFiles(DirectoryName, UnDirFileInfo.FileNames,UnDirFileInfo.RecordLengths);
UnDirFileInfo.Notes.NoteOnsets = NoteOnsets;
UnDirFileInfo.Notes.NoteOffsets = NoteOffsets;
UnDirFileInfo.Notes.NoteLabels = NoteLabels;

clear NoteOnsets NoteOffsets NoteLabels;

[Syllables, Gaps] = ProcessNoteFiles(UnDirFileInfo, Motif);

UnDirFileInfo.Syllables = Syllables;
UnDirFileInfo.Gaps = Gaps;

clear Syllables Gaps;

if (strfind(SpikeSortOrThreshold,'spikesort'))
    [SpikeTimes,SpikeAmplitudes,SpikeWaveforms,UnDirFileInfo.ClusterParameters] = LoadSpikeTimes(DirectoryName,UnDirFileInfo.FileNames,ChannelNo,UnDirFileInfo.RecordLengths, FileType);
    UnDirFileInfo.SpikeData.Times = SpikeTimes;
    UnDirFileInfo.SpikeData.Waveforms = SpikeWaveforms;
    UnDirFileInfo.SpikeData.Amplitudes = SpikeAmplitudes;
end

if (strfind(SpikeSortOrThreshold,'threshold'))
    [SpikeTimes,SpikeAmplitudes,SpikeWaveforms,UnDirFileInfo.ThresholdingParameters] = ThresholdSpikeTimes(DirectoryName,UnDirFileInfo.FileNames,ChannelNo,UnDirFileInfo.RecordLengths, FileType);
    UnDirFileInfo.SpikeData.Times = SpikeTimes;
    UnDirFileInfo.SpikeData.Waveforms = SpikeWaveforms;
    UnDirFileInfo.SpikeData.Amplitudes = SpikeAmplitudes;
end

clear SpikeTimes SpikeWaveforms SpikeAmplitudes;

disp('Finished loading up the syllable onset, offset information');

% =========================================================================
% Now, plot the average waveforms and the ISI histogram and the PST and
% raster for both directed and undirected song. Along with this, i also
% plot the spectrogram of the motif with median length
% =========================================================================
% Raster Plot Figure with raster, pst, spectrogram, average spike
% waveforms, spike amplitudes and ISI

RasterPlotFigure = figure;
set(gcf,'Color','w');
set(RasterPlotFigure,'Position',[354 35 645 910])

MotifSpectrogramPlot = axes('Position',[0.1 0.8 0.57 0.15]);
RasterPlot = axes('Position',[0.1 0.1 0.57 0.65]);
set(gca,'Box','off');

WaveFormPlot = axes('Position',[0.75 0.4 0.2 0.1]);
set(gca,'Box','off');

ISIPlot(1) = axes('Position',[0.75 0.25 0.2 0.1]);
set(gca,'Box','off');
ISIPlot(2) = axes('Position',[0.75 0.1 0.2 0.1]);
set(gca,'Box','off');

SpikeAmplitudePlot(1) = axes('Position',[0.75 0.6 0.2 0.15]);
SpikeAmplitudePlot(2) = axes('Position',[0.75 0.8 0.2 0.15]);

MedianMotif = CalculateMedianMotifLength(DirFileInfo,UnDirFileInfo,Motif,RasterPlotFigure);

[Raster, PST, Edges, WSpikeTrain, UWSpikeTrain, SpikeWaveforms] = WarpSpikeTrains(DirFileInfo, Motif, MedianMotif, BinSize, Latency);
DirFileInfo.SpikeRaster = Raster;
DirFileInfo.PST = PST;
DirFileInfo.Edges = Edges;
DirFileInfo.WSpikeTrain = WSpikeTrain;
DirFileInfo.UWSpikeTrain = UWSpikeTrain;
DirFileInfo.SpikeData.SongSpikeWaveforms = SpikeWaveforms;
clear Raster PST Edges WSpikeTrain UWSpikeTrain SpikeWaveforms;

[Raster,PST, Edges, WSpikeTrain, UWSpikeTrain, SpikeWaveforms] = WarpSpikeTrains(UnDirFileInfo, Motif, MedianMotif, BinSize, Latency);
UnDirFileInfo.SpikeRaster = Raster;
UnDirFileInfo.PST = PST;
UnDirFileInfo.Edges = Edges;
UnDirFileInfo.WSpikeTrain = WSpikeTrain;
UnDirFileInfo.UWSpikeTrain = UWSpikeTrain;
UnDirFileInfo.SpikeData.SongSpikeWaveforms = SpikeWaveforms;
clear Raster PST Edges WSpikeTrain UWSpikeTrain SpikeWaveforms;

PlotMotifSpectrogram(DirectoryName, MedianMotif, RasterPlotFigure, MotifSpectrogramPlot, FileType);
PlotRasterPST(DirFileInfo, UnDirFileInfo, MedianMotif, BinSize, RasterPlotFigure, RasterPlot);
PlotWaveformsISI(DirFileInfo.SpikeData, UnDirFileInfo.SpikeData, RasterPlotFigure, WaveFormPlot, ISIPlot);
PlotSpikeAmplitudes(DirFileInfo, UnDirFileInfo, MedianMotif, RasterPlotFigure, SpikeAmplitudePlot);
annotation('textbox',[0.375 0.97 0.25 0.025],'FontSize',16,'FontWeight','bold','String',upper(FileName),'LineStyle','none');

if (strfind(SpikeSortOrThreshold,'threshold'))
    annotation('textbox',[0.05 0.05 0.8 0.0025],'FontSize',8,'FontWeight','bold','String',['Lower Threshold is ',[DirFileInfo.ThresholdingParameters{1}],', ',[UnDirFileInfo.ThresholdingParameters{1}],' and Upper Threshold is ',[DirFileInfo.ThresholdingParameters{2}],', ',[UnDirFileInfo.ThresholdingParameters{2}]],'LineStyle','none');
else
    if (strfind(SpikeSortOrThreshold,'spikesort'))    
        annotation('textbox',[0.05 0.05 0.8 0.0025],'FontSize',8,'FontWeight','bold','String',['Cluster Nos are ',[DirFileInfo.ClusterParameters{1}],', ',[UnDirFileInfo.ClusterParameters{1}],', Total number of clusters is ',[DirFileInfo.ClusterParameters{2}],', ',[UnDirFileInfo.ClusterParameters{2}],' and Outliers include is set at ',[DirFileInfo.ClusterParameters{3}],', ',[UnDirFileInfo.ClusterParameters{3}]],'LineStyle','none');    
    end
end

annotation('textbox',[0.05 0.03 0.8 0.01],'FontSize',8,'FontWeight','bold','String',['File Nos are Directed:',num2str(FileNos.Directed),' and Undirected:',num2str(FileNos.Undirected)],'LineStyle','none');

% =========================================================================
% Raw Waveform figure - has 3 representative waveforms of directed and
% undirected song

RawWaveformFigure = figure;
set(gcf,'Color','w');
set(RawWaveformFigure,'Position',[354 35 645 910])
PlotRawWaveforms(DirectoryName, DirFileInfo, UnDirFileInfo, ChannelNo, RawWaveformFigure, FileType);
annotation('textbox',[0.375 0.97 0.25 0.025],'FontSize',16,'FontWeight','bold','String',upper(FileName),'LineStyle','none');

if (strfind(SpikeSortOrThreshold,'threshold'))
    annotation('textbox',[0.05 0.08 0.8 0.0025],'FontSize',8,'FontWeight','bold','String',['Lower Threshold is ',[DirFileInfo.ThresholdingParameters{1}],', ',[UnDirFileInfo.ThresholdingParameters{1}],' and Upper Threshold is ',[DirFileInfo.ThresholdingParameters{2}],', ',[UnDirFileInfo.ThresholdingParameters{2}]],'LineStyle','none');
else
    if (strfind(SpikeSortOrThreshold,'spikesort'))    
        annotation('textbox',[0.05 0.08 0.8 0.0025],'FontSize',8,'FontWeight','bold','String',['Cluster Nos are ',[DirFileInfo.ClusterParameters{1}],', ',[UnDirFileInfo.ClusterParameters{1}],', Total number of clusters is ',[DirFileInfo.ClusterParameters{2}],', ',[UnDirFileInfo.ClusterParameters{2}],' and Outliers include is set at ',[DirFileInfo.ClusterParameters{3}],', ',[UnDirFileInfo.ClusterParameters{3}]],'LineStyle','none');    
    end
end

annotation('textbox',[0.05 0.04 0.8 0.01],'FontSize',8,'FontWeight','bold','String',['File Nos are Directed:',num2str(FileNos.Directed),' and Undirected:',num2str(FileNos.Undirected)],'LineStyle','none');

% =========================================================================
% Detailed Analyis figure - has the length of motifs for directed and
% undirected song, length of different syllables and gaps for both directed
% and undirected song, pairwise correlations between spike trains for each
% rendition of the motif.
% This figure also has the graphs for no of spikes/event, spike
% jitter/event and the ISI statistics for each event for both directed and
% undirected song.

DetailedAnalysisFigure = figure;
set(gcf,'Color','w');
set(DetailedAnalysisFigure,'Position',[354 35 645 910])
annotation('textbox',[0.375 0.97 0.25 0.025],'FontSize',16,'FontWeight','bold','String',upper(FileName),'LineStyle','none');

if (strfind(SpikeSortOrThreshold,'threshold'))
    annotation('textbox',[0.05 0.08 0.8 0.0025],'FontSize',8,'FontWeight','bold','String',['Lower Threshold is ',[DirFileInfo.ThresholdingParameters{1}],', ',[UnDirFileInfo.ThresholdingParameters{1}],' and Upper Threshold is ',[DirFileInfo.ThresholdingParameters{2}],', ',[UnDirFileInfo.ThresholdingParameters{2}]],'LineStyle','none');
else
    if (strfind(SpikeSortOrThreshold,'spikesort'))    
        annotation('textbox',[0.05 0.08 0.8 0.0025],'FontSize',8,'FontWeight','bold','String',['Cluster Nos are ',[DirFileInfo.ClusterParameters{1}],', ',[UnDirFileInfo.ClusterParameters{1}],', Total number of clusters is ',[DirFileInfo.ClusterParameters{2}],', ',[UnDirFileInfo.ClusterParameters{2}],' and Outliers include is set at ',[DirFileInfo.ClusterParameters{3}],', ',[UnDirFileInfo.ClusterParameters{3}]],'LineStyle','none');    
    end
end

annotation('textbox',[0.05 0.04 0.8 0.01],'FontSize',8,'FontWeight','bold','String',['File Nos are Directed:',num2str(FileNos.Directed),' and Undirected:',num2str(FileNos.Undirected)],'LineStyle','none');

% First calculate pairwise correlations between spike trains and the event
% related statistics.

[Correlation] = CalculateIFRCorr(DirFileInfo.WSpikeTrain, MedianMotif, Latency);
DirFileInfo.IFRCorrelation = Correlation;
clear Correlation;

[Correlation] = CalculateIFRCorr(UnDirFileInfo.WSpikeTrain, MedianMotif, Latency);
UnDirFileInfo.IFRCorrelation = Correlation;
clear Correlation;

Index = 1;
for GaussianWidth = [1 5 10 25],
    [Correlation] = CalculateCorrGaussSmooth(DirFileInfo.WSpikeTrain, MedianMotif, Latency, GaussianWidth/1000);
    DirFileInfo.GaussianCorrelations(Index,:) = [GaussianWidth Correlation];
    clear Correlation;
    Index = Index + 1;
end

Index = 1;
for GaussianWidth = [1 5 10 25],
    [Correlation] = CalculateCorrGaussSmooth(UnDirFileInfo.WSpikeTrain, MedianMotif, Latency, GaussianWidth/1000);
    UnDirFileInfo.GaussianCorrelations(Index,:) = [GaussianWidth Correlation];
    clear Correlation;
    Index = Index + 1;
end

if (length(DirFileInfo.PST) > 0)
    [EventParameters] = CalculateEventParametersWarpedSpikeTrain(DirFileInfo, MedianMotif, BinSize, Latency, UnDirFileInfo);
    DirFileInfo.WEventParameters = EventParameters;
    clear EventParameters;
    [EventParameters] = CalculateEventParametersUnWarpedSpikeTrain(DirFileInfo, MedianMotif, BinSize, Latency, UnDirFileInfo);
    DirFileInfo.UWEventParameters = EventParameters;
    clear EventParameters;
    [FanoFactor_30, FanoFactor_100] = CalculateFanoFactor(DirFileInfo,MedianMotif,Latency);
    DirFileInfo.FanoFactor_30 = FanoFactor_30;
    DirFileInfo.FanoFactor_100 = FanoFactor_100;    
    clear FanoFactor_30 FanoFactor_100;
else
    DirFileInfo.UWEventParameters = [];
    DirFileInfo.WEventParameters = [];
end

if (length(UnDirFileInfo.PST) > 0)
    [EventParameters] = CalculateEventParametersWarpedSpikeTrain(UnDirFileInfo, MedianMotif, BinSize, Latency, DirFileInfo);
    UnDirFileInfo.WEventParameters = EventParameters;
    clear EventParameters;
    [EventParameters] = CalculateEventParametersUnWarpedSpikeTrain(UnDirFileInfo, MedianMotif, BinSize, Latency, DirFileInfo);
    UnDirFileInfo.UWEventParameters = EventParameters;
    clear EventParameters;
    [FanoFactor_30, FanoFactor_100] = CalculateFanoFactor(UnDirFileInfo,MedianMotif,Latency);
    UnDirFileInfo.FanoFactor_30 = FanoFactor_30;
    UnDirFileInfo.FanoFactor_100 = FanoFactor_100;
    clear FanoFactor_30 FanoFactor_100;
else
    UnDirFileInfo.WEventParameters = [];
    UnDirFileInfo.UWEventParameters = [];
end

if (length(DirFileInfo.Syllables) > 0)
    DirFileInfo.SongLengths = DirFileInfo.Syllables.End(:,end) - DirFileInfo.Syllables.Start(:,1);
end

if (length(UnDirFileInfo.Syllables) > 0)
    UnDirFileInfo.SongLengths = UnDirFileInfo.Syllables.End(:,end) - UnDirFileInfo.Syllables.Start(:,1);
end


% Now do the plots

PlotCorrelations(DirFileInfo, UnDirFileInfo, DetailedAnalysisFigure);

PlotFanoFactor(DirFileInfo, UnDirFileInfo, MedianMotif, DetailedAnalysisFigure);

DetailedAnalysisMotifSpectrogramPlot = axes('Position',[0.15 0.25 0.8 0.1]);
PlotMotifSpectrogram(DirectoryName, MedianMotif, DetailedAnalysisFigure, DetailedAnalysisMotifSpectrogramPlot, FileType);
axis([-Latency MedianMotif.Length 300 10000]);

if ((isfield(DirFileInfo,'WSpikeTrain')) || (isfield(UnDirFileInfo,'WSpikeTrain')))
    [DirMotifFiringRate, UnDirMotifFiringRate] = PlotMotifFiringRate(DirFileInfo, UnDirFileInfo, DetailedAnalysisFigure);
    DirFileInfo.MotifFiringRate = DirMotifFiringRate;
    UnDirFileInfo.MotifFiringRate = UnDirMotifFiringRate;
    clear DirMotifFiringRate UnDirMotifFiringRate;
end

if ((isfield(DirFileInfo,'UWEventParameters')) || (isfield(UnDirFileInfo,'UWEventParameters')))
    PlotEventParameters(DirFileInfo, UnDirFileInfo, MedianMotif, DetailedAnalysisFigure);
%    PlotEventWaveforms(DirectoryName, ChannelNo, DirFileInfo, UnDirFileInfo, MedianMotif, DetailedAnalysisFigure);
end

% =========================================================================
% Song Analyis figure - has the length of motifs for directed and
% undirected song, length of different syllables and gaps for both directed
% and undirected song

SongAnalysisFigure = figure;
set(gcf,'Color','w');
set(SongAnalysisFigure,'Position',[354 35 645 910])
annotation('textbox',[0.375 0.97 0.25 0.025],'FontSize',16,'FontWeight','bold','String',upper(FileName),'LineStyle','none');

if (strfind(SpikeSortOrThreshold,'threshold'))
    annotation('textbox',[0.05 0.08 0.8 0.0025],'FontSize',8,'FontWeight','bold','String',['Lower Threshold is ',[DirFileInfo.ThresholdingParameters{1}],', ',[UnDirFileInfo.ThresholdingParameters{1}],' and Upper Threshold is ',[DirFileInfo.ThresholdingParameters{2}],', ',[UnDirFileInfo.ThresholdingParameters{2}]],'LineStyle','none');
else
    if (strfind(SpikeSortOrThreshold,'spikesort'))    
        annotation('textbox',[0.05 0.08 0.8 0.0025],'FontSize',8,'FontWeight','bold','String',['Cluster Nos are ',[DirFileInfo.ClusterParameters{1}],', ',[UnDirFileInfo.ClusterParameters{1}],', Total number of clusters is ',[DirFileInfo.ClusterParameters{2}],', ',[UnDirFileInfo.ClusterParameters{2}],' and Outliers include is set at ',[DirFileInfo.ClusterParameters{3}],', ',[UnDirFileInfo.ClusterParameters{3}]],'LineStyle','none');    
    end
end

annotation('textbox',[0.05 0.04 0.8 0.01],'FontSize',8,'FontWeight','bold','String',['File Nos are Directed:',num2str(FileNos.Directed),' and Undirected:',num2str(FileNos.Undirected)],'LineStyle','none');

if ((length(DirFileInfo.Syllables) > 0) && (length(UnDirFileInfo.Syllables) > 0))
    PlotSongLengths(DirFileInfo, UnDirFileInfo, MedianMotif, SongAnalysisFigure);
    PlotSyllableGapStatistics(DirFileInfo, UnDirFileInfo, Motif, MedianMotif, SongAnalysisFigure);
end

% =========================================================================
% Bout figure - the spike times are plotted as a function of bouts with the
% labels for the various sounds added into the figure.

% BoutFigure = figure;
% set(gcf,'Color','w');
% set(BoutFigure,'Position',[354 35 645 910])
% annotation('textbox',[0.375 0.97 0.25 0.025],'FontSize',16,'FontWeight','bold','String',upper(FileName),'LineStyle','none');
% 
% if (strfind(SpikeSortOrThreshold,'threshold'))
%     annotation('textbox',[0.05 0.08 0.8 0.0025],'FontSize',8,'FontWeight','bold','String',['Lower Threshold is ',[DirFileInfo.ThresholdingParameters{1}],', ',[UnDirFileInfo.ThresholdingParameters{1}],' and Upper Threshold is ',[DirFileInfo.ThresholdingParameters{2}],', ',[UnDirFileInfo.ThresholdingParameters{2}]],'LineStyle','none');
% else
%     if (strfind(SpikeSortOrThreshold,'spikesort'))    
%         annotation('textbox',[0.05 0.08 0.8 0.0025],'FontSize',8,'FontWeight','bold','String',['Cluster Nos are ',[DirFileInfo.ClusterParameters{1}],', ',[UnDirFileInfo.ClusterParameters{1}],', Total number of clusters is ',[DirFileInfo.ClusterParameters{2}],', ',[UnDirFileInfo.ClusterParameters{2}],' and Outliers include is set at ',[DirFileInfo.ClusterParameters{3}],', ',[UnDirFileInfo.ClusterParameters{3}]],'LineStyle','none');    
%     end
% end
% 
% annotation('textbox',[0.05 0.04 0.8 0.01],'FontSize',8,'FontWeight','bold','String',['File Nos are Directed:',num2str(FileNos.Directed),' and Undirected:',num2str(FileNos.Undirected)],'LineStyle','none');
% 
% PlotBoutSpikeTimes(DirFileInfo,UnDirFileInfo,Motif,BoutFigure);

cd(PresentDirectory);
