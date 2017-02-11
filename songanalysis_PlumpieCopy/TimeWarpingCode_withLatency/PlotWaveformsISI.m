function [] = PlotWaveformsISI(DirSpikeData, UnDirSpikeData, MainFigure, WaveformPlot, ISIPlot)

figure(MainFigure);

Edges = 0:1:100;

DirISI = diff(DirSpikeData.Times * 1000);
UnDirISI = diff(UnDirSpikeData.Times * 1000);

if (length(DirSpikeData.Times) > 0)
    axes(ISIPlot(2));
    DirISIHIST = histc(DirISI,Edges);
    DirISIHIST = DirISIHIST * 100/length(DirISI);
    DirISIPLOT = bar(Edges,DirISIHIST,'histc');
    axis tight;
    set(gca,'Box','off');
    set(DirISIPLOT, 'EdgeColor', ' k', 'FaceColor', 'k');
    set(gca,'FontSize',8,'FontWeight','bold');
    xlabel('Time (ms)','FontSize',8,'FontWeight','bold');
    ylabel('Percentage of intervals','FontSize',8,'FontWeight','bold');
    title('ISI Histogram - Directed','FontSize',8,'FontWeight','bold');    
end

if (length(UnDirSpikeData.Times) > 0)
    axes(ISIPlot(1));
    UnDirISIHIST = histc(UnDirISI,Edges);
    UnDirISIHIST = UnDirISIHIST * 100/length(UnDirISI);
    UnDirISIPLOT = bar(Edges,UnDirISIHIST,'histc');
    axis tight;
    set(gca,'Box','off');
    set(UnDirISIPLOT, 'EdgeColor', [0.6 0.6 0.6], 'FaceColor', [0.6 0.6 0.6]);
    set(gca,'FontSize',8,'FontWeight','bold');
    xlabel('Time (ms)','FontSize',8,'FontWeight','bold');
    ylabel('Percentage of intervals','FontSize',8,'FontWeight','bold');    
    title('ISI Histogram - UnDirected','FontSize',8,'FontWeight','bold');    
end

if (length(UnDirSpikeData.SongSpikeWaveforms) > 0)
    axes(WaveformPlot);
%     errorbar((0:1/32:(1 - 1/32)),mean(UnDirSpikeData.Waveforms),std(UnDirSpikeData.Waveforms),'Color', [0.6 0.6 0.6]);
    plot((0:1/32:(1.75 - 1/32)),UnDirSpikeData.SongSpikeWaveforms,'Color', [0.6 0.6 0.6]);
    axis tight;
    set(gca,'Box','off');
    set(gca,'FontSize',8,'FontWeight','bold');
    set(gca,'XTicklabel',[]);
    hold on;
    if (length(DirSpikeData.SongSpikeWaveforms) > 0)
%         errorbar((1.5:1/32:(2.5 - 1/32)),mean(DirSpikeData.Waveforms),std(DirSpikeData.Waveforms),'k');
        plot((2.25:1/32:(4 - 1/32)),DirSpikeData.SongSpikeWaveforms,'k');
    end
    axis tight;
    xlabel('Time (msec)','FontSize',8,'FontWeight','bold');
    ylabel('Amplitude (\muV)','FontSize',8,'FontWeight','bold');    
    title('Spike Waveforms','FontSize',8,'FontWeight','bold');        
else
    if (length(DirSpikeData.SongSpikeWaveforms) > 0)
        axes(WaveformPlot);
        set(gca,'Box','off');
        set(gca,'FontSize',8,'FontWeight','bold');
        set(gca,'XTicklabel',[]);
        hold on;
%         errorbar((1.5:1/32:(2.5 - 1/32)),mean(DirSpikeData.Waveforms),std(DirSpikeData.Waveforms),'k');
        plot((0:1/32:(1.75 - 1/32)),DirSpikeData.SongSpikeWaveforms,'k');
        axis tight;
        xlabel('Time (msec)','FontSize',8,'FontWeight','bold');
        ylabel('Amplitude (\muV)','FontSize',8,'FontWeight','bold');    
        title('Spike Waveforms','FontSize',8,'FontWeight','bold');        
    end
end

disp('Finished plotting ISI histograms and average waveforms');