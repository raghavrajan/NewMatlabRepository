function [] = PlotSpikeAmplitudes(DirFileInfo, UnDirFileInfo, MedianMotif, MainFigure, SpikeAmplitudePlot)

figure(MainFigure);

if (length(DirFileInfo.SpikeRaster) > 0)
    axes(SpikeAmplitudePlot(1));
    hold on;

    plot(DirFileInfo.SpikeRaster(:,1),DirFileInfo.SpikeRaster(:,3),'k+','MarkerSize',2);
    axis tight;
    temp1 = axis;

    set(gca,'FontSize',8,'FontWeight','bold');
    set(gca,'Box','off');
    xlabel('Time (sec)','FontSize',8,'FontWeight','bold');
    ylabel('Amplitude (\muV)','FontSize',8,'FontWeight','bold');        
    title('Spike Amplitudes - Directed','FontSize',8,'FontWeight','bold');        
end

if (length(UnDirFileInfo.SpikeRaster) > 0)
    axes(SpikeAmplitudePlot(2));
    plot(UnDirFileInfo.SpikeRaster(:,1),UnDirFileInfo.SpikeRaster(:,3),'LineStyle','none','Marker','+','MarkerEdgeColor', [0.6 0.6 0.6],'MarkerSize',2);
    set(gca,'FontSize',8,'FontWeight','bold');
    set(gca,'Box','off');
    xlabel('Time (sec)','FontSize',8,'FontWeight','bold');    
    ylabel('Amplitude (\muV)','FontSize',8,'FontWeight','bold');    
    axis tight;
    temp2 = axis;
    title('Spike Amplitudes - UnDirected','FontSize',8,'FontWeight','bold');    
end

if ((length(DirFileInfo.SpikeRaster) > 0) && (length(UnDirFileInfo.SpikeRaster) > 0))
    axes(SpikeAmplitudePlot(1));
    axis([-0.2 (MedianMotif.Length + 0.2) min(temp1(3),temp2(3)) max(temp1(4),temp2(4))]);

    axes(SpikeAmplitudePlot(2));
    axis([-0.2 (MedianMotif.Length + 0.2) min(temp1(3),temp2(3)) max(temp1(4),temp2(4))]);
end
