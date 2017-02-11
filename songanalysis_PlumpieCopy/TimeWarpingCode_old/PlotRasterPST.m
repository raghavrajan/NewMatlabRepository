function [] = PlotRasterPST(DirFileInfo, UnDirFileInfo, MedianMotif, BinSize, MainFigure,RasterPlot)

figure(MainFigure);
axes(RasterPlot);
hold on;

PreDataTime = 0.05;
Edges = -PreDataTime:BinSize:(MedianMotif.Length);
Edges = Edges + BinSize/2;

plot(Edges, mean(DirFileInfo.PST), 'k');
plot(Edges, mean(UnDirFileInfo.PST),'Color',[0.6 0.6 0.6]);

MaxPST = max(max(mean(DirFileInfo.PST)),max(mean(UnDirFileInfo.PST)));

if (length(DirFileInfo.SpikeRaster) > 0)
    DirFileInfo.SpikeRaster(:,2) = DirFileInfo.SpikeRaster(:,2) * MaxPST/20;
    DirFileInfo.SpikeRaster(:,2) = DirFileInfo.SpikeRaster(:,2) + MaxPST + MaxPST/10;
end

if (length(UnDirFileInfo.SpikeRaster) > 0)
    UnDirFileInfo.SpikeRaster(:,2) = UnDirFileInfo.SpikeRaster(:,2) * MaxPST/20;
    if (length(DirFileInfo.SpikeRaster) > 0)
        UnDirFileInfo.SpikeRaster(:,2) = UnDirFileInfo.SpikeRaster(:,2) + DirFileInfo.SpikeRaster(end,2) + MaxPST/3;
    else
        UnDirFileInfo.SpikeRaster(:,2) = UnDirFileInfo.SpikeRaster(:,2) + MaxPST + MaxPST/10;
    end
end

if (length(DirFileInfo.SpikeRaster) > 0)
    plot(DirFileInfo.SpikeRaster(:,1),DirFileInfo.SpikeRaster(:,2),'w+');
    MarkerString = repmat('|',size(DirFileInfo.SpikeRaster(:,1),1),1);
    text(DirFileInfo.SpikeRaster(:,1),DirFileInfo.SpikeRaster(:,2),MarkerString,'HorizontalAlignment','center','VerticalAlignment','middle','FontSize',2,'FontName','FixedWidth','FontUnits','pixels','Color','k','FontWeight','bold');
end

if (length(UnDirFileInfo.SpikeRaster) > 0)
    plot(UnDirFileInfo.SpikeRaster(:,1),UnDirFileInfo.SpikeRaster(:,2),'w+');
    MarkerString = repmat('|',size(UnDirFileInfo.SpikeRaster(:,1),1),1);
%    text(UnDirFileInfo.SpikeRaster(:,1),UnDirFileInfo.SpikeRaster(:,2),MarkerString,'HorizontalAlignment','center','VerticalAlignment','middle','FontSize',2,'FontName','FixedWidth','FontUnits','pixels','Color',[0.6 0.6 0.6],'FontWeight','bold');
    text(UnDirFileInfo.SpikeRaster(:,1),UnDirFileInfo.SpikeRaster(:,2),MarkerString,'HorizontalAlignment','center','VerticalAlignment','middle','FontSize',2,'FontName','FixedWidth','FontUnits','pixels','Color','k','FontWeight','bold');
end

axis([-PreDataTime (MedianMotif.Length) 0 UnDirFileInfo.SpikeRaster(end,2)]);

set(gca,'FontSize',14,'FontWeight','bold');
set(gca,'ytick',[0 round(MaxPST/2) round(MaxPST)]);

xlabel('Time (sec)','FontSize',16,'FontWeight','bold');
ylabel('Firing Rate (Hz)','FontSize',16,'FontWeight','bold');
