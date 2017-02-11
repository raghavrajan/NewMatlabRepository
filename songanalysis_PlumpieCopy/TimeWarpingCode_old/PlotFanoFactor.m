function [] = PlotFanoFactor(DirFileInfo, UnDirFileInfo, MedianMotif, MainFigure)

figure(MainFigure);

axes('Position',[0.15 0.15 0.8 0.1]);
set(gca,'Box','off');
hold on;

Edges = 0:0.001:(MedianMotif.Length + DirFileInfo.Latency);
Edges = Edges - DirFileInfo.Latency;

MaxFanoFactor = 0;

if (isfield(DirFileInfo,'FanoFactor'))
    plot(Edges, DirFileInfo.FanoFactor, 'k');
    MaxFanoFactor = max(MaxFanoFactor, max(DirFileInfo.FanoFactor));
end

hold on;

if (isfield(UnDirFileInfo,'FanoFactor'))
    plot(Edges, UnDirFileInfo.FanoFactor,'Color',[0.6 0.6 0.6]);
    MaxFanoFactor = max(MaxFanoFactor, max(UnDirFileInfo.FanoFactor));    
end

axis([-DirFileInfo.Latency (MedianMotif.Length) 0 MaxFanoFactor]);

set(gca,'FontSize',14,'FontWeight','bold');

xlabel('Time (sec)','FontSize',16,'FontWeight','bold');
ylabel('Fano Factor','FontSize',16,'FontWeight','bold');
