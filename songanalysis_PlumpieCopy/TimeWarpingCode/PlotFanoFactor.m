function [] = PlotFanoFactor(DirFileInfo, UnDirFileInfo, MedianMotif, MainFigure)

figure(MainFigure);

axes('Position',[0.15 0.15 0.8 0.1]);
set(gca,'Box','off');
hold on;

Edges = -0.1:0.001:(MedianMotif.Length + DirFileInfo.Latency - 0.03);
Edges = Edges - DirFileInfo.Latency;

MaxFanoFactor = 0;

if (isfield(DirFileInfo,'FanoFactor_30'))
    plot(Edges, DirFileInfo.FanoFactor_30, 'k');
    MaxFanoFactor = max(MaxFanoFactor, max(DirFileInfo.FanoFactor_30));
end

hold on;

if (isfield(UnDirFileInfo,'FanoFactor_30'))
    plot(Edges, UnDirFileInfo.FanoFactor_30,'Color',[0.6 0.6 0.6]);
    MaxFanoFactor = max(MaxFanoFactor, max(UnDirFileInfo.FanoFactor_30));    
end

axis([(DirFileInfo.Latency -0.1) (MedianMotif.Length) 0 MaxFanoFactor]);

set(gca,'FontSize',14,'FontWeight','bold');

xlabel('Time (sec)','FontSize',16,'FontWeight','bold');
ylabel('Fano Factor','FontSize',16,'FontWeight','bold');
