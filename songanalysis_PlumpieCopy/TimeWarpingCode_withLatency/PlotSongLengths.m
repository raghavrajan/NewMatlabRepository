function [] = PlotSongLengths(DirFileInfo, UnDirFileInfo, MedianMotif, MainFigure)

figure(MainFigure);

axes('Position',[0.15 0.15 0.2 0.2]);
set(gca,'Box','off');
hold on;

MaxSongLength = 0;
MinSongLength = 1;

if (isfield(DirFileInfo,'SongLengths'))
    plot(ones(length(DirFileInfo.SongLengths),1),DirFileInfo.SongLengths,'k+','MarkerSize',2);
    MaxSongLength = max(max(DirFileInfo.SongLengths),MaxSongLength);
    MinSongLength = min(min(DirFileInfo.SongLengths),MinSongLength);    
end

if (isfield(UnDirFileInfo,'SongLengths'))
    plot((ones(length(UnDirFileInfo.SongLengths),1)*2),UnDirFileInfo.SongLengths,'Marker','+','LineStyle','none','MarkerEdgeColor',[0.6 0.6 0.6],'MarkerSize',2);
    MaxSongLength = max(max(UnDirFileInfo.SongLengths),MaxSongLength);
    MinSongLength = min(min(UnDirFileInfo.SongLengths),MinSongLength);    
end

hold on;
plot([0.5 2.5],[MedianMotif.Length MedianMotif.Length],'k--');
axis([0.5 2.5 (MinSongLength*0.95) (MaxSongLength*1.05)]);
set(gca,'FontSize',10,'FontWeight','bold');
set(gca,'Xtick',[1 2],'XTickLabel',[{'Directed'} {'Undirected'}]);
ylabel('Time (sec)','FontSize',10,'FontWeight','bold');    
title('Motif Lengths','FontSize',12,'FontWeight','bold');    