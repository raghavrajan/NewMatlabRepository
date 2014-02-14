function [] = PlotMotifFiringRate(DirFileInfo, UnDirFileInfo, MainFigure)

figure(MainFigure);

axes('Position',[0.15 0.82 0.15 0.1]);
set(gca,'Box','off');
hold on;

MaxNoofSpikes = 0;

if (isfield(DirFileInfo,'SpikeTrain'))
    if (length(DirFileInfo.SpikeTrain) > 0)
        for i = 1:length(DirFileInfo.SpikeTrain),
            DirMotifNoofSpikes(i) = length([DirFileInfo.SpikeTrain{i}]);
        end
        DirNoofSpikes = bar(1,mean(DirMotifNoofSpikes));
        set(DirNoofSpikes, 'EdgeColor', ' k', 'FaceColor', 'k');
        errorbar(1,mean(DirMotifNoofSpikes),std(DirMotifNoofSpikes),'k.');
        MaxNoofSpikes = max(MaxNoofSpikes,(mean(DirMotifNoofSpikes) + std(DirMotifNoofSpikes)));
    end
end

if (isfield(UnDirFileInfo,'SpikeTrain'))
    if (length(UnDirFileInfo.SpikeTrain) > 0)
        for i = 1:length(UnDirFileInfo.SpikeTrain),
            UnDirMotifNoofSpikes(i) = length([UnDirFileInfo.SpikeTrain{i}]);
        end
        UnDirNoofSpikes = bar(2,mean(UnDirMotifNoofSpikes));
        set(UnDirNoofSpikes, 'EdgeColor', [0.6 0.6 0.6], 'FaceColor', [0.6 0.6 0.6]);
        errorbar(2,mean(UnDirMotifNoofSpikes),std(UnDirMotifNoofSpikes),'Color',[0.6 0.6 0.6],'Marker','.');        
        MaxNoofSpikes = max(MaxNoofSpikes,(mean(UnDirMotifNoofSpikes) + std(UnDirMotifNoofSpikes)));        
    end
end

axis([0 3 0 MaxNoofSpikes]);
set(gca,'FontSize',8,'FontWeight','bold');
set(gca,'Xtick',[1 2],'XTickLabel',[{'Directed'} {'Undirected'}]);
ylabel('Mean no of spikes','FontSize',8,'FontWeight','bold');    
title('No of spikes in motif','FontSize',8,'FontWeight','bold');    