function [] = PlotCorrelations(DirFileInfo, UnDirFileInfo, MainFigure)

figure(MainFigure);

axes('Position',[0.37 0.82 0.15 0.1]);
set(gca,'Box','off');
hold on;

if (length(DirFileInfo.GaussianCorrelations) > 0)
    plot(DirFileInfo.GaussianCorrelations(:,1),DirFileInfo.GaussianCorrelations(:,2),'ks-');
end

if (length(UnDirFileInfo.GaussianCorrelations) > 0)
    plot(UnDirFileInfo.GaussianCorrelations(:,1),UnDirFileInfo.GaussianCorrelations(:,2),'Color',[0.6 0.6 0.6],'Marker','d','MarkerEdgeColor',[0.6 0.6 0.6]);
end

axis([0 30 0 1]);
set(gca,'FontSize',8,'FontWeight','bold');
xlabel('Width of Gaussian (ms)','FontSize',8,'FontWeight','bold');
ylabel('Mean Correlation','FontSize',8,'FontWeight','bold');    
title('Pairwise Correlations','FontSize',8,'FontWeight','bold');    
