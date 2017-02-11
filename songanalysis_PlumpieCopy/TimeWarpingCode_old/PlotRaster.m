function [] = PlotRaster(Raster, MedianMotif, MainFigure, RasterPlot, LabelString)

PreDataTime = 0.05;

figure(MainFigure);
axes(RasterPlot);
hold off;
%plot(Raster(:,1),Raster(:,2),'w+');
%hold on;
%MarkerString = repmat('|',size(Raster,1),1);
%text(Raster(:,1),Raster(:,2),MarkerString,'HorizontalAlignment','center','VerticalAlignment','middle','FontSize',6,'Color','k', 'FontWeight', 'bold');
%text(Raster(:,1),Raster(:,2),MarkerString,'HorizontalAlignment','center','VerticalAlignment','middle','FontSize',6,'Color','k');
for i = 1:Raster(end,2),
    RasterX_values = Raster((find(Raster(:,2) == (i-1))),1)';
    RasterX_values = repmat(RasterX_values,2,1);
    RasterY_values = ones(size(RasterX_values,1),size(RasterX_values,2)) * (i-1);
    RasterY_values(1,:) = RasterY_values(1,:) - 0.3;
    RasterY_values(2,:) = RasterY_values(2,:) + 0.3;
    line(RasterX_values, RasterY_values,'Color','k','LineWidth',0.1);
end
axis([-PreDataTime (MedianMotif.Length) -1 (Raster(end,2) + 1)]);
set(gca, 'Box', 'off');
set(gca, 'XTick', []);
set(gca, 'XColor', 'w');
ylabel(LabelString);